# -*- coding=utf-8 -*-
import os
from copy import deepcopy

from sbmlutils.modelcreator import creator
from sbmlutils.units import *
from sbmlutils.annotation.sbo import *
from sbmlutils.factory import *

from panmin.model import templates


# -----------------------------------------------------------------------------
# Pancreas Metabolism
# -----------------------------------------------------------------------------
mid = 'pancreas_min'
version = 1
# -----------------------------------------------------------------------------
notes = templates.notes(tissue="Pancreas")
creators = templates.CREATORS
model_units = templates.MODEL_UNITS

# --- units ---
units = deepcopy(templates.UNITS)

# volume for model (FIXME: scale to beta-cell)
pancreas_volume = 0.5  # [L]

# --- compartments ---
compartments = [
    Compartment('Vpa', value=pancreas_volume, unit=UNIT_KIND_LITRE, constant=False, name='pancreas tissue', port=True),
    Compartment('Vext', value=0.1, unit=UNIT_KIND_LITRE, constant=False, name='pancreas blood', port=True),
]

# --- species ---
species = [
    Species('Aext_glc', compartment="Vext", initialConcentration=5.0, substanceUnit='mmole',
               name="glucose", hasOnlySubstanceUnits=True, port=True, boundaryCondition=True),
    Species('Aext_lac', compartment="Vext", initialConcentration=0.8, substanceUnit='mmole',
               name="lactate", hasOnlySubstanceUnits=True, port=True),
    Species('Aext_ins', compartment="Vext", initialConcentration=60E-6, substanceUnit='mmole',
               name="insulin", hasOnlySubstanceUnits=True, port=True),
    Species('Aext_cpep', compartment="Vext", initialConcentration=0, substanceUnit='mmole',
               name="c-peptide", hasOnlySubstanceUnits=True, port=True),


    Species('Apa_glc', compartment="Vpa", initialConcentration=5.0, substanceUnit='mmole',
               name="glucose", hasOnlySubstanceUnits=True),
    Species('Apa_lac', compartment="Vpa", initialConcentration=0.8, substanceUnit='mmole',
               name="lactate", hasOnlySubstanceUnits=True),
    # Species('Apa_atp', compartment="Vpa", initialConcentration=3.0, substanceUnit='mmole',
    #           name="ATP", hasOnlySubstanceUnits=True),
    # Species('Apa_adp', compartment="Vpa", initialConcentration=3.0, substanceUnit='mmole',
    #           name="ADP", hasOnlySubstanceUnits=True),
]

# --- assignment rules
rules = []
for s in species:
    # concentration rules
    rules.append(
        AssignmentRule(f'C{s.sid[1:]}', f'{s.sid}/{s.compartment}',
                          'mM', name=f'{s.name} concentration ({s.compartment})'),
    )

# --- reactions ---
reactions = [

    Reaction(
        sid="GLCIM",
        name="glucose import",
        equation="Aext_glc <-> Apa_glc",
        compartment='Vpa',
        sboTerm=SBO_TRANSPORT_REACTION,
        pars=[
            Parameter('GLCIM_Vmax', 100.0, 'mmole_per_minl',
                         name='Glucose import'),
            Parameter('GLCIM_Km', 1.0, 'mM'),
        ],
        rules=[],

        formula=("Vpa * GLCIM_Vmax/GLCIM_Km * (Cext_glc-Cpa_glc)/(1 dimensionless + Cext_glc/GLCIM_Km + Cpa_glc/GLCIM_Km)", 'mmole_per_min'),
    ),

    Reaction(
        sid="LACEX",
        name="lactate export",
        equation="Apa_lac <-> Aext_lac",
        compartment='Vpa',
        sboTerm=SBO_TRANSPORT_REACTION,
        pars=[
            Parameter('LACEX_Vmax', 100.0, 'mmole_per_minl',
                         name='Lactate import'),
            Parameter('LACEX_Km', 0.5, 'mM'),
        ],
        rules=[],

        formula=(
        "Vpa * LACEX_Vmax/LACEX_Km * (Cpa_lac-Cext_lac)/(1 dimensionless + Cext_lac/LACEX_Km + Cpa_lac/LACEX_Km)", 'mmole_per_min'),
    ),

    Reaction(
        sid="GLC2LAC",
        name="glycolysis",
        equation="Apa_glc -> 2 Apa_lac",
        compartment='Vpa',
        sboTerm=SBO_BIOCHEMICAL_REACTION,
        pars=[
            Parameter('GLC2LAC_Vmax', 0.1, 'mmole_per_minl',
                         name='Glucose utilization'),
            Parameter('GLC2LAC_Km', 4.5, 'mM'),  # effective Km of glycolysis
        ],
        rules=[],

        formula=("Vpa * GLC2LAC_Vmax * (Cpa_glc/(Cpa_glc + GLC2LAC_Km))", 'mmole_per_min'),
    ),

    # Insulin secretion
    # - equimolar secretion of insulin and c-peptide
    # - hill kinetics of secretion depending on glucose
    Reaction(
        sid="IRS",
        name="IRS insulin secretion",
        equation="-> Aext_ins + Aext_cpep [Apa_glc]",
        compartment='Vpa',
        pars=[
            Parameter('IRS_Vmax', 1.6E-6, 'mmole_per_minl',  # 40/1000/60
                      name='Insulin secretion'),
            Parameter('IRS_n_glc', 4, 'dimensionless'),
            Parameter('IRS_Km_glc', 7.0, 'mM'),
        ],
        rules=[],

        formula=(
            "Vpa * IRS_Vmax * power(Cpa_glc, IRS_n_glc) / "
            "(power(Cpa_glc, IRS_n_glc) + power(IRS_Km_glc, IRS_n_glc))",
            'mmole_per_min'),
    ),
]


def create_model():
    """ Creates minimal pancreas model

    :return: path to SBML file
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))

    target_dir = os.path.abspath(
        os.path.join(base_dir, "../../models")
    )
    _, _, sbml_path = creator.create_model(
        modules=['pancreas_model'],
        target_dir=target_dir,
        annotations=None,
        create_report=True
    )

    return sbml_path


if __name__ == "__main__":
    create_model()
