from pathlib import Path

from sbmlutils.factory import *
from sbmlutils.metadata import *
from sbmlutils.examples.templates import terms_of_use
from sbmlutils.cytoscape import visualize_sbml

from panmin.model import annotations

class U(Units):
    """UnitDefinitions."""

    min = UnitDefinition("min", "min")
    s = UnitDefinition("s", "s")
    kg = UnitDefinition("kg", "kg")
    m2 = UnitDefinition("m2", "meter^2")
    mg = UnitDefinition("mg", "mg")
    ml = UnitDefinition("ml", "ml")
    mmole = UnitDefinition("mmole", "mmole")
    per_min = UnitDefinition("per_min", "1/min")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_minl = UnitDefinition("mmole_per_minl", "mmole/min/l")


# -----------------------------------------------------------------------------
# Pancreas Metabolism
# -----------------------------------------------------------------------------
mid = "pancreas_min"
version = 3

_m = Model(
    sid=f"{mid}_{version}",
    name=f"Pancreas minimal glucose model version {version}",
    notes="""
    # Pancreas minimal glucose model.
    """
    + terms_of_use,
    units=U,
    model_units=ModelUnits(
        time=U.min,
        extent=U.mmole,
        substance=U.mmole,
        length=U.meter,
        area=U.m2,
        volume=U.liter,
    ),
    creators=[
        Creator(
            familyName="Koenig",
            givenName="Matthias",
            email="koenigmx@hu-berlin.de",
            organization="Humboldt-University Berlin, Institute for Theoretical Biology",
            site="https://livermetabolism.com",
        ),
    ],
)


# volume for model (FIXME: scale to beta-cell)
pancreas_volume = 0.5  # [L]
species_in_amounts = False

_m.compartments = [
    Compartment(
        "Vpa",
        value=pancreas_volume,
        unit=U.liter,
        constant=True,
        name="pancreas tissue",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["pa"],
    ),
    Compartment(
        "Vext", value=5.0, unit=U.liter, constant=True, name="pancreas blood", port=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["blood"],
    ),
    Compartment(
        "Vmem",
        value=1.0,
        unit=U.m2,
        constant=True,
        name="pancreas plasma membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma membrane"],
    ),
]

_m.species = [
    Species(
        "Aext_glc",
        compartment="Vext",
        initialConcentration=5.0,
        substanceUnit=U.mmole,
        name="glucose",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        port=True,
        boundaryCondition=True,
        annotations=annotations.species["glc"],
    ),
    Species(
        "Aext_lac",
        compartment="Vext",
        initialConcentration=0.8,
        substanceUnit=U.mmole,
        name="lactate",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        port=True,
        annotations=annotations.species["lac"],
    ),
    Species(
        "Aext_ins",
        compartment="Vext",
        initialConcentration=60e-9,
        substanceUnit=U.mmole,
        name="insulin",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.MACROMOLECULE,
        port=True,
        annotations=annotations.species["ins"],
    ),
    Species(
        "Aext_cpep",
        compartment="Vext",
        initialConcentration=0,
        substanceUnit=U.mmole,
        name="c-peptide",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.MACROMOLECULE,
        port=True,
        annotations=annotations.species["cpep"],
    ),
    Species(
        "Apa_glc",
        compartment="Vpa",
        initialConcentration=5.0,
        substanceUnit=U.mmole,
        name="glucose",
        sboTerm=SBO.SIMPLE_CHEMICAL,
        hasOnlySubstanceUnits=True,
        annotations=annotations.species["glc"],
    ),
    Species(
        "Apa_lac",
        compartment="Vpa",
        initialConcentration=0.8,
        substanceUnit=U.mmole,
        name="lactate",
        sboTerm=SBO.SIMPLE_CHEMICAL,
        hasOnlySubstanceUnits=True,
    annotations=annotations.species["lac"],
    ),
]

rules = []
if species_in_amounts:
    # concentration rules
    for s in _m.species:
        rules.append(
            AssignmentRule(
                f"C{s.sid[1:]}",
                f"{s.sid}/{s.compartment}",
                U.mM,
                name=f"{s.name} concentration ({s.compartment})",
            ),
        )

_m.reactions = [
    Reaction(
        sid="GLCIM",
        name="glucose import",
        equation="Aext_glc <-> Apa_glc" if species_in_amounts else "Cext_glc <-> Cpa_glc",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter("GLCIM_Vmax", 100.0, U.mmole_per_minl,
                      name="Vmax glucose import", sboTerm=SBO.MAXIMAL_VELOCITY),
            Parameter("GLCIM_Km", 1.0, U.mM,
            name="Km glucose import", sboTerm=SBO.MICHAELIS_CONSTANT
                      ),
        ],
        rules=[],
        formula=(
            "Vpa * GLCIM_Vmax/GLCIM_Km * (Cext_glc-Cpa_glc)/(1 dimensionless + Cext_glc/GLCIM_Km + Cpa_glc/GLCIM_Km)",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="LACEX",
        name="lactate export",
        equation="Apa_lac <-> Aext_lac"
        if species_in_amounts
        else "Cpa_lac <-> Cext_lac",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter("LACEX_Vmax", 100.0, U.mmole_per_minl, name="Vmax lactate export",
                      sboTerm=SBO.MAXIMAL_VELOCITY),
            Parameter("LACEX_Km", 0.5, U.mM,
                      name="Km lactate export", sboTerm=SBO.MICHAELIS_CONSTANT),
        ],
        rules=[],
        formula=(
            "Vpa * LACEX_Vmax/LACEX_Km * (Cpa_lac-Cext_lac)/(1 dimensionless + Cext_lac/LACEX_Km + Cpa_lac/LACEX_Km)",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="GLC2LAC",
        name="glycolysis",
        equation="Apa_glc -> 2 Apa_lac"
        if species_in_amounts
        else "Cpa_glc -> 2 Cpa_lac",
        compartment="Vpa",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "GLC2LAC_Vmax", 0.1, U.mmole_per_minl, name="Vmax glucose utilization (glycolysis)",
                sboTerm=SBO.MAXIMAL_VELOCITY
            ),
            Parameter("GLC2LAC_Km", 4.5, U.mM,
                      name="Km effective glucose utilization (glycolysis)", sboTerm=SBO.MICHAELIS_CONSTANT
                      ),
        ],
        rules=[],
        formula=(
            "Vpa * GLC2LAC_Vmax * (Cpa_glc/(Cpa_glc + GLC2LAC_Km))",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="IRS",
        name="IRS insulin secretion",
        equation="-> Aext_ins + Aext_cpep [Apa_glc]"
        if species_in_amounts
        else "-> Cext_ins + Cext_cpep [Cpa_glc]",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "IRS_Vmax",
                1.6e-6,
                U.mmole_per_minl,  # 40/1000/60
                name="Vmax insulin secretion",
                sboTerm=SBO.MAXIMAL_VELOCITY
            ),
            Parameter("IRS_n_glc", 4, U.dimensionless,
                      name="Hill coeffient glucose in insulin secretion", sboTerm=SBO.HILL_COEFFICIENT),
            Parameter("IRS_Km_glc", 7.0, U.mM,
                      name="Km glucose insulin secretion", sboTerm=SBO.MICHAELIS_CONSTANT),
        ],
        rules=[],
        formula=(
            "Vpa * IRS_Vmax * power(Cpa_glc, IRS_n_glc) / "
            "(power(Cpa_glc, IRS_n_glc) + power(IRS_Km_glc, IRS_n_glc))",
            U.mmole_per_min,
        ),
        notes="""
        # Insulin secretion
        - equimolar secretion of insulin and c-peptide
        - hill kinetics of secretion depending on glucose
        """
    ),
]

if not species_in_amounts:
    # replace species ids
    replacements = {"Aext": "Cext", "Apa": "Cpa"}
    for s in _m.species:
        sid_new = s.sid
        for key, value in replacements.items():
            sid_new = sid_new.replace(key, value)
        s.sid = sid_new
        s.hasOnlySubstanceUnits = False

pancreas_model = _m


def create_pancreas_model() -> FactoryResult:
    """Creates minimal pancreas model

    :return: path to SBML file
    """
    output_dir = Path(__file__).parent.parent.parent / "models" / f"v{version}"
    output_dir.mkdir(parents=True, exist_ok=True)

    results: FactoryResult = create_model(
        models=_m,
        output_dir=output_dir,
        annotations=None,
        create_report=True,
    )

    return results


if __name__ == "__main__":
    results = create_pancreas_model()
    print(results.sbml_path)
    visualize_sbml(results.sbml_path)
