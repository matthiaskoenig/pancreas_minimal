"""
Basic simulations with pancreas model.
"""
# TODO: basic simulations and parameter scan
import os
from collections import OrderedDict
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import sbmlsim
from sbmlsim import load_model, timecourse, Timecourse

from panmin.model.pancreas_model import create_model, species_in_amounts


# TODO: simulation experiment for this
def simulate_glucose_dependency(r, model_name):
    """ Scans the external glucose. """
    glc_vec = np.linspace(0.5, 20, num=20)  # [mM]

    # simulation
    results = []
    for glc in glc_vec:
        if species_in_amounts:
            changes = {'[Aext_glc]': glc}
        else:
            changes = {'[Cext_glc]': glc}

        s = timecourse(r, Timecourse(start=0, end=10, steps=200, changes=changes))
        results.append(s)

    results_path = "../../models"

    # --------------------------
    # Time course
    # --------------------------
    f1 = analysis_plot1(results, xid="time", xlabel="time [min]")
    f1.savefig(f"{results_path}/{model_name}_timecourse.png", bbox_inches="tight")
    plt.show()

    # --------------------------
    # Dose response
    # --------------------------
    s_scan = collect_scan(results, dose_vec=glc_vec)
    f2 = analysis_plot1(s_scan, xid="index", xlabel="Glucose [mM]")
    f2.savefig(f"{results_path}/{model_name}_glcscan.png", bbox_inches="tight")

    plt.show()


def collect_scan(results, dose_vec):
    """ Collect scan data.

    :param dose_vec:
    :param results:
    :return:
    """

    columns = results[0].columns
    scan_data = OrderedDict()
    for col_name in columns:
        vec = np.empty_like(dose_vec)
        for k, s in enumerate(results):
            # last data point (often steady state)
            vec[k] = s[col_name].values[-1]
        scan_data[col_name] = vec

    df_scan = pd.DataFrame(data=scan_data, index=dose_vec)

    return df_scan


def analysis_plot1(results, xid, xlabel):
    """ Plots the timecourses.

    :return:
    """
    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, figsize=(14, 12))
    f.subplots_adjust(wspace=.3, hspace=.3)
    axes = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12)

    kwargs = {
        "linestyle": '-',
        "marker": 'None'
    }

    if not isinstance(results, list):
        # scan or single simulation
        results = [results]
        kwargs["marker"] = 'o'

    for s in results:
        if xid == "index":
            x = s.index
        else:
            x = s[xid]

        if species_in_amounts:
            # glucose
            ax1.plot(x, s.Cpa_glc, color="darkblue", label="Cpa_glc", **kwargs)
            ax1.plot(x, s.Cext_glc, color="black", label="Cext_glc", **kwargs)
            # lactate
            ax2.plot(x, s.Cpa_lac, color="darkblue", label="Cpa_lac", **kwargs)
            ax2.plot(x, s.Cext_lac, color="black", label="Cext_lac", **kwargs)
            # insulin
            ax3.plot(x, s.Cext_ins*1E9, color="black", label="Cext_ins", **kwargs)
            # c-peptide
            ax4.plot(x, s.Cext_cpep*1E9, color="black", label="Cext_cpep", **kwargs)

        else:
            # glucose
            ax1.plot(x, s["[Cpa_glc]"], color="darkblue", label="Cpa_glc", **kwargs)
            ax1.plot(x, s["[Cext_glc]"], color="black", label="Cext_glc", **kwargs)
            # lactate
            ax2.plot(x, s["[Cpa_lac]"], color="darkblue", label="Cpa_lac", **kwargs)
            ax2.plot(x, s["[Cext_lac]"], color="black", label="Cext_lac", **kwargs)
            # insulin
            ax3.plot(x, s["[Cext_ins]"]*1E9, color="black", label="Cext_ins", **kwargs)
            # c-peptide
            ax4.plot(x, s["[Cext_cpep]"]*1E9, color="black", label="Cext_cpep", **kwargs)

        # rates
        ax5.plot(x, s.GLCIM, color="black", **kwargs)
        ax5.plot(x, s.GLC2LAC, color="red", **kwargs)

        ax6.plot(x, s.LACEX, color="black", **kwargs)
        ax7.plot(x, s.IRS*1E9, color="black", **kwargs)
        ax8.plot(x, s.IRS*1E9, color="black", **kwargs)

        # rates scaled (volume)
        Vpa = s.Vpa.values[0]*1000  # [ml] pancreas model volume
        ax9.plot(x, s.GLCIM/Vpa*1E6, color="black", **kwargs)
        ax10.plot(x, s.LACEX/Vpa*1E6, color="black", **kwargs)
        ax11.plot(x, s.IRS/Vpa*1E9, color="black", **kwargs)
        ax12.plot(x, s.IRS/Vpa*1E9, color="black", **kwargs)

    # rates concentration
    ax1.set_title("Glucose")
    ax1.set_ylabel("Glucose [mM]")
    ax2.set_title("Lactate")
    ax2.set_ylabel("Lactate [mM]")
    ax3.set_title("Insulin")
    ax3.set_ylabel("Insulin [pmole/l]")
    ax4.set_title("C-peptide")
    ax4.set_ylabel("C-peptide [pmole/l]")

    # rates
    ax5.set_title("Glucose import (model)")
    ax5.set_ylabel("GLCIM [mmole/min]")
    ax6.set_title("Lactate export (model)")
    ax6.set_ylabel("LACEX [mmole/min]")
    ax7.set_title("Insulin secretion (model)")
    ax7.set_ylabel("IRS [pmole/min]")
    ax8.set_title("C-peptide secretion (model)")
    ax8.set_ylabel("IRS [pmole/min]")

    # rates (volume)
    ax9.set_title("Glucose import (volume)")
    ax9.set_ylabel("GLCIM [nmole/min/ml]")
    ax10.set_title("Lactate export (volume)")
    ax10.set_ylabel("LACEX [nmole/min/ml]")
    ax11.set_title("Insulin secretion (volume)")
    ax11.set_ylabel("IRS [pmole/min/ml]")
    ax12.set_title("C-peptide secretion (volume)")
    ax12.set_ylabel("IRS [pmole/min/ml]")

    for ax in axes:
        ax.set_xlabel(xlabel)

    # only positive rates and concentration
    for ax in axes:
        ax.set_ylim(bottom=0)
        ylim = ax.get_ylim()
        ax.set_ylim(top=ylim[1]*1.02)

    return f


if __name__ == "__main__":
    # create latest model version
    sbml_path = create_model()
    model_name = os.path.splitext(os.path.basename(sbml_path))[0]
    print(model_name)


    # load model
    r = load_model(sbml_path)
    integrator_kwargs = {
        "relative_tolerance": 1E-8,
        "absolute_tolerance": 1E-8,
    }
    integrator = sbmlsim.model.set_integrator_settings(r, **integrator_kwargs)
    print(integrator)

    # run simulation
    simulate_glucose_dependency(r, model_name=model_name)
