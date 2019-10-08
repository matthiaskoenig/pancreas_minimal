"""
Basic simulations with pancreas model.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sbmlsim import load_model, simulate

from panmin.model.pancreas_model import create_model


def simulate_glucose_dependency(r):
    """ Scans the external glucose. """
    r.resetToOrigin()
    glc_vec = np.linspace(3.0, 20, num=10)  # [mM]

    results = []
    for glc in glc_vec:
        r.resetToOrigin()
        r.setValue('[Aext_glc]', glc)
        s = simulate(r, start=0, end=100, steps=200)  # [min]
        results.append(s)

    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(14, 8))
    f.subplots_adjust(wspace=.3, hspace=.3)
    axes = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)

    for ax in axes:
        ax.set_xlabel("time [min]")

    for s in results:
        ax1.plot(s.time, s.Cpa_glc, '-', color="darkblue", label="Cpa_glc")
        ax1.plot(s.time, s.Cext_glc, '-', color="black", label="Cext_glc")

        # ax2.plot(s.time, s.IRS*1000, '-', color="k")
        # ax2.plot(s.time, s.Cext_cpep*1E6, '-', color="k")
        # ax3.plot(s.time, s.Cext_ins*1E6, '-', color="k")

    for ax in (ax1, ax2, ax3):
        ax.set_ylim(bottom=0)

    ax1.set_ylabel("Glucose [mM]")
    ax2.set_ylabel("Lactate [mM]")

    ax3.set_ylabel("IRS (insulin secretion) [nmole/min]")


    ax2.set_ylabel("C-peptide concentration [pM]")
    ax3.set_ylabel("Insulin concentration [pM]")
    plt.show()

    print(s)


if __name__ == "__main__":
    # create latest model version
    sbml_path = create_model()

    # load model
    r = load_model(sbml_path)

    # run simulation
    simulate_glucose_dependency(r)
