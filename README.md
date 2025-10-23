# CharmFitter

CharmFitter is a frequentist implementation of a fitter to charm measurements of mixing and CP violation, to compute 1D
and 2D confidence regions on the underlying physics parameters.
The fit implements both the phenomenological and theoretical parametrisations, see the following references for details,

  * T. Pajero, *Search for time-dependent CP violation in D0 -> K+ K- and D0 -> pi+ pi- decays*, PhD thesis,
    Scuola Normale Superiore, 2021, [CERN-THESIS-2020-231](https://cds.cern.ch/record/2747731);
  * T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*,
    [JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950);

which are based on the very useful reference

  * A. L. Kagan and L. Silvestrini, *Dispersive and absorptive CP violation in D0/anti-D0 mixing*,
    [Phys. Rev. D **103** 053008, 2021](https://inspirehep.net/literature/1776611).

The fitter is a submodule of the GammaCombo package; see [gammacombo](https://gammacombo.github.io) for details on the
gammacombo framework and on how to run the fitter.
Simple example plots can be produced with the `make-plots.py` script.

Additionally, the folder `blue` contains CPP scripts to perform single-observable averages (e.g. for DeltaY and for CP
violation in D0 -> KS KS decays), based on the [BLUE](https://blue.hepforge.org/) combiner, and Python scripts to plot
the results.

# Citation

Please consider citing the following article if you use this software:

T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*,
[JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950).


# Installation

To build CharmCombo cmake is needed in version 3.14 or higher.

Get the source code from GitHub:

    git clone https://github.com/gammacombo/gammacombo.git
    cd gammacombo
    git checkout for-charm-fitter
    cd gammacombo
    git clone https://github.com/tpajero/charm-fitter.git

Then follow the instructions in the `README.md` of the main gammacombo project for building and installing.


# Credits

We would like to thank Nico Kleijne for a careful review of the code.
