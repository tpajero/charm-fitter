CharmFitter
============

CharmFitter is a frequentist implementation of a fitter to charm measurements of mixing and CP violation, to compute 1D and 2D confidence regions on the underlying physics parameters.
The fit implements both the phenomenological and theoretical parametrisations, see the following references for details:

  * T. Pajero, *Search for time-dependent CP violation in D0 -> K+ K- and D0 -> pi+ pi- decays*, PhD thesis, Scuola Normale Superiore, 2021, [CERN-THESIS-2020-231](https://cds.cern.ch/record/2747731);
  * T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*, [JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950).

Another useful reference is

  * A. L. Kagan and L. Silvestrini, *Dispersive and absorptive CP violation in D0/anti-D0 mixing*, [Phys. Rev. D **103** 053008, 2021](https://inspirehep.net/literature/1776611).

The fitter is based on the GammaCombo package; see [gammacombo](https://gammacombo.github.io) for details and for how to run the fitter.
The custom code for the charm fit resides in the `charmcombo` directory.
Simple example plots can be produced with the `charmcombo/make_plots.py` script.

Citation
========

If you use this software, please cite the following article:

T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*, [JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950).

Installation
============

To build CharmCombo cmake is needed in version 3.14 or higher.

Get the source code from GitHub:

    git clone https://github.com/tpajero/charm-fitter.git
    cd charm-fitter

Build the project in a `<build>` directory and install it:

    cmake -B <build>
    cmake --build <build>
    cmake --install <build>

Credits
=======

We would like to thank Nico Kleijne for a careful review of the code.
