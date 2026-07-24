# CharmFitter

A global analysis of charm measurements of mixing and CP violation, interpreting the experimental observables in terms
the the relevant theory parameters.
Both the theoretical or phenomenological parametrisations can be used, adopting the formulas in the following
references,

  * T. Pajero, *Search for time-dependent CP violation in D0 -> K+ K- and D0 -> pi+ pi- decays*, PhD thesis,
    Scuola Normale Superiore, 2021, Appendix B [CERN-THESIS-2020-231](https://cds.cern.ch/record/2747731);
  * T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*,
    [JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950);

which are mostly based on the article

  * A. L. Kagan and L. Silvestrini, *Dispersive and absorptive CP violation in D0/anti-D0 mixing*,
    [Phys. Rev. D **103** 053008, 2021](https://inspirehep.net/literature/1776611).

The impact of using different assumptions on the universality of DeltaY(D0 -> h- h+) can be tested, too.

The statistical treatment is frequentist and relies on the GammaCombo package (see
[gammacombo](https://gammacombo.github.io) for details).

Additionally, the folder `scripts/BLUE` contains C++ scripts to perform single-observable averages (e.g. for DeltaY, CP
violation in D0 -> KS KS decays, and D(s)+ -> eta(')h+ branching fractions), based on the
[BLUE](https://blue.hepforge.org/) combiner, and Python scripts to plot the results.

## Install

Get the source code from GitHub, add the charm-fitter submodule and build the project:

    git clone https://github.com/gammacombo/gammacombo.git core
    cd core
    git checkout charm-fitter
    git submodule init
    git submodule update
    # The following command is optional if you need to setup a software environment on a machine with CVMFS installed
    source scripts/setup-env-cvmfs.sh
    cmake -B <build-dir>
    cmake --build <build-dir> [-j <n-cores>]
    cmake --install <build-dir>
    python -m pip install -e . -e charm-fitter

If you have sourced the script to setup the software environment, you will need to do that in each new shell session
before running any executables.

## Usage

### Main combiner

One- and two-dimensional confidence intervals for the theory parameters of interest can be plotted with
[Matplotlib](https://matplotlib.org/) using the script [scripts/mpl-plots.py](scripts/mpl-plots.py), while the list of
parameters to plot, their ranges and the axes titles are set in the configuration files in [config](config).
For example, you can plot the results of the world average of 2025 with

    python scripts/mpl-plots.py --config config/2025.py -a all --rescan --extra-opts " --pr --ps 1 "

To see more options, please run `python scripts/mpl-plots.py -h`.

Please see the [GammaCombo manual](https://gammacombo.github.io/manual.pdf) for instructions on how to add new
measurements to the combination.

### BLUE combinations

Please refer to the anotations at the beginning of C++ and Python files in the [scripts/BLUE](scripts/BLUE) folder.

## Maintainers

[@tpajero](tommaso.pajero@cern.ch)

## Contributors

We would like to thank Nico Kleijne for a careful review of the code.

## Citation

Please consider citing the following article if you use this software:

T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*,
[JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950).
