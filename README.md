# CharmFitter

Global combination of charm measurements of mixing and CP violation, interpreting the experimental observables in terms
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

The impact of using different assumptions for the final-state corrections for DeltaY(D0 -> h- h+) can be tested, too.

The statistical treatment is frequentist and relies on the GammaCombo package (see
[gammacombo](https://gammacombo.github.io) for details).

Additionally, the folder `BLUE/main` contains executables to perform single-observable averages of charm quantities
(e.g. for DeltaY, CP violation in D0 -> KS KS decays, D(s)+ -> eta(')h+ branching fractions and CP asymmetries, yCP
and yCP - yCP(RS)).
The executables are based on the [BLUE](https://blue.hepforge.org/) package, and are direcly linked against a local
copy of its library.

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
[Matplotlib](https://matplotlib.org/) using the script [scripts/charm-combo.py](scripts/charm-combo.py),
while the list of parameters to plot, their ranges and the axes titles are set in the configuration files
in [config/plotting](config/plotting).
For example, you can plot the results of the world average of 2025 with

    python scripts/charm-combo.py --config config/plotting/charm-2025.py -a all --rescan

To see all available options, run `python scripts/charm-combo.py -h`.
Analogous plots for the subset of WS/RS D0 -> Kpi measurements can be obtained through a sibling script,
see `python scripts/ws-combo.py -h`.

Please refer to the [GammaCombo manual](https://gammacombo.github.io/manual.pdf) for instructions on how to add new
measurements to the combination.

### BLUE combinations

The BLUE combinations can be performed using the executables built from [BLUE/main](BLUE/main), e.g.

    bin/BLUE/<exec> <combination-id>

The executables support the `-h` option, to list which combinations are supported.
The results of the combinations can be plotted using the Python scripts in [BLUE/scripts](BLUE/scripts)
(run them with `-h` to explore available options).

## Maintainers

[@tpajero](tommaso.pajero@cern.ch)

## Contributors

We would like to thank Nico Kleijne for a careful review of the code.

## Citation

Please cite the following article if you use this software:

T. Pajero and M. J. Morello, *Mixing and CP violation in D0 -> K- pi+ decays*,
[JHEP **03** (2022) 162](https://inspirehep.net/literature/1866950).
