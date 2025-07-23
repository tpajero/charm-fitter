# 1D scans ---------------------------------------------------------------------

# python3 make-plots.py -c 55 --1d -p theoretical --parfile plots/start/theoretical-wa2024.dat
# python3 make-plots.py -c 56 --1d -p theoretical --parfile plots/start/theoretical-wa2024.dat

# 2D scans ---------------------------------------------------------------------

python3 make-plots.py -c 55 --2d -p theoretical --parfile plots/start/theoretical-wa2024.dat
python3 make-plots.py -c 56 --2d -p theoretical --parfile plots/start/theoretical-wa2024.dat

# Plots ------------------------------------------------------------------------

python3 make-plots.py -c 55 56 --1d --2d -a -p theoretical -t "World average" "World average + this measurement (blinded)"
