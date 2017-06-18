#!/bin/bash

#g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
g++ -fopenmp -o uniandino_neutrino_parallel_exe uniandino_neutrino_parallel.cpp `gsl-config --cflags --libs`

./uniandino_neutrino_parallel_exe > probsTest.csv

mpiexec -n 4 python spectra_sampling_MCMC.py | sed 's/[0-9]\./\n&/g' | sed '/^\s*$/d' > allowed_energies.csv
python plotTheThing.py
