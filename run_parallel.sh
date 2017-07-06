#!/bin/bash


g++ -fopenmp -o uniandino_neutrino_parallel.o uniandino_neutrino_parallel.cpp `gsl-config --cflags --libs`

./uniandino_neutrino_parallel.o > probsTest.csv

#mpiexec -n 4 python spectra_sampling_MCMC.py | sed 's/[0-9]\./\n&/g' | sed '/^\s*$/d' | sed '/-/g' > allowed_energies.csv
g++ -o spectra_sampling_MCMC.o spectra_sampling_MCMC.cpp `gsl-config --cflags --libs`

./spectra_sampling_MCMC.o > allowed_energies.csv

python plotTheThing.py
