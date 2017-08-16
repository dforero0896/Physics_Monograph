#!/bin/bash


g++ -fopenmp -o uniandino_neutrino_parallel.o uniandino_neutrino_parallel.cpp `gsl-config --cflags --libs`

./uniandino_neutrino_parallel.o > probsTest.csv

mpic++ -o spectra_sampling_MCMC.o spectra_sampling_MCMC.cpp `gsl-config --cflags --libs`

mpirun -n 2 ./spectra_sampling_MCMC.o

python plotTheThing.py
