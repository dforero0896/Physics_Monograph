#!/bin/bash

#g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
g++ -fopenmp -o uniandino_neutrino_parallel_exe uniandino_neutrino_parallel.cpp `gsl-config --cflags --libs`

./uniandino_neutrino_parallel_exe > probsTest.csv

python plotTheThing.py
