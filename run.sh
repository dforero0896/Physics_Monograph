#!/bin/bash

#g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
g++ -o uniandino_neutrino_exe uniandino_neutrino.cpp `gsl-config --cflags --libs`

./uniandino_neutrino_exe > probsTest.csv

python plotTheThing.py
