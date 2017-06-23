#!/bin/bash

#g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
g++ -o uniandino_neutrino.o uniandino_neutrino.cpp `gsl-config --cflags --libs`

./uniandino_neutrino.o > probsTest.csv

python plotTheThing.py
