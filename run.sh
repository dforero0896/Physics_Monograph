#!/bin/bash

#g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
g++ -o oscillations uniandino_neutrino.cpp `gsl-config --cflags --libs`

./oscillations > probsTest.csv

python plotTheThing.py


