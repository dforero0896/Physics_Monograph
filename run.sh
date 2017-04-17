#!/bin/bash

g++ -fopenmp -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`

echo Compilation ready

./oscillations


