#!/usr/bin/env python

import numpy as np
data_W=np.loadtxt('prob_weight.dat', dtype=float)
data_P=np.loadtxt('raw_probs.csv', delimiter=',', dtype=float)
data_Pee = data_P[:,1]
avg = np.average(data_Pee[np.logical_not(np.isnan(data_Pee))])
print avg
for nan in np.argwhere(np.isnan(data_Pee)):
    data_Pee[nan]= avg
print data_Pee
avg_prob=sum(data_W*data_P[:,1])

print avg_prob
