#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


probData=np.loadtxt('probsTest.csv', delimiter=',', dtype=float)
fig=plt.figure()
plt.plot(probData[:,0], probData[:,1])
plt.gcf()
plt.savefig('probPlot.png')

