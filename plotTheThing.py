#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


probData=np.loadtxt('probsTest.csv', delimiter=',', dtype=float)
fig=plt.figure()
plt.plot(probData[:,0], probData[:,1], 'or')
ax=plt.gca()
ax.set_xscale('log')
plt.xlabel('$E_{\\nu}$(eV)', fontsize=15)
plt.ylabel('$P_{\mu e}$', fontsize=15)
plt.gcf()
plt.savefig('probPlot.png')


energyData=np.loadtxt('allowed_energies.csv', dtype=float)
enerfig=plt.figure()
plt.hist(energyData, bins=100, normed=True)
ax=plt.gca()
ax.set_yscale('log')
plt.gcf()
plt.savefig('allowed_energies_hist.png')
