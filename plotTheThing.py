#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


probData=np.loadtxt('probsTest.csv', delimiter=',', dtype=float)
probabilities, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].plot(probData[:,0], probData[:,1])
ax[1].plot(probData[:,0], probData[:,2])
ax[2].plot(probData[:,0], probData[:,3])
ax[0].set_ylabel('$P_{e e}$', fontsize=15)
ax[1].set_ylabel('$P_{\mu e}$', fontsize=15)
ax[2].set_ylabel('$P_{\\tau e}$', fontsize=15)

for i in range(3):
    ax[i].set_xscale('log')
    ax[i].set_xlabel('$E_{\\nu}$(eV)', fontsize=15)
    ax[i].set_xlim(1e5, 1e11)
plt.gcf()
plt.savefig('probPlot.png')


energyData=np.loadtxt('allowed_energies.csv', dtype=float)
enerfig=plt.figure()
plt.hist(energyData, bins=100, normed=True)
ax=plt.gca()
ax.set_yscale('log')
plt.gcf()
plt.savefig('allowed_energies_hist.png')
