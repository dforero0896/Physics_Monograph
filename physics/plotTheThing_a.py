#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
'''
test_path=np.loadtxt('test_path.csv')
prob_test=np.loadtxt('prob_path.csv', delimiter=',')
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
ax[0].plot(test_path)
ax[0].set_title('Path tested', fontsize=20)
ax[0].set_xlabel('$abs(R_0-r)$ ($\\times 50km$)', fontsize=20)
ax[0].set_ylabel('Potential ($eV$)', fontsize=20)
ax[1].plot(prob_test[:,0], prob_test[:,1])
ax[1].set_xlabel('$E_\\nu (eV)$', fontsize=20)
ax[1].set_ylabel('$P_{e\\rightarrow e}$', fontsize=25)
ax[1].set_title('Survival probability for the path', fontsize=20)
ax[1].set_xscale('log')
#ax[1].set_xlim(4400e3, 4500e3)
fig.tight_layout()
plt.gcf()
plt.savefig('test_path_fig.png', dpi=500)



probData=np.loadtxt('probsTest.csv', delimiter=',', dtype=float)
probabilities, ax = plt.subplots(2, 2, figsize=(10, 5))
ax[0, 0].plot(probData[:,0], probData[:,1])
ax[0, 1].plot(probData[:,0], probData[:,2])
ax[1, 0].plot(probData[:,0], probData[:,3])
ax[1, 1].plot(probData[:,0], probData[:,1]+ probData[:,2]+ probData[:,3] )
ax[0, 0].set_ylabel('$P_{e e}$', fontsize=15)
ax[0, 1].set_ylabel('$P_{\mu e}$', fontsize=15)
ax[1, 0].set_ylabel('$P_{\\tau e}$', fontsize=15)
ax[1, 1].set_ylabel('$P_{\\tau e}+P_{\mu e}+ P_{e e}$', fontsize=15)
ax[1, 1].set_ylim(1-0.001, 1+0.001)

for i in range(2):
    for k in range(2):
        ax[i,k].set_xscale('log')
        ax[i,k].set_xlabel('$E_{\\nu}$(eV)', fontsize=15)
        ax[i,k].set_xlim(1e9, 1e10)
plt.tight_layout()
plt.gcf()
plt.savefig('probPlot.png', dpi=500)


comparison, ax = plt.subplots(1, 1)
ax.plot(probData[:,0], probData[:,2])
ax.set_ylabel('$P_{\mu e}$', fontsize=15)
ax.set_xscale('log')
ax.set_xlabel('$E_{\\nu}$(eV)', fontsize=15)
ax.set_ylim(0, 0.6)
ax.set_xlim(1e9, 1e10+3e9)
plt.gcf()
plt.savefig('document/probComp.png', dpi=500)

'''
'''
energyData=np.loadtxt('energy_repo_238U.knt', dtype=float)
energyData2=np.loadtxt('energy_repo_232Th.knt', dtype=float)
enerfig=plt.figure(figsize=(15, 10))
plt.hist(energyData, bins=50, normed=True, histtype="step", label='$^{238}U$', lw=2)
plt.hist(energyData2, bins=50, normed=True, histtype="step", label='$^{232}Th$', lw=2)
plt.legend(loc=0, fontsize=15)
plt.xlabel('Neutrino energy $E_\\nu$ (MeV)', fontsize=15)
plt.ylabel('Neutrinos', fontsize=15)
plt.xlim(0, 5)
ax=plt.gca()
ax.set_yscale('log')
plt.gcf()
plt.savefig('allowed_energies_hist.png', dpi=500)
'''

earth_model=np.loadtxt('earth_simul_plots.csv', delimiter=',', dtype=float)
earthplot=plt.figure(figsize=(12, 15))
plt.imshow(earth_model, interpolation="None")
ax=plt.gca()
ax.set_xlabel('x', fontsize=25)
ax.set_ylabel('z', fontsize=25)
cb=earthplot.colorbar(ax.imshow(earth_model))
plt.scatter(0, len(earth_model[:,0]), s=100, color='k', label='Detector')
plt.legend(loc=3, fontsize=25)
ax.tick_params(axis='both', which='major', labelsize=25)
cb.ax.tick_params(labelsize=25)
cb.ax.set_title('Radius ($km$)', fontsize=25)
ax.invert_yaxis()
plt.gcf()
plt.savefig('earth_test.png', transparecy=True, dpi=500)

