#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

probData=np.loadtxt('raw_probs.csv', delimiter=',', dtype=float)
def plot(data):
	probabilities, ax = plt.subplots(2, 2, figsize=(10, 5))
	ax[0, 0].plot(probData[:,0], probData[:,1])
	ax[0, 1].plot(probData[:,0], probData[:,2])
	ax[1, 0].plot(probData[:,0], probData[:,3])
	ax[1, 1].plot(probData[:,0], probData[:,1]+ probData[:,2]+ probData[:,3] )
	ax[0, 0].set_ylabel('$P_{e e}$', fontsize=15)
	ax[0, 1].set_ylabel('$P_{\mu e}$', fontsize=15)
	ax[1, 0].set_ylabel('$P_{\\tau e}$', fontsize=15)
	ax[1, 1].set_ylabel('$P_{\\tau e}+P_{\mu e}+ P_{e e}$', fontsize=15)
	#ax[1, 1].set_ylim(1-0.01, 1+0.01)

	for i in range(2):
	    for k in range(2):
		#ax[i,k].set_xscale('log')
		ax[i,k].set_xlabel('$E_{\\nu}$(eV)', fontsize=15)
		ax[i,k].set_xlim(500, 4.5e6)
		ax[i,k].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
		ax[i,k].locator_params(axis='x', nbins=5)
	#ax[1,1].set_ylim(1-0.00001, 1+0.00001)
	ax[0,1].set_ylim(0,0.5)
	ax[1,0].set_ylim(0,0.5)
	ax[0,0].set_ylim(0,1)
	plt.tight_layout()
	return probabilities
def normalize(data):
	tot = data[:,1]+data[:,2]+data[:,3]
	data[:,1]/=tot
	data[:,2]/=tot
	data[:,3]/=tot
p1=plot(probData)
normalize(probData)
p2=plot(probData)
#plt.gcf()
print type(p1)
plt.savefig('probPlot.png', dpi=300)

newfig = plt.figure()
pathData = np.loadtxt('potentialTest.csv')
plt.plot(pathData)
plt.gcf()
#plt.gca().set_yscale('log')
plt.savefig('potentialPlot.png', dpi=300)
