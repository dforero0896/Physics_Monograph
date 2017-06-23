#!/usr/bin/env python
#Run with mpiexec -n 4 python spectra_sampling_MCMC.py

import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_size=comm.Get_size()

#Import the spectra.
data40K=np.loadtxt('../AntineutrinoSpectrum_all/AntineutrinoSpectrum_40K.knt')
data232Th=np.loadtxt('../AntineutrinoSpectrum_all/AntineutrinoSpectrum_232Th.knt')
data235U=np.loadtxt('../AntineutrinoSpectrum_all/AntineutrinoSpectrum_235U.knt')
data238U=np.loadtxt('../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt')

#Define necessary functions for Metropolis-Hastings sampling algorithm.
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def evaluate_spectrum(spectrum, x):
    index=find_nearest(spectrum[:,0]/1000, x)
    return spectrum[:,1][index]
#This function can either print or store the markov chain. For parallel computing it is easier to print.
def MH_Spectrum_Sampling(spectrum, niter):
    markov_chain=[]
    max_interval=spectrum[:,0][spectrum[:,1]!=0.][-1]/1000
    initial_value=max_interval*np.random.random()
    print initial_value
    markov_chain.append(initial_value)
    for i in range(niter):
        possible_jump=np.random.normal(markov_chain[i], 0.1)
        criteria = evaluate_spectrum(spectrum, possible_jump)/evaluate_spectrum(spectrum, markov_chain[i])
        if(criteria>=1.):
            print abs(possible_jump)
            markov_chain.append(abs(possible_jump))
        else:
            other_random = np.random.random()
            if(other_random<=criteria):
                print possible_jump
                markov_chain.append(possible_jump)
            else:
                print markov_chain[i]
                markov_chain.append(markov_chain[i])
    return np.array(markov_chain)

#Define the number of geoneutrinos to be generated.
if rank==0:
#    n_geoneutrinos=int(input("Please enter the number of Geoneutrinos to generate: "))
     n_geoneutrinos=5000

else:
    n_geoneutrinos=None
n_geoneutrinos = comm.bcast(n_geoneutrinos, root=0)

energy=MH_Spectrum_Sampling(data238U, int(n_geoneutrinos/world_size))
'''
#Define function that calculates how many geoneutrinos are produced by isotope of a given concentration, half_life, geoneutrinos per decy and in time T.
def calculate_geoneutrino_number(concentration, half_life, geo_nu_per_decay, T=1.):
    return geo_nu_per_decay*(1-2**-(T/half_life))*concentration

#Define concentration of Uranium in the Earth given a BSE model.
U_conc=1.
#List to hold the spectra
isotopes=[data238U, data235U, data232Th, data40K]
#Function that calculates how many of n_geoneutrinos correspond to each decay chain.
def calculate_weights(n_geoneutrinos):
    produced_geonu=np.array([calculate_geoneutrino_number(0.99274*U_conc, 4.47, 6.),calculate_geoneutrino_number(0.720e-2*U_conc, 0.7, 4.),calculate_geoneutrino_number(4.*U_conc, 14, 4.),calculate_geoneutrino_number(77350.*U_conc, 1.28, 2.)])
    return (n_geoneutrinos*produced_geonu/sum(produced_geonu)).astype(int)
weights= calculate_weights(n_geoneutrinos)

if rank<3:
    e_r02=MH_Spectrum_Sampling(isotopes[rank], weights[rank])


new_weights=int(weights[3]/world_size)
e_r3N=MH_Spectrum_Sampling(isotopes[3], new_weights)
'''
