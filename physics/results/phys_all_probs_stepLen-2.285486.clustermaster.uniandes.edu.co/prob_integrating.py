#!/usr/bin/env python

import numpy as np
import subprocess
import os




nodes_to_consider = np.loadtxt('planet_coords.csv', delimiter=' ', dtype=int)
iteration = 0

Emin = str(0.0006)
Emax = str(4.5004)
Steps_in_energy = str(100)
distribs = ['unif', 'two_layer']
bse_models = ['geoch', 'cosmo', 'geodyn']
for hpe_dist in distribs:
	for bse_model in bse_models:
		print "distrib: ", distrib, ", model: ", bse_model
		filename = "probability_planet_%s_%s.csv"%(hpe_dist,bse_model)
		final_file=open(filename, "w")
		for node in nodes_to_consider:
		    i_c = node[0]
		    k_c = node[1]
		    if(i_c==0):
			continue
		    command_w=[Emin,Emax,Steps_in_energy,hpe_dist,bse_model,str(i_c),str(k_c)]
		    command_p=list(command_w)
		    command_p.append(str(iteration))
		    if(iteration==0):
			os.system('g++ -fopenmp -o prob_weight.o earth_simul.cpp prob_weight.cpp `gsl-config --cflags --libs`')
			os.system('g++ -fopenmp -o raw_probs.o uandino.cpp earth_simul.cpp raw_probs.cpp `gsl-config --cflags --libs`')
		    print command_w
		    command_w.insert(0,'./prob_weight.o')
		    command_p.insert(0,'./raw_probs.o')
		    subprocess.call(command_w)
		#Run raw probabilities once since it doesn't depend on distribution nor bse model.
		    if(hpe_dist == distribs[0] and bse_model==bse_models[0]):
		    	subprocess.call(command_p)
		    data_W=np.loadtxt('prob_weight.dat', dtype=float)
		    data_P=np.loadtxt('raw_probs.csv', delimiter=',', dtype=float)
		    avg_prob=sum(data_W*data_P[:,1])
		    print avg_prob
		    final_file.write("%i %i %f\n"%(i_c, k_c, avg_prob))
		    iteration+=1

		final_file.close()
