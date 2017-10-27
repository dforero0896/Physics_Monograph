/*
g++ -fopenmp -o simulate_fluxes.o earth_simul.cpp simulate_fluxes.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>




string distributions[2]={"unif", "two_layer"};
string bse_models[3] = {"geoch", "geodyn", "cosmo"};
ofstream outfile;
int main(int argc, char const *argv[]) {
  outfile.open("flux_simulation_results.csv");
  outfile << "#Columns correspond to geoch, geodyn, cosmo BSE models\n#Rows correspond to the uniform and two_layer HPE mantle distributions.\n#There are tuples in each component that are [Uflux, Thflux]" << endl;
  for(int n=0;n<2;n++){
    for(int m=0;m<3;m++){
      string this_dist = distributions[n];
      string this_bse = bse_models[m];
      Planet *earth = new Planet();
      earth->Planet::initialize(this_dist, this_bse);
      outfile << "["<<earth->totalUFlux<<","<<earth->totalThFlux<<"]" << ";";
      delete earth;
    }
    outfile << "[0,0]" << endl;
  }
  outfile.close();
  return 0;
}
