/*
g++ -fopenmp -o simulate_fluxes.o earth_simul.cpp simulate_fluxes.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>




string distributions[2]={"unif", "two_layer"};
string bse_models[3] = {"geoch", "cosmo", "geodyn"};
ofstream outfile;
int main(int argc, char const *argv[]) {
  outfile.open("flux_simulation_results.csv");
  outfile << "#Columns correspond to geoch, cosmo, geodyn BSE models\n#Rows correspond to the uniform and two_layer HPE mantle distributions.\n#There are tuples in each component that are [Uflux, Thflux]" << endl;
  for(int n=0;n<2;n++){
    outfile << distributions[n] << endl;
    for(int m=0;m<3;m++){
      string this_dist = distributions[n];
      string this_bse = bse_models[m];
      Planet *earth = new Planet();
      earth->Planet::initialize(this_dist, this_bse);
      earth->Planet::initializeFluxes(1, this_dist, this_bse);
      //outfile << "["<<earth->totalUFlux*4.07+ earth->totalThFlux*12.8<<"]" << ";";
      //outfile << "["<<earth->totalUFlux*4.07<< "&" << earth->totalThFlux*12.8<<"]" << ";";
      outfile << bse_models[m] << " ["<<earth->totalUFlux<< "," << earth->totalThFlux<<"]" << "\n";

      delete earth;
    }
    outfile << "\n\n\n" << endl;
  }
  outfile.close();
  return 0;
}
