/*
g++ -fopenmp -o test_simul.o earth_simul.cpp test_simul.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>
#include <cmath>
using namespace std;

int main(int argc, char const *argv[]) {
  int t_i=1;int t_k = 990;
  Planet *earth = new Planet();
  string hpe_dist = "unif";
  string bse_model = "geoch";
  earth->Planet::initialize(hpe_dist, bse_model);
  earth->Planet::initializePaths(0,t_i,t_k);
  //cout << "distance " << earth->asArray[t_i][t_k].path[1] << endl;
  cout << "dstance to detector " << earth->asArray[t_i][t_k].distanceToDetector << endl;
  cout << "path len " << earth->asArray[t_i][t_k].pathLen << endl;
  cout << "total flux " << earth->totalFlux << endl;
  //cout << "crust mass " << earth->crustMass << endl;
  //cout << "mantle mass " << earth->mantleMass << endl;
  earth->Planet::initializeFluxes(1, hpe_dist, bse_model);
  cout << "total oscillated flux " << earth->totalFlux << endl;
  cout << "total oscillated Th flux " << earth->totalThFlux << endl;
  cout << "total oscillated U flux " << earth->totalUFlux << endl;
  ofstream outfile;
  outfile.open("earth_simul_plots.csv");
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      float TFlux = earth->asArray[i][k].neutrinoThFlux + earth->asArray[i][k].neutrinoUFlux;
      outfile << earth->asArray[i][k].isEarth*earth->asArray[i][k].neutrinoUFlux << ',' ;
      }
      outfile <<0<< endl;
    }
  outfile.close();

int xtest, ytest;
xtest=3;
ytest=998;
  ofstream test_path_file;
  test_path_file.open("test_path.csv");
  int test_N = int(earth->asArray[xtest][ytest].pathLen);
  for(int step=0;step<test_N;step++){
    test_path_file << earth->asArray[xtest][ytest].path[step] << endl;
  }
  test_path_file.close();

  long double sum_test1, sum_test2;
  sum_test2=0;
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      if(earth->asArray[i][k].isEarth){//integrate over the whole Earth
        float x= earth->asArray[i][k].x;
        float z = earth->asArray[i][k].z;
        sum_test2+=earth->asArray[i][k].volume*earth->asArray[i][k].solidAngle;
      }
      else{;}
    }
  }
  cout << "integral over whole Earth " << sum_test2 << endl;
  cout << "fraction of Earth " << sum_test2/6371 << endl;


  //cout << "an energy " << (earth->asArray[xtest][ytest].allowedEnergiesU[30]) << endl;
  //cout << "the prob " << calculateProbability(test_N, earth->asArray[xtest][ytest].path, (earth->asArray[xtest][ytest].allowedEnergiesU[30])) << endl;
  //calculateProbabilitiesFunctionEnergy(test_N, earth->asArray[xtest][ytest].path);
  //cout << "total mass " << earth->totalMass << endl;

  delete earth;
  return 0;
}
