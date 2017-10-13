/*
g++ -fopenmp -o test_simul.o earth_simul.cpp test_simul.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>
using namespace std;
int N =1000;

int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->Planet::initialize("two_layer", "cosmo");
  cout << "total flux " << earth->totalFlux << endl;
  //cout << "crust mass " << earth->crustMass << endl;
  //cout << "mantle mass " << earth->mantleMass << endl;
  earth->Planet::initializeFluxes(1);
  cout << "total oscillated flux " << earth->totalFlux << endl;
  ofstream outfile;
  outfile.open("earth_simul_plots.csv");
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      outfile << earth->asArray[i][k].isEarth*earth->asArray[i][k].neutrinoFlux << ',' ;
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
      if(earth->asArray[i][k].isMantle){
        float x= earth->asArray[i][k].x;
        float z = earth->asArray[i][k].z;
        sum_test2+=earth->asArray[i][k].volume*earth->asArray[i][k].solidAngle;
      }
      else{;}
    }
  }
  cout << "variables test " << sum_test2 << endl;
  sum_test1=0;

for(int n=0;n<10000;n++){
  sum_test1+=(1.-pow(0.45, 2.*n+3.))/((2.*n+1.)*(2.*n+3.));
}
cout << "legendre pols test  " << sum_test1 << endl;

  //cout << "an energy " << (earth->asArray[xtest][ytest].allowedEnergiesU[30]) << endl;
  //cout << "the prob " << calculateProbability(test_N, earth->asArray[xtest][ytest].path, (earth->asArray[xtest][ytest].allowedEnergiesU[30])) << endl;
  //calculateProbabilitiesFunctionEnergy(test_N, earth->asArray[xtest][ytest].path);
  //cout << "total mass " << earth->totalMass << endl;

  delete earth;
  return 0;
}