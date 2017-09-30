/*
g++ -fopenmp -o test_simul.o earth_simul.cpp test_simul.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>
using namespace std;
int N =1000;
int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->Planet::initialize("two_layer", "geodyn");
  cout << "total flux " << earth->totalFlux << endl;
  ofstream outfile;
  outfile.open("earth_simul_plots.csv");
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      outfile << earth->asArray[i][k].neutrinoFlux << ',' ;
      }
      outfile << 0 << endl;
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
  cout << "an energy " << (earth->asArray[xtest][ytest].allowedEnergiesU[30]) << endl;
  //cout << "the prob " << calculateProbability(test_N, earth->asArray[xtest][ytest].path, (earth->asArray[xtest][ytest].allowedEnergiesU[30])) << endl;
  //calculateProbabilitiesFunctionEnergy(test_N, earth->asArray[xtest][ytest].path);
  cout << "total mass " << earth->totalMass << endl;

  delete earth;
  return 0;
}
