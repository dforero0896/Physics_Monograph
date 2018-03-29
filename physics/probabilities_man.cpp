/*
g++ -fopenmp -o probabilities_man.o earth_simul.cpp probabilities_man.cpp `gsl-config --cflags --libs`

*/
#include "earth_simul.h"


int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->initialize("two_layer", "geodyn");
  cout << "total flux " << earth->totalFlux << endl;
  /*In File outfile, saves a matrix corresponding to the neutrino flux in each node, according to the parameters introduced.*/
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
/*In file test_path_file, saves the density (potential) profile for the test node (x,y) given.*/
  ofstream test_path_file;
  test_path_file.open("test_path.csv");
  int test_N = int(earth->asArray[xtest][ytest].pathLen);
  for(int step=0;step<test_N;step++){
    test_path_file << earth->asArray[xtest][ytest].path[step] << endl;
  }
  test_path_file.close();

  delete earth;
  return 0;
}
