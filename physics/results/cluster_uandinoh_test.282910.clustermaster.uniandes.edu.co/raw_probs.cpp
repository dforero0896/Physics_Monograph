/*
g++ -fopenmp -o raw_probs.o uandino.cpp earth_simul.cpp raw_probs.cpp `gsl-config --cflags --libs`
*/

#include "uandino.h"
//This program expects Emin, Emax, Steps_in_energy, dist, bse_model, i, k, run_num

int main(int argc, char const *argv[]) {
  if(argc!=9){
    cout << "ERROR\nInvalid number of arguments passed." << endl;
    return 1;
  }
  string dist, bse_model;
  dist=argv[4];
  bse_model=argv[5];
  float E_min, E_max;
  int Steps;
  E_min = atof(argv[1]);
  E_max = atof(argv[2]);
  Steps = atoi(argv[3]);
  if(E_min<0.0005||E_max>4.5005){
    cout << "ERROR\nChoose energy values in [0.0005, 4.5005] MeV" << endl;
    return 1;
  }
  int i,k;
  i=atoi(argv[6]);
  k=atoi(argv[7]);
  Planet *earth = new Planet();
  if(atoi(argv[8])==1){
    earth->Planet::initializeCoords(1);
  }
  else{
    earth->Planet::initializeCoords(0);
  }
  earth->Planet::initializeDensity();
  earth->Planet::initializeAbundanceCrust();
  earth->Planet::initializeAbundanceMantle(dist, bse_model);
  earth->Planet::initializePaths(0, i, k);
  int path_len = int(earth->asArray[i][k].pathLen);
  cout << path_len << endl;

  calculateProbabilities(earth->asArray[i][k].path, Steps, path_len, E_min, E_max );

  

  return 0;
}
