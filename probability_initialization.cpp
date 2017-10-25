/*
g++ -fopenmp -o probability_initialization.o earth_simul.cpp probability_initialization.cpp `gsl-config --cflags --libs`
*/
#include "earth_simul.h"
#include <iostream>
using namespace std;
int N =1000;
int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();
  earth -> initializeCoords();
  
  for(int i =0 ; i<500;i++){
    for(int k = 0;k<N;k++){
      if(earth->asArray[i][k].isSE){
        earth->asArray[i][k].meanSurvProb = 0.55;}
    }
  }
  ofstream outfile;
  outfile.open("probability_planet.csv");
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      outfile << i<< ' ' << k <<' ' << earth->asArray[i][k].meanSurvProb << endl;
      }
    }
  outfile.close();
  return 0;
}
