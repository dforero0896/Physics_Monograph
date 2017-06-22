/*
g++ -o array_importing.o array_importing.cpp `gsl-config --cflags --libs`
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>
//"../../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt"
void spectrum_array2D(string filename, float to_fill[4500][2]){
  ifstream infile(filename.c_str());

  string x, y;
  int i=0;
  //vector<float> x_vec;
  //vector<float> y_vec;
  while(infile >> y >> x){
    i++; //i corresponds to iteration number
    if(i>=12){
      istringstream x_str(x);
      istringstream y_str(y);
      double x_num, y_num;
      x_str>>x_num ;
      y_str>>y_num;
      //y_vec.push_back(y_num);
      to_fill[i-12][0]=x_num;
      to_fill[i-12][1]=y_num;
      //x_vec.push_back(x_num);
      //cout << y  << endl;
      x_str.clear();
      y_str.clear();
      //cout << y_num << endl;

    }

  }
}


int main(int argc, char const *argv[]) {
float arr[4500][2];
spectrum_array2D("../../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt", arr);
cout << arr[30][0]<< endl;

  return 0;
}
