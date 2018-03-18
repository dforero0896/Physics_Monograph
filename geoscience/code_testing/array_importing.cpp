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



void import_model(string filename, double model_matrix[199][10]){
  float rad, depth, density, Vpv, Vph, Vsv, Vsh, eta, Q_mu, Q_kappa;
  string line;
  ifstream infile(filename.c_str());
  int i = 0;
  while(getline(infile, line) && i<=199){
    istringstream splittable_line(line);
    string field;
    vector<double> row;
    row.reserve(10);
    while(getline(splittable_line, field, ',')){
      istringstream field_ss(field);
      double field_num;
      field_ss >> field_num;
      row.push_back(field_num);
    }
    //splittable_line >> depth;
    //splittable_line >> Vpv;

/*
    model_matrix[i][0]=rad;
    model_matrix[i][1]=depth;
    model_matrix[i][2]=density;
    model_matrix[i][3]=Vpv;
    model_matrix[i][4]=Vph;
    model_matrix[i][5]=Vsv;
    model_matrix[i][6]=Vsh;
    model_matrix[i][7]=eta;
    model_matrix[i][8]=Q_mu;
    model_matrix[i][9]=Q_kappa;
    i++;
*/
    cout << row[4] << endl;
    splittable_line.clear();
  }
}
vector < vector <float> > import_probability(string filename){
  string line;
  ifstream infile(filename.c_str());
  int i = 0; //iteration/line
  vector < vector <float> > prob_matrix;
  prob_matrix.reserve(1000);
  while(getline(infile, line) && i<1000){
    istringstream splittable_line(line);
    string field;
    vector<float> row;
    row.reserve(500);
    while(getline(splittable_line, field, ',')){
      istringstream field_ss(field);
      float field_num;
      field_ss >> field_num;
      row.push_back(field_num);
    }
    prob_matrix.push_back(row);
    splittable_line.clear();
    i++;
  }
  return prob_matrix;
}
int main(int argc, char const *argv[]) {
  vector< vector<float> > prob_matrix = import_probability("../probability_planet.csv");
  for(int i = 0;i<1000;i++){
    for(int k = 0;k<500;k++){
      cout << prob_matrix[i][k] << ',';
    }
    cout << '0' << endl;
  }

  return 0;
}
