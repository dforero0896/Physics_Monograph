/*
g++ -o PREM_testing.o PREM_testing.cpp `gsl-config --cflags --libs`
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>

int counter=0;
vector<double> copy_vector(vector<double> to_copy){
  vector<double> copy;
  copy.reserve(10);
  for(int n=0;n<10;n++){
    copy.push_back(to_copy[n]);
  }
  return copy;
}

vector< vector<double> > import_model(string filename){
  float rad, depth, density, Vpv, Vph, Vsv, Vsh, eta, Q_mu, Q_kappa;
  string line;
  ifstream infile(filename.c_str());
  int i = 0;
  vector< vector<double> > model_matrix;
  model_matrix.reserve(199);
  vector<double> last_row;
  last_row.reserve(10);
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


    if(i==0){
      //cout << row[0] << ',' << model_matrix[i-1][0] << endl;
      model_matrix.push_back(row);
      last_row = copy_vector(row);
      counter++;
    }

    else if(row[0]<last_row[0]){
      model_matrix.push_back(row);
      last_row = copy_vector(row);
      counter++;
    }
    splittable_line.clear();
    i++;

  }
  return model_matrix;
}
int main(int argc, char const *argv[]) {
  vector< vector<double> > PREM_complete;

  PREM_complete = import_model("../../Models/PREM_1s.csv");

  for(int i=0; i<counter;i++){
    cout << PREM_complete[i][0]<< "," << PREM_complete[i][2]<< endl; //"," << PREM_complete[i][2]<< "," << PREM_complete[i][3]<< "," << PREM_complete[i][4]<< "," << PREM_complete[i][5]<< "," << PREM_complete[i][6]<< "," << PREM_complete[i][7]<< "," << PREM_complete[i][8]<< "," << PREM_complete[i][9] << endl;
  }

  return 0;
}
