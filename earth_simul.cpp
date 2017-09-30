  /*
g++ -fopenmp -o earth_simul.o earth_simul.cpp `gsl-config --cflags --libs`

 */
#include "earth_simul.h"
/*
#include<iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>
using namespace std;
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include "omp.h"
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/*
#include<iostream>
#include <vector>
using namespace std;
#include<gsl/gsl_sf_trig.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<string>
#include<sstream>
#include <omp.h>
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

const int N = 1000;
const int N_x=500;
const int N_z=1000;
const int PREM_len = 187;
const int totalNeutrinos=10000000;
const int path_resolution = 50; //km
int counter =0;
const float R_earth = 6371.;
float dz = 2*R_earth/N;
float dx = dz;
const float PI = 3.1415962589793238;
float abundanceU_BSE[3]={12., 20., 35.};
float abundanceTh_BSE[3]={43., 80., 140.};
//format = {cosm, geoch, geodyn}
const gsl_rng_type *Gen_Type; //Type of random number generator to use.
gsl_rng *Gen; //The actual generator object.
gsl_spline *spectrum_spline;
gsl_interp_accel *acc;

//RingNode asArray[N_x][N_z];

void read_file_into_2D_array(string filename, double to_fill[4500][2]){
  //Fills the array to_fill with 4500 elements for energy and neutrino/energy.

  //Open the file filename
  ifstream spectrum_file(filename.c_str());
  //Strings to save the input as there's header.
  string x, y;
  int it=0; //Count the lines with header.
  //Read file line by line.
  while(spectrum_file >> y >> x){
    it++; //corresponds to iteration number.
    //When header has been read.
    if(it>=12 && it<4512){ //Avoid stack smashing!
      //Convert read strings in stringstreams.
      istringstream x_str(x);
      istringstream y_str(y);
      double x_num, y_num;
      //Cast stringstream into doubles.
      x_str >> x_num;
      y_str >> y_num;
      //Fill array.
      to_fill[it-12][0]=x_num/1e3; //MeV
      to_fill[it-12][1]=y_num*1e3; //1/MeV
      //Clear the stringstreams to avoid problems in next iteration.
      x_str.clear();
      y_str.clear();
    }
  }
}
void split_array(double to_split[4500][2], double to_return[4500], int comp){
  for(int k=0;k<4500;k++){
    to_return[k]=to_split[k][comp];
  }
}

float *retrieve_energies(string filename){
  ifstream energy_repo(filename.c_str());
  float *energy_repo_arr = new float[10000000];
  string value_str;
  int ind=0;
  while (energy_repo >> value_str && ind <10000000){
    istringstream value_ss(value_str);
    float value_num;
    value_ss >> value_num;
    energy_repo_arr[ind]=value_num;
    ind++;
  }
  return energy_repo_arr;
}
vector<float> linspace(float min, float max, int arraylen){
  vector<float> array;
  array.reserve(arraylen);
	float step=(max-min)/arraylen;
	for(int i=0; i<=arraylen; i++){
		array.push_back(min+i*step);
	}
  return array;
}
float the_r(float x, float z, float t, char component){
    float the_r_z=6371-z;
    float the_r_x = x;
    float the_r_mag = sqrt(the_r_x*the_r_x + the_r_z*the_r_z);
    if (component=='x'){
      return x-t*300000*the_r_x/the_r_mag;
    }
    else if(component=='z'){
      return z+t*300000*the_r_z/the_r_mag;
    }
}

float get_r(float x, float y){
  return sqrt(x*x + y*y);
}


float density_polynomials(float radius){
    float x = radius/6371.;
    //inner core
    if( radius<= 1221.5){
        return 13.0885-8.8381*x*x;}
    //outer core
    else if (radius<=3480){
        return 12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x;}
    //lower mantle
    else if (radius <= 5701){
        return 7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x;}
    //transition zone
    else if (radius <= 5771){
        return 5.3197-1.4836*x;}
    else if (radius <= 5971){
        return 11.2494-8.0298*x;}
    else if (radius <= 6151){
        return 7.1089-3.8045*x;}
    //lvz + lid
    else if (radius <= 6346.6){
        return 2.6910+0.6924*x;}
    //crust
    else if (radius <= 6356){
        return 2.9;}
    else if (radius <= 6368){
        return 2.6;}
    //ocean
    else if (radius <= 6371){
        return 1.020;}
}
vector<float> calculateMantleAbundances(float c_mass, float m_mass, float t_mass, float abundance_isot[3], float crust_abundance_isot){
  vector <float> abundances;
  abundances.reserve(3);
  for(int i =0;i<3;i++){
    abundances.push_back((t_mass*abundance_isot[i] - c_mass*crust_abundance_isot)/m_mass);
  }
  return abundances;
}
vector<float> copy_vector(vector<float> to_copy){
  vector<float> copy;
  copy.reserve(10);
  for(int n=0;n<10;n++){
    copy.push_back(to_copy[n]);
  }
  return copy;
}
float density_to_potential(float dty, bool antineutrino){
  float to_return = (1/sqrt(2))*dty*1e-3*8.96189e-47*1e9   /1.672e-27;
  if(antineutrino){
    return -1*to_return;
  }
  else{
    return to_return;
  }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*

//Constants
//Mass differences
double dM32 = 3.2E-3; //eV^2
double dm21 = 0.0; //eV^2
//Vacuum mixing angles
double thetaA = 45.; //Degrees
double thetaB = 5.; //Degrees
double thetaC = 45.; //Degrees
//CKM Elements
double Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3;
gsl_matrix *CKM;

//Functions
double longitude_units_conversion(double lon_in_km){
	return lon_in_km*1e3/(1.972e-7);
}
double deg2rad(double deg){
	return deg*PI/180.;
}
void fill_real_matrix(gsl_matrix *empty, double elem_11, double elem_12, double elem_13, double elem_21, double elem_22, double elem_23, double elem_31, double elem_32, double elem_33){
  gsl_matrix_set(empty, 0, 0, elem_11);
  gsl_matrix_set(empty, 0, 1, elem_12);
  gsl_matrix_set(empty, 0, 2, elem_13);
  gsl_matrix_set(empty, 1, 0, elem_21);
  gsl_matrix_set(empty, 1, 1, elem_22);
  gsl_matrix_set(empty, 1, 2, elem_23);
  gsl_matrix_set(empty, 2, 0, elem_31);
  gsl_matrix_set(empty, 2, 1, elem_32);
  gsl_matrix_set(empty, 2, 2, elem_33);
}

void toFlavor (const gsl_matrix *toTransform, gsl_matrix *destiny, const gsl_matrix *CKM){
	int alpha, beta, a, b;
	long double sum;
	for (alpha=0;alpha<3;alpha++){
		for (beta=0;beta<3; beta++){
			sum=0.;
			for (a=0;a<3;a++){
				for (b=0;b<3;b++){
					sum += gsl_matrix_get(CKM, alpha, a)*gsl_matrix_get(CKM, beta, b)*gsl_matrix_get(toTransform, a, b);

				}
			}
			gsl_matrix_set(destiny, alpha, beta, sum);
		}
	}
}
gsl_matrix generate_real_identity(gsl_matrix *matrix){
	int i, k;
	for(i=0;i<3;i++){
		for(k=0;k<3;k++){
			if(i==k){
				gsl_matrix_set(matrix, i, k, 1);
			}
			else{
				gsl_matrix_set(matrix, i, k,0);
			}
		}
	}
	return *matrix;
}
gsl_matrix scale_real_matrix(gsl_matrix *to_scale, double factor){
  gsl_matrix *result =gsl_matrix_alloc(3,3);
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      double element=gsl_matrix_get(to_scale, i, k);
      element*=factor;
      gsl_matrix_set(result, i,k,element);
    }
  }
  return *result;
}
gsl_matrix add_real_matrices(gsl_matrix *term_1, gsl_matrix *term_2){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      double element_1=gsl_matrix_get(term_1, i, k);
      double element_2=gsl_matrix_get(term_2, i, k);
      element_1+=element_2;
      gsl_matrix_set(term_1, i,k,element_1);
    }
  }
  return *term_1;
}
gsl_matrix_complex copy_to_complex_from_real(gsl_matrix *real, gsl_matrix_complex *container){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_matrix_complex_set(container, i,k, gsl_complex_rect(gsl_matrix_get(real, i, k), 0));
    }
  }
}
gsl_matrix_complex scale_complex_matrix(gsl_matrix_complex *to_scale, gsl_complex complex_factor, double real_factor){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_complex element=gsl_matrix_complex_get(to_scale, i, k);
      element = gsl_complex_mul(element, complex_factor);
      element = gsl_complex_mul_real(element, real_factor);
      gsl_matrix_complex_set(to_scale, i,k,element);
    }
  }
  return *to_scale;}

void fill_complex_matrix(gsl_matrix_complex *empty, gsl_complex elem_11, gsl_complex elem_12, gsl_complex elem_13, gsl_complex elem_21, gsl_complex elem_22, gsl_complex elem_23, gsl_complex elem_31, gsl_complex elem_32, gsl_complex elem_33){
  gsl_matrix_complex_set(empty, 0, 0, elem_11);
  gsl_matrix_complex_set(empty, 0, 1, elem_12);
  gsl_matrix_complex_set(empty, 0, 2, elem_13);
  gsl_matrix_complex_set(empty, 1, 0, elem_21);
  gsl_matrix_complex_set(empty, 1, 1, elem_22);
  gsl_matrix_complex_set(empty, 1, 2, elem_23);
  gsl_matrix_complex_set(empty, 2, 0, elem_31);
  gsl_matrix_complex_set(empty, 2, 1, elem_32);
  gsl_matrix_complex_set(empty, 2, 2, elem_33);
}
gsl_matrix_complex copy_to_complex_from_complex(gsl_matrix_complex *complex, gsl_matrix_complex *container){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_matrix_complex_set(container, i,k, gsl_matrix_complex_get(complex, i, k));
    }
  }
}
gsl_matrix_complex calculateOperator(double neutrinoEnergy, double A, double L){
  double E21=dm21/(2*neutrinoEnergy);
  double E32=dM32/(2*neutrinoEnergy);
  double E12=-E21;
  double E23=-E32;
  double E31=E12-E23;
  double E13=-E31;
  //Elements of the Tmatrix in mass basis
  double T_11=A*Ue1*Ue1-(1./3)*A+(1./3)*(E12+E13);
  double T_12=A*Ue1*Ue2;
  double T_13=A*Ue1*Ue3;
  double T_21=T_12;
  double T_22=A*Ue2*Ue2-(1./3)*A+(1./3)*(E21+E23);
  double T_23=A*Ue2*Ue3;
  double T_31=T_13;
  double T_32=T_23;
  double T_33=A*Ue3*Ue3-(1./3)*A+(1./3)*(E31+E32);
  gsl_matrix *T_mass_mat = gsl_matrix_alloc(3,3);
  fill_real_matrix(T_mass_mat, T_11, T_12, T_13, T_21, T_22, T_23, T_31, T_32, T_33);
  //cout << "T matrix is:"<< endl;
  //print_real_matrix(T_mass_mat);
  //Elements of the T**2 matrix in mass basis
  double T_sq_11=(1./3)*(A*A*(Ue1*Ue1+(1./3))+2*A*(Ue1*Ue1-(1./3))*(E21+E13)+(1./3)*(E12+E13)*(E12+E13));
  double T_sq_12=(1./3)*Ue1*Ue2*A*(A+E13+E23);
  double T_sq_13=(1./3)*Ue1*Ue3*A*(A+E21+E31);
  double T_sq_21=T_sq_12;
  double T_sq_22=(1./3)*(A*A*(Ue2*Ue2+(1./3))+2*A*(Ue2*Ue2-(1./3))*(E21+E23)+(1./3)*(E21+E23)*(E21+E23));
  double T_sq_23=(1./3)*Ue2*Ue3*A*(A+E21+E31);
  double T_sq_31=T_sq_13;
  double T_sq_32=T_sq_23;
  double T_sq_33=(1./3)*(A*A*(Ue3*Ue3+(1./3))+2*A*(Ue3*Ue3-(1./3))*(E31+E32)+(1./3)*(E31+E32)*(E31+E32));
  gsl_matrix *T_sq_mass_mat = gsl_matrix_alloc(3, 3);
  fill_real_matrix(T_sq_mass_mat, T_sq_11, T_sq_12, T_sq_13, T_sq_21, T_sq_22, T_sq_23, T_sq_31, T_sq_32, T_sq_33);
  //cout << "T**2 matrix is:"<< endl;
  //print_real_matrix(T_sq_mass_mat);
  //T matrix in flavor basis
  gsl_matrix *T_flav_mat = gsl_matrix_alloc(3, 3);
  toFlavor(T_mass_mat, T_flav_mat, CKM);
  //cout << "T matrix in flavor basis is:"<< endl;
  //print_real_matrix(T_flav_mat);
  //T**2 matrix in flavor basis
  gsl_matrix *T_sq_flav_mat = gsl_matrix_alloc(3, 3);
  toFlavor(T_sq_mass_mat, T_sq_flav_mat, CKM);
  //cout << "T**2 matrix in flavor basis is:"<< endl;
  //print_real_matrix(T_sq_flav_mat);
  //Get rid of T and T**2 in mass basis as they are no longer useful
  gsl_matrix_free(T_mass_mat);
  gsl_matrix_free(T_sq_mass_mat);
  //Calculate c's
  double c1=T_11*T_22-T_12*T_21+T_11*T_33-T_13*T_31+T_22*T_33-T_23*T_32;
  double c0=-(T_11*T_22*T_33-T_11*T_23*T_32-T_12*T_21*T_33+T_12*T_31*T_23+T_13*T_21*T_32-T_13*T_31*T_22);
  //Calculate eigenvalues
  long double q=c1/3;
  long double r=-0.5*c0;
  //Calculate eigenvalues.
  gsl_complex atanArg = gsl_complex_rect((1./c0)*sqrt(-c0*c0-(4./27.)*c1*c1*c1), 0);
  gsl_complex atanVal=gsl_complex_mul_real(gsl_complex_arctan(atanArg), 1./3.);
  gsl_complex half = gsl_complex_rect(2*sqrt((-1./3.)*c1), 0);
  gsl_complex s1Ps2 = gsl_complex_mul(half, gsl_complex_cos(atanVal));
  gsl_complex dummy_s1Ms2 = gsl_complex_mul(half, gsl_complex_sin(atanVal));
  gsl_complex s1Ms2 = gsl_complex_mul(dummy_s1Ms2, gsl_complex_rect(0., 1.));

  gsl_complex lam1 = gsl_complex_add(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));

  gsl_complex lam2 = gsl_complex_sub(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));

  gsl_complex lam3 = s1Ps2;

  //cout << "Eigenvals" << endl;
  //print_complex_number(lam1);
  //print_complex_number(lam2);
  //print_complex_number(lam3);
  //print_complex_number(gsl_complex_mul_real(gsl_complex_sub(lam3,lam2), 2*neutrinoEnergy));
  //Calculate Operator
  double trace_hamiltonian=0.5*E21+E32+3*neutrinoEnergy+A;
  gsl_complex phi_phase = gsl_complex_polar(1., -L*trace_hamiltonian/3);
  //print_complex_number(phi_phase);


  gsl_complex eigenvalues[3]={lam1, lam2, lam3};
  gsl_matrix_complex *evol_operator = gsl_matrix_complex_alloc(3, 3);
  gsl_complex complex_zero=gsl_complex_rect(0,0);
  fill_complex_matrix(evol_operator, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero);
  for(int n=0;n<3;n++){
    gsl_matrix *I_term = gsl_matrix_alloc(3, 3);
    generate_real_identity(I_term);
    //cout << "Created Identity"<< endl;
    //print_real_matrix(I_term);
    gsl_matrix_scale(I_term, GSL_REAL(eigenvalues[n])*GSL_REAL(eigenvalues[n])+c1);
    //cout << "Scaled Identity"<< endl;
    //print_real_matrix(I_term);
    gsl_matrix *T_term = gsl_matrix_alloc(3, 3);
    *T_term = scale_real_matrix(T_flav_mat, GSL_REAL(eigenvalues[n]));
    gsl_matrix_add(I_term, T_term);
    gsl_matrix_add(I_term, T_sq_flav_mat);
    //cout <<"The three terms"<<endl;
    gsl_matrix_scale(I_term, 1./(3*GSL_REAL(eigenvalues[n])*GSL_REAL(eigenvalues[n])+c1));
    //print_real_matrix(I_term);
    gsl_matrix_complex *sum_term = gsl_matrix_complex_alloc(3, 3);
    copy_to_complex_from_real(I_term, sum_term);
    gsl_matrix_complex_scale(sum_term, gsl_complex_polar(1., -L*GSL_REAL(eigenvalues[n])));
    gsl_matrix_complex_scale(sum_term, phi_phase);
    //print_complex_matrix(sum_term);
    gsl_matrix_complex_add(evol_operator, sum_term);
    gsl_matrix_free(I_term);
    gsl_matrix_free(T_term);
    gsl_matrix_complex_free(sum_term);
  }
  //print_complex_matrix(evol_operator);
  return *evol_operator;
}
void calculateProbabilitiesFunctionEnergy(int steps, vector<float> path){
	int threads =4;
  //CKM matrix elements calculated just once.
	double theta1=deg2rad(thetaA);
	double theta2=deg2rad(thetaB);
	double theta3=deg2rad(thetaC);
  Ue1 = gsl_sf_cos(theta2)*gsl_sf_cos(theta3);
  Ue2 = gsl_sf_sin(theta3)*gsl_sf_cos(theta2);
  Ue3 = gsl_sf_sin(theta2);
  Umu1=-gsl_sf_sin(theta3)*gsl_sf_cos(theta1)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3);
  Umu2=gsl_sf_cos(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_sin(theta3);
  Umu3=gsl_sf_sin(theta1)*gsl_sf_cos(theta2);
  Ut1=gsl_sf_sin(theta1)*gsl_sf_sin(theta3)-gsl_sf_sin(theta2)*gsl_sf_cos(theta1)*gsl_sf_cos(theta3);
  Ut2=-gsl_sf_sin(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta2)*gsl_sf_sin(theta3)*gsl_sf_cos(theta1);
  Ut3=gsl_sf_cos(theta1)*gsl_sf_cos(theta2);
  CKM=gsl_matrix_alloc(3, 3);
  fill_real_matrix(CKM, Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3);

  int NN=10000;
	//vector<float> EnergyLins=linspace(500, 5000000, NN);
  vector<float> EnergyLins=linspace(2, 12, NN);
  omp_set_num_threads(threads);
	int i,k;
	double Probabilities[NN][3];
	#pragma omp parallel for private(i)
	for(i=0;i<NN;i++){
	  double energy=5*pow(10,EnergyLins[i]);
    //double energy = EnergyLins[i];
	  gsl_matrix *Id =gsl_matrix_alloc(3, 3);
	  gsl_matrix_complex *operator_product = gsl_matrix_complex_alloc(3, 3);
	  generate_real_identity(Id);
	  copy_to_complex_from_real(Id, operator_product);
		#pragma omp parallel for private(k)
	  for(k=0;k<steps;k++){
	    double density=double(path[k]);
			double len = double(path_resolution);
	    gsl_matrix_complex *iter_operator = gsl_matrix_complex_alloc(3,3);
	    *iter_operator=calculateOperator(energy, density, longitude_units_conversion(len));
	    gsl_matrix_complex *operator_product_copy = gsl_matrix_complex_alloc(3,3);
	    copy_to_complex_from_complex(operator_product, operator_product_copy);
	    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1., 0), iter_operator, operator_product_copy, gsl_complex_rect(0., 0.),operator_product);
	    gsl_matrix_complex_free(operator_product_copy);
	    gsl_matrix_complex_free(iter_operator);
	  }
		//Probabilities[i] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,1));
		Probabilities[i][0] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,0));
		Probabilities[i][1] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,1));
		Probabilities[i][2] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,2));
    gsl_matrix_complex_free(operator_product);
		//cout << EnergyLins[i] << "," << Probabilities[i][0] << "," << Probabilities[i][1] << "," << Probabilities[i][2] << endl;

	}
  ofstream prob_path_file;
  prob_path_file.open("prob_path.csv");

  for(i=0;i<NN;i++){
		prob_path_file <<5*pow(10,EnergyLins[i]) << "," << Probabilities[i][0] <<  endl;
		//cout << EnergyLins[i] << "," << Probabilities[i] << endl;
	}
  prob_path_file.close();



}
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


/*
vector< vector<float> > import_model(string filename){
  float rad, depth, density, Vpv, Vph, Vsv, Vsh, eta, Q_mu, Q_kappa;
  string line;
  ifstream infile(filename.c_str());
  int i = 0;
  vector< vector<float> > model_matrix;
  model_matrix.reserve(199);
  vector<float> last_row;
  last_row.reserve(10);
  while(getline(infile, line) && i<=199){
    istringstream splittable_line(line);
    string field;
    vector<float> row;
    row.reserve(10);
    while(getline(splittable_line, field, ',')){
      istringstream field_ss(field);
      float field_num;
      field_ss >> field_num;
      row.push_back(field_num);
    }
    if(i==0){
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

void split_array(vector< vector<float> > to_split, float container[PREM_len], int comp, bool invert){
  if(invert){
    for(int i=0;i<PREM_len;i++){
      container[PREM_len-1-i]=to_split[i][comp];
    }
  }
  else{
    for(int i=0;i<PREM_len;i++){
      container[i]=to_split[i][comp];
    }  }
}
*/
/*
float calculateProbability(int steps, vector<float> path, float energy_param){
  int threads =4;
  //CKM matrix elements calculated just once.
  double theta1=deg2rad(thetaA);
  double theta2=deg2rad(thetaB);
  double theta3=deg2rad(thetaC);
  Ue1 = gsl_sf_cos(theta2)*gsl_sf_cos(theta3);
  Ue2 = gsl_sf_sin(theta3)*gsl_sf_cos(theta2);
  Ue3 = gsl_sf_sin(theta2);
  Umu1=-gsl_sf_sin(theta3)*gsl_sf_cos(theta1)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3);
  Umu2=gsl_sf_cos(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_sin(theta3);
  Umu3=gsl_sf_sin(theta1)*gsl_sf_cos(theta2);
  Ut1=gsl_sf_sin(theta1)*gsl_sf_sin(theta3)-gsl_sf_sin(theta2)*gsl_sf_cos(theta1)*gsl_sf_cos(theta3);
  Ut2=-gsl_sf_sin(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta2)*gsl_sf_sin(theta3)*gsl_sf_cos(theta1);
  Ut3=gsl_sf_cos(theta1)*gsl_sf_cos(theta2);
  CKM=gsl_matrix_alloc(3, 3);
  fill_real_matrix(CKM, Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3);

  int NN=10000;
  //vector<float> EnergyLins=linspace(500, 5000000, NN);
  vector<float> EnergyLins=linspace(2, 12, NN);
  omp_set_num_threads(threads);
  int i,k;
  double Probabilities[NN][3];
i=0;
    double energy=double(energy_param);
    //double energy = EnergyLins[i];
    gsl_matrix *Id =gsl_matrix_alloc(3, 3);
    gsl_matrix_complex *operator_product = gsl_matrix_complex_alloc(3, 3);
    generate_real_identity(Id);
    copy_to_complex_from_real(Id, operator_product);
    for(k=0;k<steps;k++){
      double density=double(path[k]);
      double len = double(path_resolution);
      gsl_matrix_complex *iter_operator = gsl_matrix_complex_alloc(3,3);
      *iter_operator=calculateOperator(energy, density, longitude_units_conversion(len));
      gsl_matrix_complex *operator_product_copy = gsl_matrix_complex_alloc(3,3);
      copy_to_complex_from_complex(operator_product, operator_product_copy);
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1., 0), iter_operator, operator_product_copy, gsl_complex_rect(0., 0.),operator_product);
      gsl_matrix_complex_free(operator_product_copy);
      gsl_matrix_complex_free(iter_operator);
    }
    //Probabilities[i] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,1));
    Probabilities[i][0] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,0));
    Probabilities[i][1] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,1));
    Probabilities[i][2] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,2));
    gsl_matrix_complex_free(operator_product);
    //cout << EnergyLins[i] << "," << Probabilities[i][0] << "," << Probabilities[i][1] << "," << Probabilities[i][2] << endl;

  return float(Probabilities[0][0]);
}
*/
    float RingNode::getRadius(){
      r = sqrt(x*x + z*z);
      return r;
    }
    float RingNode::getSolidAngle(){
      solidAngle = (1./((R_earth-z)*(R_earth-z) + x*x));
      return solidAngle;
    }
    float volume;
    float RingNode::getVolume(){
      volume = 2*PI*x*dx*dz;
      return volume;
    }
    void Planet::Planet::Planet::initializeCoords(){
      cout << "Initializing Coordinates" << endl;
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          asArray[i][k].x=i*dx;
          asArray[i][k].z=-6371. + k*dz;
          float r = asArray[i][k].getRadius();
          if(r<6371){
            asArray[i][k].isEarth=1;
            float dk=1e3-k;
            float di=-i;
            if(r>3480){
              asArray[i][k].isSE=1;
              if(r>R_earth-32.5){
                asArray[i][k].isCrust=1;
                asArray[i][k].isMantle=0;
              }
              else{
                asArray[i][k].isCrust=0;
                asArray[i][k].isMantle=1;
              }
            }
            else{
              asArray[i][k].isSE=0;
            }
          }
          else{asArray[i][k].isEarth=0;}
          float dummy_sa = asArray[i][k].getSolidAngle();
          float dummy_vol = asArray[i][k].getVolume();
        }
      }
    }
    void Planet::initializeDensity(){
      cout << "Initializing Density" << endl;
      /*
      vector< vector<float> > PREM_complete;
      PREM_complete = ::import_model("../Models/PREM_1s.csv");
      float radiusArray[PREM_len], densityArray[PREM_len];
      ::split_array(PREM_complete, radiusArray, 0, 1);
      ::split_array(PREM_complete, densityArray, 2, 1);
      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, PREM_len);
      gsl_spline_init(spline, radiusArray, densityArray, PREM_len);
      */

      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isEarth){
            //asArray[i][k].massDensity=gsl_spline_eval(spline, asArray[i][k].r, acc);
            asArray[i][k].massDensity=density_polynomials(asArray[i][k].r);
            asArray[i][k].mass=asArray[i][k].massDensity*1e3*asArray[i][k].volume*1e9;
            totalMass+=asArray[i][k].mass;
            if(asArray[i][k].isCrust){
              crustMass+=asArray[i][k].mass;
            }
            else if(asArray[i][k].isMantle){
              mantleMass+=asArray[i][k].mass;
            }
          }
          else{asArray[i][k].massDensity=-1;}
        }
      }
      //gsl_spline_free(spline);
      //gsl_interp_accel_free(acc);
    }
    void Planet::initializeAbundanceCrust(){
      cout << "Initializing Abundances in the Crust" << endl;
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isCrust){
            asArray[i][k].abundanceU=1.31; //ppb
            asArray[i][k].abundanceTh=5.61; //ppb
          }
          else{
            asArray[i][k].abundanceU=0; //%
            asArray[i][k].abundanceTh=0; //%
          }
        }
      }
    }
    void Planet::initializeAbundanceMantle(string key, string bse_model){
      cout << "Initializing Abundances in the Mantle" << endl;
      int model;
      if(bse_model=="cosmo"){model=0;}
      else if(bse_model=="geoch"){model=1;}
      else if(bse_model=="geodyn"){model=2;}
      if(key=="unif"){
        for(int i =0 ; i<N/2;i++){
          for(int k = 0;k<N;k++){
            if(asArray[i][k].isMantle){
              asArray[i][k].abundanceU=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE, 1.31)[model]; //ppb
              asArray[i][k].abundanceTh=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 5.61)[model]; //ppb
            }
          }
        }
      }
      else if(key=="two_layer"){
        float bulk_mantle_U = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE, 1.31)[model];
        float bulk_mantle_Th = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 5.61)[model];
        float mantleMass_fraction = 0.1*mantleMass;
        float mass_count=0;
        float limit_rad;
        int n=0;
        do {
          mass_count+=4*PI*(asArray[n][500].r)*(asArray[n][500].r)*dx*1e3*1e9*asArray[n][500].massDensity;
          limit_rad=asArray[n][500].r;
          n++;
        } while(mass_count<=mantleMass_fraction);
        for(int i =0 ; i<N/2;i++){
          for(int k = 0;k<N;k++){
            if(asArray[i][k].isMantle && asArray[i][k].r>limit_rad){
              asArray[i][k].abundanceU=4.7; //ppb
              asArray[i][k].abundanceTh=13.7; //ppb
            }
            else if(asArray[i][k].isMantle && asArray[i][k].r<=limit_rad){
              asArray[i][k].abundanceU=(-4.7*(mantleMass-mass_count)+bulk_mantle_U*mantleMass)/mass_count; //ppb
              asArray[i][k].abundanceTh=(-13.7*(mantleMass-mass_count)+bulk_mantle_Th*mantleMass)/mass_count; //ppb
            }
          }
        }

      }
    }
    void Planet::initializeFluxes(){
      cout << "Initializing Fluxes" << endl;
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          if(asArray[i][k].isSE){
            asArray[i][k].neutrinoUFlux=(6.)*(asArray[i][k].abundanceU*1e-9)*(0.9927)*(4.916*1e-18*1e-6)*(asArray[i][k].massDensity*1e-3)*asArray[i][k].volume*asArray[i][k].solidAngle*1e5/(238.051*1.661e-27);
            asArray[i][k].neutrinoThFlux=(4.)*(asArray[i][k].abundanceTh*1e-9)*(1.)*(1.563*1e-18*1e-6)*(asArray[i][k].massDensity*1e-3)*asArray[i][k].volume*asArray[i][k].solidAngle*1e5/(235.044*1.661e-27);
            asArray[i][k].neutrinoFlux=  asArray[i][k].neutrinoUFlux+asArray[i][k].neutrinoThFlux;
            asArray[i][k].relativeNeutrinoU=asArray[i][k].neutrinoUFlux/asArray[i][k].neutrinoFlux;
            asArray[i][k].relativeNeutrinoTh=asArray[i][k].neutrinoThFlux/asArray[i][k].neutrinoFlux;
            totalFlux+=asArray[i][k].neutrinoFlux;
          }
          else{
            asArray[i][k].neutrinoUFlux=0;
            asArray[i][k].neutrinoThFlux=0;
            asArray[i][k].neutrinoFlux= 0;
            asArray[i][k].relativeNeutrinoU=0;
            asArray[i][k].relativeNeutrinoTh=0;
          }
        }
      }
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          asArray[i][k].relativeNeutrino=asArray[i][k].neutrinoFlux/totalFlux;
          asArray[i][k].neutrinosProduced=  roundf(asArray[i][k].relativeNeutrino*totalNeutrinos);
          asArray[i][k].neutrinosProducedU=roundf(asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoU);
          asArray[i][k].neutrinosProducedTh=roundf(asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoTh);
          totalNeut+=asArray[i][k].neutrinosProduced;
        }
      }
    }
    void Planet::initializePaths(){
      cout << "Initializing Potential (density) paths" << endl;
      omp_set_num_threads(4);
      int i, k;
      for(i=0;i<N/2;i++){
        for(k=0;k<N;k++){
          if(asArray[i][k].isSE){
            float z, x;
            z=asArray[i][k].z;
            x=asArray[i][k].x;
            vector<float> path;
            float path_len = get_r(R_earth-z,x );
            asArray[i][k].distanceToDetector=path_len;
            float element_num = roundf(path_len/path_resolution);
            asArray[i][k].pathLen=element_num;
            path.reserve(int(element_num));
            vector<float> times;
            times=linspace(0, path_len/300000, int(element_num));
            for(int n=0;n<int(element_num);n++){
              float t = times[n];
              float dty = density_polynomials(get_r(the_r(x, z, t, 'x'), the_r(x, z, t, 'z')));
              path.push_back(density_to_potential(dty, 1)); //eV
              }
            asArray[i][k].path=path;
          }
        }
      }
    }
    void Planet::initializeEnergySamples(){
      cout << "Initializing Energies for Neutrinos" << endl;
      float *uran_energy_repo;
      float *thor_energy_repo;
      uran_energy_repo = retrieve_energies("energy_repo_238U.knt");
      thor_energy_repo = retrieve_energies("energy_repo_232Th.knt");
      //cout << thor_energy_repo[9999239] << ',' << uran_energy_repo[9935793] << endl;
      int uran_i=0;
      int thor_i=0;
      omp_set_num_threads(4);
      int i, k, n, m;
      for(i=0;i<N/2;i++){
        for(k=0;k<N;k++){
          if(asArray[i][k].isSE){
            vector <float> path_U;
            vector <float> path_Th;
            int U_len = int(asArray[i][k].neutrinosProducedU);
            if (U_len>0){
              path_U.reserve(U_len);
              asArray[i][k].allowedEnergiesU.reserve(U_len);
              for(n=0;n<U_len;n++){
                path_U.push_back(1e6*uran_energy_repo[uran_i]);
                uran_i++;
              }
            }
            int Th_len = int(asArray[i][k].neutrinosProducedTh);
            if(Th_len>0){
              asArray[i][k].allowedEnergiesTh.reserve(Th_len);
              path_Th.reserve(Th_len);
              for(m=0;m<Th_len;m++){
                path_Th.push_back(1e6*thor_energy_repo[thor_i]);
                thor_i++;
              }
            }
            asArray[i][k].allowedEnergiesU=path_U;
            asArray[i][k].allowedEnergiesTh=path_Th;
          }
        }
      }
    }
    /*
    void simulateProbabilities(){
      cout << "Simulating Oscillations" << endl;
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          int U_len = int(asArray[i][k].neutrinosProducedU);
          int Th_len = int(asArray[i][k].neutrinosProducedTh);
    //      vector <double> probabilitiesU;
  //        vector <double> probabilitiesTh;
          int n, m;
          if (U_len>0){
//            probabilitiesU.reserve(U_len);
            asArray[i][k].probabilitiesU.reserve(U_len);
            for(n=0;n<U_len;n++){
              asArray[i][k].probabilitiesU.push_back(calculateProbability(int(asArray[i][k].pathLen), asArray[i][k].path, asArray[i][k].allowedEnergiesU[n]));
            }

          }
          if(Th_len>0){
            asArray[i][k].probabilitiesTh.reserve(Th_len);
      //      probabilitiesTh.reserve(Th_len);
            for(m=0;m<Th_len;m++){
              asArray[i][k].probabilitiesTh.push_back(calculateProbability(int(asArray[i][k].pathLen), asArray[i][k].path, asArray[i][k].allowedEnergiesTh[n]));
            }

          }
        }
      }
    }
*/

    void Planet::initialize(string key, string bse_model){
      cout << "Building Planet" << endl;
      Planet::initializeCoords();
      Planet::initializeDensity();
      Planet::initializeAbundanceCrust();
      Planet::initializeAbundanceMantle(key, bse_model);
      Planet::initializeFluxes();
      Planet::initializePaths();
      Planet::initializeEnergySamples();
      //simulateProbabilities();
      cout << "Done" << endl;
    }



/*

int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->initialize("two_layer", "geodyn");
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
*/
