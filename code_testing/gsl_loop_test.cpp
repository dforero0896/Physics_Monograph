#include<iostream>
using namespace std;
#include<gsl/gsl_sf_trig.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<string>
#include<sstream>


//Vacuum mixing angles
double thetaA = 45.; //Degrees
double thetaB = 5.; //Degrees
double thetaC = 45.; //Degrees
double PI = 3.1415926589793238;
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
void linspace(double array[], double max, double min, int arraylen){
	int i;
	double step=(max-min)/arraylen;
	for(i=0; i<=arraylen; i++){
		array[i]=min+i*step;
	}
}
void densityStep ( double *fill , double *dist, int arraysize){

	int i;
	for(i=0;i<arraysize;i++){
		if(dist[i] < 2885.){fill[i]=1.7e-13;}
		else {fill[i] = 4.4e-13;}
	}

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
void print_real_matrix(gsl_matrix *to_print){
  for(int i=0;i<3;i++){
    cout << gsl_matrix_get(to_print, i,0) << "," << gsl_matrix_get(to_print, i, 1) << "," << gsl_matrix_get(to_print, i, 2) << endl;

  }
}
string print_complex_number(gsl_complex to_print){
  stringstream real;
  stringstream imag;
  real << GSL_REAL(to_print);
  imag << GSL_IMAG(to_print);
  string result= real.str() + "+ j" + imag.str();
  //cout << result << endl;
  return result;

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
void print_complex_matrix(gsl_matrix_complex *to_print){
  for(int i=0;i<3;i++){
    cout << print_complex_number(gsl_matrix_complex_get(to_print, i,0)) << "," << print_complex_number(gsl_matrix_complex_get(to_print, i,1)) << "," << print_complex_number(gsl_matrix_complex_get(to_print, i,2)) << endl;

  }
}
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
int main(int argc, char const *argv[]) {
  //Lets test the unitarity os the CKM matrix
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
//  cout << "CKM Matrix is:" << endl;
//  cout << Ue1<< ","<< Ue2<<"," << Ue3 <<endl;
//  cout << Umu1<<"," << Umu2<<"," << Umu3 <<endl;
//  cout << Ut1<<","<< Ut2<<"," << Ut3 <<endl;
  CKM=gsl_matrix_alloc(3, 3);
  fill_real_matrix(CKM, Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3);
  gsl_matrix *trans_CKM = gsl_matrix_alloc(3, 3);
  fill_real_matrix(trans_CKM, Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3);
  gsl_matrix *result = gsl_matrix_alloc(3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., CKM, trans_CKM, 0., result);
  cout << "This is the product of the CKM matrix and it's transpose" << endl;
  print_real_matrix(result);
  cout << "It is Unitary enough" << endl;

  //Now, let us test looping on matrix products:
  gsl_matrix *Id =gsl_matrix_alloc(3, 3);
  gsl_matrix_complex *operator_product = gsl_matrix_complex_alloc(3, 3);
  generate_real_identity(Id);
  copy_to_complex_from_real(Id, operator_product);
  for(int n=1;n<5;n++){
    gsl_complex element=gsl_complex_rect(n, n-1);
    gsl_matrix_complex *iter_operator =gsl_matrix_complex_alloc(3, 3);
    fill_complex_matrix(iter_operator, element, element, element, element, element, element, element, element, element);
    gsl_matrix_complex *operator_product_copy = gsl_matrix_complex_alloc(3,3);
    copy_to_complex_from_complex(operator_product, operator_product_copy);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1., 0), iter_operator, operator_product_copy, gsl_complex_rect(0., 0.),operator_product);
    gsl_matrix_complex_free(operator_product_copy);
    gsl_matrix_complex_free(iter_operator);
  }
  cout << "This is a complex matrix productorial test" << endl;
  print_complex_matrix(operator_product);
  cout << "Passed test, for comparison look ipynb code." << endl;
  return 0;
}
