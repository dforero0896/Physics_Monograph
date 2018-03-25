#ifndef UANDINO_H
#define UANDINO_H
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
#include "earth_simul.h"
//Constants
//Mass differences

//Functions
double longitude_units_conversion(double lon_in_km);
double deg2rad(double deg);
void fill_real_matrix(gsl_matrix *empty, double elem_11, double elem_12, double elem_13, double elem_21, double elem_22, double elem_23, double elem_31, double elem_32, double elem_33);
void toFlavor (const gsl_matrix *toTransform, gsl_matrix *destiny, const gsl_matrix *CKM);
gsl_matrix generate_real_identity(gsl_matrix *matrix);
gsl_matrix scale_real_matrix(gsl_matrix *to_scale, double factor);
gsl_matrix add_real_matrices(gsl_matrix *term_1, gsl_matrix *term_2);
void copy_to_complex_from_real(gsl_matrix *real, gsl_matrix_complex *container);
void copy_to_real_from_real(gsl_matrix *real, gsl_matrix *container);
gsl_matrix_complex scale_complex_matrix(gsl_matrix_complex *to_scale, gsl_complex complex_factor, double real_factor);
void fill_complex_matrix(gsl_matrix_complex *empty, gsl_complex elem_11, gsl_complex elem_12, gsl_complex elem_13, gsl_complex elem_21, gsl_complex elem_22, gsl_complex elem_23, gsl_complex elem_31, gsl_complex elem_32, gsl_complex elem_33);
gsl_matrix_complex copy_to_complex_from_complex(gsl_matrix_complex *complex, gsl_matrix_complex *container);
void calculateOperator(double neutrinoEnergy, double A, double L, gsl_matrix_complex* matrix);
void calculateProbabilities(vector<float> path, int N, int Steps, float E_min, float E_max);
#endif
