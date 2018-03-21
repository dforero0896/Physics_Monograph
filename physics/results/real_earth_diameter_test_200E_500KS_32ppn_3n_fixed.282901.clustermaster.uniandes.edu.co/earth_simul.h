#ifndef EARTH_SIMUL_H
#define EARTH_SIMUL_H
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
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
const int PREM_len = 187;
const int totalNeutrinos=10000000;
const float path_resolution = 1e-2;
const int N = 1000;


extern const int N_x;
extern const int N_z;
void read_file_into_2D_array(string filename, double to_fill[4500][2]);
void split_array(double to_split[4500][2], double to_return[4500], int comp);
float *retrieve_energies(string filename);
vector<float> linspace(float min, float max, int arraylen);
float the_r(float x, float z, float t, char component);
float get_r(float x, float y);
float density_polynomials(float radius);
vector<float> calculateMantleAbundances(float c_mass, float m_mass, float t_mass, float abundance_isot[3], float crust_abundance_isot);
vector<float> copy_vector(vector<float> to_copy);
float density_to_potential(float dty, bool antineutrino);

class RingNode{
  public:
    float x;
    float z;
    float r;
    float getRadius();
    bool isEarth;
    bool isSE;
    bool isCrust;
    bool isMantle;
    float mass;
    float massDensity;
    float solidAngle;
    float getSolidAngle();
    float volume;
    float getVolume();
    float abundanceU;
    float abundanceTh;
    float neutrinoFlux;
    float neutrinoThFlux;
    float neutrinoUFlux;
    float neutrinoFluxMeasure;
    float neutrinoThFluxMeasure;
    float neutrinoUFluxMeasure;
    float relativeNeutrinoTh;
    float relativeNeutrinoU;
    float relativeNeutrino;
    float neutrinosProduced;
    float neutrinosProducedU;
    float neutrinosProducedTh;
    vector<float> path;
    float pathLen;
    float distanceToDetector;
    vector <float> allowedEnergiesTh;
    vector <float> allowedEnergiesU;
    vector <float> probabilitiesTh;
    vector <float> probabilitiesU;
    float meanSurvProb;


};

class Planet{
  public:
    RingNode asArray[500][1000];
    float totalMass;
    float crustMass;
    float mantleMass;
    float totalFlux;
    float totalUFlux;
    float totalThFlux;
    int totalNeut;
    void initializeCoords(bool expo);
    void initializeDensity();
    void initializeAbundanceCrust();
    void initializeAbundanceMantle(string key, string bse_model);
    void initializeFluxes(bool oscillated, string hpe_dist, string bse_model);
    void initializePaths(bool all, int i, int k);
    void initializeEnergySamples();
    void initializeProbabilities();
    void initialize(string key, string bse_model);

};
#endif
