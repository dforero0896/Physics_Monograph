//g++ -o random_numbers.o random_numbers.cpp `gsl-config --cflags --libs`
#include <iostream>
using namespace std;
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<cmath>
#include <sys/time.h>

struct timeval t;
const gsl_rng_type *Gen_Type;
gsl_rng *Gen;
int i, N;
void gen_random_uniform(){
  cout << gsl_rng_uniform(Gen) << endl;
}

 int main(int argc, char const *argv[]) {

   //Define a variable Gen_Type for the type of random number generator

   //cout << "How many numbers do you want?" << endl;
   //cin >> N;
   N=10000;
   //Setup environment variables for generator
   gsl_rng_env_setup();
   //Seed the generator
   gettimeofday(&t,NULL);
   unsigned long int seed = (t.tv_sec * 1000) + (t.tv_usec / 1000);
   cout <<seed<< endl;
   //Set type as taus because it is the fastest
   Gen_Type = gsl_rng_taus;
   //Allocate necessary memory space
   Gen = gsl_rng_alloc(Gen_Type);
   gsl_rng_set(Gen, seed);

   for(i=0;i<N;i++){
     //Returns a positive number uniformly distributed in the range [0,1)
     //double rand_num = 3.5*gen_random_uniform();
     //Return values that follow a normal distribution.
     double rand_num = gsl_ran_gaussian(Gen, 0.1)+1.;
     cout << rand_num << endl;
   }
   //Free allocated memory
   gsl_rng_free(Gen);
   return 0;
 }
