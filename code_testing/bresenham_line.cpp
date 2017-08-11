#include<iostream>
using namespace std;
#include<cmath>
/*
g++ -o bresenham_line.o bresenham_line.cpp
*/



int main(int argc, char const *argv[]) {
  bool testing[500][1000];
  for(int i=0;i<500;i++){
    for(int k=0;k<1000;k++){
      testing[i][k]=0;
    }
  }
  for(int n = 2; n<3; n+=150){

    int i_init = 0;
    int i_end = 25;
    int k_init = 999;
    int k_end = n;
    int deltaX = i_end-i_init;
    int deltaY = k_end - k_init;
    float delta_error = abs(float(deltaY)/deltaX);
    int k, i;
    if(abs(deltaX)>=abs(deltaY)){
      float error = delta_error-0.5;
      k = k_init;
      for(i=i_init;i<i_end;i++){
        testing[i][k]=10;
        error+=delta_error;
        if(error>=0.5){
          k--;
          error-=1.;
        }
      }
    }
    else{
      float error = 1./delta_error+0.5;
      i = i_init;
      for(k=k_init;k>k_end;k--){
        testing[i][k]=10;
        error+=1./delta_error;
        if(error>=0.5){
          i++;
          error-=1.;
        }
      }
    }
}




    for(int i=0;i<50;i++){
      for(int k=0;k<100;k++){
        cout << testing[i][k]  << ',' ;
        }
        cout << 0 << endl;
    }
  return 0;
}
