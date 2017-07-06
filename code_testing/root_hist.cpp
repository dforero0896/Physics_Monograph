#include <iostream>
using namespace std;
#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
void histogram_test(){
  TF1 *any_function = new TF1("any_function", "function_name", 0, 3);
  any_function->SetParameters(10., 1.5, 0.5);
  any_function->Draw();
  (any_function->GetHistogram())->Draw();
  TH1F *myhist = new TH1F("myhist", "Histogram", 50, 0, 3);
  myhist ->Draw();
}
# ifndef __CINT__
int main(){
  histogram_test();
  return 0;
}
# endif
