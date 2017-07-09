//Compile with
//g++ -L/usr/local/lib/root -lPhysics -lPostscript -lGraf3d -lImt -L/home/daniel/anaconda2/lib/ -lpcre -o hello.o hello.cpp `root-config --cflags --glibs`
/*
# include <iostream>

using namespace std;

void hello()
{
  cout << "hello world!" << endl;
}

# ifndef __CINT__
int main()
{
  hello();
  return 0;
}
# endif
*/
# include <iostream>
# include "TRandom.h"
# include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
using namespace std;
/*
void test()
{
  TRandom *rnd = new TRandom(time(0));
  Double_t x = rnd->Rndm();
  cout << "x = " << x << endl;
}
*/
void test(){
  Float_t xmin = 0.0005;
  Float_t xmax=4.5;
  Int_t nbins = 1000;
  Float_t tofill = 3.5;
  TH1F *hist = new TH1F("hist", "AntineutrinoSpectrum", nbins, xmin, xmax);
  hist->Fill(tofill);
  TFile *outfile = new TFile("hist.root", "recreate");
  hist->SetDirectory(outfile);
  //hist->GetXaxis()->SetTitle("Energy (MeV)");
  //Double_t scale = norm/hist->Integral();
  //hist->Scale(scale);
  //hist->Draw();
  hist->Write();
  outfile->Close();
}
# ifndef __CINT__
int main()
{
  test();
  return 0;
}
# endif
