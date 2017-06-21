#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <vector>
int main(int argc, char const *argv[]) {

  ifstream fs("../../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt");
  string line;
  while (std::getline(fs, line)) {
    cout<<line<<endl;
  }
  // ...

  return 0;
}
