#include "myClass.C"
#include "TROOT.h"


int main(){
  gROOT->ProcessLine("#include <vector>");
  //gROOT->ProcessLine(".L loader.C+");
  gROOT->ProcessLine("#include <map>");
  //myClass m(23);
  myClass m(11);
}
