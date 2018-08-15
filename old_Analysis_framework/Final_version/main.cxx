#include "TApplication.h"
#include "TROOT.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "makePlots.h"


int main (int argc, char* argv[])
{

  char *cvalue = NULL;
  char *cvalue2= NULL;
  int c;
  int pidflag=1;
  unsigned int truthflag=1;//MC is the default
  while ((c = getopt (argc, argv, "c:i:t:")) != -1)
    switch (c) {
    case 't':
      cvalue2 = optarg;
      break;
    case 'c':
      cvalue = optarg;
      break;
    case '?':
      if (optopt == 'c')
        std::cout<<"Option -%c requires an argument"<<std::endl;
      else if (isprint (optopt))
        std::cout<<"Unknown option `-%c'.\n"<<std::endl;
      else
        std::cout<<"Unknown option character" <<optopt<<std::endl;
      return 1;
    default:
      abort();
    }
  
  if(cvalue){
    pidflag=atoi(cvalue);
  }
  if(cvalue2){
    truthflag=atoi(cvalue2);
  }

  std::ifstream intxtFile;
  //intxtFile.open("data.txt"); 
  intxtFile.open("MC.txt"); 

  std::string rootdir;
  if(truthflag) {
    rootdir = "HGCalTBAnalyzer/HGCTB";
  }
  else {
    rootdir = "hgcaltbntuple/HGC_Events";
  }
  TChain *chain=new TChain(rootdir.c_str());


  std::string filename;
  while(!intxtFile.eof())
    {
      std::getline(intxtFile, filename);
      if (filename.size()<2)
	{
	  break;
	}
      chain->Add(filename.c_str());
    }



  std::cout << "To run it you type: " << std::endl;
  std::cout << "./makePlots -t <truthflag>"<< std::endl;
  //std::cout << "pidflag: Run-type "<< std::endl;
  std::cout << "<truthflag>: 0 or 1  ...  (1 is true) "<< std::endl;


  bool doTruth=truthflag;
  if(!doTruth) {
    std::cout << "doTruth: False (ie do Data)"<<std::endl;
  }
  else {
    std::cout << "doTruth: True (ie do MC)"<<std::endl;
  }

  makePlots makePlots(chain);

  makePlots.doTruth=doTruth;
  makePlots.Loop();

  cout << "Last line in main ..." << endl;
  //return 0;

}
