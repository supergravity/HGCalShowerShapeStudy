#include "TROOT.h"
#include "TChain.h"
#include "TProof.h"

void runAll(string fileToOpen, string outfilename){
  

  //gROOT->LoadMacro("HistManager.cc+");
  TChain *fChain = new TChain("HGCalTBAnalyzer/HGCTB");
  /*fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_1.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_1.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_10.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_100.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_101.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_102.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_103.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_104.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_105.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_106.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_107.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_108.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_109.root");
  fChain->Add("root://eoscms.cern.ch//eos/cms/store/group/upgrade/HGCAL/simulation/1moduleIv1/mc/CRAB_PrivateMC/crab_Ele100GeV/160926_131407/0000/TBGenSim_110.root");
  */


  ifstream file;
  file.open(fileToOpen.c_str(), ifstream::in );
  char filename[200];
  while(!file.eof()){
    file >> filename;
    if(strncmp(filename,"#",1)==0)
      {
        continue;
      }
    cout<<"filename = "<<filename<<endl;
    
    fChain->Add(filename);   

  }
  
    
  TProof *plite = TProof::Open("workers=6");
  fChain->SetProof();
  fChain->Process("HGCTB.C+", outfilename.c_str());
  fChain->SetProof(0);

  






}
