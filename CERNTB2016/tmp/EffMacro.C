

void EffMacro() {
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TDirectory.h"
  
using namespace std;
  
  TFile *_file0;
  _file0 = new TFile::Open("rfio:/castor/cern.ch/user/s/sudha/PYTHIA7TeV/PhotonJet/MCEff/TestPhJt_20_30_looseID_1.root");
  TDirectory *t =  _file0->Get("PhotonMCEfficiencyStudies/mcPhotonIDCutEfficiencyIsLoose");
  t->ls();
  TH1F *AllLooseReco=(TH1F *)histAllRecoPhotonsSCLooseCutDist->Project3D("z")->Clone("AllLooseReco");
  TH1F *AllReco=(TH1F *)histAllRecoPhotonsSuperClusterDist->Project3D("z")->Clone("AllReco");
  TH1F *LooseEff = new TH1F("LooseEff","ratio",200,0,100);
  LooseEff->Divide(AllLooseReco,AllReco,1,1,"");
  TH1F *Pt_photon=(TH1F *)histAllRecoPhotonsSuperClusterDist->Project3D("x")->Clone("Pt_photon");

  TH2F *EffvsPt = new TH2F("EffvsPt","eff vs Pt_photon",50,0,5,500,0,1000);
  EffvsPt->Fill(Pt_photon,LooseEff);
  EffvsPt->Draw();
  
}
