#ifndef makePlots_cxx
#define makePlots_cxx

#include "makePlots.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TRandom1.h"
#include  <cmath>
#include <iostream>
#include <algorithm>
#include "TProfile.h"


#define ZMASS 91187.0
//Constructor-initializer
makePlots::makePlots(TChain* inchain ):fChain(inchain)
{

  //Constants:

  //Logicals:

  //Output:
  //eve = new std::vector<int>();
  //run = new std::vector<int>();
  //m4l = new std::vector<double>();
  
}

makePlots::~makePlots()
{
  cout << "Calling Destructor" << endl;

}


void makePlots::Init()
{


  // Set object pointer (if there are vectors etc... need to be initialized to 0)
   simHitLayEn1E = 0;
   simHitLayEn2E = 0;
   simHitLayEn1H = 0;
   simHitLayEn2H = 0;
   simHitCellIdE = 0;
   simHitCellEnE = 0;
   simHitCellIdH = 0;
   simHitCellEnH = 0;

   cellID = 0;
   x = 0;
   y = 0;
   z = 0;
   energy = 0;
   cluster_x = 0;
   cluster_y = 0;
   cluster_z = 0;
   cluster_energy = 0;
   cluster_size = 0;

   // Set branch addresses and branch pointers
   fCurrent = -1;
   fChain->SetMakeClass(1);

   if(doTruth) {
     fChain->SetBranchAddress("simHitLayEn1E", &simHitLayEn1E, &b_simHitLayEn1E);
     fChain->SetBranchAddress("simHitLayEn2E", &simHitLayEn2E, &b_simHitLayEn2E);
     fChain->SetBranchAddress("simHitLayEn1H", &simHitLayEn1H, &b_simHitLayEn1H);
     fChain->SetBranchAddress("simHitLayEn2H", &simHitLayEn2H, &b_simHitLayEn2H);
     fChain->SetBranchAddress("xBeam", &xBeam, &b_xBeam);
     fChain->SetBranchAddress("yBeam", &yBeam, &b_yBeam);
     fChain->SetBranchAddress("zBeam", &zBeam, &b_zBeam);
     fChain->SetBranchAddress("pBeam", &pBeam, &b_pBeam);
     fChain->SetBranchAddress("simHitCellIdE", &simHitCellIdE, &b_simHitCellIdE);
     fChain->SetBranchAddress("simHitCellEnE", &simHitCellEnE, &b_simHitCellEnE);
     fChain->SetBranchAddress("simHitCellIdH", &simHitCellIdH, &b_simHitCellIdH);
     fChain->SetBranchAddress("simHitCellEnH", &simHitCellEnH, &b_simHitCellEnH);
   } else {
     fChain->SetBranchAddress("evtID", &evtID, &b_evtID);
     fChain->SetBranchAddress("nhit", &nhit, &b_nhit);
     fChain->SetBranchAddress("cellID", &cellID, &b_cellID);
     //fChain->SetBranchAddress("cellID", &simHitCellIdE, &b_simHitCellIdE);
     fChain->SetBranchAddress("x", &x, &b_x);
     fChain->SetBranchAddress("y", &y, &b_y);
     fChain->SetBranchAddress("z", &z, &b_z);
     fChain->SetBranchAddress("energy", &energy, &b_energy);
     //fChain->SetBranchAddress("energy", &simHitCellIdE, &b_simHitCellIdE);
     fChain->SetBranchAddress("thrustX0", &thrustX0, &b_thrustX0);
     fChain->SetBranchAddress("thrustX", &thrustX, &b_thrustX);
     fChain->SetBranchAddress("thrustY0", &thrustY0, &b_thrustY0);
     fChain->SetBranchAddress("thrustY", &thrustY, &b_thrustY);
     fChain->SetBranchAddress("cluster_x", &cluster_x, &b_cluster_x);
     fChain->SetBranchAddress("cluster_y", &cluster_y, &b_cluster_y);
     fChain->SetBranchAddress("cluster_z", &cluster_z, &b_cluster_z);
     fChain->SetBranchAddress("cluster_energy", &cluster_energy, &b_cluster_energy);
     fChain->SetBranchAddress("cluster_size", &cluster_size, &b_cluster_size);
   }

   Notify();



}
Int_t makePlots::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
void makePlots::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makePlots::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Long64_t makePlots::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
Bool_t makePlots::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makePlots::Loop()
{

  Init();


  TFile* outfile;
  
  if(doTruth)
  {
    //outfile=new TFile("outputMC.root","RECREATE");
    outfile=new TFile("output.root","RECREATE");
  }
  else
  {
    //outfile=new TFile("outputData.root","RECREATE");
    outfile=new TFile("output.root","RECREATE");
  }
  double upbound =0;
  if(doTruth)
    {
      upbound = 10;
    }
    else
    {
      upbound = 10;
    }
    
  //  for (int iL = 0;iL < NLAYERS;iL++)
  //  {
  //      TwoPad_corr100[iL] = new TH2D(Form("TwoPad_corr100Layer%d",iL),"",9,0.5,9.5,10000000,0.,upbound);
  //  } 
   //TH2D* TwoPad_corr100=new TH2D("TwoPad_corr100","",9,0.5,9.5,10000000,0.,upbound);
   TH2D* TwoPad_corr100=new TH2D("TwoPad_corr100","",9,0.5,9.5,1000000,0.,upbound); // Number of Bins: 10^6
  //
  //* Counters
  //
  int EventsPassed=0;

  //
  //* Loop over Events (!)
  //
  Long64_t nentries = fChain->GetEntries();
  //nentries = 93997;
  std::cout << "Entries: " << nentries << std::endl;
  //for (Long64_t entry=0;entry<nentries;++entry){
  for (Long64_t entry=0;entry<nentries;++entry){
    if(entry%1000==0) std::cout << "Processed ... " << entry 
                               << "/" << nentries << " events" << std::endl;

    fChain->GetEntry(entry); 

    //
    //* How to "see" the contents of an event:
    //
    //Show(entry);

    //
    //* Build Layer objects
    //
    int ibuild = buildCompObjects();
        
    
    // ===== Preselection =====
        if ( !Layer.size() ) continue;
        double xmax_pre = Layer.at(0)->GetXMax();
        double ymax_pre = Layer.at(0)->GetYMax();
   
        //cout <<"The xmax (seed) is at "<< ymax_pre << "mm"<<endl;
        //cout <<"The ymax (seed) is at "<< xmax_pre << "mm"<<endl;
        //getchar();

    
        //if(!(fabs(xmax_pre)< 20.0) && !(fabs(ymax_pre)< 20.0))
        //if(!(fabs(xmax_pre)< 2.0 && fabs(ymax_pre)< 2.0))
        if(!(fabs(xmax_pre)< 2.0 && fabs(ymax_pre)< 2.0))
        {   
            Clear();
            continue; 
        }
   // ====== Elayer 2 cut ======     
        double E2vis = 0.0;
        int Nhits = Layer.at(1) -> GetNhits();
        for (int ih = 0;ih < Nhits;ih++)
        {
            double Ehit = Layer.at(1) -> GetErawHit(ih);
            E2vis = E2vis + Ehit;
        }
        E2vis = E2vis*GEVTOMEV;
        if (!( E2vis > 3))//E2vis > 3 MeV
        {
            Clear();
            continue;
        }


        //
        //* Calculate Shower Depth
        //
        double shDepth = 0.0;
        double sumEmax = 0.0;
        for(int iL=0; iL<NLAYERS; iL++) {
          double lX0s=Layer.at(iL)->GetX0depth(iL);
          double emax=Layer.at(iL)->GetEneMax();
          shDepth += (emax*lX0s);
          sumEmax += emax;
          //cout << "L=" << iL << " X0:" << lX0s << " Emax:" << emax << endl;
        }
        //getchar();
        shDepth = shDepth/sumEmax;//shower depth using emax in layer.
    
        //Two Point Correlation study
        //
        for (int iL = 0;iL < NLAYERS;iL++)
        {
            if(iL == 1)
            {
                int nhits = Layer.at(iL) -> GetNhits();
                for (int ih = 0;ih < nhits;ih++)//outer hit loop
                {
                    for (int ij = ih+1; ij  < nhits;ij++)
                    {
                        double xhit1=Layer.at(iL)->GetXHit(ih);//mm
                        double yhit1=Layer.at(iL)->GetYHit(ih);
                        double xhit2=Layer.at(iL)->GetXHit(ij);
                        double yhit2=Layer.at(iL)->GetYHit(ij);
                        double dRij =sqrt((xhit1-xhit2)*(xhit1-xhit2)+
                              (yhit1-yhit2)*(yhit1-yhit2))*MMtoCM;//cm
                        double EMeV1=Layer.at(iL)->GetErawHit(ih)*GEVTOMEV;
                        double EMeV2=Layer.at(iL)->GetErawHit(ij)*GEVTOMEV;
                        double emax=Layer.at(iL)->GetEneMax();
                        double EmaxMeV=emax*GEVTOMEV;
                        double Ecorr=10.*EMeV1*EMeV2/(EmaxMeV*EmaxMeV);
                        double dRing=0;
                        //if (EMeV1 < 1.5*ENEPERMIP) cout <<"Bingo"<<endl;
                        if(dRij>0.2 && (EMeV1*MEVTOGEV>1.*ENEPERMIP || EMeV2*MEVTOGEV>1.*ENEPERMIP)) {
                          //cout << "dRij= " << dRij << " Ecorr= " <<Ecorr << endl;
                          if(dRij<1.3) { 
                            dRing=0.; 
                          } 
                          else if(dRij<2.35) {dRing=1.;}
                          else if(dRij<3.40) {dRing=2.;}
                          else if(dRij<4.20) {dRing=3.;}
                          else if(dRij<5.30) {dRing=4.;}
                          else if(dRij<6.30) {dRing=5.;}
                          else if(dRij<7.30) {dRing=6.;}
                          else if(dRij<8.40) {dRing=7.;}
                          else if(dRij<9.40) {dRing=8.;}
                          else {dRing=9.;}
                          
                          TwoPad_corr100 ->Fill(dRing,Ecorr);
                          //TwoPad_corr100 [iL]->Fill(dRing,Ecorr);
                        }//IF ends 
                
                    }//inner loop ends
                }//outer loop ends
        
            }//layer if ends 
        }//layer loop ends
        EventsPassed++;
        Clear();
  }//event loop ends
  cout << "PassedEvents :" << EventsPassed << endl;
  cout << "Entries      :" << nentries << endl;
  cout << "Efficiency is :" << ((float) EventsPassed)/nentries << endl;

  outfile->Write("",TObject::kOverwrite);
  outfile->Close();

}//loop function ends
#endif // #ifdef makePlots_cxx
