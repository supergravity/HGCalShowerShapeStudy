#define HGCTB_cxx
// The class definition in HGCTB.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("HGCTB.C")
// root> T->Process("HGCTB.C","some options")
// root> T->Process("HGCTB.C+")
//


#include "HGCTB.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>

#include "TMath.h"
//#include "HGCalTBPlots.C"

//double ADCperMIP = 16; 
//double ADCperMIP = 17/1.06; 
//double ADCperMIP = 17; 
const double ALLCELLS_THRESHOLD = 50.;
const double SEVENCELLS_THRESHOLD = 50.;
const double NINETEENCELLS_THRESHOLD = 50.;

HexGeometry geomc(false);

//bool debug = true;
bool debug = false;

void HGCTB::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void HGCTB::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).


   TString option = GetOption();

   //HexTopology ht1(false);
   //fout = new TFile("histoSIM.root","RECREATE");

   //EperMIP = 55.91e-6; //taken from pin run at 125 GeV
   string OUTFILENAME = "";
   cout<<"option = "<<option.Data()<<endl;
   TObjArray *args = (TObjArray*)option.Tokenize(" ");
   if (args->GetSize()>1)
     OUTFILENAME = (string)((TObjString*)args->At(0))->GetString();


   if(debug) cout<<"Inside SlaveBegin"<<endl;

   //EperMIP = 62.671e-6; //taken from EUDET ... from https://www.eudet.org/e26/e28/e182/e308/eudet-memo-2007-17.pdf page 3
   //EperMIP = 64.19e-06; //taken from Arabella which makes the comparison better for her
   //EperMIP = 52.8e-6; //taken from muon at 500 MeV

   // from G4 simulation EMM consistent with data
   EperMIP = 56.31e-06;  // 125pion MPV response
   double MIP2ParticleCalib = 1.33;  // muon 500MeV mean / pion 125GeV MPV (74.76e-06 / 56.31e-06)
   EperMIP = EperMIP * MIP2ParticleCalib;

   fProofFile = new TProofOutputFile(OUTFILENAME.c_str()); 
   fOutput->Add(fProofFile);  
   fout = fProofFile->OpenFile("RECREATE"); 

   ///x y position of the cells for TH2Poly
   std::pair<double,double> xy[200];


   double xPos[300] = {0.};
   double yPos[300] = {0.};

   for (int k = 0; k < 133; ++k) {
      xy[k] = geomc.position(k);
      //std::cout << "Coarse Cell[" << k << "] " << xy[k].first << ":" << xy[k].second<< std::endl;
      xPos[k] = xy[k].first;
      yPos[k] = xy[k].second;
      
   }

   /*
   int nbins = 200;
   double xmin = 0;
   double xmax = 2000;
   */

   int nbins = 40010;
   double xmin = -10;
   double xmax = 40000;


   hmap["h_e19_all"] = new TH1F("h_e19_all","h_e19_all",nbins, xmin, xmax);
   hmap["h_e19_all"]->Sumw2();

   hmap["h_e7_all"] = new TH1F("h_e7_all","h_e7_all",nbins, xmin, xmax);
   hmap["h_e7_all"]->Sumw2();
   
   hmap[Form("layerE_all")] = new TH1F(Form("layerE_all"),Form("layerE_all"),nbins ,xmin,xmax);
   hmap[Form("layerE_all")]->Sumw2();
   
   
   hmap[Form("allcellE_all")] = new TH1F(Form("allcellE_all"),Form("allcellE_all"),nbins ,xmin,xmax);
   hmap[Form("allcellE_all")]->Sumw2();
   
   hmap[Form("allcellE_2mip_all")] = new TH1F(Form("allcellE_2mip_all"),Form("allcellE_2mip_all"),nbins ,xmin,xmax);
   hmap[Form("allcellE_2mip_all")]->Sumw2();

   hmap[Form("h_eAll_all")] = new TH1F(Form("h_eAll_all"),Form("h_eAll_all"),nbins ,xmin,xmax);
   hmap[Form("h_eAll_all")]->Sumw2();

   for(int il=1; il<=nlayers; il++){
     
          nbins = 40010;
     xmin = -10;
     xmax = 40000;
     
     hmap[Form("h_eMax_L%d",il)] = new TH1F(Form("h_eMax_L%d",il),Form("h_eMax_L%d",il),nbins, xmin, xmax);
     hmap[Form("h_eMax_L%d",il)]->Sumw2();

     hmap[Form("h_e7_L%d",il)] = new TH1F(Form("h_e7_L%d",il),Form("h_e7_L%d",il),nbins, xmin, xmax);
     hmap[Form("h_e7_L%d",il)]->Sumw2();
     
     hmap[Form("h_e19_L%d",il)] = new TH1F(Form("h_e19_L%d",il),Form("h_e19_L%d",il),nbins, xmin, xmax);
     hmap[Form("h_e19_L%d",il)]->Sumw2();

     
    hmap[Form("h_eAll_L%d",il)] = new TH1F(Form("h_eAll_L%d",il),Form("h_eAll_L%d",il),nbins ,xmin,xmax);
    hmap[Form("h_eAll_L%d",il)]->Sumw2();
    

     
     /*
     nbins = 100;
     xmin = -1;
     xmax = 1.5;
     */

     nbins = 5000;
     xmin = -5;
     xmax = 5;
     
     hmap[Form("h_E1oE7_L%d",il)] = new TH1F(Form("h_E1oE7_L%d",il),Form("h_e1by7_L%d",il),nbins, xmin, xmax);
     hmap[Form("h_E1oE7_L%d",il)]->Sumw2();
     
     hmap[Form("h_E1oE19_L%d",il)] = new TH1F(Form("h_E1oE19_L%d",il),Form("h_e1by19_L%d",il),nbins, xmin, xmax);
     hmap[Form("h_E1oE19_L%d",il)]->Sumw2();

     int  nbinsx = 200;
     xmin = -5;
     xmax = 1;
     int nbinsy = 200;
     double ymin = -5;
     double ymax = 1.;
     
     h2Dmap[Form("h2_f7_f19_L%d",il)] = new TH2D(Form("h2_f7_f19_L%d",il),Form("h2_f7_f19_L%d",il),nbinsx, xmin,xmax, nbinsy,ymin,ymax);
     h2Dmap[Form("h2_f7_f19_L%d",il)]->Sumw2();
     
     hmap[Form("h_logF7_L%d",il)] = new TH1F(Form("h_logF7_L%d",il),Form("h_logF7_L%d",il),nbinsx, xmin, xmax);
     hmap[Form("h_logF7_L%d",il)]->Sumw2();
     
     hmap[Form("h_logF19_L%d",il)] = new TH1F(Form("h_logF19_L%d",il),Form("h_logF19_L%d",il),nbinsy, ymin, ymax);
     hmap[Form("h_logF19_L%d",il)]->Sumw2();

     for(int iev=0; iev<=20; iev++)
       {
	 hmap_poly[Form("logF19_L%d_ev%d",il,iev)] = new TH2Poly();
	 hmap_poly[Form("logF19_L%d_ev%d",il,iev)]->AddBin(133, xPos, yPos);
       }


   }//for(int il=0; il<=nlayers; il++)


   ///Energy VS nMIPS
   //hmap[Form("h_EnVSX0")] = new TH1F("h_EnVSX0","",);
   

   if(debug) cout<<"Inside SlaveBegin - All Initialization done!!!"<<endl;

}

Bool_t HGCTB::Process(Long64_t entry)
{

  using namespace std;
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


  if(debug) cout<<"Inside Process!!!"<<endl;
   
   fReader.SetEntry(entry);

   //std::cout<<"Number of cells is "<<simHitCellEnE.GetSize()<<std::endl;


   if(debug) cout<<"Event number is "<<entry<<endl;
   //cout<<"Event number "<<entry<<endl;
   if(TMath::Abs(*xBeam)<2 && TMath::Abs(*yBeam)<2){
     std::vector<unsigned int> CellId;
     std::vector<float> CellE;
     for(int ii=0; ii<simHitCellIdE.GetSize(); ii++){
       
       CellId.push_back(simHitCellIdE[ii]);
       CellE.push_back(simHitCellEnE[ii]);
     }
     
     double mipthresh = 2.;

     double totE19 = 0;
     double totE7 = 0;

     double totLE = 0;
     double totalcellE_nomip = 0;
     double totalcellE = 0;
     
     if(debug) cout<<"MIP threshold is "<<mipthresh<<endl;

     if(debug) cout<<"simHitLayEn2E.GetSize() is "<<simHitLayEn2E.GetSize()<<endl;
     for (int ilayer=0; ilayer<simHitLayEn2E.GetSize(); ilayer++){
       //unsigned int locMaxId = ht1.localMax( (simHitCellIdE), (simHitCellEnE), ilayer+1 );
       //double clusterE2 = ht1.cluster( (simHitCellIdE), (simHitCellEnE), locMaxId, 2 );
       
       if(debug) cout<<"ilayer  is "<<ilayer<<endl;


       //hmap[Form("layerE_L%d",ilayer)]->Fill(simHitLayEn2E[ilayer]);


       unsigned int locMaxId = ht1.localMax( (CellId), (CellE), ilayer+1 );
       double clusterEMax = ht1.cluster( (CellId), (CellE), locMaxId, 0, EperMIP, mipthresh );
       double clusterE7 = ht1.cluster( (CellId), (CellE), locMaxId, 1, EperMIP, mipthresh );
       double clusterE19 = ht1.cluster( (CellId), (CellE), locMaxId, 2, EperMIP, mipthresh );


       double allcell_nomip = ht1.cluster( (CellId), (CellE), locMaxId, 7, EperMIP, 0 );
       double allcell = ht1.cluster( (CellId), (CellE), locMaxId, 7, EperMIP, mipthresh );

       double allL = simHitLayEn2E[ilayer];
       //     clusterEn2E[ilayer]=clusterE2;
       //cout<<"cluster energy in 19 cells is "<<clusterE19/EperMIP<<endl;
     
       totE19 += clusterE19;
       totE7 += clusterE7;
       
       
       totLE += simHitLayEn2E[ilayer];
       totalcellE_nomip += allcell_nomip;
       totalcellE += allcell;


       hmap[Form("h_eAll_L%d",ilayer+1)]->Fill(allcell);

       double e1Overe7 = clusterEMax/clusterE7;
       double e1Overe19 = clusterEMax/clusterE19;
       double F7 = log( (clusterE7-clusterEMax)/clusterE7 );
       double F19 = log( (clusterE19-clusterE7)/clusterE19 );

       if(debug) cout<<"called HGCALplots.cc "<<endl;
       
       hmap[Form("h_eMax_L%d",ilayer+1)]->Fill(clusterEMax/EperMIP);

       if(debug) cout<<"done with filling eMax "<<endl;

       //       if(clusterE7/EperMIP > SEVENCELLS_THRESHOLD )
       hmap[Form("h_e7_L%d",ilayer+1)]->Fill(clusterE7/EperMIP);

       if(debug) cout<<"done with filling e7 "<<endl;

       //       if(clusterE19/EperMIP > NINETEENCELLS_THRESHOLD )
       hmap[Form("h_e19_L%d",ilayer+1)]->Fill(clusterE19/EperMIP);
    
       if(debug) cout<<"done with filling e19 "<<endl;

       hmap[Form("h_E1oE7_L%d",ilayer+1)]->Fill(e1Overe7);
       if(debug) cout<<"done with filling e1/e7 "<<endl;

       hmap[Form("h_E1oE19_L%d",ilayer+1)]->Fill(e1Overe19);
       if(debug) cout<<"done with filling e1/e19 "<<endl;

       hmap[Form("h_logF7_L%d",ilayer+1)]->Fill(F7);
       if(debug) cout<<"done with filling Log(e7) "<<endl;

       hmap[Form("h_logF19_L%d",ilayer+1)]->Fill(F19);
       if(debug) cout<<"done with filling Log(e19) "<<endl;

       h2Dmap[Form("h2_f7_f19_L%d",ilayer+1)]->Fill(F7,F19);
       if(debug) cout<<"done with filling Log(e7)VSLog(e19) "<<endl;

       /*
       if(entry<20){
	 ///try a TH2PolyMap
	 for(int icell=0; icell<CellId.size(); icell++){
	   std::pair<double,double> xy = geomc.position(CellId[icell]);
	   hmap_poly[Form("logF19_L%d_ev%d",ilayer,(int)entry)]->Fill(xy.first, xy.second, CellE[icell]);
	   //hmap_poly[Form("logF19_L%d_ev%d",ilayer,0)]->Fill(xy.first, xy.second, CellE[icell]);
	 }//for(int icell=0; icell<CellId.size(); icell++)
       }//if(entry<20)
       */

       //cout<<"Layer 1 energy "<<simHitLayEn2E[ilayer]<<endl;
     }//for (int ilayer=0; ilayer<simHitLayEn2E.GetSize(); ilayer++)

     hmap["h_e19_all"]->Fill(totE19/EperMIP);
     hmap["h_e7_all"]->Fill(totE7/EperMIP);
       

     hmap[Form("layerE_all")]->Fill(totLE/EperMIP);
     hmap[Form("allcellE_all")]->Fill(totalcellE_nomip/EperMIP);
     hmap[Form("allcellE_2mip_all")]->Fill(totalcellE/EperMIP);

     hmap[Form("h_eAll_all")]->Fill(totalcellE/EperMIP);
     
   }//if(fabs(xBeam)<2 && fabs(yBeam)<2)
     
   return kTRUE;
}
   
   void HGCTB::SlaveTerminate()
   {
     // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  
  fout->cd();
  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    hmap[it->first]->Write();
  }


  for(map<string,TH2D*>::iterator it = h2Dmap.begin(); it != h2Dmap.end(); ++it) {
    h2Dmap[it->first]->Write();
  }


  for(map<string,TH2Poly*>::iterator it = hmap_poly.begin(); it != hmap_poly.end(); ++it) {
    hmap_poly[it->first]->Write();
  }

  fout->Write();
  fout->Close();


}

void HGCTB::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  cout<<"Back to master - Everything is done now ...."<<endl;
}
