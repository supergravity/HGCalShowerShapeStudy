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
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

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
  outfile=new TFile("output.root","RECREATE");
  outTree = new TTree("DataTree","Data");
  outTree->Branch("Mgg",&m_Mgg,"Mgg/F");
  outTree->Branch("Mjj",&m_Mjj,"Mjj/F");

  TCanvas *c1  = new TCanvas("c1","two_point_correlation",100,10,700,500);
  //
  //* Define some histograms
  //
  int   binshErec =50  ; 
  float starthErec=90. ;
  float stophErec =110.;
  TH1F* Ereco =new TH1F("Ereco","",binshErec,starthErec,stophErec);
  TH1F* hLayerEnergySum =new TH1F("hLayerEnergySum","",100,0,0.1);
  TH1D* dRdist =new TH1D("dRdist","",100,0,12.);
  for(int iL=0; iL<NLAYERS; iL++) {
    char histoName[20];
    sprintf(histoName,"hetaphiprof_%d",iL);
    hetaphiprof[iL] =new TH2D(histoName,"",17,-8.5,8.5,17,-8.5,8.5);//cm
    sprintf(histoName,"dRprof_%d",iL);
    dRprof[iL] =new TH1D(histoName,"",15,-0.5,14.5);
  }
  TH2D* dRVlayer =new TH2D("dRVlayer","",10,-0.5,9.5,31,0,31);     //cm
  TH1D* dRprofAll =new TH1D("dRprofAll","",15,-0.5,14.5);//cm
  TH1D* shDepthAbs=new TH1D("shDepthAbs","",60,0.,30.);

  const int bin_size = 50;
  //========== Two Point Correlation ==========
    for (int iL = 0; iL < NLAYERS ;iL++)
    {
         char histoName[20];
         sprintf(histoName,"hEiEj_Multiplication_of_layer%d",iL);
         hEiEj_Multiplication [iL] = new TH1D (histoName,"",100,0,10);
         sprintf(histoName,"hDist_Two_hits_of_layer%d",iL);
         hdist_Two_hits [iL] = new TH1D (histoName,"",100,0,15);
         sprintf(histoName,"hTwo_point_correlation_of_layer%d",iL);
         Two_point_correlation [iL] = new TH2D (histoName,"",bin_size,0,15,100,0,70);
         //sprintf(histoName,"Two_point_correlation_of_layer%d",iL);//Don't know what is it for 
         Two_point_correl [iL] = new TGraph ();//
         p1               [iL] = new TProfile ("hProf","Profile of Energy Multiplication versus Distance btw Hits",bin_size,0,15,0,70); // rename p1 
         //r[iL] = new TFitResultPtr();
    }
  //========= Fitting the Profile ========== 
   TF1 *f1_profFit = new TF1("f1_profFit","[0]+[1]*1/x");
   TF1 *f2_profFit = new TF1("f2_profFit","[0]*exp([1]*x+[2]*(x**2))");
   TF1 *f3_profFit = new TF1("f3_profFit","[0]*exp([1]*x)+[2]*exp([3]*x)",0.5,15);  
   TF1 *f4_profFit = new TF1("f4_profFit","expo");  
   TF1 *f5_profFit = new TF1("f5_profFit","[0]*exp([1]*x+[2]*(x**2)+[3]*(x**3))",0.5,15);
   TF1 *f6_profFit = new TF1("f6_profFit","[0]*(1/x)*exp([1]*x+[2]*(x**2)+[3]*(x**3))");
   TF1 *f7_profFit = new TF1("f7_profFit","[0]*exp([1]*x)+[2]*exp([3]*x)+[4]*exp([5]*x)",0.5,15);
  
   TFitResultPtr r ;

  //
  //* Counters
  //
  int EventsPassed=0;

  //
  //* Loop over Events (!)
  //
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Entries: " << nentries << std::endl;
  for (Long64_t entry=0;entry< nentries ;++entry){
    if(entry%100==0) std::cout << "Processed ... " << entry 
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

    //
    //* Calculate Shower Depth
    //
    double shDepth = 0.0;
    double sumEmax = 0.0;
    for(int iL=0; iL<NLAYERS; iL++) {
      double lX0s=Layer.at(iL)->GetX0depth(iL);
      double emax=Layer.at(iL)->GetEneMax();
      //printf("iL:%u\n",iL);
      shDepth += (emax*lX0s);
      sumEmax += emax;
      //cout << "Layer number = " << iL+1 << " X0:" << lX0s << " Emax:" << emax << endl;
    }
    //getchar();
    shDepth = shDepth/sumEmax;//shower depth using emax in layer.
    shDepthAbs->Fill(shDepth);

    //
    //* Calculate lateral shapes:
    //
    //cout << "lay   Xmax    Ymax    dx    dy  EGeV  Emax cell" << endl;
    for(int iL=0; iL<NLAYERS; iL++) {
      int nhits=Layer.at(iL)->GetNhits();
      
      for(int ih=0; ih<nhits; ih++) {
        int    celln=Layer.at(iL)->GetCellNo(ih);
        double xhit=Layer.at(iL)->GetXHit(ih);
        double yhit=Layer.at(iL)->GetYHit(ih);
        double xmax=Layer.at(iL)->GetXMax();
        double ymax=Layer.at(iL)->GetYMax();
        if(!(fabs(xmax)<1.0 && fabs(ymax)<1.0)) continue;//look only in the center.
        double EGeV=Layer.at(iL)->GetErawHit(ih);
        //cout << "EGeV: " << EGeV << endl;
        double lX0s=Layer.at(iL)->GetX0depth(iL);
        double emax=Layer.at(iL)->GetEneMax();
        //double dR  =sqrt(deta*deta+dphi*dphi);
        double dx  =Layer.at(iL)->GetDX(ih)*MMtoCM; 
        double dy  =Layer.at(iL)->GetDY(ih)*MMtoCM;
        double dRxy=Layer.at(iL)->GetDRxy(ih)*MMtoCM;
        double areaInPadUnits = GetAreaInPadUnits(iL,dRxy);
        double tfactor = (lX0s/shDepth);//norm shower depth weight
        double t       = 1.0;
        double cterm   = 0.1;
        double efactor = 1.+log(emax);//what is efactor?
        
        efactor=0;
        // cout << "The following parameters are the contents: " << endl << 
        //         "Layer number: " << iL+1 <<"   "<< "number of hits of this event in the layer: "<< nhits << endl <<
        //         "Cell number: "  << celln<<"   "<< "the x coordinate of the hit?: "<< xhit <<"    "<< "the y coordinate of the hit?: "<< yhit << endl <<
        //         "the X coordinate of the maximum hit of the layer?: "<< xmax <<"     "<< "the Y coordinate of the maximum hit of the layer?: "<< ymax << endl <<
        //         "the Raw Energy of the layer of the cell of the hit?: "<< EGeV << "     "<< "the thickness of the layer in unit of X0: " <<  lX0s << endl<<
        //          "    " << "the dx of the layer of the hit?: "<< dx << endl <<
        //         "the dy of the layer of the hit?: "<< dy << "   "<< "the distance between the hit and the seed which has the maximum energy" << dRxy <<endl<< 
        //         "the area in Pad Units?: " << areaInPadUnits << "    "<<"the normalized shower depth weighting?" <<tfactor <<endl;  
        //Layer.at(iL)->GetXMax() << " " 
        //     << Layer.at(iL)->GetYMax() << " " << dx << " " << dy
        //    << " " << EGeV << " " << emax << " " << celln << endl; 
         
        double cucuWX0 = Layer.at(iL)->GetCuX0()+Layer.at(iL)->GetCuWX0();
        double mfactor = ((0.55>cucuWX0) ? 0.5-cucuWX0 : 0.) ;
        if(iL==0) t=pow(tfactor,0.12+efactor);//
        if(iL==1) t=pow(tfactor,0.125+efactor);
        if(iL==2) t=pow(tfactor,0.1+efactor); //matches with 1
        if(iL==3) t=pow(tfactor,0.1+efactor); //matches with 1
        if(iL==4) t=pow(tfactor,0.7+efactor); //too high (but 0.7 worksOK)
        if(iL==5) t=pow(tfactor,0.1+efactor); //also matches with 1
        if(iL==6) t=pow(tfactor,0.7+efactor); //too high! (but 0.5 worksOK)
        if(iL==7) t=pow(tfactor,0.1+efactor); //also matches with 1
        //t=pow(tfactor,cterm+mfactor);
        t=1.0;
        //cout <<"l:"<<iL<<" lX0s:" << lX0s << " tfactor:"
        //     <<tfactor<< " t:"<<t<<" shDepth:"<<shDepth<<endl;
        //cout << "dR=" << dR << " dRxy="<<dRxy << endl;
        //hetaphiprof[iL]->Fill(dx,dy,mips);
        hetaphiprof[iL]->Fill(dx,dy,EGeV/ENEPERMIP);
        dRprof[iL]->Fill(dRxy/t,EGeV/areaInPadUnits);
        dRVlayer->Fill(dRxy,iL,EGeV/areaInPadUnits);
        dRprofAll->Fill(dRxy/t,EGeV/areaInPadUnits);
        dRdist->Fill(dRxy);
        //cout << "dRxy= " << dRxy << endl;
      }
    }
    //cout << "Press Return to continue ..." << endl;
    //getchar();

    double LayerEnergySum=0.0;
    //for (int ilayer=0;ilayer<simHitLayEn2E->size();++ilayer){
    //      LayerEnergySum+=simHitLayEn2E->at(ilayer);
    //cout << simHitLayEn2E->at(ilayer) << "  unit?" << endl;
    //}
    //cout << "Total layer energy: " << LayerEnergySum << endl;
    //getchar();

// ================= Two Point Correlation =================
    //Initialization 
    double E_two_point[8][10000] = {0};
    double dist_two_point[8][10000] = {0};
    int    Layer_counter[8] = {0};
    //double *E_two_point = new double [NLAYERS*10000]; 
    //double *dist_two_point = new double [NLAYERS*10000];
    
    double Ehit1    [NLAYERS][100] = {0};
    double Ehit2    [NLAYERS][100] = {0};
    double dist_1_X [NLAYERS][100] = {0};
    double dist_1_Y [NLAYERS][100] = {0};
    double dist_2_X [NLAYERS][100] = {0};
    double dist_2_Y [NLAYERS][100] = {0};
    
    //Filling the arrays
    for (int iL = 0;iL < NLAYERS; iL++)
    {
        int nhits = Layer.at(iL)->GetNhits();
        //printf("the number of hit in layer %d is: %d\n",iL,nhits);
        for (int ih = 0;ih < nhits; ih++)
        {
            Ehit1[iL][ih] = Layer.at(iL)->GetErawHit(ih);
            Ehit2[iL][ih] = Layer.at(iL)->GetErawHit(ih);
            dist_1_X[iL][ih] = Layer.at(iL)->GetXHit(ih);
            dist_1_Y[iL][ih] = Layer.at(iL)->GetYHit(ih);
            dist_2_X[iL][ih] = Layer.at(iL)->GetXHit(ih);
            dist_2_Y[iL][ih] = Layer.at(iL)->GetYHit(ih);
            
        }
    }

    // Computing the correlation    
    for (int iL = 0;iL < NLAYERS;iL++)   // the for loop for layers 
    {
        int nhits=Layer.at(iL)->GetNhits();
        int counter = 0;
        for (int ih1 = 0; ih1 < nhits;ih1++) // for loop for hits(index number 1)
        {
            for (int ih2 = 0; ih2 < ih1;ih2++) //for loop for hits(index number 2)
            {
                E_two_point[iL][counter] = Ehit1[iL][ih1]*Ehit2[iL][ih2];
                dist_two_point[iL][counter] = sqrt(pow(dist_1_X[iL][ih1]-dist_2_X[iL][ih2],2)+pow(dist_1_Y[iL][ih1]-dist_2_Y[iL][ih2],2));
                //cout << "The Multiplication of Energies of the hits "<<ih1<<" and "<<ih2<<" : "<<E_two_point[iL][counter]*pow(GEVTOMEV,2)<<"(MeV)^2"<<endl<<"The distance between the hits "<<ih1<<" and "<<ih2<<" : "<<dist_two_point[iL][counter]*MMtoCM<<" (CM)"<<endl<<"The x components of hit1 and hit2 are: "<<dist_1_X[iL][ih1]<<"(MM)"<<" and "<<dist_2_X[iL][ih2]<<"(MM)"<<endl<<"The y components of hit1 and hit2 are: "<<dist_1_Y[iL][ih1]<<"(MM)"<<" and "<<dist_2_Y[iL][ih2]<<"(MM)"<<endl;
                hEiEj_Multiplication   [iL]  -> Fill(E_two_point[iL][counter]*pow(GEVTOMEV,2));
                hdist_Two_hits         [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM);
                Two_point_correlation  [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM,E_two_point[iL][counter]*pow(GEVTOMEV,2));
                Two_point_correl       [iL]  -> SetPoint(counter,dist_two_point[iL][counter]*MMtoCM,E_two_point[iL][counter]*pow(GEVTOMEV,2));
                p1                     [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM,E_two_point[iL][counter]*pow(GEVTOMEV,2),1);
                counter = counter + 1;
                //cout << "counter: "<<counter<<endl<< "ih2: " <<ih2<<"\n"<<endl; 
            }
        }
        Layer_counter[iL] = counter;
    }// layer loop ends

    // Draw the plots //per event
    //TGraph  *gr[NLAYERS];
    //for (int iL = 0 ;iL < NLAYERS; iL++)
    //{
    //    gr[iL] = new TGraph();
    //    for (int idx = 0;idx < Layer_counter[iL]; idx++)
    //    {
    //        gr[iL] -> SetPoint(idx,dist_two_point[iL][idx],E_two_point[iL][idx]);
    //        cout <<idx<<endl<<dist_two_point[iL][idx]<<endl<<E_two_point[iL][idx]<<endl<<"\n"<<endl;
    //    }
    //    gr[iL] -> Draw("AP*");
    //    //c1 -> Draw("Same");
    //    c1 -> Print(Form("Two_point_correlation_Layer_%d_event_%d.png",iL,entry));
    //}
    
    //TGraph *.....


// ================= Two Point Correlation analysis end ================
    //
    //* Fill histograms
    //
    //Ereco->Fill(Etotal);
    hLayerEnergySum->Fill(LayerEnergySum);

    //
    //* Fill Ntuple
    //
    m_Mgg=1.0;
    m_Mjj=1.0;

    outTree->Fill();
    EventsPassed++;
    //
    //* Delete Layers (note that we cannot 'continue' after we build them)
    //
    Clear();

  }

  //* Normalize profile histograms per event
  for(int iL=0; iL<NLAYERS; iL++) {
    hetaphiprof[iL]->Scale(1./EventsPassed);
  }
  dRVlayer->Scale(1./EventsPassed);
  dRprofAll->Sumw2();
  dRprofAll->Scale(1./EventsPassed);
  //========== Make Plots ==========
  
     double chi_sqr = 0.0;

    for(int iL = 0;iL < 1; iL++)
    {      
        // hEiEj_Multiplication [iL] -> SetXTitle("(MeV)^2");
        // hEiEj_Multiplication [iL] -> SetYTitle("Number_of_Hits");
        // hEiEj_Multiplication [iL] -> Draw();
        // c1 -> Print(Form("EiEj_Multiplication_Layer_%d_BinSize_%d.png",iL,bin_size));
        // hdist_Two_hits       [iL] -> SetXTitle("CM");
        // hdist_Two_hits       [iL] -> SetYTitle("Number_of_Hits");
        // hdist_Two_hits       [iL] -> Draw();
        // c1 -> Print(Form("Distance_of_two_Hits_Layer_%d_BinSize_%d.png",iL,bin_size));
        // Two_point_correlation[iL] -> SetXTitle("CM");
        // Two_point_correlation[iL] -> SetYTitle("(MeV)^2");
        // Two_point_correlation[iL] -> Draw();
        // c1 -> Print(Form("Two_point_correlation_Layer_%d_BinSize_%d.png",iL,bin_size));   
        // Two_point_correl     [iL] -> GetXaxis() -> SetTitle("CM");
        // Two_point_correl     [iL] -> GetYaxis() -> SetTitle("(MeV)^2");
        // Two_point_correl     [iL] -> Draw("AP*");
         ////c1 -> SetLogy();
         //c1 -> Print(Form("Graph_Two_point_correlation_yLogscale_Layer_%d_BinSize_%d.png",iL,bin_size));
         //c1 -> Print(Form("Graph_Two_point_correlation_Layer_%d_BinSize_%d.png",iL,bin_size));
         p1                   [iL] -> Draw("SAME");
         //c1 -> SetLogy();
         if (iL == 0)
         {
             
           //r  =   p1                   [iL] -> Fit("f5_profFit","RS");
           r  =   p1                   [iL] -> Fit("f7_profFit","RS");
           f7_profFit                       -> SetLineColor(iL+3);    
         }
         else
         {
           //r =  p1                   [iL] -> Fit("f3_profFit","RS");
           r =  p1                   [iL] -> Fit("f7_profFit","RS");
           f7_profFit                     -> SetLineColor(iL+3);    
         }
         r                         -> Print("V");  
         p1                   [iL] -> SetXTitle("CM");
         p1                   [iL] -> SetYTitle("(MeV)^2");
         p1                   [iL] -> SetAxisRange(0.1,5,"Y");
         p1                   [iL] -> SetMarkerStyle(21+iL);
         p1                   [iL] -> SetStats(kFALSE);//
         //p1                   [iL] -> SetMarkerColor(iL+1)
         //p1                   [iL] -> //Legend adding chi_sqr value 
         //c1 -> Print(Form("Energy_Multiplication_yLogscale_Fit_Exp_versus_HitDis_Layer%d_BinSize_%d.png",iL,bin_size));
         //c1 -> Print(Form("Energy_Multiplication_versus_HitDis_Fit_Event_28000_Exp_Layer%d__BinSize_%d.png",iL,bin_size));
    }
    c1 -> SetLogy();
    c1 -> Print(Form("Energy_Multiplication_versus_HitDis_Fit_Event_28000_Exp_Total_Layer__BinSize_%d.png",bin_size));

  cout << "PassedEvents :" << EventsPassed << endl;
  cout << "Entries      :" << nentries << endl;
  cout << "Efficiency is :" << ((float) EventsPassed)/nentries << endl;

  outfile->Write("",TObject::kOverwrite);
  outfile->Close();
}




#endif // #ifdef makePlots_cxx
 

