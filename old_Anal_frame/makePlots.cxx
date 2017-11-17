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
#include "TH2D.h"
#include  <cmath>
#include <iostream>
#include <algorithm>
#include "TProfile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "makePlots_function.h"

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

    // Tell input file MC/Data
    printf("doTruth: %d",doTruth); 

    // Set branch addresses and branch pointers
    fCurrent = -1;
    fChain->SetMakeClass(1);

    if(doTruth) 
    {
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
    } 
    else 
    {
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

//Analysis Starts

void makePlots::Loop()
{


    Init();

    const int input_energy = 125;
    const int bin_size = 200;
    const int bin_number = 50;
    const int bin_number_2 = 15;

    // Open the File

    TFile* outfile;
    if (doTruth)
    {
        //outfile=new TFile(Form("output_MC_%dGeV.root",input_energy),"RECREATE");
        outfile=new TFile(Form("output_electron_MC_%dGeV.root",input_energy),"RECREATE");
    }
    else
    {
        outfile=new TFile(Form("output_Data_%dGeV.root",input_energy),"RECREATE");
    }
    outTree = new TTree("DataTree","Data");
    outTree->Branch("Mgg",&m_Mgg,"Mgg/F");
    outTree->Branch("Mjj",&m_Mjj,"Mjj/F");

    //
    //* Define some histograms
    //

    double Elayertotal [NLAYERS]      = {0};
    /*
    TH1D* dRdist =new TH1D("dRdist","",100,0,12.);
    TH2D* dRVlayer =   new TH2D("dRVlayer","",10,-0.5,9.5,31,0,31);     //cm
    TH1D* dRprofAll =  new TH1D("dRprofAll","",15,-0.5,14.5);//cm
    */
    TH1D* shDepthAbs = new TH1D("shDepthAbs","",60,0.,30.);





    //========== Two Point Correlation ==========
    if(doTruth)
    {
        helayerRawtotal  = new TH1D ("The energy distribution of total layer per Event of MC","",bin_number,0,200);
        helayerRawtotal -> Sumw2();
        hene_total_layer = new TH1D ("The total layer_energy distribution of MC","",8,0,8); 
        hene_total_layer -> Sumw2();
    }
    else
    {
        helayerRawtotal  = new TH1D ("The energy distribution of total layer per Event of Data","",bin_number,0,200);
        helayerRawtotal -> Sumw2();
        hene_total_layer = new TH1D ("The total layer_energy distribution of Data","",8,0,8); 
        hene_total_layer -> Sumw2();
    }
    for (int iL = 0; iL < NLAYERS ;iL++)
    {
         if (doTruth)
         {
              char histoName[20];
              char histoName_2[20];
              sprintf(histoName,"htotal_ene_of_layer%d_MC",iL);
              //helayerRaw [iL] = new TH1D (histoName,"",bin_number,0,50);
              helayerRaw [iL] = new TH1D (histoName,"",25,0,50);
              helayerRaw [iL] -> Sumw2 ();
              sprintf(histoName,"hEiEj_Multi_normal_of_layer%d_MC",iL);
              hEiEj_Multi_normal [iL] = new TH1D (histoName,"",bin_number,0,0.04);
              hEiEj_Multi_normal [iL] -> Sumw2();
              sprintf(histoName,"hDist_Two_hits_of_layer%d_MC",iL);
              hdist_Two_hits [iL] = new TH1D (histoName,"",bin_number_2,0,15);
              hdist_Two_hits [iL] -> Sumw2();
              //sprintf(histoName,"Two_point_correlation_of_layer%d",iL);//Don't know what is it for 
              //EiEj_nor_rel_dis [iL] = new TGraph ();//
              sprintf(histoName,"hTwo_point_correlation_of_layer%d_MC",iL);
              EiEj_nor_rel_dis[iL] = new TH2D (histoName,"",bin_size,0,15,100,0,1);
              sprintf(histoName,"hTwo_point_correlation_pad_of_layer%d_MC",iL);
              EiEj_nor_dis_pad[iL] = new TH2D (histoName,"",bin_size,0,8,100,0,1);
              sprintf(histoName,"hProf_layer_pad_layer%d_MC",iL);
              EiEj_nor_dis_avg[iL] = new TH2D (histoName,"",bin_size,0,8,100,0,1);
              //sprintf(histoName,"Profile of Energy Multiplication versus Distance btw Hits_of_layer%d_MC",iL);
              sprintf(histoName_2,"hProf_layer%d_MC",iL);
              p1               [iL] = new TProfile (histoName_2,"",bin_size,0,15,0,70); // rename p1 
              EiEj_prof        [iL] = new TProfile (histoName,"",bin_size,0,8,0,1);
              EiEj_nor_dis_graph [iL] = new TGraphErrors();
              char histoName_6[40];
              for (int r = 0; r < Ring;r++)
              {
                  sprintf(histoName_6,"hist_layer_%d_ring_%d_MC",iL+1,r+1); 
                  Pad_ring_hist [iL][r] = new TH1D (histoName_6,histoName_6,bin_size,0,0.1);
              }
                  //sprintf(histoName_2,"hProf_layer_pad%d_MC",iL);
              //p2               [iL] = new TProfile (histoName_2,"",bin_size,0,8,0,70); // rename p1 
              //r[iL] = new TFitResultPtr();
        }
        else
        {
              char histoName[20];
              char histoName_2[20];
              sprintf(histoName,"htotal_ene_of_layer%d_Data",iL);
              //helayerRaw [iL] = new TH1D (histoName,"",bin_number,0,50);
              helayerRaw [iL] = new TH1D (histoName,"",25,0,50);
              helayerRaw [iL] -> Sumw2 ();
              sprintf(histoName,"hEiEj_Multi_normal_of_layer%d_Data",iL);
              hEiEj_Multi_normal [iL] = new TH1D (histoName,"",bin_number,0,0.04);
              hEiEj_Multi_normal [iL] -> Sumw2();
              sprintf(histoName,"hDist_Two_hits_of_layer%d_Data",iL);
              hdist_Two_hits [iL] = new TH1D (histoName,"",bin_number_2,0,15);
              hdist_Two_hits [iL] -> Sumw2();
              //sprintf(histoName,"Two_point_correlation_of_layer%d",iL);//Don't know what is it for 
              //EiEj_nor_rel_dis [iL] = new TGraph ();//
              sprintf(histoName,"hTwo_point_correlation_of_layer%d_Data",iL);
              EiEj_nor_rel_dis[iL] = new TH2D (histoName,"",bin_size,0,15,100,0,1);
              sprintf(histoName,"hTwo_point_correlation_pad_of_layer%d_Data",iL);
              EiEj_nor_dis_pad[iL] = new TH2D (histoName,"",bin_size,0,8,100,0,1);
              sprintf(histoName,"hProf_layer_pad_layer%d_Data",iL);
              EiEj_nor_dis_avg[iL] = new TH2D (histoName,"",bin_size,0,8,100,0,1);
               //sprintf(histoName,"Profile of Energy Multiplication versus Distance btw Hits_of_layer%d_Data",iL);
              sprintf(histoName_2,"hProf_layer%d_Data",iL);
              p1               [iL] = new TProfile (histoName_2,"",bin_size,0,15,0,70); // rename p1 
              EiEj_prof        [iL] = new TProfile (histoName,"",bin_size,0,8,0,1);
              EiEj_nor_dis_graph [iL] = new TGraphErrors();
              //sprintf(histoName_2,"hProf_layer_pad%d_Data",iL);
              //p2               [iL] = new TProfile (histoName_2,"",bin_size,0,8,0,70); // rename p1 
              char histoName_6[40];
              for (int r = 0; r < Ring;r++)
              {
                  sprintf(histoName_6,"hist_layer_%d_ring_%d_Data",iL+1,r+1); 
                  Pad_ring_hist [iL][r] = new TH1D (histoName_6,histoName_6,bin_size,0,0.1);
              }
              //r[iL] = new TFitResultPtr();
         
        }
    }
    //======== Setting TH2D plots  ==========

    
    //======== Declaring the plots ========== 
        
    
    c1 = new TCanvas("c1","",100,10,700,500);


    //========= Fitting the Profile ========== 
    TF1 *f1_profFit = new TF1("f1_profFit","[0]+[1]*1/x");
    TF1 *f2_profFit = new TF1("f2_profFit","[0]*exp([1]*x+[2]*(x**2))");
    TF1 *f3_profFit = new TF1("f3_profFit","[0]*exp([1]*x)+[2]*exp([3]*x)",0.5,15);  
    TF1 *f4_profFit = new TF1("f4_profFit","expo");  
    TF1 *f5_profFit = new TF1("f5_profFit","[0]*exp([1]*x+[2]*(x**2)+[3]*(x**3))",0.5,15);
    TF1 *f6_profFit = new TF1("f6_profFit","[0]*(1/x)*exp([1]*x+[2]*(x**2)+[3]*(x**3))");
    TF1 *f7_profFit = new TF1("f7_profFit","[0]*exp([1]*x)+[2]*exp([3]*x)+[4]*exp([5]*x)",0.5,15);

    TFitResultPtr r ;

    //std::vector< std::pair<double,double> >* test[7];
    //for ( int i=0;i<7;++i )
    //{
    //    test[i] = new std::vector< std::pair<double,double> >();
    //}
    //double* test1[7];
    //double* test2[7];
    //int MAX_NUMBER = nhits * nhits /2;
    //for ( int i=0; i<7; ++i )
    //{
    //    test1[i] = new double[MAX_NUMBER];
    //    test2[i] = new double[MAX_NUMBER];
    //}
    std::vector<std::pair<double, double> >*  prof_pad_1[7];
    std::vector<std::pair<double, double> >*  prof_pad_2[7];
    std::vector<std::pair<double, double> >*  prof_pad_3[7];
    std::vector<std::pair<double, double> >*  prof_pad_4[7];
    std::vector<std::pair<double, double> >*  prof_pad_5[7];
    std::vector<std::pair<double, double> >*  prof_pad_6[7];
    std::vector<std::pair<double, double> >*  prof_pad_7[7];
    std::vector<std::pair<double, double> >*  prof_pad_8[7];
    for (int ring =  0; ring < 7;ring++)
    {
        prof_pad_1[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_2[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_3[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_4[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_5[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_6[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_7[ring] = new std::vector<std::pair<double,double > > ();
        prof_pad_8[ring] = new std::vector<std::pair<double,double > > ();
    }
    // for (int ring = 0;ring < 7; ring++)
   // {
   //     
   // }
    //std::vector<std::pair<double, double> > prof_pad_2[7];
    //std::vector<std::pair<double, double> > prof_pad_3[7];
    //std::vector<std::pair<double, double> > prof_pad_4[7];
    //std::vector<std::pair<double, double> > prof_pad_5[7];
    //std::vector<std::pair<double, double> > prof_pad_6[7];
    //std::vector<std::pair<double, double> > prof_pad_7[7];
    //std::vector<std::pair<double, double> > prof_pad_8[7];

    std::vector<double> avg_pad_1;
    std::vector<double> avg_pad_2;
    std::vector<double> avg_pad_3;
    std::vector<double> avg_pad_4;
    std::vector<double> avg_pad_5;
    std::vector<double> avg_pad_6;
    std::vector<double> avg_pad_7;
    std::vector<double> avg_pad_8;

    std::vector<double> err_pad_1;
    std::vector<double> err_pad_2;
    std::vector<double> err_pad_3;
    std::vector<double> err_pad_4;
    std::vector<double> err_pad_5;
    std::vector<double> err_pad_6;
    std::vector<double> err_pad_7;
    std::vector<double> err_pad_8;
    //
    //* Counters
    //
    int EventsPassed = 0;

    //
    //* Loop over Events (!)
    //
    Long64_t nentries = fChain->GetEntries();

    std::cout << "Entries: " << nentries << std::endl;
    for (Long64_t entry = 0;entry < nentries ;++entry)
    {
        if(entry%100==0) std::cout << "Processed ... " << entry << "/" << nentries << " events" << std::endl;

        //std::cout << "entry : " << entry << std::endl;
        fChain->GetEntry(entry); 

        //
        //* How to "see" the contents of an event:
        //
        //Show(entry);

        //
        //* Build Layer objects
        //
        int ibuild = buildCompObjects();



    

        // ==========  Event Cuts Application ========== 
    
        // ===== Preselection =====
        if ( !Layer.size() ) continue;
        double xmax_pre = Layer.at(0)->GetXMax();
        double ymax_pre = Layer.at(0)->GetYMax();
   
        //cout <<"The xmax (seed) is at "<< ymax_pre << "mm"<<endl;
        //cout <<"The ymax (seed) is at "<< xmax_pre << "mm"<<endl;
        //getchar();

    
        if(!(fabs(xmax_pre)< 20.0) && !(fabs(ymax_pre)< 20.0))
        {   
            Clear();
            continue; 
        }
        
        // ==========  Event Cuts Application Ends ==========      




        //
        //* Calculate Shower Depth
        //
        double shDepth = 0.0;
        double sumEmax = 0.0;
        for(int iL=0; iL<NLAYERS; iL++) 
        {
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
        // Define energy distribution of total layer per event 
        double elayertot = 0.0;

        //cout << "lay   Xmax    Ymax    dx    dy  EGeV  Emax cell" << endl;
        for(int iL=0; iL<NLAYERS; iL++) 
        {
            int     nhits=Layer.at(iL)->GetNhits();
            double  xmax=Layer.at(iL)->GetXMax();
            double  ymax=Layer.at(iL)->GetYMax();
            double  emax=Layer.at(iL)->GetEneMax();
            double  elayer=0;//initialize
           //     std::cout << "reserve start\n";
           //     std::cout << "nhits = " << nhits << std::endl;
           //     std::cout << "element number : " << prof_pad_1[0].size() << std::endl;
           //     std::cout << "element number6 : " << prof_pad_7[0].size() << std::endl;
           // for ( int i=0;i<7; ++i )
           // {
           //     if ( iL == 0 )
           //     prof_pad_1[i].reserve( nhits*nhits/2 );
           //     if ( iL == 1 )
           //     prof_pad_2[i].reserve( nhits*nhits/2 );
           //     if ( iL == 2 )
           //     prof_pad_3[i].reserve( nhits*nhits/2 );
           //     if ( iL == 3 )
           //     prof_pad_4[i].reserve( nhits*nhits/2 );
           //     if ( iL == 4 )
           //     prof_pad_5[i].reserve( nhits*nhits/2 );
           //     if ( iL == 5 )
           //     prof_pad_6[i].reserve( nhits*nhits/2 );
           //     if ( iL == 6 )
           //     prof_pad_7[i].reserve( nhits*nhits/2 );
           // }
           //     std::cout << "element number2 : " << prof_pad_1[0].size() << std::endl;
           //     std::cout << "reserve end\n";
            
            for(int ih=0; ih<nhits; ih++) 
            {
                int    celln=Layer.at(iL)->GetCellNo(ih);
                double xhit=Layer.at(iL)->GetXHit(ih);
                double yhit=Layer.at(iL)->GetYHit(ih);
                double EGeV=Layer.at(iL)->GetErawHit(ih);
                elayer=elayer+EGeV;
                //cout << "EGeV: " << EGeV << endl;
                //double dR  =sqrt(deta*deta+dphi*dphi);
                double lX0s=Layer.at(iL)->GetX0depth(iL);
                double dx  =Layer.at(iL)->GetDX(ih)*MMtoCM; 
                double dy  =Layer.at(iL)->GetDY(ih)*MMtoCM;
                double dRxy=Layer.at(iL)->GetDRxy(ih)*MMtoCM;
                double areaInPadUnits = GetAreaInPadUnits(iL,dRxy);
        
                //cout << "dR=" << dR << " dRxy="<<dRxy << endl;
                //dRprof[iL]->Fill(dRxy,EGeV/areaInPadUnits);
                //dRVlayer->Fill(dRxy,iL,EGeV/areaInPadUnits);
                //dRprofAll->Fill(dRxy,EGeV/areaInPadUnits);
                //dRdist->Fill(dRxy);
                //cout << "dRxy= " << dRxy << endl;
            }
            elayertot = elayertot + elayer;
            if((iL==0) && (elayer >= 2.*ENEPERMIP)) helayerRaw[0] -> Fill(elayer*GEVTOMEV);
            if((iL==1) && (elayer >= 2.*ENEPERMIP)) helayerRaw[1] -> Fill(elayer*GEVTOMEV);
            if((iL==2) && (elayer >= 2.*ENEPERMIP)) helayerRaw[2] -> Fill(elayer*GEVTOMEV);
            if((iL==3) && (elayer >= 2.*ENEPERMIP)) helayerRaw[3] -> Fill(elayer*GEVTOMEV);
            if((iL==4) && (elayer >= 2.*ENEPERMIP)) helayerRaw[4] -> Fill(elayer*GEVTOMEV);
            if((iL==5) && (elayer >= 2.*ENEPERMIP)) helayerRaw[5] -> Fill(elayer*GEVTOMEV);
            if((iL==6) && (elayer >= 2.*ENEPERMIP)) helayerRaw[6] -> Fill(elayer*GEVTOMEV);
            if((iL==7) && (elayer >= 2.*ENEPERMIP)) helayerRaw[7] -> Fill(elayer*GEVTOMEV);

        }
        //cout << "Press Return to continue ..." << endl;
        //cout << "energy distribution of total layer Event "<<entry<<" is: "<<elayertot*GEVTOMEV<<" MeV"<<endl;
        helayerRawtotal -> Fill(elayertot*GEVTOMEV);
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
        double dist_padunit [8][8][10000] = {0};
        int    Layer_counter[8] = {0};
        //double *E_two_point = new double [NLAYERS*10000]; 
        //double *dist_two_point = new double [NLAYERS*10000];

        double Ehit1       [NLAYERS][100] = {0};
        double Ehit2       [NLAYERS][100] = {0};
        double dist_1_X    [NLAYERS][100] = {0};
        double dist_1_Y    [NLAYERS][100] = {0};
        double dist_2_X    [NLAYERS][100] = {0};
        double dist_2_Y    [NLAYERS][100] = {0};


        //Filling the arrays
        for (Int_t iL = 0;iL < NLAYERS; iL++)
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
                Elayertotal [iL] += Ehit1[iL][ih]; 
                //cout<<"the_x_position_hit:   "<<dist_1_X[iL][ih]<<endl;
                //cout<<"the_y_position_hit:   "<<dist_1_Y[iL][ih]<<endl;
            }
            
        }
        //std::cout<<"checkpoint_01"<<endl;
        // Computing the correlation    
        for (int iL = 0;iL < NLAYERS;iL++)   // the for loop for layers 
        {
            int nhits=Layer.at(iL)->GetNhits();
            //std::cout << "Nhits: " << nhits << std::endl;
            int counter = 0;
            double pid_nor_variable = 0;
            for (int ih1 = 0; ih1 < nhits;ih1++) // for loop for hits(index number 1)
            {
                for (int ih2 = 0; ih2 < ih1;ih2++) //for loop for hits(index number 2)
                {
                    E_two_point[iL][counter] = Ehit1[iL][ih1]*Ehit2[iL][ih2];
                    dist_two_point[iL][counter] = sqrt(pow(dist_1_X[iL][ih1]-dist_2_X[iL][ih2],2)+pow(dist_1_Y[iL][ih1]-dist_2_Y[iL][ih2],2));
                    //cout << "The Multiplication of Energies of the hits "<<ih1<<" and "<<ih2<<" : "<<E_two_point[iL][counter]*pow(GEVTOMEV,2)<<"(MeV)^2"<<endl<<"The distance between the hits "<<ih1<<" and "<<ih2<<" : "<<dist_two_point[iL][counter]*MMtoCM<<" (CM)"<<endl<<"The x components of hit1 and hit2 are: "<<dist_1_X[iL][ih1]<<"(MM)"<<" and "<<dist_2_X[iL][ih2]<<"(MM)"<<endl<<"The y components of hit1 and hit2 are: "<<dist_1_Y[iL][ih1]<<"(MM)"<<" and "<<dist_2_Y[iL][ih2]<<"(MM)"<<endl;
                    pid_nor_variable    =  E_two_point[iL][counter]/pow(Elayertotal[iL],2);
                    if (1000*pid_nor_variable >= 0.002)
                    {
                        hEiEj_Multi_normal     [iL]  -> Fill(1000*pid_nor_variable);
                        hdist_Two_hits         [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM);
                        EiEj_nor_rel_dis       [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM,1000*pid_nor_variable);
                        p1                     [iL]  -> Fill(dist_two_point[iL][counter]*MMtoCM,1000*pid_nor_variable,1);// Getting the profile (the mean value of )
                        
                        //std::cout<<"checkpoint_02"<<endl;
                        //std::cout << " iL = " << iL << ", counter = " << counter << std::endl;
                        //std::cout << "test value = " << dist_two_point[iL][counter] << std::endl;
                        //std::cout << "test end\n";
                        if((10.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 19.))
                        {
                            //std::cout << "hi\n";
                            if (iL == 0)
                            {
                                prof_pad_1[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            //std::cout << "hi_2\n";
                            if (iL == 1)
                            {
                                prof_pad_2[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[0] -> push_back(std::pair<double,double>(1,1000*pid_nor_variable));
                                Pad_ring_hist[iL][0] -> Fill(1000*pid_nor_variable);
                            }
                            
                            EiEj_nor_dis_pad[iL] -> Fill(1,1000*pid_nor_variable);
                            //EiEj_prof       [iL] -> Fill(1,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(1,1000*pid_nor_variable); 
                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                        }
                        //std::cout<<"checkpoint_03"<<endl;
                        if((19.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 29.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[1] -> push_back(std::pair<double,double>(2,1000*pid_nor_variable));
                                Pad_ring_hist[iL][1] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(2,1000*pid_nor_variable);
                            //EiEj_prof       [iL] -> Fill(2,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(2,1000*pid_nor_variable); 
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                       }
//                        std::cout<<"checkpoint_04"<<endl;
//                        //{dis_pad[iL][1]-> push_back(std::pair<double,double>(dist_two_point[iL][counter],1000*pid_nor_variable));}
                        if((29.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 34.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[2] -> push_back(std::pair<double,double>(3,1000*pid_nor_variable));
                                Pad_ring_hist[iL][2] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(3,1000*pid_nor_variable);
                            //EiEj_prof       [iL] -> Fill(3,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(3,1000*pid_nor_variable); 
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                        }
//        std::cout<<"checkpoint_05"<<endl;
                        if((38.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 50.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[3] -> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[3] -> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[3] -> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[3] -> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[3]-> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[3]-> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[3]-> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[3]-> push_back(std::pair<double,double>(4,1000*pid_nor_variable));
                                Pad_ring_hist[iL][3] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(4,1000*pid_nor_variable);
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                            //EiEj_prof       [iL] -> Fill(4,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(4,1000*pid_nor_variable); 
                        }
//        std::cout<<"checkpoint_06"<<endl;
                        if((51.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 58.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[4] -> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[4]-> push_back(std::pair<double,double>(5,1000*pid_nor_variable));
                                Pad_ring_hist[iL][4] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(5,1000*pid_nor_variable);
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                            //EiEj_prof       [iL] -> Fill(5,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(5,1000*pid_nor_variable); 
                        }
//        std::cout<<"checkpoint_07"<<endl;
                        if((58.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 68.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[5]->push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[5]-> push_back(std::pair<double,double>(6,1000*pid_nor_variable));
                                Pad_ring_hist[iL][5] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(6,1000*pid_nor_variable);
                            //EiEj_prof       [iL] -> Fill(6,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(6,1000*pid_nor_variable); 
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                        }
                        if((68.<dist_two_point[iL][counter]) && (dist_two_point[iL][counter]< 77.))
                        {
                            if (iL == 0)
                            {
                                prof_pad_1[6] -> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 1)
                            {
                                prof_pad_2[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 2)
                            {
                                prof_pad_3[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 3)
                            {
                                prof_pad_4[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 4)
                            {
                                prof_pad_5[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 5)
                            {
                                prof_pad_6[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 6)
                            {
                                prof_pad_7[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            if (iL == 7)
                            {
                                prof_pad_8[6]-> push_back(std::pair<double,double>(7,1000*pid_nor_variable));
                                Pad_ring_hist[iL][6] -> Fill(1000*pid_nor_variable);
                            }
                            EiEj_nor_dis_pad[iL] -> Fill(7,1000*pid_nor_variable);
                            //EiEj_prof       [iL] -> Fill(7,1000*pid_nor_variable,1); 
                            EiEj_prof       [iL] -> Fill(7,1000*pid_nor_variable); 
//                            //p2              [iL] -> Fill(1,1000*pid_nor_variable,1);
                        }
                        //
                    }
                    counter = counter + 1;
                    //cout << "counter: "<<counter<<endl<< "ih2: " <<ih2<<"\n"<<endl; 
                }
            }
            Layer_counter[iL] = counter;// number of points  
        }// layer loop ends

    //    std::cout<<"checkpoint_08"<<endl;

        // ================= Two Point Correlation analysis end ================
        //
        //* Fill histograms
         
    
    
        //

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
        //std::cout << "chkpoint_081\n";
        Clear();

    }
        //std::cout<<"checkpoint_09"<<endl;
    Average_Calculation(prof_pad_1,avg_pad_1); 
    Average_Calculation(prof_pad_2,avg_pad_2); 
    Average_Calculation(prof_pad_3,avg_pad_3); 
    Average_Calculation(prof_pad_4,avg_pad_4); 
    Average_Calculation(prof_pad_5,avg_pad_5); 
    Average_Calculation(prof_pad_6,avg_pad_6); 
    Average_Calculation(prof_pad_7,avg_pad_7); 
    Average_Calculation(prof_pad_8,avg_pad_8); 
   
    error_on_mean(prof_pad_1,avg_pad_1,err_pad_1);
    error_on_mean(prof_pad_2,avg_pad_2,err_pad_2);
    error_on_mean(prof_pad_3,avg_pad_3,err_pad_3);
    error_on_mean(prof_pad_4,avg_pad_4,err_pad_4);
    error_on_mean(prof_pad_5,avg_pad_5,err_pad_5);
    error_on_mean(prof_pad_6,avg_pad_6,err_pad_6);
    error_on_mean(prof_pad_7,avg_pad_7,err_pad_7);
    error_on_mean(prof_pad_8,avg_pad_8,err_pad_8);


    for (int i = 0;i < 7;i++)
    {
        EiEj_nor_dis_graph [0] -> SetPoint(i,i+1,avg_pad_1[i]);
        EiEj_nor_dis_graph [0] -> SetPointError(i,0,err_pad_1[i]); 
        EiEj_nor_dis_graph [1] -> SetPoint(i,i+1,avg_pad_2[i]);
        EiEj_nor_dis_graph [1] -> SetPointError(i,0,err_pad_2[i]); 
        EiEj_nor_dis_graph [2] -> SetPoint(i,i+1,avg_pad_3[i]);
        EiEj_nor_dis_graph [2] -> SetPointError(i,0,err_pad_3[i]); 
        EiEj_nor_dis_graph [3] -> SetPoint(i,i+1,avg_pad_4[i]);
        EiEj_nor_dis_graph [3] -> SetPointError(i,0,err_pad_4[i]); 
        EiEj_nor_dis_graph [4] -> SetPoint(i,i+1,avg_pad_5[i]);
        EiEj_nor_dis_graph [4] -> SetPointError(i,0,err_pad_5[i]); 
        EiEj_nor_dis_graph [5] -> SetPoint(i,i+1,avg_pad_6[i]);
        EiEj_nor_dis_graph [5] -> SetPointError(i,0,err_pad_6[i]); 
        EiEj_nor_dis_graph [6] -> SetPoint(i,i+1,avg_pad_7[i]);
        EiEj_nor_dis_graph [6] -> SetPointError(i,0,err_pad_7[i]); 
        EiEj_nor_dis_graph [7] -> SetPoint(i,i+1,avg_pad_8[i]);
        EiEj_nor_dis_graph [7] -> SetPointError(i,0,err_pad_8[i]); 
    }
    vector<double>::iterator iter_1 = avg_pad_1.begin();
    vector<double>::iterator iend_1 = avg_pad_1.end();
    
    vector<double>::iterator iter_2 = avg_pad_2.begin();
    vector<double>::iterator iend_2 = avg_pad_2.end();
    vector<double>::iterator iter_3 = avg_pad_3.begin();
    vector<double>::iterator iend_3 = avg_pad_3.end();
    vector<double>::iterator iter_4 = avg_pad_4.begin();
    vector<double>::iterator iend_4 = avg_pad_4.end();
    vector<double>::iterator iter_5 = avg_pad_5.begin();
    vector<double>::iterator iend_5 = avg_pad_5.end();
    vector<double>::iterator iter_6 = avg_pad_6.begin();
    vector<double>::iterator iend_6 = avg_pad_6.end();
    vector<double>::iterator iter_7 = avg_pad_7.begin();
    vector<double>::iterator iend_7 = avg_pad_7.end();
    vector<double>::iterator iter_8 = avg_pad_8.begin();
    vector<double>::iterator iend_8 = avg_pad_8.end();
        
    vector<double>::iterator kter_1 = err_pad_1.begin();
    vector<double>::iterator kend_1 = err_pad_1.end();
    
    vector<double>::iterator kter_2 = err_pad_2.begin();
    vector<double>::iterator kend_2 = err_pad_2.end();
    vector<double>::iterator kter_3 = err_pad_3.begin();
    vector<double>::iterator kend_3 = err_pad_3.end();
    vector<double>::iterator kter_4 = err_pad_4.begin();
    vector<double>::iterator kend_4 = err_pad_4.end();
    vector<double>::iterator kter_5 = err_pad_5.begin();
    vector<double>::iterator kend_5 = err_pad_5.end();
    vector<double>::iterator kter_6 = err_pad_6.begin();
    vector<double>::iterator kend_6 = err_pad_6.end();
    vector<double>::iterator kter_7 = err_pad_7.begin();
    vector<double>::iterator kend_7 = err_pad_7.end();
    vector<double>::iterator kter_8 = err_pad_8.begin();
    vector<double>::iterator kend_8 = err_pad_8.end();
        /* 1
        vector<double>::iterator iter[8];
        vector<double>::iterator iend[8];
        for ( int i=0;i<8;++i)
        {
            iter[i] = avg_pag[i].begin();
            iend[i] = avg_pag[i].end  ();
        }
        // 2
#define aaa(x) avg_pad_##x
        aaa(3) -> avg_pad_3
        vector<double>::iterator iter = aaa(x).begin();

#define MAXEVENT 999
        double aaa[MAXEVENT];
        double aaa[999];
        */
        
        std::cout << "checkpoint_10"<<endl;
        //std::cout << "aaa   " <<  iter_1 << std::endl;
        std::cout << "aaa   " << *iter_1 << std::endl;
        std::cout << "avg_pad_1 size:"<<avg_pad_1.size()<<endl;
        std::cout << "avg1"<<avg_pad_1[0]<<endl;
       // while (iter_1 != iend_1)
       // {
       //     cout<<"avg1 element:"<<*iter_1<<endl; 
       //     iter_1++;
       // }
        
        int ring = 0;
        for (iter_1;iter_1 !=iend_1;iter_1++)
        {
            EiEj_nor_dis_avg    [0] -> Fill(ring+1,*iter_1); 
            cout<<"avg1 element:"<<*iter_1<<endl; 
            ring++;
        }
        
        ring = 0;
        for (iter_2;iter_2 !=iend_2;iter_2++)
        {
            EiEj_nor_dis_avg    [1] -> Fill(ring+1,*iter_2); 
            cout<<"avg2 element:"<<*iter_2<<endl; 
            ring++;
        }
         ring = 0;
        for (iter_3;iter_3 !=iend_3;iter_3++)
        {
            EiEj_nor_dis_avg    [2] -> Fill(ring+1,*iter_3); 
            cout<<"avg3 element:"<<*iter_3<<endl; 
            ring++;
        }
         ring = 0;
        for (iter_4;iter_4 !=iend_4;iter_4++)
        {
            EiEj_nor_dis_avg    [3] -> Fill(ring+1,*iter_4); 
            cout<<"avg4 element:"<<*iter_4<<endl; 
            ring++;
        }
        ring = 0;
        for (iter_5;iter_5 !=iend_5;iter_5++)
        {
            EiEj_nor_dis_avg    [4] -> Fill(ring+1,*iter_5); 
            cout<<"avg5 element:"<<*iter_5<<endl; 
            ring++;
        }
        ring = 0;
        for (iter_6;iter_6 !=iend_6;iter_6++)
        {
            EiEj_nor_dis_avg    [5] -> Fill(ring+1,*iter_6); 
            cout<<"avg6 element:"<<*iter_6<<endl; 
            ring++;
        }
        ring = 0;
        for (iter_7;iter_7 !=iend_7;iter_7++)
        {
            EiEj_nor_dis_avg    [6] -> Fill(ring+1,*iter_7); 
            cout<<"avg7 element:"<<*iter_7<<endl; 
            ring++;
        }
        ring = 0;
        for (iter_8;iter_8 !=iend_8;iter_8++)
        {
            EiEj_nor_dis_avg    [7] -> Fill(ring+1,*iter_8); 
            cout<<"avg8 element:"<<*iter_8<<endl; 
            ring++;
        }
        
        ring = 0;
        for (kter_1;kter_1 !=kend_1;kter_1++)
        {
            EiEj_nor_dis_avg    [0] -> SetBinError(ring+1,*kter_1);
            cout<<"err1 element:"<<*kter_1<<endl;
            ring++;
        }
        ring = 0;
        for (kter_2;kter_2 !=kend_2;kter_2++)
        {
            EiEj_nor_dis_avg    [1] -> SetBinError(ring+1,*kter_2);
            cout<<"err2 element:"<<*kter_2<<endl;
            ring++;
        }
        ring = 0;
        for (kter_3;kter_3 !=kend_3;kter_3++)
        {
            EiEj_nor_dis_avg    [2] -> SetBinError(ring+1,*kter_3);
            cout<<"err3 element:"<<*kter_3<<endl;
            ring++;
        }
        ring = 0;
        for (kter_4;kter_4 !=kend_4;kter_4++)
        {
            EiEj_nor_dis_avg    [3] -> SetBinError(ring+1,*kter_4);
            cout<<"err4 element:"<<*kter_4<<endl;
            ring++;
        }
        ring = 0;
        for (kter_5;kter_5 !=kend_5;kter_5++)
        {
            EiEj_nor_dis_avg    [4] -> SetBinError(ring+1,*kter_5);
            cout<<"err5 element:"<<*kter_5<<endl;
            ring++;
        }
        ring = 0;
        for (kter_6;kter_6 !=kend_6;kter_6++)
        {
            EiEj_nor_dis_avg    [5] -> SetBinError(ring+1,*kter_6);
            cout<<"err6 element:"<<*kter_6<<endl;
            ring++;
        }
        ring = 0;
        for (kter_7;kter_7 !=kend_7;kter_7++)
        {
            EiEj_nor_dis_avg    [6] -> SetBinError(ring+1,*kter_7);
            cout<<"err7 element:"<<*kter_7<<endl;
            ring++;
        }
        ring = 0;
        for (kter_8;kter_8 !=kend_8;kter_8++)
        {
            EiEj_nor_dis_avg    [7] -> SetBinError(ring+1,*kter_8);
            cout<<"err8 element:"<<*kter_8<<endl;
            ring++;
        }

    //Fill Total energy of each layer
    for (int iL = 0;iL < NLAYERS;iL++)
    {
        cout << "The Normalized by event total energy of layer "<<iL+1<<" is "<<Elayertotal[iL]*GEVTOMEV/EventsPassed<<"MeV"<<endl;
        hene_total_layer -> SetBinContent(iL+1,Elayertotal[iL]*GEVTOMEV/EventsPassed);
    }
    
    // Fill TGraph2D 
    double profile_entry[NLAYERS][bin_size];
    double profile_content_X[NLAYERS][bin_size];
    double profile_content_Y[NLAYERS][bin_size];
    char histoName_3[100];
    sprintf(histoName_3,"Interlayer_Distribution_of_Two_point_correlation_with_TGraph2D_MC");
    char histoName_4[100];
    sprintf(histoName_4,"Interlayer_Distribution_of_Two_point_correlation_with_TGraph2D_Data");
    h2_prof_mc = new TGraph2D();
    h2_prof_data = new TGraph2D();
    for(int iL = 0;iL < NLAYERS;iL++)
    {
        double total_entry = 0;
        for(int bin = 0;bin < bin_size;bin++)
        {

            profile_entry[iL][bin] = p1[iL] -> GetBinEntries(bin);
            //cout <<"The entry number in bin "<<bin<<"is:"<<profile_entry[iL][bin]<<endl;
            //total_entry = total_entry + profile_entry[iL][bin];
            if(profile_entry[iL][bin] != 0)
            {
                profile_content_X[iL][bin] = p1[iL] -> GetBinCenter(bin);
                profile_content_Y[iL][bin] = p1[iL] -> GetBinContent(bin);
                
                if(doTruth)
                {
                    h2_prof_mc -> SetPoint(bin_size*iL+bin,profile_content_X[iL][bin],profile_content_Y[iL][bin],iL);
                }
                else
                {
                    h2_prof_data -> SetPoint(bin_size*iL+bin,profile_content_X[iL][bin],profile_content_Y[iL][bin],iL);
                }
                
            }
        }
            //cout <<"The Total_entry_number in layer "<<iL<<" is "<<total_entry<<endl;
    }
    //* Normalize profile histograms per event
    //dRVlayer->Scale(1./EventsPassed);
    //dRprofAll->Sumw2();
    //dRprofAll->Scale(1./EventsPassed);
    
    
    //========== Make Plots ==========

    double chi_sqr = 0.0;
    int logy = 0; // 1: plots with log scale 0:normal scale 

    //Fill the Histogram of the total energy of each layer       
    hene_total_layer -> SetXTitle("Number_of_Layer");
    hene_total_layer -> SetYTitle("MeV");
    hene_total_layer -> Draw();

    //Fill the Histogram of the total energy of the detector(8 layer)
    helayerRawtotal -> Scale(1./EventsPassed);
    helayerRawtotal -> SetXTitle("MeV");
    helayerRawtotal -> SetYTitle("");
    helayerRawtotal -> Draw ();

    
    for(int iL = 0;iL < NLAYERS; iL++)
    {  
        //Fill the Histogram of the energy of each layer per event
        helayerRaw [iL] -> Scale(1./EventsPassed);
        helayerRaw [iL] -> SetXTitle("MeV");
        helayerRaw [iL] -> SetYTitle("");
        helayerRaw [iL] -> Draw();


        // Fill the Histogram of the_normalized_Multiplication of the energies of two hits 
        hEiEj_Multi_normal   [iL] -> SetXTitle("1000 x (Ei*Ej)/(sum Ei)^2");
        hEiEj_Multi_normal   [iL] -> SetYTitle("Number_of_Hits");
        hEiEj_Multi_normal   [iL] -> Draw();


        // Fill the Histogram of the Distance of the two hits
        hdist_Two_hits       [iL] -> SetXTitle("cm");
        hdist_Two_hits       [iL] -> SetYTitle("Number_of_Hits");
        hdist_Two_hits       [iL] -> Draw();
            
        //1000 x (Ei*Ej)/(sum Ei)^21000 x (Ei*Ej)/(sum Ei)^2
        TH1D *hbin_number_x  [iL]; 
        hbin_number_x        [iL]   = new TH1D("n1","",bin_number_2,0,15);
            
        //c1 -> cd ();    
        EiEj_nor_rel_dis     [iL] -> GetXaxis() -> SetTitle("The relative distance between two hits (cm)");
        EiEj_nor_rel_dis     [iL] -> GetYaxis() -> SetTitle("1000*(EiEj/(sum Ei)^2");
        //EiEj_nor_rel_dis     [iL] -> Draw("AP*");
        EiEj_nor_rel_dis     [iL] -> Draw("");
        //c1 -> Update();
        EiEj_nor_dis_pad     [iL] -> GetXaxis() -> SetTitle("The relative distance between two hits (in ring)");
        EiEj_nor_dis_pad     [iL] -> GetYaxis() -> SetTitle("1000*(EiEj/(sum Ei)^2");
        //EiEj_nor_dis_pad     [iL] -> Draw("AP*");
        EiEj_nor_dis_pad   [iL] -> Draw("");
        
        EiEj_nor_dis_avg      [iL] -> GetXaxis() -> SetTitle("The relative distance between two hits (in ring)");
        EiEj_nor_dis_avg      [iL] -> GetYaxis() -> SetTitle("1000*(EiEj/(sum Ei)^2");
        EiEj_nor_dis_avg      [iL] -> SetMarkerStyle(20);
        EiEj_nor_dis_avg      [iL] -> SetMarkerColor();
        //EiEj_nor_dis_avg      [iL] -> SetMarkerColor(kRed);
        EiEj_nor_dis_avg      [iL] -> Draw("EP");
        //EiEj_nor_dis_avg    [iL] -> Draw("");
       
        char histoName_5 [100];
        sprintf(histoName_5,"graph%d",iL);
        EiEj_nor_dis_graph [iL] -> Write(histoName_5);
        EiEj_nor_dis_graph [iL] -> SetTitle("");
        EiEj_nor_dis_graph [iL] -> SetMarkerColor(4);
        EiEj_nor_dis_graph [iL] -> SetMarkerStyle();
        EiEj_nor_dis_graph [iL] -> Draw("AP");


        /* 
        double* xaxis = EiEj_nor_rel_dis  [iL]->GetX();
        double* yaxis = EiEj_nor_rel_dis  [iL]->GetY();
        cout<<"xaxis num "<<endl;
        for(int i=0;i<EiEj_nor_rel_dis[iL]->GetBinEntries() ;i++)
        {
            if (yaxis[i] >= 0.08 )
                hbin_number_x [iL] -> Fill(xaxis[i]);
            //cout<<xaxis[i]<<endl;
        }
        
        hbin_number_x [iL] -> SetXTitle("CM");
        hbin_number_x [iL] -> SetYTitle("#");
        hbin_number_x [iL] -> Draw();
        */
        for(int r = 0; r < Ring;r++)
        {    
            Pad_ring_hist[iL][ring] -> SetXTitle("1000 x (Ei*Ej)/(sum Ei)^2");
            Pad_ring_hist[iL][ring] -> SetYTitle("Number_of_Hits"); 
            Pad_ring_hist[iL][ring] -> Draw();
        }    
            
        // Draw the Profile of the variables
        p1                   [iL] -> SetXTitle("cm");
        p1                   [iL] -> SetYTitle("1000*(EiEj/(sum Ei)^2");
        //p1                   [iL] -> SetAxisRange(0,0.1,"Y");
        
       EiEj_prof             [iL] -> SetXTitle("in ring");
       EiEj_prof             [iL] -> SetYTitle("1000*(EiEj/(sum Ei)^2");
       //EiEj_prof             [iL] -> SetAxisRange(0,0.1,"Y");
       EiEj_prof             [iL] -> Draw();
       

       // p2                   [iL] -> SetXTitle("in ring");
       // p2                   [iL] -> SetYTitle("1000*(EiEj/(sum Ei)^2");
       // p2                   [iL] -> SetAxisRange(0,0.1,"Y");
            
        if(doTruth)
        {
            p1                   [iL] -> Draw("HIST");
            //p2                   [iL] -> Draw("HIST");
        }
        else 
        {
            p1                   [iL] -> Draw();
            //p2                   [iL] -> Draw();
        }
    
                 
    }        
    
    //========== TGraph2D part ==========
    if(doTruth)
    {
        h2_prof_mc -> SetTitle(histoName_3);
            
            
        gStyle->SetPalette(1);
        h2_prof_mc->SetMarkerStyle(20);
        h2_prof_mc ->Draw("pcol");
    }
    else
    {
        h2_prof_data -> SetTitle(histoName_4);
        gStyle->SetPalette(1);
        h2_prof_mc->SetMarkerStyle(20);
        h2_prof_mc ->Draw("pcol");
    }
    

    cout << "PassedEvents :" << EventsPassed << endl;
    cout << "Entries      :" << nentries << endl;
    cout << "Efficiency is :" << ((float) EventsPassed)/nentries << endl;

    //std::cout << "a1\n";
    h2_prof_mc -> SetName(histoName_3);
    //std::cout << "a2\n";
    h2_prof_data -> SetName(histoName_4);
    //std::cout << "a3\n";
    h2_prof_mc -> Write();
    //std::cout << "a4\n";
    h2_prof_data -> Write();
    //std::cout << "a5\n";
    outfile->Write("",TObject::kOverwrite);
    //std::cout << "a6\n";
    outfile->Close();
    //std::cout << "a7\n";
   
    for (int i = 0;i < 7;i++  )
    {
        delete   prof_pad_1[i];
        delete   prof_pad_2[i];
        delete   prof_pad_3[i];
        delete   prof_pad_4[i];
        delete   prof_pad_5[i];
        delete   prof_pad_6[i];
        delete   prof_pad_7[i];
        delete   prof_pad_8[i];
     
    }
    //for ( int i=0; i<7 ; ++i )
    //{
    //    delete[] test1[i];
    //    delete[] test2[i];
    //}
}

// Remember to Delete the vector memory


#endif // #ifdef makePlots_cxx
 

