#define myClass_cxx
#include "myClass.h"
#include "math.h"
#include <fstream>
#include <iostream>
#include <string>
#include "utilities.cc"
#include "TRFIOFile.h"
//#include "TCastorFile.h"

#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector< vector<double> >+;
#endif

/*#include <map>
#ifdef __CINT__
#pragma link C++ class map<string,int>+;
#pragma link C++ class map<string,int>::iterator;
#pragma link C++ class pair<string,int>;
//#pragma link C++ class map<string,int,less<string>,allocator<pair<const string,int> > >
#endif
*/

//template class std::vector< std::vector <double> >;

using namespace std;

void myClass::setweight()
{
  using namespace std;
  for( int ii = 0; ii < nFiles; ii++)
    weight.push_back(1.0);
  
  double wantat_lumi = 50;  //pb-1
  double ne          = 10000;
  
  //all are in pb
  double cs_zgam               = 7.3/3.0 ;
  double fe_zgam               = 1;
  int    ne_zgam               = 33536;  // Since zgamma sample is summed over all three leptons, so here only the number of events decaying to ee should come which I got by doing : myEvent->Draw("gen_ZdaughterE[0]","gen_ZdaughterE[0]!=-99.0")                                                                                         
  double genlumi_zgam          = ne_zgam/cs_zgam;

  double cs_zjet_15to20        = 1.45e+2 ;
  double fe_zjet_15to20        = 1;
  int    ne_zjet_15to20        = 209740;
  double genlumi_zjet_15to20   = ne_zjet_15to20/cs_zjet_15to20;

  double cs_zjet_20to30        = 1.318e+2 ;
  double fe_zjet_20to30        = 1;
  int    ne_zjet_20to30        = 175590;
  double genlumi_zjet_20to30   = ne_zjet_20to30/cs_zjet_20to30;

  double cs_zjet_30to50        = 8.438e+1;
  double fe_zjet_30to50        = 1;
  int    ne_zjet_30to50        = 189024;
  double genlumi_zjet_30to50   = ne_zjet_30to50/cs_zjet_30to50;

  double cs_zjet_50to80        = 3.235e+1;
  double fe_zjet_50to80        = 1;
  int    ne_zjet_50to80        = 118383;
  double genlumi_zjet_50to80   = ne_zjet_50to80/cs_zjet_50to80;

  double cs_zjet_80to120       = 9.981;
  double fe_zjet_80to120       = 1;
  int    ne_zjet_80to120       = 150600;
  double genlumi_zjet_80to120  = ne_zjet_80to120/cs_zjet_80to120;

  double cs_zjet_120to170      = 2.76;
  double fe_zjet_120to170      = 1;
  int    ne_zjet_120to170      = 127820;
  double genlumi_zjet_120to170 = ne_zjet_120to170/cs_zjet_120to170;

  double cs_zjet_170to230      = 7.241e-1;
  double fe_zjet_170to230      = 1;
  int    ne_zjet_170to230      = 153000;
  double genlumi_zjet_170to230 = ne_zjet_170to230/cs_zjet_170to230;

  double cs_zjet_230to300      = 1.946e-1;
  double fe_zjet_230to300      = 1;
  int    ne_zjet_230to300      = 110720;
  double genlumi_zjet_230to300 = ne_zjet_230to300/cs_zjet_230to300;

  double cs_zjet_300toinf      = 7.627e-2;
  double fe_zjet_300toinf      = 1;
  int    ne_zjet_300toinf      = 113675;
  double genlumi_zjet_300toinf = ne_zjet_300toinf/cs_zjet_300toinf;


  double cs_estar_1tev      = 5.267e-3;
  double fe_estar_1tev      = 1;
  int    ne_estar_1tev      = 22000;
  double genlumi_estar_1tev = ne_estar_1tev/cs_estar_1tev;


  weight[0] = wantat_lumi/genlumi_zgam;

  weight[1] = wantat_lumi/genlumi_zjet_15to20;
  weight[2] = wantat_lumi/genlumi_zjet_20to30;
  weight[3] = wantat_lumi/genlumi_zjet_30to50;

  weight[4] = wantat_lumi/genlumi_zjet_50to80;
  weight[5] = wantat_lumi/genlumi_zjet_80to120;

  weight[6] = wantat_lumi/genlumi_zjet_120to170;
  weight[7] = wantat_lumi/genlumi_zjet_170to230;

  weight[8] = wantat_lumi/genlumi_zjet_230to300;
  weight[9] = wantat_lumi/genlumi_zjet_300toinf;


  weight[10] = wantat_lumi/genlumi_estar_1tev; 


   if( weight.size() != nFiles)
   {
     cout<<"Stopping the processing further!"<<endl;
     cout<<"number of files provided is "<<nFiles;
     cout<<"But weight is available only for "<<weight.size()<<" files"<<endl;
   }
     
   return;
}



int myClass::eID( char filename[200], int eventno, int ielec )
{
  
  using namespace std;

  
  int flag_e1_EB = 0;
  int flag_e1_EE = 0;

  int flag_e2_EB = 0;
  int flag_e2_EE = 0;

  int selectE    = 0;


  //e1
  double e1_et               = pat_electronet[ielec];
  double e1_pt               = pat_electronpt[ielec];
  double e1_SCeta            = pat_electronSCeta[ielec];
  double e1_dphiin           = pat_electrondeltaPhiSuperClusterTrackAtVtx[ielec];
  double e1_detain           = pat_electrondeltaEtaSuperClusterTrackAtVtx[ielec];
  double e1_HoverE           = pat_electronhadronicOverEm[ielec];
  double e1_sigmaiEtaiEta    = pat_electronscsigmalEtalEta[ielec];
  double e1_ratio_E1x5toE5x5 = (pat_electronscE1X5[ielec])/(pat_electronscE5X5[ielec]);
  double e1_ratio_E2x5toE5x5 = (pat_electronscE2X5Max[ielec])/(pat_electronscE5X5[ielec]);
  double e1_trackIso         = pat_electrondr03dr03TkSumPt[ielec];
  double e1_isoEM            = pat_electrondr03EcalRecHitSumEt[ielec];
  double e1_hadDepth1        = pat_electrondr03HcalDepth1TowerSumEt[ielec];
  double e1_hadDepth2        = pat_electrondr03HcalDepth2TowerSumEt[ielec];
  Char_t e1_isEcalDriven     = pat_electronisEcalDriven[ielec];

  

  //some checks
  
  //cout<<""<<endl;
  //cout<<"====Inside eID======"<<endl;
  //cout<<"filename = "<<filename<<endl;
  //cout<<"event no = "<<eventno<<endl;
  //cout<<"p.pat_electronet[ielec] = "<<pat_electronet[ielec]<<endl;
  //cout<<"p.pat_electronSCeta[ielec] = "<<pat_electronSCeta[ielec]<<endl;
  //cout<<"========================"<<endl;

  //e2
  double e2_et               = pat_electronet[ielec+1];
  double e2_pt               = pat_electronpt[ielec+1];
  double e2_SCeta            = pat_electronSCeta[ielec+1];
  double e2_dphiin           = pat_electrondeltaPhiSuperClusterTrackAtVtx[ielec+1];
  double e2_detain           = pat_electrondeltaEtaSuperClusterTrackAtVtx[ielec+1];
  double e2_HoverE           = pat_electronhadronicOverEm[ielec+1];
  double e2_sigmaiEtaiEta    = pat_electronscsigmalEtalEta[ielec+1];
  double e2_ratio_E1x5toE5x5 = (pat_electronscE1X5[ielec+1])/(pat_electronscE5X5[ielec+1]);
  double e2_ratio_E2x5toE5x5 = (pat_electronscE2X5Max[ielec+1])/(pat_electronscE5X5[ielec+1]);
  double e2_trackIso         = pat_electrondr03dr03TkSumPt[ielec+1];
  double e2_isoEM            = pat_electrondr03EcalRecHitSumEt[ielec+1];
  double e2_hadDepth1        = pat_electrondr03HcalDepth1TowerSumEt[ielec+1];
  double e2_hadDepth2        = pat_electrondr03HcalDepth2TowerSumEt[ielec+1];
  Char_t e2_isEcalDriven     = pat_electronisEcalDriven[ielec+1];


 

  //isolation values for EB
  double etmin_EB              = 25.0;
  char isEcalDriven_EB[10]     = "true";
  double HoverE_min_EB         = 0.05;
  double ratio_E2x5toE5x5min_EB= 0.94; 
  double ratio_E1x5toE5x5min_EB= 0.83;  
  double detain_min_EB         = 0.005;
  double dphiin_min_EB         = 0.09;
  double trackIsomin_EB        = 7.5;

  //pt dependent cut...will be different for e1 and e2
  double e1_isoEm_HadDepth1min_EB    = 2.0 + 0.03*e1_pt;
  double e2_isoEm_HadDepth1min_EB    = 2.0 + 0.03*e2_pt;

  
  //isolation values for EE
  double etmin_EE              = 25.0;
  char isEcalDriven_EE[10]     = "true";
  double HoverE_min_EE         = 0.05;
  double sigma_iEta_iEtamin    = 0.03; 
  double detain_min_EE         = 0.007;
  double dphiin_min_EE         = 0.09;
  double isoHadDepth2_min_EE   = 0.5;
  double trackIsomin_EE        = 15.0;
  
  //pt dependent cut...will be different for e1 and e2
  double e1_isoEm_HadDepth1min_EE;
  double e2_isoEm_HadDepth1min_EE;
  if( e1_pt<50.0 ) e1_isoEm_HadDepth1min_EE = 2.5;
  if( e1_pt>50.0 ) e1_isoEm_HadDepth1min_EE = 2.5 + 0.03*(e1_pt - 50.0);

  if( e2_pt<50.0 ) e2_isoEm_HadDepth1min_EE = 2.5;
  if( e2_pt>50.0 ) e2_isoEm_HadDepth1min_EE = 2.5 + 0.03*(e2_pt - 50.0);

  
  



  
  //==============e1 =====================
  
  //EB
  if( fabs(e1_SCeta)<1.442 )
    {
      if( e1_et>etmin_EB && fabs(e1_detain)<detain_min_EB && fabs(e1_dphiin)<dphiin_min_EB && e1_HoverE<HoverE_min_EB && (e1_ratio_E2x5toE5x5>ratio_E2x5toE5x5min_EB || e1_ratio_E1x5toE5x5>ratio_E1x5toE5x5min_EB) && e1_trackIso<trackIsomin_EB && (e1_isEcalDriven) && ((e1_isoEM+e1_hadDepth1)<e1_isoEm_HadDepth1min_EB) )
	flag_e1_EB = 1; 
    }

   //EE
  if( fabs(e1_SCeta)<2.5 && fabs(e1_SCeta)>1.560  )
    {
      if( e1_et>etmin_EE && fabs(e1_detain)<detain_min_EE && fabs(e1_dphiin)<dphiin_min_EE && e1_HoverE<HoverE_min_EE && e1_sigmaiEtaiEta<sigma_iEta_iEtamin  && e1_trackIso<trackIsomin_EE && (e1_isEcalDriven) && ((e1_isoEM+e1_hadDepth1)<e1_isoEm_HadDepth1min_EE) && e1_hadDepth2<isoHadDepth2_min_EE )
	flag_e1_EE = 1; 
    }
  

  //==============e2 =====================
  //EB                                                                                                                                                                  
  if( fabs(e2_SCeta)<1.442 )
    {
      if( e2_et>etmin_EB && fabs(e2_detain)<detain_min_EB && fabs(e2_dphiin)<dphiin_min_EB && e2_HoverE<HoverE_min_EB && (e2_ratio_E2x5toE5x5>ratio_E2x5toE5x5min_EB ||e2_ratio_E2x5toE5x5>ratio_E2x5toE5x5min_EB) && e2_trackIso<trackIsomin_EB && (e2_isEcalDriven) && ((e2_isoEM+e2_hadDepth1)<e2_isoEm_HadDepth1min_EB) )
	flag_e2_EB = 1;
    }

  //EE                                                                                                                                                                 
  if( fabs(e2_SCeta)<2.5 && fabs(e2_SCeta)>1.560  )
    {
      if( e2_et>etmin_EE && fabs(e2_detain)<detain_min_EE && fabs(e2_dphiin)<dphiin_min_EE && e2_HoverE<HoverE_min_EE && e2_sigmaiEtaiEta<sigma_iEta_iEtamin  && e2_trackIso<trackIsomin_EE && (e2_isEcalDriven) && ((e2_isoEM+e2_hadDepth1)<e2_isoEm_HadDepth1min_EE) && e2_hadDepth2<isoHadDepth2_min_EE )
	flag_e2_EE = 1;
    }


  
  if( (flag_e1_EB==1 || flag_e1_EE==1) && (flag_e2_EB==1 || flag_e2_EE==1) ) 
    selectE = 1;

  return selectE;
}


int myClass::pID( char filename[200], int eventno, int ipho )
{

  using namespace std;

  //postAnalyzer p;

  //TFile *file = TFile::Open(filename);
  //TTree *tree = (TTree*)file->Get("myEvent");
  //p.Init(tree);
  //tree->GetEntry(eventno,1);

  //cout<<""<<endl;
  //cout<<"==========Inside pID========="<<endl;
  //cout<<"pat_photonisoecalRecHit[ipho] = "<<pat_photonisoecalRecHit[ipho]<<endl;
  //cout<<"pat_photonisohcalRecHit[ipho] = "<<pat_photonisohcalRecHit[ipho]<<endl;
  //cout<<"=============================="<<endl;

  int selectP = 0;

  double isolation_Ecal = pat_photonecalRecHitSumEtConeDR04[ipho];
  double isolation_Hcal = pat_photonhcalTowerSumEtConeDR04[ipho];
  double HoverE         = pat_phohadronicOverEm[ipho];
  double hollowtrackIso = pat_photontrkSumPtHollowConeDR04[ipho];
  double etaWidth       = pat_photonsigmaIetaIeta[ipho];
  double pt             = pat_photonpt[ipho];
  double eta            = pat_photoneta[ipho];

  //isTight criteria
  //EB
  double hollowtrackIsoMin_EB = 9.0;
  double isolation_EcalMin_EB = 5.0 + 0.004*pt;
  double isolation_HcalMin_EB = 5.0;
  double HoverEMin_EB         = 0.15;
  double etaWidthMin_EB       = 0.013;

  //EE
  double hollowtrackIsoMin_EE = 9.0;
  double isolation_EcalMin_EE = 5.0 + 0.0021*pt;
  double isolation_HcalMin_EE = 5.0;
  double HoverEMin_EE         = 0.15;
  //No sigmaIeta cut in the endcap

  //In case of photon detector eta is same as the photon momentum eta

  //EB
  if( fabs(eta)<=1.442 )
    {
      if( isolation_Ecal<isolation_EcalMin_EB && isolation_Hcal<isolation_HcalMin_EB && HoverE<HoverEMin_EB && hollowtrackIso<hollowtrackIsoMin_EB && etaWidth<etaWidthMin_EB )
	selectP = 1;
    }

  //EE
  if( fabs(eta)<=2.5 && fabs(eta)>1.560 )
    {
      if( isolation_Ecal<isolation_EcalMin_EE && isolation_Hcal<isolation_HcalMin_EE && HoverE<HoverEMin_EE && hollowtrackIso<hollowtrackIsoMin_EE )
	selectP = 1;
    }
  return selectP;
}




myClass::myClass(int filenum)
{
  //using namespace std;
  
  nFiles = filenum;
  //myTree = new TTree("myTree","a tree with leaves");
  
  //TFile *rootfile[nFiles];
  //now the weights will be set
  setweight();
  
  //start opening the files
  ifstream infile;
  ifstream outfile;
  infile.open("inputFilenames.list", ifstream::in );
  outfile.open("outputFilenames.list",ifstream::in);

  ntimes = 0;

  fileindex = 0;
  initializeTreeVar();
  while(!infile.eof()){
    infile >> infilename;
    outfile >> outfilename;
    if(strncmp(infilename,"#",1)==0)
      {
	fileindex++;
	continue; 
      }
   
    cout<<""<<endl;
    cout<<"-----some useful info------"<<endl;
    cout<<"started with the "<<fileindex<<"th file(this is fileindex)"<<endl;
    cout<<"open file "<<infilename<<endl;
    cout<<"final output filename = "<<outfilename<<endl;
    cout<<"weight of this file = "<<weight[fileindex]<<endl;
    //rootfile[fileindex] = new TFile(outfilename,"RECREATE");

    //actual function is called here
    //cout<<"Start opening above file in root"<<endl;
    //myfile = TFile::Open(infilename);
    //myfile = TFile::Open(infilename);
    myfile = new TRFIOFile(infilename,"READ");
    //cout<<"opened the above file in root"<<endl;
    bool iszombie = myfile->IsZombie();
    //cout<<"zombie status of file = "<<iszombie<<endl;
    if(myfile->IsZombie()) cout<<"could not access file"<<endl;
    tree = (TTree*)myfile->Get("myEvent");
    if(tree == 0) cout<<"empty tree"<<endl;
    
    //cout<<"got hold of tree"<<endl;
    Init(tree);
    //cout<<"Initialized tree"<<endl;
    //setbranch();
    //cout<<"Calling loop"<<endl;
    Loop();
    //cout<<"back from loop"<<endl;
    fileindex++;
    //cout<<"fileindex now is "<<fileindex<<endl;
    //cout<<"going to directory of root file"<<endl;
    //rootfile[fileindex]->cd();
    //myTree->Write();
    //rootfile[fileindex]->Close();
    //delete myTree;
    
  }//end of while(!infile.eof())

}//end of constructor

myClass::~myClass()
{}

void myClass::Loop()
{
  using namespace std;
  //   In a ROOT session, you can do:
  //      Root > .L postAnalyzer.C
  //      Root > postAnalyzer t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   cout<<"No of entries inside it = "<<nentries<<endl;

   Long64_t nbytes = 0, nb = 0;

   int flagfill_ieta_iphi_elecm = 1; 
   int flagfill_ieta_iphi_phot  = 1;
   //cout<<"opening a root file"<<endl;
   
   //TFile *rootfile = new TFile(outfilename,"RECREATE");
   TRFIOFile *rootfile = new TRFIOFile(outfilename,"RECREATE");

   //cout<<"opening a Tree"<<endl;
   TTree *myTree = new TTree("myTree","a tree with leaves"); // tree where everything will be written
   cout<<"After opening now setting the weights"<<endl;
   myTree->SetWeight(weight[fileindex]);
   myTree->AutoSave();
   cout<<"Autosave"<<endl;
   cout<<"now setting the branch"<<endl;
   setbranch(myTree);
   cout<<"set the branch"<<endl;

   cout<<"booking the hist now"<<endl;
   bookhisto();
   
   //initialize every tree variable for each event 
   //initializeTreeVar();

   int ievent = 0;

   //nentries = 1;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      //set some flags first for every event
      int flag_e = 0;
      int flag_p = 0; 
      int flag_j = 0; 

      //if(is_HLT_DoubleEle10_SW_L1R_event == 1)
      //cout<<" is_HLT_DoubleEle10_SW_L1R_event = "<<is_HLT_DoubleEle10_SW_L1R_event<<endl;

      //if(is_HLT_DoubleEle10_SW_L1R_event == 0)
      //cout<<" is_HLT_DoubleEle10_SW_L1R_event is zero"<<endl;

      //cout<<"setting all flags"<<endl;
      if( pat_electronsize>1 ) 
	 {
	   flag_e = 1;
	 }


       if(pat_photonsize>0)
	 {
	   flag_p = 1;
	 }
       

       if(pat_jetsize>0)
 	 {
	   flag_j = 1;
	 } 
      
       if( flag_e == 1 && flag_p == 1 && (gen_ZdaughterE[0]!=-99.0 || gen_ZbosonE[0]==-99.0) ) ////gen_ZdaughterE[0]!=-99.0 ensures that I am not choosing electrons which are coming from tau decay. Because in zgamma sample, all the three decays are present and I am storing gen_Zdaughter info only when electron is found.Whenever electron is not found then the value stored is -99.0. Since in my signal sample I there is no Z and hence gen_ZdaughterE is always = -99.0 so to make my signal events enter this loop, I am putting an extra condition of gen_Zboson==-99.0 which is always true for my signal sample. Signal is e*e production   
	 {
	   cout<<"Now inside flags conditions"<<endl;

	   triggermap->clear();
	   
	   //fill trigger info here
	   triggermap->insert(pair<string,int>("is_HLT_Ele10_SW_L1R_event",is_HLT_Ele10_SW_L1R_event));  
	   triggermap->insert(pair<string,int>("is_HLT_Ele15_SW_L1R_event",is_HLT_Ele15_SW_L1R_event));  
	   triggermap->insert(pair<string,int>("is_HLT_Ele15_SW_EleId_L1R_event",is_HLT_Ele15_SW_EleId_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele15_SW_LooseTrackIso_L1R_event",is_HLT_Ele15_SW_LooseTrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event",is_HLT_Ele15_SC15_SW_LooseTrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele15_SC15_SW_EleId_L1R_event",is_HLT_Ele15_SC15_SW_EleId_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele20_SW_L1R_event",is_HLT_Ele20_SW_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele20_SC15_SW_L1R_event",is_HLT_Ele20_SC15_SW_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele25_SW_L1R_event",is_HLT_Ele25_SW_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event",is_HLT_Ele25_SW_EleId_LooseTrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoubleEle5_SW_Jpsi_L1R_event",is_HLT_DoubleEle5_SW_Jpsi_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoubleEle5_SW_Upsilon_L1R_event",is_HLT_DoubleEle5_SW_Upsilon_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoubleEle10_SW_L1R_event",is_HLT_DoubleEle10_SW_L1R_event));
	   
	   triggermap->insert(pair<string,int>("is_HLT_Photon10_L1R_event",is_HLT_Photon10_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event",is_HLT_Photon10_LooseEcalIso_TrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon15_L1R_event",is_HLT_Photon15_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event",is_HLT_Photon20_LooseEcalIso_TrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon25_L1R_event",is_HLT_Photon25_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event",is_HLT_Photon25_LooseEcalIso_TrackIso_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon30_L1R_1E31_event",is_HLT_Photon30_L1R_1E31_event));
	   triggermap->insert(pair<string,int>("is_HLT_Photon70_L1R_event",is_HLT_Photon70_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoublePhoton10_L1R_event",is_HLT_DoublePhoton10_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoublePhoton15_L1R_event",is_HLT_DoublePhoton15_L1R_event));
	   triggermap->insert(pair<string,int>("is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event",is_HLT_DoublePhoton15_VeryLooseEcalIso_L1R_event));
	   
	   
	   //clear the vectors
	   tree_eID->clear();
	   tree_pID->clear();
	   //for every event get eID and pID for every e and for every photon
	   
	   int eid = eID( infilename,jentry,0 );
	   tree_eID->push_back(eid); 
	   
	   int pid = pID( infilename,jentry,0 );
	   tree_pID->push_back(pid);
	   	   
	   tree_genZdaughterE->clear();
	   tree_genZdaughterE->push_back(gen_ZdaughterE[0]);
	   
	   

	   //============================begin photon tree variables=====================================================

	   //first clear all these arrays 
	   tree_patphotonpt->clear();
	   tree_patphotonpx->clear();
	   tree_patphotonpy->clear();
	   tree_patphotonpz->clear();
           tree_patphotoneta->clear();
	   tree_patphotonphi->clear();
           tree_patphotonE->clear();
	   tree_patphotonHoverE->clear();

           tree_patphotonEcalRecHit03->clear();
           tree_patphotonHcalRecHit03->clear();
           tree_patphotonEcalRecHit04->clear();
           tree_patphotonHcalRecHit04->clear();
	   tree_patphotonHcalDepth104->clear();
	   tree_patphotonHcalDepth204->clear();
           tree_patphotonIsoHollowtrkcone03->clear();
           tree_patphotonIsoHollowtrkcone04->clear();
           tree_patphotonSigmaIeta->clear();
	   tree_patphotonisConverted->clear();
	   //tree_patphotonSCeta->clear();  //somehow not filled in cmssw314
	   //tree_patphotonSCphi->clear();  //somehow not filled in cmssw314

	   //now fill
	   tree_patphotonpt->push_back(pat_photonpt[0]);
	   tree_patphotonpx->push_back(pat_photon_px[0]);
	   tree_patphotonpy->push_back(pat_photon_py[0]);
	   tree_patphotonpz->push_back(pat_photonpz[0]);
	   tree_patphotoneta->push_back(pat_photoneta[0]);
	   tree_patphotonphi->push_back(pat_photonphi[0]);
	   tree_patphotonE->push_back(pat_photonE[0]);
	   tree_patphotonHoverE->push_back(pat_phohadronicOverEm[0]);

	   tree_patphotonEcalRecHit03->push_back(pat_photonisoecalRecHit[0]);
	   tree_patphotonHcalRecHit03->push_back(pat_photonisohcalRecHit[0]);
	   tree_patphotonIsoHollowtrkcone03->push_back(pat_photonisohollowtrkcone[0]);
	   tree_patphotonEcalRecHit04->push_back(pat_photonecalRecHitSumEtConeDR04[0]);
	   tree_patphotonHcalRecHit04->push_back(pat_photonhcalTowerSumEtConeDR04[0]);
	   tree_patphotonHcalDepth104->push_back(pat_photonhcalDepth1TowerSumEtConeDR04[0]);
	   tree_patphotonHcalDepth204->push_back(pat_photonhcalDepth2TowerSumEtConeDR04[0]);
	   tree_patphotonIsoHollowtrkcone04->push_back(pat_photontrkSumPtHollowConeDR04[0]);
	   tree_patphotonSigmaIeta->push_back(pat_photonsigmaIetaIeta[0]);
	   tree_patphotonisConverted->push_back(pat_photonisConvertedPhoton[0]);
	  
	   //cout<<""<<endl;
	   //cout<<"pat_photonpt[0]  = "<<pat_photonpt[0]<<endl;
	   //cout<<"tree_patphotonpt = "<<(*tree_patphotonpt)[0]<<endl;
	   //==========================end of photon tree variables===================================================
	   

	   //=============================begin with e- tree variables===================================================
	   //e-
	   //clear the vectors

	   //cout<<"clearing e's variables"<<endl;
	   tree_pateeta->clear();
	   tree_patephi->clear();

           tree_pateE->clear();
           tree_patept->clear();
	   tree_patepx->clear();
           tree_patepy->clear();
           tree_patepz->clear();
           tree_pateEt->clear();
           tree_pateSCeta->clear();
	   //cout<<"cleared till SCeta"<<endl;
           //tree_pateSCphi->clear();
           tree_patedeta_in->clear();
	   //cout<<"cleared till detain"<<endl;
           tree_patedphi_in->clear();
	   //cout<<"cleared till dephiin"<<endl;
           tree_pateHoverE->clear();
	   //cout<<"cleared till hOverE"<<endl;
           tree_patesigIetaIeta->clear();
	   //cout<<"clearing sigIeta  "<<endl;
	   tree_patesigetaeta->clear();
	   //cout<<"cleared till before energy ratio1"<<endl;
	   tree_pateratio_E1x5OverE5x5->clear();
	   //cout<<"cleared till before energy ratio1"<<endl;
	   tree_pateratio_E2x5OverE5x5->clear();
	   //cout<<"cleared till ratios"<<endl;
	   
	   tree_patedr03EcalRecHitSumEt->clear();
	   tree_patedr03HcalDepth1TowerSumEt->clear();
	   tree_patedr03HcalDepth2TowerSumEt->clear();
	   tree_patedr03HcalTowerSumEt->clear();
	   tree_patedr03TkSumPt->clear();
	   //cout<<"cleared till dr03cone"<<endl;

	   tree_patedr04EcalRecHitSumEt->clear();
	   tree_patedr04HcalDepth1TowerSumEt->clear();
	   tree_patedr04HcalDepth2TowerSumEt->clear();
	   tree_patedr04HcalTowerSumEt->clear();
	   tree_patedr04TkSumPt->clear();
	   //cout<<"cleared till dr04"<<endl;

	   tree_patecharge->clear();
	   tree_pateisEcalDriven->clear();
	   
	   //cout<<"filling e's variables"<<endl;

	   //now fill
	   for(int ielec=0; ielec<2; ielec++)
	     {
	       //cout<<"ielec = "<<ielec<<endl;
	       tree_pateeta->push_back(pat_electroneta[ielec]);
	       tree_patephi->push_back(pat_electronphi[ielec]);
	       tree_pateE->push_back(pat_electronE[ielec]);
	       tree_patept->push_back(pat_electronpt[ielec]);
	       tree_patepx->push_back(pat_electron_px[ielec]);
	       tree_patepy->push_back(pat_electron_py[ielec]);
	       tree_patepz->push_back(pat_electronpz[ielec]);
	       tree_pateEt->push_back(pat_electronet[ielec]);
	       tree_pateSCeta->push_back(pat_electronSCeta[ielec]);
	       //tree_pateSCphi->push_back(pat_electronSCphi[ielec]);
	       tree_patedeta_in->push_back(pat_electrondeltaEtaSuperClusterTrackAtVtx[ielec]);
	       tree_patedphi_in->push_back(pat_electrondeltaPhiSuperClusterTrackAtVtx[ielec]);
	       tree_pateHoverE->push_back(pat_electronhadronicOverEm[ielec]);
	       tree_patesigIetaIeta->push_back(pat_electronscsigmalEtalEta[ielec]);
	       tree_patesigetaeta->push_back(pat_electronscsigmaEtaEta[ielec]);
	       
	       double eratio_E1x5OverE5x5 = pat_electronscE1X5[ielec]/pat_electronscE5X5[ielec];
	       double eratio_E2x5OverE5x5 = pat_electronscE2X5Max[ielec]/pat_electronscE5X5[ielec];
	       
	       tree_pateratio_E1x5OverE5x5->push_back(eratio_E1x5OverE5x5);
	       tree_pateratio_E2x5OverE5x5->push_back(eratio_E2x5OverE5x5);
	       
	       //cout<<"filling the cone isolations now"<<endl;

	       tree_patedr03EcalRecHitSumEt->push_back(pat_electrondr03EcalRecHitSumEt[ielec]);
	       tree_patedr03HcalDepth1TowerSumEt->push_back(pat_electrondr03HcalDepth1TowerSumEt[ielec]);
	       tree_patedr03HcalDepth2TowerSumEt->push_back(pat_electrondr03HcalDepth2TowerSumEt[ielec]);
	       tree_patedr03HcalTowerSumEt->push_back(pat_electrondr03HcalTowerSumEt[ielec]);
	       tree_patedr03TkSumPt->push_back(pat_electrondr03dr03TkSumPt[ielec]);
	       
	       tree_patedr04EcalRecHitSumEt->push_back(pat_electrondr04EcalRecHitSumEt[ielec]);
	       tree_patedr04HcalDepth1TowerSumEt->push_back(pat_electrondr04HcalDepth1TowerSumEt[ielec]);
	       tree_patedr04HcalDepth2TowerSumEt->push_back(pat_electrondr04HcalDepth2TowerSumEt[ielec]);
	       tree_patedr04HcalTowerSumEt->push_back(pat_electrondr04HcalTowerSumEt[ielec]);
	       tree_patedr04TkSumPt->push_back(pat_electrondr04dr04TkSumPt[ielec]);
	       
	       tree_patecharge->push_back(pat_electroncharge[ielec]);
	       tree_pateisEcalDriven->push_back(pat_electronisEcalDriven[ielec]);
	       //cout<<"filled up"<<endl;
	     }
	   //cout<<"Filled e's variables"<<endl;

	   //=============================end of e- tree variables===================================================
	   
	   
	   //===========================get the rechit map of sc from e and gamma for one event======================
	   
           //rechit map of e-
           /*if(flagfill_ieta_iphi_elecm==1)
             {
	     char file_num[100];
	       sprintf(file_num, "%d", fileindex);
	       std::string i_file = file_num;
               //cout<<"inside if(flagfill_ieta_iphi_elecm==1)"<<endl;
               if(pat_electronxtalsize[0]==0){ flagfill_ieta_iphi_elecm = 1; }//cout<<"but pat_electronxtalsize[0]==0"<<endl;
               if(pat_electronxtalsize[0]!=0)
                 {
                   flagfill_ieta_iphi_elecm=0;

                   for(int xtalno=0; xtalno<pat_electronxtalsize[0]; xtalno++)
                     {
                       if(pat_electronxtalE[0][xtalno]>0)
                         {
			   
			   hist2D_patelectron_ieta_iphi["emIetaIphi_"+i_file]->Fill(pat_electronieta[0][xtalno],pat_electroniphi[0][xtalno],log(pat_electronxtalE[0][xtalno]));

                           //cout<<"filled hist_patelectron_ieta_iphi[fileindex]"<<endl;
                         }
                     } //for(int xtalno=0; xtalno<pat_elecm_xtalsize; xtalno++)
                 } //if(pat_elecm_xtalsize!=0)
             } //if(flagfill_ieta_iphi_elecm==1)

	   //rechit map of gamma
           if(flagfill_ieta_iphi_phot==1)
             {
	       char file_num[100];
               sprintf(file_num, "%d", fileindex);
	       std::string i_file = file_num;
               if(pat_photonxtalsize[0]==0) flagfill_ieta_iphi_phot = 1;
               if(pat_photonxtalsize[0]!=0)
                 {
                   flagfill_ieta_iphi_phot=0;
                   for(int xtalno=0; xtalno<pat_photonxtalsize[0]; xtalno++)
                     {
                       if(pat_photonxtalE[0][xtalno]>0)
                         hist2D_patphoton_ieta_iphi["photonIetaIphi_"+i_file]->Fill(pat_photonieta[0][xtalno],pat_photoniphi[0][xtalno],log(pat_photonxtalE[0][xtalno]));
                     } //for(int xtalno=0; xtalno<pat_elecm_xtalsize; xtalno++)
                 } //if(pat_elecm_xtalsize!=0)
             } //if(flagfill_ieta_iphi_elecm==1)

	   //get the rechit of e- when E2x5/E5x5>1
	   if( ievent<10 && ((pat_electronscE2X5Max[0]/pat_electronscE5X5[0]) > 1.0) )
             {
               char file_num[100];
               sprintf(file_num, "%d", fileindex);
	       std::string i_file = file_num;
               
               char i_event[100];
	       sprintf(i_event, "%d", ievent);
	       std::string eventnum = i_event;
               if(pat_electronxtalsize[0]!=0)
                 {
		   ievent = ievent + 1 ;
		   
                   for(int xtalno=0; xtalno<pat_electronxtalsize[0]; xtalno++)
                     {
		       if(ievent<3)
			 {
			   std::cout<<""<<endl;
			   std::cout<<" ievent = "<<ievent<<endl;
			   std::cout<<"xtalE   = "<<pat_electronxtalE[0][xtalno]<<endl;
			   std::cout<<"ieta    = "<<pat_electronieta[0][xtalno]<<endl;
			   std::cout<<"iphi    = "<<pat_electroniphi[0][xtalno]<<endl;
			 }

                       if(pat_electronxtalE[0][xtalno]>0)
                         {
			   
                           hist2D_patelectron_ieta_iphi_E2x5Ratio["emIetaIphi_"+i_file+"_"+eventnum]->Fill(pat_electronieta[0][xtalno],pat_electroniphi[0][xtalno],log(pat_electronxtalE[0][xtalno]));

                           //cout<<"filled hist_patelectron_ieta_iphi[fileindex]"<<endl;
                         }
                     } //for(int xtalno=0; xtalno<pat_elecm_xtalsize; xtalno++)
                 } //if(pat_elecm_xtalsize!=0)
             }//if( ievent<10 && ((pat_electronscE2X5Max[0]/pat_electronscE5X5[0]) > 1.0) )



	   //====================construct some quantities==================================
	 
           //===========================for e- ============================================
           double em_energy_prod_ratio[10];
           double em_energy_prod[10];
           double em_energy_prod8_ratio[10];
           double em_energy_prod8[10];
           double em_energy_prod25_ratio[10];
           double em_energy_prod25[10];
           for( int ii=0;ii<10;ii++ )
             {
               em_energy_prod_ratio[ii]=1;
               em_energy_prod[ii]=0;
               em_energy_prod8_ratio[ii]=1;
               em_energy_prod8[ii]=0;
               em_energy_prod25_ratio[ii]=1;
               em_energy_prod25[ii]=0;
	       
             }

	   //clear the vectors
	   tree_em_ratio_E1x5OverE5x5->clear();
	   tree_em_ratio_E2x5OverE5x5->clear();
	   
	   tree_em_energy_prod_ratio->clear();
	   //cout<<"tree_em_energy_prod_ratio->size() = "<<tree_em_energy_prod_ratio->size()<<endl;  //debug
	   tree_em_energy_prod->clear();
	   
	   tree_em_energy_prod8_ratio->clear();
	   tree_em_energy_prod8->clear();

	   tree_em_energy_prod25_ratio->clear();
	   tree_em_energy_prod25->clear();
	   
	   tree_em_index->clear();
	   tree_em_8xtalindex->clear();
	   tree_em_25xtalindex->clear();
	   
	   tree_em0_energy_prod_ratio->clear();
	   tree_em0_energy_prod->clear();
	   tree_em1_energy_prod_ratio->clear();
	   tree_em1_energy_prod->clear();

	   tree_em0_energy_prod8_ratio->clear();
           tree_em0_energy_prod8->clear();
           tree_em1_energy_prod8_ratio->clear();
           tree_em1_energy_prod8->clear();

	   tree_em0_energy_prod25_ratio->clear();
           tree_em0_energy_prod25->clear();
           tree_em1_energy_prod25_ratio->clear();
           tree_em1_energy_prod25->clear();


	   */

	   //cout<<"entry no = "<<jentry<<endl;
	   //if( pat_electronsize !=0 )
	   //{
	       //resize various vectors
	       /*tree_em_energy_prod_ratio->resize(pat_electronsize);
	       tree_em_energy_prod->resize(pat_electronsize);

	       tree_em_energy_prod8_ratio->resize(pat_electronsize);
	       tree_em_energy_prod8->resize(pat_electronsize);

	       tree_em_energy_prod25_ratio->resize(pat_electronsize);
	       tree_em_energy_prod25->resize(pat_electronsize);
	       */
      
	       /*  for( int ielec=0; ielec<pat_electronsize;ielec++ )
		 {
		   // fill the vectors for the scatter plot of E1x5 and E2x5 vs the quantities below
		   tree_em_ratio_E1x5OverE5x5->push_back(pat_electronscE1X5[ielec]/pat_electronscE5X5[ielec]);
                   tree_em_ratio_E2x5OverE5x5->push_back(pat_electronscE2X5Max[ielec]/pat_electronscE5X5[ielec]);
		   
		   int flag_fill   = 0;
		   int flag_fill8  = 0;
		   int flag_fill25 = 0;
		   int index[1000];
		   if( pat_electronxtalsize[ielec] != 0 )
		     {
		       //tree_em_energy_prod_ratio->resize(pat_electronsize);
		       TMath::Sort(pat_electronxtalsize[ielec],pat_electronxtalE[ielec],index);
		       
		       //construct that qty for all crystals
		       for ( int jcrys = 1; jcrys < pat_electronxtalsize[ielec]; jcrys++ )
			 {
			   //if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index[jcrys]]) && (pat_electronxtalE[0][index[0]] == pat_electronxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			   if( (!isnan(pat_electronxtalE[ielec][index[jcrys]]) && !isnan(pat_electronxtalE[ielec][index[0]])) &&  (!isinf(pat_electronxtalE[ielec][index[jcrys]]) && !isinf(pat_electronxtalE[ielec][index[0]]) ) )
			     {
			       if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec][index[jcrys]] > 0.09  && pat_electronxtalE[ielec][index[0]] < 20000 && pat_electronxtalE[ielec][index[0]] > 0.09 )
				 {
				   flag_fill = 1;
				   em_energy_prod_ratio[ielec] *= (pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]])/pat_electronxtalE[ielec][index[0]];
				   em_energy_prod[ielec] += log(pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]]);
				 } //if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec][index...])
			     }//if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index..])
			 }//for (jcrys = 1; jcrys < pat_electronxtalsize[ielec]; jcrys++ )
		       
		       
		       
		       if( flag_fill == 1 )
			 {
			   //fill the vector 
			   //vector<double> em;
			   //em.push_back(em_energy_prod_ratio[ielec]);
			   //cout<<"filling tree_em_energy_prod_ratio"<<endl;
			   tree_em_energy_prod_ratio->push_back(em_energy_prod_ratio[ielec]);
			   tree_em_energy_prod->push_back(em_energy_prod[ielec]);
			   tree_em_index->push_back(ielec);

			   if(ielec==0)
			     {
			       tree_em0_energy_prod_ratio->push_back(em_energy_prod_ratio[ielec]);
			       tree_em0_energy_prod->push_back(em_energy_prod[ielec]);
			     }
			   
			   if(ielec==1)
			     {
			       tree_em1_energy_prod_ratio->push_back(em_energy_prod_ratio[ielec]);
			       tree_em1_energy_prod->push_back(em_energy_prod[ielec]);
			     }

			   //(*tree_em_energy_prod_ratio)[ielec] = (em_energy_prod_ratio[ielec]);
			   //(*tree_em_energy_prod_ratio)[ielec].push_back(em_energy_prod_ratio[ielec]);
			   //try[ielec] = em_energy_prod_ratio[ielec];
			   //cout<<"filled tree_em_energy_prod_ratio"<<endl;
			   //(*tree_em_energy_prod)[ielec] = (em_energy_prod[ielec]);
			   //cout<<"filled tree_em_energy_prod"<<endl;

			   
			   if(ielec == 0) 
			     {
			       ntimes++;
			       //std::cout<<"One more entry inside flag_fill with ntimes = "<<ntimes<<endl; 
			       //cout<<"jentry = "<<jentry<<endl;
			       //cout<<""<<endl;
			     }

			   //i have got E1x5/E5x5 and E2x5/E5x5 for every electron to get some scatter plots 
			   //get the same when eId and pID is true
			 }//if( flag_fill == 1)
		       
		       //construct that qty for 8 highest energy crstals
		       if(pat_electronxtalsize[ielec]>=8)
			 {
			   for ( int jcrys = 1; jcrys < 8; jcrys++ )
			     {
			       //if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index[jcrys]]) && (pat_electronxtalE[0][index[0]] == pat_electronxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_electronxtalE[ielec][index[jcrys]]) && !isnan(pat_electronxtalE[ielec][index[0]])) && (!isinf(pat_electronxtalE[ielec][index[jcrys]]) && !isinf(pat_electronxtalE[ielec][index[0]]) ) )
				 {
				   
				   if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec][index[jcrys]] > 0.09 && pat_electronxtalE[ielec][index[0]] < 20000 && pat_electronxtalE[ielec][index[0]] > 0.09 )
				   {
				     flag_fill8 = 1;
				     em_energy_prod8_ratio[ielec] *= (pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]])/pat_electronxtalE[ielec][index[0]];
				       em_energy_prod8[ielec] += log(pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]]);
				   } //if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec]..)
				 }// if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index...])
			     }//for (jcrys = 1; jcrys < 8; jcrys++)
			   
			   if( flag_fill8 == 1)
			     {
			       //(*tree_em_energy_prod8_ratio)[ielec] = (em_energy_prod8_ratio[ielec]);
			       //cout<<"filled tree_em_energy_prod8_ratio"<<endl;
			       //(*tree_em_energy_prod8)[ielec] = (em_energy_prod8[ielec]);
			       //cout<<"filled tree_em_energy_prod8:em"<<endl;
			       //tree_em_8xtalindex->push_back(ielec);

			       tree_em_energy_prod8_ratio->push_back(em_energy_prod8_ratio[ielec]);
			       tree_em_energy_prod8->push_back(em_energy_prod8[ielec]);
			       tree_em_8xtalindex->push_back(ielec);

			       if(ielec==0)
				 {
				   tree_em0_energy_prod8_ratio->push_back(em_energy_prod8_ratio[ielec]);
				   tree_em0_energy_prod8->push_back(em_energy_prod8[ielec]);
				 }

			       if(ielec==1)
				 {
				   tree_em1_energy_prod8_ratio->push_back(em_energy_prod8_ratio[ielec]);
				   tree_em1_energy_prod8->push_back(em_energy_prod8[ielec]);
				 }

			     }//if( flag_fill8 == 1)
			 }//if(pat_electronxtalsize[ielec]>=8)
		       
		       //construct that qty for 25 highest energy crstals
		       if(pat_electronxtalsize[ielec]>=25)
			 {
			   for ( int jcrys = 1; jcrys < 25; jcrys++ )
			     {
			       //if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index[jcrys]]) && (pat_electronxtalE[0][index[0]] == pat_electronxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_electronxtalE[ielec][index[jcrys]]) && !isnan(pat_electronxtalE[ielec][index[0]])) && (!isinf(pat_electronxtalE[ielec][index[jcrys]]) && !isinf(pat_electronxtalE[ielec][index[0]])) )
				 {
				   if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec][index[jcrys]] > 0.09 && pat_electronxtalE[ielec][index[0]] < 20000 && pat_electronxtalE[ielec][index[0]] > 0.09 )
				     {
				       flag_fill25 =1;
				       em_energy_prod25_ratio[ielec] *= (pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]])/pat_electronxtalE[ielec][index[0]];
				       em_energy_prod25[ielec] += log(pat_electronxtalE[ielec][index[0]] - pat_electronxtalE[ielec][index[jcrys]]);
				     }//if( pat_electronxtalE[ielec][index[jcrys]] < 20000 && pat_electronxtalE[ielec][index])
				 }//if( (pat_electronxtalE[ielec][index[jcrys]] == pat_electronxtalE[ielec][index..])
			     }//for (jcrys = 1; jcrys < 25; jcrys++)

			   if( flag_fill25 == 1)
			     {
			       //(*tree_em_energy_prod25_ratio)[ielec] = (em_energy_prod25_ratio[ielec]);
			       //cout<<"filled tree_em_energy_prod25_ratio"<<endl;
			       //(*tree_em_energy_prod25)[ielec] = (em_energy_prod25[ielec]);
			       //cout<<"filled tree_em_energy_prod25:em"<<endl;
			       //tree_em_25xtalindex->push_back(ielec);
			       
			       tree_em_energy_prod25_ratio->push_back(em_energy_prod25_ratio[ielec]);
                               tree_em_energy_prod25->push_back(em_energy_prod25[ielec]);
                               tree_em_25xtalindex->push_back(ielec);

                               if(ielec==0)
                                 {
                                   tree_em0_energy_prod25_ratio->push_back(em_energy_prod25_ratio[ielec]);
                                   tree_em0_energy_prod25->push_back(em_energy_prod25[ielec]);
                                 }

                               if(ielec==1)
                                 {
                                   tree_em1_energy_prod25_ratio->push_back(em_energy_prod25_ratio[ielec]);
                                   tree_em1_energy_prod25->push_back(em_energy_prod25[ielec]);
                                 }


			     }//if( flag_fill25 == 1)
			 }//if(pat_electronxtalsize[ielec]>=25)
		     }// if( pat_electronxtalsize[ielec] != 0 )
		 } //for(int ielec=0; ielec<pat_electronsize;ielec++)
	     }//if( pat_electronsize!=0 )
           //=======================for e- ends here===========================



           //===========================for e+ ============================================
           double ep_energy_prod_ratio[10];
           double ep_energy_prod[10];
           double ep_energy_prod8_ratio[10];
           double ep_energy_prod8[10];
           double ep_energy_prod25_ratio[10];
           double ep_energy_prod25[10];
           for( int ii=0;ii<10;ii++ )
             {
               ep_energy_prod_ratio[ii]=1;
               ep_energy_prod[ii]=0;
               ep_energy_prod8_ratio[ii]=1;
               ep_energy_prod8[ii]=0;
               ep_energy_prod25_ratio[ii]=1;
               ep_energy_prod25[ii]=0;
	       
             }

	   //clear the vectors
	   tree_ep_ratio_E1x5OverE5x5->clear();
	   tree_ep_ratio_E2x5OverE5x5->clear();
	   
	   tree_ep_energy_prod_ratio->clear();
	   tree_ep_energy_prod->clear();
	   
	   tree_ep_energy_prod8_ratio->clear();
	   tree_ep_energy_prod8->clear();
	   
	   tree_ep_energy_prod25_ratio->clear();
	   tree_ep_energy_prod25->clear();

	   tree_ep_index->clear();
           tree_ep_8xtalindex->clear();
           tree_ep_25xtalindex->clear();

	   tree_ep0_energy_prod_ratio->clear();
           tree_ep0_energy_prod->clear();
           tree_ep1_energy_prod_ratio->clear();
           tree_ep1_energy_prod->clear();

           tree_ep0_energy_prod8_ratio->clear();
           tree_ep0_energy_prod8->clear();
           tree_ep1_energy_prod8_ratio->clear();
           tree_ep1_energy_prod8->clear();

           tree_ep0_energy_prod25_ratio->clear();
           tree_ep0_energy_prod25->clear();
           tree_ep1_energy_prod25_ratio->clear();
           tree_ep1_energy_prod25->clear();

	   if( pat_electronpsize !=0 )
	     {
	       //resize various vectors
               tree_ep_energy_prod_ratio->resize(pat_electronpsize);
               tree_ep_energy_prod->resize(pat_electronpsize);

               tree_ep_energy_prod8_ratio->resize(pat_electronpsize);
               tree_ep_energy_prod8->resize(pat_electronpsize);

               tree_ep_energy_prod25_ratio->resize(pat_electronpsize);
               tree_ep_energy_prod25->resize(pat_electronpsize);

	       for( int ielec=0; ielec<pat_electronpsize;ielec++ )
		 {
		   // for the scatter plot of E1x5 and E2x5 vs the quantities below
		   tree_ep_ratio_E1x5OverE5x5->push_back(pat_electronpscE1X5[ielec]/pat_electronpscE5X5[ielec]);
		   tree_ep_ratio_E2x5OverE5x5->push_back(pat_electronpscE2X5Max[ielec]/pat_electronpscE5X5[ielec]);
		   
		   int flag_fill   = 0;
		   int flag_fill8  = 0;
		   int flag_fill25 = 0;
		   int index[1000];
		   if( pat_electronptalsize[ielec] != 0 )
		     {
		       TMath::Sort(pat_electronptalsize[ielec],pat_electronpxtalE[ielec],index);
		       
		       //construct that qty for all crystals
		       for ( int jcrys = 1; jcrys < pat_electronptalsize[ielec]; jcrys++ )
			 {
			   //if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index[jcrys]]) && (pat_electronpxtalE[0][index[0]] == pat_electronpxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			   if( (!isnan(pat_electronpxtalE[ielec][index[jcrys]]) && !isnan(pat_electronpxtalE[ielec][index[0]])) &&  (!isinf(pat_electronpxtalE[ielec][index[jcrys]]) && !isinf(pat_electronpxtalE[ielec][index[0]]) ) )
			     {
			       if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec][index[jcrys]] > 0.09  && pat_electronpxtalE[ielec][index[0]] < 20000 && pat_electronpxtalE[ielec][index[0]] > 0.09 )
				 {
				   flag_fill = 1;
				   ep_energy_prod_ratio[ielec] *= (pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]])/pat_electronpxtalE[ielec][index[0]];
				   ep_energy_prod[ielec] += log(pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]]);
				 } //if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec][index...])
			     }//if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index..])
			 }//for (jcrys = 1; jcrys < pat_electronptalsize[ielec]; jcrys++ )
		       
		       if( flag_fill == 1 )
			 {
			   //(*tree_ep_energy_prod_ratio)[ielec] = (ep_energy_prod_ratio[ielec]);
			   //(*tree_ep_energy_prod)[ielec] = (ep_energy_prod[ielec]);
			   //cout<<"filled tree_em_energy_prod_ratio & prod only:ep"<<endl;
			   //tree_ep_index->push_back(ielec);
			   
			   tree_ep_energy_prod_ratio->push_back(ep_energy_prod_ratio[ielec]);
                           tree_ep_energy_prod->push_back(ep_energy_prod[ielec]);
                           tree_ep_index->push_back(ielec);

                           if(ielec==0)
                             {
                               tree_ep0_energy_prod_ratio->push_back(ep_energy_prod_ratio[ielec]);
                               tree_ep0_energy_prod->push_back(ep_energy_prod[ielec]);
                             }

                           if(ielec==1)
                             {
                               tree_ep1_energy_prod_ratio->push_back(ep_energy_prod_ratio[ielec]);
                               tree_ep1_energy_prod->push_back(ep_energy_prod[ielec]);
                             }


			   //i have got E1x5/E5x5 and E2x5/E5x5 for every electron to get some scatter plots 
			   //get the same when eId and pID is true
			 }//if( flag_fill == 1)
		       
		       //construct that qty for 8 highest energy crstals
		       if(pat_electronptalsize[ielec]>=8)
			 {
			   for ( int jcrys = 1; jcrys < 8; jcrys++ )
			     {
			       //if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index[jcrys]]) && (pat_electronpxtalE[0][index[0]] == pat_electronpxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_electronpxtalE[ielec][index[jcrys]]) && !isnan(pat_electronpxtalE[ielec][index[0]])) && (!isinf(pat_electronpxtalE[ielec][index[jcrys]]) && !isinf(pat_electronpxtalE[ielec][index[0]]) ) )
				 {
				   
				   if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec][index[jcrys]] > 0.09 && pat_electronpxtalE[ielec][index[0]] < 20000 && pat_electronpxtalE[ielec][index[0]] > 0.09 )
				   {
				     flag_fill8 = 1;
				     ep_energy_prod8_ratio[ielec] *= (pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]])/pat_electronpxtalE[ielec][index[0]];
				       ep_energy_prod8[ielec] += log(pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]]);
				   } //if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec]..)
				 }// if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index...])
			     }//for (jcrys = 1; jcrys < 8; jcrys++)
			   
			   if( flag_fill8 == 1)
			     {
			       //(*tree_ep_energy_prod8_ratio)[ielec] = (ep_energy_prod8_ratio[ielec]);
			       //(*tree_ep_energy_prod8)[ielec] = (ep_energy_prod8[ielec]);
			       //cout<<"filled tree_em_energy_prod8_ratio & prod only:ep"<<endl;
			       //tree_ep_8xtalindex->push_back(ielec);
			       
			       tree_ep_energy_prod8_ratio->push_back(ep_energy_prod8_ratio[ielec]);
                               tree_ep_energy_prod8->push_back(ep_energy_prod8[ielec]);
                               tree_ep_8xtalindex->push_back(ielec);

                               if(ielec==0)
                                 {
                                   tree_ep0_energy_prod8_ratio->push_back(ep_energy_prod8_ratio[ielec]);
                                   tree_ep0_energy_prod8->push_back(ep_energy_prod8[ielec]);
                                 }

                               if(ielec==1)
                                 {
                                   tree_ep1_energy_prod8_ratio->push_back(ep_energy_prod8_ratio[ielec]);
                                   tree_ep1_energy_prod8->push_back(ep_energy_prod8[ielec]);
                                 }
			       
			     }//if( flag_fill8 == 1)
			 }//if(pat_electronptalsize[ielec]>=8)
		       
		       //construct that qty for 25 highest energy crstals
		       if(pat_electronptalsize[ielec]>=25)
			 {
			   for ( int jcrys = 1; jcrys < 25; jcrys++ )
			     {
			       //if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index[jcrys]]) && (pat_electronpxtalE[0][index[0]] == pat_electronpxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_electronpxtalE[ielec][index[jcrys]]) && !isnan(pat_electronpxtalE[ielec][index[0]])) && (!isinf(pat_electronpxtalE[ielec][index[jcrys]]) && !isinf(pat_electronpxtalE[ielec][index[0]])) )
				 {
				   if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec][index[jcrys]] > 0.09 && pat_electronpxtalE[ielec][index[0]] < 20000 && pat_electronpxtalE[ielec][index[0]] > 0.09 )
				     {
				       flag_fill25 =1;
				       ep_energy_prod25_ratio[ielec] *= (pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]])/pat_electronpxtalE[ielec][index[0]];
				       ep_energy_prod25[ielec] += log(pat_electronpxtalE[ielec][index[0]] - pat_electronpxtalE[ielec][index[jcrys]]);
				     }//if( pat_electronpxtalE[ielec][index[jcrys]] < 20000 && pat_electronpxtalE[ielec][index])
				 }//if( (pat_electronpxtalE[ielec][index[jcrys]] == pat_electronpxtalE[ielec][index..])
			     }//for (jcrys = 1; jcrys < 25; jcrys++)

			   if( flag_fill25 == 1)
			     {
			       //(*tree_ep_energy_prod25_ratio)[ielec] = (ep_energy_prod25_ratio[ielec]);
			       //(*tree_ep_energy_prod25)[ielec] = (ep_energy_prod25[ielec]);
			       //cout<<"filled tree_em_energy_prod25_ratio & prod only:ep"<<endl;
			       //tree_ep_25xtalindex->push_back(ielec);

			       tree_ep_energy_prod25_ratio->push_back(ep_energy_prod25_ratio[ielec]);
                               tree_ep_energy_prod25->push_back(ep_energy_prod25[ielec]);
                               tree_ep_25xtalindex->push_back(ielec);

                               if(ielec==0)
                                 {
                                   tree_ep0_energy_prod25_ratio->push_back(ep_energy_prod25_ratio[ielec]);
                                   tree_ep0_energy_prod25->push_back(ep_energy_prod25[ielec]);
                                 }

                               if(ielec==1)
                                 {
                                   tree_ep1_energy_prod25_ratio->push_back(ep_energy_prod25_ratio[ielec]);
                                   tree_ep1_energy_prod25->push_back(ep_energy_prod25[ielec]);
                                 }

			     }//if( flag_fill25 == 1)
			 }//if(pat_electronptalsize[ielec]>=25)
		     }// if( pat_electronptalsize[ielec] != 0 )
		 } //for(int ielec=0; ielec<pat_electronpsize;ielec++)
	     }//if( pat_electronpsize!=0 )
           //=======================for e+ ends here===========================




           //===========================for photon ============================================
           double pho_energy_prod_ratio[10];
           double pho_energy_prod[10];
           double pho_energy_prod8_ratio[10];
           double pho_energy_prod8[10];
           double pho_energy_prod25_ratio[10];
           double pho_energy_prod25[10];
           for( int ii=0;ii<10;ii++ )
             {
               pho_energy_prod_ratio[ii]=1;
               pho_energy_prod[ii]=0;
               pho_energy_prod8_ratio[ii]=1;
               pho_energy_prod8[ii]=0;
               pho_energy_prod25_ratio[ii]=1;
               pho_energy_prod25[ii]=0;
	       
             }

	   //clear the vectors
	   tree_pho_energy_prod_ratio->clear();
	   tree_pho_energy_prod->clear();
	   
	   tree_pho_energy_prod8_ratio->clear();
	   tree_pho_energy_prod8->clear();
	   
	   tree_pho_energy_prod25_ratio->clear();
	   tree_pho_energy_prod25->clear();
	   
	   tree_pho_index->clear();
           tree_pho_8xtalindex->clear();
           tree_pho_25xtalindex->clear();


	   tree_pho0_energy_prod_ratio->clear();
           tree_pho0_energy_prod->clear();
           tree_pho1_energy_prod_ratio->clear();
           tree_pho1_energy_prod->clear();

           tree_pho0_energy_prod8_ratio->clear();
           tree_pho0_energy_prod8->clear();
           tree_pho1_energy_prod8_ratio->clear();
           tree_pho1_energy_prod8->clear();

           tree_pho0_energy_prod25_ratio->clear();
           tree_pho0_energy_prod25->clear();
           tree_pho1_energy_prod25_ratio->clear();
           tree_pho1_energy_prod25->clear();

	   if( pat_photonsize !=0 )
	     {
	       //resize various vectors
               tree_pho_energy_prod_ratio->resize(pat_photonsize);
               tree_pho_energy_prod->resize(pat_photonsize);

               tree_pho_energy_prod8_ratio->resize(pat_photonsize);
               tree_pho_energy_prod8->resize(pat_photonsize);

               tree_pho_energy_prod25_ratio->resize(pat_photonsize);
               tree_pho_energy_prod25->resize(pat_photonsize);
	       
	       for( int ipho=0; ipho<pat_photonsize;ipho++ )
		 {
		   int flag_fill   = 0;
		   int flag_fill8  = 0;
		   int flag_fill25 = 0;
		   int index[1000];
		   if( pat_photonxtalsize[ipho] != 0 )
		     {
		       TMath::Sort(pat_photonxtalsize[ipho],pat_photonxtalE[ipho],index);
		       
		       //construct that qty for all crystals
		       for ( int jcrys = 1; jcrys < pat_photonxtalsize[ipho]; jcrys++ )
			 {
			   //if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index[jcrys]]) && (pat_photonxtalE[0][index[0]] == pat_photonxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			   if( (!isnan(pat_photonxtalE[ipho][index[jcrys]]) && !isnan(pat_photonxtalE[ipho][index[0]])) &&  (!isinf(pat_photonxtalE[ipho][index[jcrys]]) && !isinf(pat_photonxtalE[ipho][index[0]]) ) )
			     {
			       if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho][index[jcrys]] > 0.09  && pat_photonxtalE[ipho][index[0]] < 20000 && pat_photonxtalE[ipho][index[0]] > 0.09 )
				 {
				   flag_fill = 1;
				   pho_energy_prod_ratio[ipho] *= (pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]])/pat_photonxtalE[ipho][index[0]];
				   pho_energy_prod[ipho] += log(pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]]);
				 } //if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho][index...])
			     }//if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index..])
			 }//for (jcrys = 1; jcrys < pat_photonxtalsize[ipho]; jcrys++ )
		       
		       if( flag_fill == 1 )
			 {
			   //(*tree_pho_energy_prod_ratio)[ipho] = (pho_energy_prod_ratio[ipho]);
			   //(*tree_pho_energy_prod)[ipho] = (pho_energy_prod[ipho]);
			   //cout<<"filled tree_em_energy_prod_ratio & prod only:pho"<<endl;
			   //tree_pho_index->push_back(ipho);

			   tree_pho_energy_prod_ratio->push_back(pho_energy_prod_ratio[ipho]);
                           tree_pho_energy_prod->push_back(pho_energy_prod[ipho]);
                           tree_pho_index->push_back(ipho);

                           if(ipho==0)
                             {
                               tree_pho0_energy_prod_ratio->push_back(pho_energy_prod_ratio[ipho]);
                               tree_pho0_energy_prod->push_back(pho_energy_prod[ipho]);
                             }

                           if(ipho==1)
                             {
                               tree_pho1_energy_prod_ratio->push_back(pho_energy_prod_ratio[ipho]);
                               tree_pho1_energy_prod->push_back(pho_energy_prod[ipho]);
                             }

			   //i have got E1x5/E5x5 and E2x5/E5x5 for every electron to get some scatter plots 
			   //get the same when eId and pID is true
			 }//if( flag_fill == 1)
		       
		       //construct that qty for 8 highest energy crstals
		       if(pat_photonxtalsize[ipho]>=8)
			 {
			   for ( int jcrys = 1; jcrys < 8; jcrys++ )
			     {
			       //if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index[jcrys]]) && (pat_photonxtalE[0][index[0]] == pat_photonxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_photonxtalE[ipho][index[jcrys]]) && !isnan(pat_photonxtalE[ipho][index[0]])) && (!isinf(pat_photonxtalE[ipho][index[jcrys]]) && !isinf(pat_photonxtalE[ipho][index[0]]) ) )
				 {
				   
				   if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho][index[jcrys]] > 0.09 && pat_photonxtalE[ipho][index[0]] < 20000 && pat_photonxtalE[ipho][index[0]] > 0.09 )
				   {
				     flag_fill8 = 1;
				     pho_energy_prod8_ratio[ipho] *= (pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]])/pat_photonxtalE[ipho][index[0]];
				       pho_energy_prod8[ipho] += log(pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]]);
				   } //if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho]..)
				 }// if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index...])
			     }//for (jcrys = 1; jcrys < 8; jcrys++)
			   
			   if( flag_fill8 == 1)
			     {
			       //(*tree_pho_energy_prod8_ratio)[ipho] = (pho_energy_prod8_ratio[ipho]);
			       //(*tree_pho_energy_prod8)[ipho] = (pho_energy_prod8[ipho]);
			       //cout<<"filled tree_em_energy_prod8_ratio & prod only:pho"<<endl;
			       //tree_pho_8xtalindex->push_back(ipho);

			       tree_pho_energy_prod8_ratio->push_back(pho_energy_prod8_ratio[ipho]);
			       tree_pho_energy_prod8->push_back(pho_energy_prod8[ipho]);
			       tree_pho_8xtalindex->push_back(ipho);

			       if(ipho==0)
				 {
				   tree_pho0_energy_prod8_ratio->push_back(pho_energy_prod8_ratio[ipho]);
				   tree_pho0_energy_prod8->push_back(pho_energy_prod8[ipho]);
				 }

			       if(ipho==1)
				 {
				   tree_pho1_energy_prod8_ratio->push_back(pho_energy_prod8_ratio[ipho]);
				   tree_pho1_energy_prod8->push_back(pho_energy_prod8[ipho]);
				 }

			     }//if( flag_fill8 == 1)
			 }//if(pat_photonxtalsize[ipho]>=8)
		       
		       //construct that qty for 25 highest energy crstals
		       if(pat_photonxtalsize[ipho]>=25)
			 {
			   for ( int jcrys = 1; jcrys < 25; jcrys++ )
			     {
			       //if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index[jcrys]]) && (pat_photonxtalE[0][index[0]] == pat_photonxtalE[0][index[0]]) ) //if the xtalE is nan then nan==nan retrns false
			       if( (!isnan(pat_photonxtalE[ipho][index[jcrys]]) && !isnan(pat_photonxtalE[ipho][index[0]])) && (!isinf(pat_photonxtalE[ipho][index[jcrys]]) && !isinf(pat_photonxtalE[ipho][index[0]])) )
				 {
				   if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho][index[jcrys]] > 0.09 && pat_photonxtalE[ipho][index[0]] < 20000 && pat_photonxtalE[ipho][index[0]] > 0.09 )
				     {
				       flag_fill25 =1;
				       pho_energy_prod25_ratio[ipho] *= (pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]])/pat_photonxtalE[ipho][index[0]];
				       pho_energy_prod25[ipho] += log(pat_photonxtalE[ipho][index[0]] - pat_photonxtalE[ipho][index[jcrys]]);
				     }//if( pat_photonxtalE[ipho][index[jcrys]] < 20000 && pat_photonxtalE[ipho][index])
				 }//if( (pat_photonxtalE[ipho][index[jcrys]] == pat_photonxtalE[ipho][index..])
			     }//for (jcrys = 1; jcrys < 25; jcrys++)

			   if( flag_fill25 == 1)
			     {
			       //(*tree_pho_energy_prod25_ratio)[ipho] = (pho_energy_prod25_ratio[ipho]);
			       //(*tree_pho_energy_prod25)[ipho] = (pho_energy_prod25[ipho]);
			       //cout<<"filled tree_em_energy_prod25_ratio & prod only:pho"<<endl;
			       //tree_pho_25xtalindex->push_back(ipho);

			       tree_pho_energy_prod25_ratio->push_back(pho_energy_prod25_ratio[ipho]);
                               tree_pho_energy_prod25->push_back(pho_energy_prod25[ipho]);
                               tree_pho_25xtalindex->push_back(ipho);

                               if(ipho==0)
                                 {
                                   tree_pho0_energy_prod25_ratio->push_back(pho_energy_prod25_ratio[ipho]);
                                   tree_pho0_energy_prod25->push_back(pho_energy_prod25[ipho]);
                                 }

                               if(ipho==1)
                                 {
                                   tree_pho1_energy_prod25_ratio->push_back(pho_energy_prod25_ratio[ipho]);
                                   tree_pho1_energy_prod25->push_back(pho_energy_prod25[ipho]);
                                 }

			     }//if( flag_fill25 == 1)
			 }//if(pat_photonxtalsize[ipho]>=25)
		     }// if( pat_photonxtalsize[ipho] != 0 )
		 } //for(int ipho=0; ipho<pat_photonsize;ipho++)
	     }//if( pat_photonsize!=0 )
           //=======================for photon ends here===========================

	   //cout<<"filled all the reco variables"<<endl;
	   //cout<<""<<endl;
	   */
	   //===================for z =================================
	   
	   //(here z denotes any particle whose all kinematic qtys are reconstructed by two opppositely charged e's as given below) 

	   //clear the vectors
	   tree_zE->clear();
           tree_zpt->clear();
           tree_zpx->clear();
           tree_zpy->clear();
           tree_zpz->clear();
           tree_zp->clear();
           tree_zmass->clear();
           tree_zeta->clear();
           tree_zphi->clear();

	   //fill the vectors now
	   double zE = pat_electronE[0] + pat_electronE[1];
	   tree_zE->push_back(zE);

	   double zpt = sqrt( pow((pat_electron_px[0]+pat_electron_px[1]),2) + pow((pat_electron_py[0]+pat_electron_py[1]),2));
           tree_zpt->push_back(zpt);

	   double zpx = pat_electron_px[0] + pat_electron_px[1];
	   tree_zpx->push_back(zpx);

	   double zpy = pat_electron_py[0] + pat_electron_py[1];
           tree_zpy->push_back(zpy);

	   double zpz = pat_electronpz[0] + pat_electronpz[1];
           tree_zpz->push_back(zpz);

	   double zp =  sqrt( pow((pat_electron_px[0]+pat_electron_px[1]),2) + pow((pat_electron_py[0]+pat_electron_py[1]),2) + pow((pat_electronpz[0]+pat_electronpz[1]),2) );
           tree_zp->push_back(zp);

	   double zmass = sqrt( pow(zE,2) - pow(zpx,2) - pow(zpy,2) - pow(zpz,2));
           tree_zmass->push_back(zmass);

	   double ztheta = acos(zpz/zp);
	   double zeta   = -log(tan(ztheta/2.0));
           tree_zeta->push_back(zeta);

	   double zphi = acos(zpx/zpt);
           zphi        = correct_acosphi(zpy,zphi);
           tree_zphi->push_back(zphi);
	   
	   //====================z ends here==================

	   //===========================begin of e- gamma ===========================================================

	   //my convention : estarm - highest pt electron
	   //                estarp - nxt to highest pt electron

	   double estarm_E         = pat_electronE[0] + pat_photonE[0];
           double estarm_px        = pat_electron_px[0] + pat_photon_px[0];
           double estarm_py        = pat_electron_py[0] + pat_photon_py[0];
           double estarm_pz        = pat_electronpz[0] + pat_photonpz[0];
           double estarm_pt        = sqrt( pow(estarm_px,2) + pow(estarm_py,2) );
           double estarm_p         = sqrt( pow(estarm_px,2) + pow(estarm_py,2) + pow(estarm_pz,2) );
           double estarm_mass      = sqrt( pow(estarm_E,2) - pow(estarm_p,2) );

	   //clear the vectors
	   tree_patestarmE->clear();
	   tree_patestarmpx->clear();
           tree_patestarmpy->clear();
           tree_patestarmpz->clear();
           tree_patestarmpt->clear();
           tree_patestarmp->clear();
           tree_patestarmmass->clear();

	   //now fill the vectors
           tree_patestarmE->push_back(estarm_E);
           tree_patestarmpx->push_back(estarm_px);
           tree_patestarmpy->push_back(estarm_py);
           tree_patestarmpz->push_back(estarm_pz);
           tree_patestarmpt->push_back(estarm_pt);
           tree_patestarmp->push_back(estarm_p);
           tree_patestarmmass->push_back(estarm_mass);

	   //clear the vector
	   tree_patestarmeta->clear();

	   Double_t theta_estarm   = acos(estarm_pz/estarm_p);
	   double   estarmeta      = -log(tan(theta_estarm/2.0));  
           tree_patestarmeta->push_back(estarmeta);

	   //==========================end of e- gamma==============================================================

	   //===========================begin of e+ gamma ===========================================================

           double estarp_E         = pat_electronE[1] + pat_photonE[0];
           double estarp_px        = pat_electron_px[1] + pat_photon_px[0];
           double estarp_py        = pat_electron_py[1] + pat_photon_py[0];
           double estarp_pz        = pat_electronpz[1] + pat_photonpz[0];
           double estarp_pt        = sqrt( pow(estarp_px,2) + pow(estarp_py,2) );
           double estarp_p         = sqrt( pow(estarp_px,2) + pow(estarp_py,2) + pow(estarp_pz,2) );
           double estarp_mass      = sqrt( pow(estarp_E,2) - pow(estarp_p,2) );

	   //clear the vectors
           tree_patestarpE->clear();
           tree_patestarppx->clear();
           tree_patestarppy->clear();
           tree_patestarppz->clear();
           tree_patestarppt->clear();
           tree_patestarpp->clear();
           tree_patestarpmass->clear();

           //now fill the vectors
           tree_patestarpE->push_back(estarp_E);
           tree_patestarppx->push_back(estarp_px);
           tree_patestarppy->push_back(estarp_py);
           tree_patestarppz->push_back(estarp_pz);
           tree_patestarppt->push_back(estarp_pt);
           tree_patestarpp->push_back(estarp_p);
           tree_patestarpmass->push_back(estarp_mass);

	   //clear the vector
	   tree_patestarpeta->clear();

	   Double_t theta_estarp   = acos(estarp_pz/estarp_p);
	   double   estarpeta      = -log(tan(theta_estarp/2.0));  
           tree_patestarpeta->push_back(estarpeta);

	   

           //==========================end of e- gamma==============================================================
	   
	   //===========================begin of jet, zjet, jetgam =============================================

           if(flag_j==1)
             {

	       Double_t phijet         = correct_phi(pat_jetphi[0]);
	       Double_t phi_zjet       = deltaphi(phijet,zphi);
	       Double_t detazjet       = zeta - pat_jeteta[0];
	       Double_t delta_r_zjet   = sqrt( pow(phi_zjet,2) + pow(detazjet,2));
	       Double_t phig           = correct_phi(pat_photonphi[0]);
	       Double_t phi_gjet       = deltaphi(phijet,phig);
	       Double_t deta           = pat_photoneta[0] - pat_jeteta[0];
	       Double_t deltar_jetgam  = sqrt( pow(phi_gjet,2) + pow(deta,2) );

               //jet 
	       //clear the vector
	       tree_patjetphi->clear();
               tree_patjetpt->clear();
               tree_patjeteta->clear();
	       
	       //fill the vector
               tree_patjetphi->push_back(phijet);
               tree_patjetpt->push_back(pat_jetpt[0]);
               tree_patjeteta->push_back(pat_jeteta[0]);

               //zjet 
	       //clear the vectors
	       tree_patzjet_phi->clear();
               tree_patzjet_deltar->clear();

	       //now fill the vector
               tree_patzjet_phi->push_back(phi_zjet);
               tree_patzjet_deltar->push_back(delta_r_zjet);

               //jetgam 
	       //clear the vector
	       tree_patjetgam_phi->clear();
               tree_patjetgam_deltar->clear();

               tree_patjetgam_phi->push_back(phi_gjet);
               tree_patjetgam_deltar->push_back(deltar_jetgam);
             }

           //==========================end of jet, zjet, jetgam ================================================
	   
	   
	   //cout<<"now filling the root tree"<<endl;
	   myTree->Fill();
	   cout<<"now filled the root tree"<<endl;
	 }//if(flag_e == 1 && flag_p == 1)

   } //for (Long64_t jentry=0; jentry<nentries;jentry++)

   //rootfile[fileindex]->cd();
   cout<<"going to rootfile directory"<<endl;
   rootfile->cd();
   cout<<"writing myTree"<<endl;
   myTree->Write();
   writehisto();
   rootfile->Close();
   cout<<"rootfile closed"<<endl;
   //rootfile[fileindex]->Close();

   return;
} //void myClass::Loop()


void myClass::setbranch( TTree*& myTree )
{
  using namespace std;

  //store eID and pID
  myTree->Branch("tree_eID","vector<int>",&tree_eID);
  myTree->Branch("tree_pID","vector<int>",&tree_pID);

  //photon
  myTree->Branch("tree_genZdaughterE","vector<double>",&tree_genZdaughterE);
  myTree->Branch("tree_patphotonpt","vector<double>",&tree_patphotonpt);
  myTree->Branch("tree_patphotonpx","vector<double>",&tree_patphotonpx);
  myTree->Branch("tree_patphotonpy","vector<double>",&tree_patphotonpy);
  myTree->Branch("tree_patphotonpz","vector<double>",&tree_patphotonpz);
  myTree->Branch("tree_patphotoneta","vector<double>",&tree_patphotoneta);
  myTree->Branch("tree_patphotonphi","vector<double>",&tree_patphotonphi);
  //myTree->Branch("tree_patphotonSCeta","vector<double>",&tree_patphotonSCeta);
  //myTree->Branch("tree_patphotonSCphi","vector<double>",&tree_patphotonSCphi);
  myTree->Branch("tree_patphotonE","vector<double>",&tree_patphotonE);
  myTree->Branch("tree_patphotonHoverE","vector<double>",&tree_patphotonHoverE);

  myTree->Branch("tree_patphotonEcalRecHit03","vector<double>",&tree_patphotonEcalRecHit03);
  myTree->Branch("tree_patphotonHcalRecHit03","vector<double>",&tree_patphotonHcalRecHit03);
  myTree->Branch("tree_patphotonIsoHollowtrkcone03","vector<double>",&tree_patphotonIsoHollowtrkcone03);
  myTree->Branch("tree_patphotonEcalRecHit04","vector<double>",&tree_patphotonEcalRecHit04);
  myTree->Branch("tree_patphotonHcalRecHit04","vector<double>",&tree_patphotonHcalRecHit04);
  myTree->Branch("tree_patphotonHcalDepth104","vector<double>",&tree_patphotonHcalDepth104);
  myTree->Branch("tree_patphotonHcalDepth204","vector<double>",&tree_patphotonHcalDepth204);
  
  myTree->Branch("tree_patphotonIsoHollowtrkcone04","vector<double>",&tree_patphotonIsoHollowtrkcone04);
  myTree->Branch("tree_patphotonSigmaIeta","vector<double>",&tree_patphotonSigmaIeta);
  myTree->Branch("tree_patphotonisConverted","vector<Char_t>",&tree_patphotonisConverted);
  
  //e-
  myTree->Branch("tree_pateeta","vector<double>",&tree_pateeta);
  myTree->Branch("tree_patephi","vector<double>",&tree_patephi);
  myTree->Branch("tree_pateE","vector<double>",&tree_pateE);
  myTree->Branch("tree_patept","vector<double>",&tree_patept);
  myTree->Branch("tree_patepx","vector<double>",&tree_patepx);
  myTree->Branch("tree_patepy","vector<double>",&tree_patepy);
  myTree->Branch("tree_patepz","vector<double>",&tree_patepz);
  myTree->Branch("tree_pateEt","vector<double>",&tree_pateEt);
  myTree->Branch("tree_pateSCeta","vector<double>",&tree_pateSCeta);
  //myTree->Branch("tree_pateSCphi","vector<double>",&tree_pateSCphi);
  myTree->Branch("tree_patedeta_in","vector<double>",&tree_patedeta_in);
  myTree->Branch("tree_patedphi_in","vector<double>",&tree_patedphi_in);
  myTree->Branch("tree_pateHoverE","vector<double>",&tree_pateHoverE);
  myTree->Branch("tree_patesigIetaIeta","vector<double>",&tree_patesigIetaIeta);
  myTree->Branch("tree_patesigetaeta","vector<double>",&tree_patesigetaeta);
  myTree->Branch("tree_pateratio_E1x5OverE5x5","vector<double>",&tree_pateratio_E1x5OverE5x5);
  myTree->Branch("tree_pateratio_E2x5OverE5x5","vector<double>",&tree_pateratio_E2x5OverE5x5);

  myTree->Branch("tree_patedr03EcalRecHitSumEt","vector<double>",&tree_patedr03EcalRecHitSumEt);
  myTree->Branch("tree_patedr03HcalDepth1TowerSumEt","vector<double>",&tree_patedr03HcalDepth1TowerSumEt);
  myTree->Branch("tree_patedr03HcalDepth2TowerSumEt","vector<double>",&tree_patedr03HcalDepth2TowerSumEt);
  myTree->Branch("tree_patedr03HcalTowerSumEt","vector<double>",&tree_patedr03HcalTowerSumEt);
  myTree->Branch("tree_patedr03TkSumPt","vector<double>",&tree_patedr03TkSumPt);

  myTree->Branch("tree_patedr04EcalRecHitSumEt","vector<double>",&tree_patedr04EcalRecHitSumEt);
  myTree->Branch("tree_patedr04HcalDepth1TowerSumEt","vector<double>",&tree_patedr04HcalDepth1TowerSumEt);
  myTree->Branch("tree_patedr04HcalDepth2TowerSumEt","vector<double>",&tree_patedr04HcalDepth2TowerSumEt);
  myTree->Branch("tree_patedr04HcalTowerSumEt","vector<double>",&tree_patedr04HcalTowerSumEt);
  myTree->Branch("tree_patedr04TkSumPt","vector<double>",&tree_patedr04TkSumPt);

  myTree->Branch("tree_patecharge","vector<int>",&tree_patecharge);
  myTree->Branch("tree_pateisEcalDriven","vector<Char_t>",&tree_pateisEcalDriven);

  //z
  myTree->Branch("tree_zE","vector<double>",&tree_zE);
  myTree->Branch("tree_zpt","vector<double>",&tree_zpt);
  myTree->Branch("tree_zpx","vector<double>",&tree_zpx);
  myTree->Branch("tree_zpy","vector<double>",&tree_zpy);
  myTree->Branch("tree_zpz","vector<double>",&tree_zpz);
  myTree->Branch("tree_zp","vector<double>",&tree_zp);
  myTree->Branch("tree_zmass","vector<double>",&tree_zmass);
  myTree->Branch("tree_zeta","vector<double>",&tree_zeta);
  myTree->Branch("tree_zphi","vector<double>",&tree_zphi);

  //estarm
  myTree->Branch("tree_patestarmE","vector<double>",&tree_patestarmE);
  myTree->Branch("tree_patestarmpx","vector<double>",&tree_patestarmpx);
  myTree->Branch("tree_patestarmpy","vector<double>",&tree_patestarmpy);
  myTree->Branch("tree_patestarmpz","vector<double>",&tree_patestarmpz);
  myTree->Branch("tree_patestarmpt","vector<double>",&tree_patestarmpt);
  myTree->Branch("tree_patestarmp","vector<double>",&tree_patestarmp);
  myTree->Branch("tree_patestarmeta","vector<double>",&tree_patestarmeta);
  myTree->Branch("tree_patestarmmass","vector<double>",&tree_patestarmmass);
 //estarp
  myTree->Branch("tree_patestarpE","vector<double>",&tree_patestarpE);
  myTree->Branch("tree_patestarppx","vector<double>",&tree_patestarppx);
  myTree->Branch("tree_patestarppy","vector<double>",&tree_patestarppy);
  myTree->Branch("tree_patestarppz","vector<double>",&tree_patestarppz);
  myTree->Branch("tree_patestarppt","vector<double>",&tree_patestarppt);
  myTree->Branch("tree_patestarpp","vector<double>",&tree_patestarpp);
  myTree->Branch("tree_patestarpeta","vector<double>",&tree_patestarpeta);
  myTree->Branch("tree_patestarpmass","vector<double>",&tree_patestarpmass);

  //jet,zjet, jetgam
  //jet
  myTree->Branch("tree_patjetphi","vector<double>",&tree_patjetphi);
  myTree->Branch("tree_patjetpt","vector<double>",&tree_patjetpt);
  myTree->Branch("tree_patjeteta","vector<double>",&tree_patjeteta);

  //zjet
  myTree->Branch("tree_patzjet_phi","vector<double>",&tree_patzjet_phi);
  myTree->Branch("tree_patzjet_deltar","vector<double>",&tree_patzjet_deltar);

  //jetgam
  myTree->Branch("tree_patjetgam_phi","vector<double>",&tree_patjetgam_phi);
  myTree->Branch("tree_patjetgam_deltar","vector<double>",&tree_patjetgam_deltar);
  
  //cout<<"set all the vector<double> branches"<<endl;

  //our constructed variables
  //e-
  /*myTree->Branch("tree_em_ratio_E1x5OverE5x5","vector<double>",&tree_em_ratio_E1x5OverE5x5);
  myTree->Branch("tree_em_ratio_E2x5OverE5x5","vector<double>",&tree_em_ratio_E2x5OverE5x5);

  myTree->Branch("tree_em_energy_prod_ratio","vector<double>",&tree_em_energy_prod_ratio);
  myTree->Branch("tree_em_energy_prod","vector<double>",&tree_em_energy_prod);
  myTree->Branch("tree_em_energy_prod8_ratio","vector<double>",&tree_em_energy_prod8_ratio);
  myTree->Branch("tree_em_energy_prod8","vector<double>",&tree_em_energy_prod8);
  myTree->Branch("tree_em_energy_prod25_ratio","vector<double>",&tree_em_energy_prod25_ratio);
  myTree->Branch("tree_em_energy_prod25","vector<double>",&tree_em_energy_prod25);

  myTree->Branch("tree_em0_energy_prod_ratio","vector<double>",&tree_em0_energy_prod_ratio);
  myTree->Branch("tree_em0_energy_prod","vector<double>",&tree_em0_energy_prod);
  myTree->Branch("tree_em1_energy_prod_ratio","vector<double>",&tree_em1_energy_prod_ratio);
  myTree->Branch("tree_em1_energy_prod","vector<double>",&tree_em1_energy_prod);


  myTree->Branch("tree_em0_energy_prod8_ratio","vector<double>",&tree_em0_energy_prod8_ratio);
  myTree->Branch("tree_em0_energy_prod8","vector<double>",&tree_em0_energy_prod8);
  myTree->Branch("tree_em1_energy_prod8_ratio","vector<double>",&tree_em1_energy_prod8_ratio);
  myTree->Branch("tree_em1_energy_prod8","vector<double>",&tree_em1_energy_prod8);

  myTree->Branch("tree_em0_energy_prod25_ratio","vector<double>",&tree_em0_energy_prod25_ratio);
  myTree->Branch("tree_em0_energy_prod25","vector<double>",&tree_em0_energy_prod25);
  myTree->Branch("tree_em25_energy_prod8_ratio","vector<double>",&tree_em1_energy_prod25_ratio);
  myTree->Branch("tree_em25_energy_prod8","vector<double>",&tree_em1_energy_prod25);


  myTree->Branch("tree_em_index","vector<int>",&tree_em_index);
  myTree->Branch("tree_em_8xtalindex","vector<int>",&tree_em_8xtalindex);
  myTree->Branch("tree_em_25xtalindex","vector<int>",&tree_em_25xtalindex);

  


  //e+
  myTree->Branch("tree_ep_ratio_E1x5OverE5x5","vector<double>",&tree_ep_ratio_E1x5OverE5x5);
  myTree->Branch("tree_ep_ratio_E2x5OverE5x5","vector<double>",&tree_ep_ratio_E2x5OverE5x5);

  myTree->Branch("tree_ep_energy_prod_ratio","vector<double>",&tree_ep_energy_prod_ratio);
  myTree->Branch("tree_ep_energy_prod","vector<double>",&tree_ep_energy_prod);

  myTree->Branch("tree_ep_energy_prod8_ratio","vector<double>",&tree_ep_energy_prod8_ratio);
  myTree->Branch("tree_ep_energy_prod8","vector<double>",&tree_ep_energy_prod8);

  myTree->Branch("tree_ep_energy_prod25_ratio","vector<double>",&tree_ep_energy_prod25_ratio);
  myTree->Branch("tree_ep_energy_prod25","vector<double>",&tree_ep_energy_prod25);

  myTree->Branch("tree_ep0_energy_prod_ratio","vector<double>",&tree_ep0_energy_prod_ratio);
  myTree->Branch("tree_ep0_energy_prod","vector<double>",&tree_ep0_energy_prod);
  myTree->Branch("tree_ep1_energy_prod_ratio","vector<double>",&tree_ep1_energy_prod_ratio);
  myTree->Branch("tree_ep1_energy_prod","vector<double>",&tree_ep1_energy_prod);


  myTree->Branch("tree_ep0_energy_prod8_ratio","vector<double>",&tree_ep0_energy_prod8_ratio);
  myTree->Branch("tree_ep0_energy_prod8","vector<double>",&tree_ep0_energy_prod8);
  myTree->Branch("tree_ep1_energy_prod8_ratio","vector<double>",&tree_ep1_energy_prod8_ratio);
  myTree->Branch("tree_ep1_energy_prod8","vector<double>",&tree_ep1_energy_prod8);

  myTree->Branch("tree_ep0_energy_prod25_ratio","vector<double>",&tree_ep0_energy_prod25_ratio);
  myTree->Branch("tree_ep0_energy_prod25","vector<double>",&tree_ep0_energy_prod25);
  myTree->Branch("tree_ep25_energy_prod8_ratio","vector<double>",&tree_ep1_energy_prod25_ratio);
  myTree->Branch("tree_ep25_energy_prod8","vector<double>",&tree_ep1_energy_prod25);

  myTree->Branch("tree_ep_index","vector<int>",&tree_ep_index);
  myTree->Branch("tree_ep_8xtalindex","vector<int>",&tree_ep_8xtalindex);
  myTree->Branch("tree_ep_25xtalindex","vector<int>",&tree_ep_25xtalindex);

  //photon
  myTree->Branch("tree_pho_energy_prod_ratio","vector<double>",&tree_pho_energy_prod_ratio);
  myTree->Branch("tree_pho_energy_prod","vector<double>",&tree_pho_energy_prod);

  myTree->Branch("tree_pho_energy_prod8_ratio","vector<double>",&tree_pho_energy_prod8_ratio);
  myTree->Branch("tree_pho_energy_prod8","vector<double>",&tree_pho_energy_prod8);

  myTree->Branch("tree_pho_energy_prod25_ratio","vector<double>",&tree_pho_energy_prod25_ratio);
  myTree->Branch("tree_pho_energy_prod25","vector<double>",&tree_pho_energy_prod25);

  myTree->Branch("tree_pho0_energy_prod_ratio","vector<double>",&tree_pho0_energy_prod_ratio);
  myTree->Branch("tree_pho0_energy_prod","vector<double>",&tree_pho0_energy_prod);
  myTree->Branch("tree_pho1_energy_prod_ratio","vector<double>",&tree_pho1_energy_prod_ratio);
  myTree->Branch("tree_pho1_energy_prod","vector<double>",&tree_pho1_energy_prod);


  myTree->Branch("tree_pho0_energy_prod8_ratio","vector<double>",&tree_pho0_energy_prod8_ratio);
  myTree->Branch("tree_pho0_energy_prod8","vector<double>",&tree_pho0_energy_prod8);
  myTree->Branch("tree_pho1_energy_prod8_ratio","vector<double>",&tree_pho1_energy_prod8_ratio);
  myTree->Branch("tree_pho1_energy_prod8","vector<double>",&tree_pho1_energy_prod8);

  myTree->Branch("tree_pho0_energy_prod25_ratio","vector<double>",&tree_pho0_energy_prod25_ratio);
  myTree->Branch("tree_pho0_energy_prod25","vector<double>",&tree_pho0_energy_prod25);
  myTree->Branch("tree_pho25_energy_prod8_ratio","vector<double>",&tree_pho1_energy_prod25_ratio);
  myTree->Branch("tree_pho25_energy_prod8","vector<double>",&tree_pho1_energy_prod25);

  myTree->Branch("tree_pho_index","vector<int>",&tree_pho_index);
  myTree->Branch("tree_pho_8xtalindex","vector<int>",&tree_pho_8xtalindex);
  myTree->Branch("tree_pho_25xtalindex","vector<int>",&tree_pho_25xtalindex);
  */
  //cout<<"setting the map branch"<<endl;
  myTree->Branch("triggermap","map<string,int>",&triggermap); 
  //cout<<"set the map branch"<<endl;
  
  //myTree->Branch("try",try,"try[10]/D");
  return;
}





void myClass::initializeTreeVar()
{

  //eID and pID

      tree_eID                          = new std::vector< int >;
      tree_pID                          = new std::vector< int >;


  //photon
  tree_genZdaughterE                    = new std::vector< double >;
  tree_patphotonpt                      = new std::vector< double >;
  tree_patphotonpx                      = new std::vector< double >;
  tree_patphotonpy                      = new std::vector< double >;
  tree_patphotonpz                      = new std::vector< double >;
  tree_patphotoneta                     = new std::vector< double >;
  tree_patphotonphi                     = new std::vector< double >;
  //tree_patphotonSCeta                   = new std::vector< double >;
  //tree_patphotonSCphi                   = new std::vector< double >;
  tree_patphotonE                       = new std::vector< double >;
  tree_patphotonHoverE                  = new std::vector< double >;
  
  tree_patphotonEcalRecHit03             = new std::vector< double >;
  tree_patphotonHcalRecHit03             = new std::vector< double >;
  tree_patphotonIsoHollowtrkcone03       = new std::vector< double >;
  tree_patphotonEcalRecHit04             = new std::vector< double >;
  tree_patphotonHcalRecHit04             = new std::vector< double >;
  tree_patphotonHcalDepth104             = new std::vector< double >;
  tree_patphotonHcalDepth204             = new std::vector< double >;
  tree_patphotonIsoHollowtrkcone04       = new std::vector< double >;
  tree_patphotonSigmaIeta                = new std::vector< double >;
  tree_patphotonisConverted              = new std::vector<Char_t>;

  //e-
  tree_pateeta                         = new std::vector< double >;
  tree_patephi                         = new std::vector< double >;
  tree_pateE                           = new std::vector< double >;
  tree_patept                          = new std::vector< double >;
  tree_patepx                          = new std::vector< double >;
  tree_patepy                          = new std::vector< double >;
  tree_patepz                          = new std::vector< double >;
  tree_pateEt                          = new std::vector< double >;
  tree_pateSCeta                       = new std::vector< double >;
  //tree_pateSCphi                       = new std::vector< double >;
  tree_patedeta_in                     = new std::vector< double >;
  tree_patedphi_in                     = new std::vector< double >;
  tree_pateHoverE                      = new std::vector< double >;
  tree_patesigIetaIeta                 = new std::vector< double >;
  tree_patesigetaeta                   = new std::vector< double >;
  tree_pateratio_E1x5OverE5x5          = new std::vector< double >;
  tree_pateratio_E2x5OverE5x5          = new std::vector< double >;

  tree_patedr03EcalRecHitSumEt         = new std::vector< double >;
  tree_patedr03HcalDepth1TowerSumEt    = new std::vector< double >;
  tree_patedr03HcalDepth2TowerSumEt    = new std::vector< double >;
  tree_patedr03HcalTowerSumEt          = new std::vector< double >;
  tree_patedr03TkSumPt                 = new std::vector< double >; 
  
  tree_patedr04EcalRecHitSumEt         = new std::vector< double >;
  tree_patedr04HcalDepth1TowerSumEt    = new std::vector< double >;
  tree_patedr04HcalDepth2TowerSumEt    = new std::vector< double >;
  tree_patedr04HcalTowerSumEt          = new std::vector< double >;
  tree_patedr04TkSumPt                 = new std::vector< double >; 

  tree_patecharge                      = new std::vector< int >;
  tree_pateisEcalDriven                = new std::vector< Char_t >;


  //z
  tree_zE                               = new std::vector< double >;
  tree_zpt                              = new std::vector< double >;
  tree_zpx                              = new std::vector< double >;
  tree_zpy                              = new std::vector< double >;
  tree_zpz                              = new std::vector< double >;
  tree_zp                               = new std::vector< double >;
  tree_zmass                            = new std::vector< double >;
  tree_zeta                             = new std::vector< double >;
  tree_zphi                             = new std::vector< double >;

  //estarm
  tree_patestarmE                       = new std::vector< double >;
  tree_patestarmpx                      = new std::vector< double >;
  tree_patestarmpy                      = new std::vector< double >;
  tree_patestarmpz                      = new std::vector< double >;
  tree_patestarmpt                      = new std::vector< double >;
  tree_patestarmp                       = new std::vector< double >;
  tree_patestarmmass                    = new std::vector< double >;
  tree_patestarmeta                     = new std::vector< double >;

  //estarp
  tree_patestarpE                       = new std::vector< double >;
  tree_patestarppx                      = new std::vector< double >;
  tree_patestarppy                      = new std::vector< double >;
  tree_patestarppz                      = new std::vector< double >;
  tree_patestarppt                      = new std::vector< double >;
  tree_patestarpp                       = new std::vector< double >;
  tree_patestarpmass                    = new std::vector< double >;
  tree_patestarpeta                     = new std::vector< double >;


  //jet,zjet,jetgam
  tree_patjetphi                        = new std::vector< double >;
  tree_patjetpt                         = new std::vector< double >;
  tree_patjeteta                        = new std::vector< double >;

  //zjet
  tree_patzjet_phi                      = new std::vector< double >;
  tree_patzjet_deltar                   = new std::vector< double >;

  //jetgam
  tree_patjetgam_phi                    = new std::vector< double >;
  tree_patjetgam_deltar                 = new std::vector< double >;
  
  /*//our constructed variables

      //e-
      tree_em_ratio_E1x5OverE5x5               = new std::vector< double >;
      tree_em_ratio_E2x5OverE5x5               = new std::vector< double >;
      tree_em_energy_prod_ratio          = new std::vector<double>;
      tree_em_energy_prod               = new std::vector<double>;
      tree_em_energy_prod8_ratio        = new std::vector<double>;
      tree_em_energy_prod8              = new std::vector<double>;
      tree_em_energy_prod25_ratio       = new std::vector<double>;
      tree_em_energy_prod25             = new std::vector<double>;

      tree_em0_energy_prod_ratio          = new std::vector<double>;
      tree_em0_energy_prod               = new std::vector<double>;
      tree_em1_energy_prod_ratio          = new std::vector<double>;
      tree_em1_energy_prod               = new std::vector<double>;

      tree_em0_energy_prod8_ratio        = new std::vector<double>;
      tree_em0_energy_prod8              = new std::vector<double>;
      tree_em1_energy_prod8_ratio          = new std::vector<double>;
      tree_em1_energy_prod8               = new std::vector<double>;

      tree_em0_energy_prod25_ratio       = new std::vector<double>;
      tree_em0_energy_prod25             = new std::vector<double>;
      tree_em1_energy_prod25_ratio          = new std::vector<double>;
      tree_em1_energy_prod25               = new std::vector<double>;

      
      tree_em_index                     = new std::vector< int >;
      tree_em_8xtalindex                = new std::vector< int >;
      tree_em_25xtalindex               = new std::vector< int >;

      //e+
      tree_ep_ratio_E1x5OverE5x5        = new std::vector<double>;
      tree_ep_ratio_E2x5OverE5x5        = new std::vector<double>;
      tree_ep_energy_prod_ratio         = new std::vector<double>;
      tree_ep_energy_prod               = new std::vector<double>;
      tree_ep_energy_prod8_ratio        = new std::vector<double>;
      tree_ep_energy_prod8              = new std::vector<double>;
      tree_ep_energy_prod25_ratio       = new std::vector<double>;
      tree_ep_energy_prod25             = new std::vector<double>;

      tree_ep0_energy_prod_ratio          = new std::vector<double>;
      tree_ep0_energy_prod               = new std::vector<double>;
      tree_ep1_energy_prod_ratio          = new std::vector<double>;
      tree_ep1_energy_prod               = new std::vector<double>;

      tree_ep0_energy_prod8_ratio        = new std::vector<double>;
      tree_ep0_energy_prod8              = new std::vector<double>;
      tree_ep1_energy_prod8_ratio          = new std::vector<double>;
      tree_ep1_energy_prod8               = new std::vector<double>;

      tree_ep0_energy_prod25_ratio       = new std::vector<double>;
      tree_ep0_energy_prod25             = new std::vector<double>;
      tree_ep1_energy_prod25_ratio          = new std::vector<double>;
      tree_ep1_energy_prod25               = new std::vector<double>;

      tree_ep_index                     = new std::vector< int >;
      tree_ep_8xtalindex                = new std::vector< int >;
      tree_ep_25xtalindex               = new std::vector< int >;

      //photon
      tree_pho_energy_prod_ratio         = new std::vector<double>;
      tree_pho_energy_prod               = new std::vector<double>;
      tree_pho_energy_prod8_ratio        = new std::vector<double>;
      tree_pho_energy_prod8              = new std::vector<double>;
      tree_pho_energy_prod25_ratio       = new std::vector<double>;
      tree_pho_energy_prod25             = new std::vector<double>;
      
      tree_pho0_energy_prod_ratio          = new std::vector<double>;
      tree_pho0_energy_prod               = new std::vector<double>;
      tree_pho1_energy_prod_ratio          = new std::vector<double>;
      tree_pho1_energy_prod               = new std::vector<double>;

      tree_pho0_energy_prod8_ratio        = new std::vector<double>;
      tree_pho0_energy_prod8              = new std::vector<double>;
      tree_pho1_energy_prod8_ratio          = new std::vector<double>;
      tree_pho1_energy_prod8               = new std::vector<double>;

      tree_pho0_energy_prod25_ratio       = new std::vector<double>;
      tree_pho0_energy_prod25             = new std::vector<double>;
      tree_pho1_energy_prod25_ratio          = new std::vector<double>;
      tree_pho1_energy_prod25               = new std::vector<double>;


      tree_pho_index                     = new std::vector< int >;
      tree_pho_8xtalindex                = new std::vector< int >;
      tree_pho_25xtalindex               = new std::vector< int >;

  */
  triggermap                             = new std::map<string,int>;
}



void myClass::bookhisto()
{
  //cout<<"inside the bookhisto"<<endl;
  char file_num[100];
  sprintf(file_num, "%d", fileindex);
  std::string i_file = file_num;
  hist2D_patelectronm_ieta_iphi["emIetaIphi_"+i_file] = new TH2D("hist2D_patelectronm_ieta_iphi", "rechit map of SC of e- for one event (iphi vs ieta) ",172,-86.0,86.0,361,0.0,361.0);
  //cout<<"emIetaIphi_+i_file = "<<"emIetaIphi_"+i_file<<endl;
  hist2D_patphoton_ieta_iphi["photonIetaIphi_"+i_file] = new TH2D("hist2D_patphoton_ieta_iphi", "rechit map of SC of photon for one event (iphi vs ieta) ",172,-86.0,86.0,361,0.0,361.0);  
  //cout<<"photonIetaIphi_+i_file = "<<"photonIetaIphi_"+i_file<<endl;
  
  for(int ii=0; ii<10; ii++)
    {
      char ievent[100];
      sprintf(ievent,"%d",ii);
      std::string i_event = ievent;
      TString ieventnum = ievent;
      hist2D_patelectronm_ieta_iphi_E2x5Ratio["emIetaIphi_"+i_file+"_"+i_event] = new TH2D("hist2D_patelectronm_ieta_iphi"+ieventnum, "rechit map of SC of e- for  event (iphi vs ieta) no. :  "+ieventnum,172,-86.0,86.0,361,0.0,361.0);
    }
  return;
}



void myClass::writehisto()
{
  char file_num[100];
  sprintf(file_num, "%d", fileindex);
  std::string i_file = file_num;
  hist2D_patelectronm_ieta_iphi["emIetaIphi_"+i_file]->Write();
  hist2D_patphoton_ieta_iphi["photonIetaIphi_"+i_file]->Write();
  
  for(int ii=0; ii<10; ii++)
    {
      char ievent[100];
      sprintf(ievent,"%d",ii);
      std::string i_event = ievent;
      TString ieventnum = ievent;
      hist2D_patelectronm_ieta_iphi_E2x5Ratio["emIetaIphi_"+i_file+"_"+i_event]->Write();
    }
  return;
}



