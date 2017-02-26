#define TwoTrees_cxx
#include "TwoTrees.h"
#include "commonHeaders.h"
#include "math.h"
void TwoTrees::Loop()
{

  //double gluonPtcut = 4.;
  //double gluonPtcut = 3.5;
  //double gluonPtcut = 5.0;
  //double gluonPtcut = 0.0;

  string region = "aaaa";

  //region = "EB_EB_EB";
  //region = "EB_EE_EB";
  region = "all";
  bool flag_region = 1;
  
  const int nbins = 30;

  double isobin[31];

  for(int ibin=0; ibin<30;ibin++)
    {
      isobin[ibin] = -5.0+ibin*2.0;
    }

  isobin[30] = 100;
  

  //data
  TH1F *dataphoPt = new TH1F("dataphoPt","dataphoPt", 18,20.,200.);
  TH1F *dataelePt = new TH1F("dataelePt","dataelePt", 18,20.,200.);
  TH1F *dataMaxMass = new TH1F("dataMaxMass","dataMaxMass", 48,20.,500.);
  TH1F *dataMee = new TH1F("dataMee","dataMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *dataeIso = new TH1F("dataeIso","dataeIso", 50,0.,100.);
  //TH1F *dataeIso = new TH1F("dataeIso","dataeIso", 50,-5.,100.);
  TH1F *dataeIso = new TH1F("dataeIso","dataeIso", nbins,isobin);

  TH1F *datahIso = new TH1F("datahIso","datahIso", 50,0.,100.);
  TH1F *datatIso = new TH1F("datatIso","datatIso", 50,0.,100.);
  TH1F *datahOe = new TH1F("datahOe","datahOe", 25,0.,1.);
  //TH1F *datasie = new TH1F("datasie","datasie", 100,0.,0.07);
  TH1F *datasie = new TH1F("datasie","datasie", 25,0.,0.03);

  //Zee
  TH1F *zeephoPt = new TH1F("zeephoPt","zeephoPt", 18,20.,200.);
  TH1F *zeeelePt = new TH1F("zeeelePt","zeeelePt", 18,20.,200.);
  TH1F *zeeMaxMass = new TH1F("zeeMaxMass","zeeMaxMass", 48,20.,500.);
  TH1F *zeeMee = new TH1F("zeeMee","zeeMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *zeeeIso = new TH1F("zeeeIso","zeeeIso", 50,0.,100.);
  //TH1F *zeeeIso = new TH1F("zeeeIso","zeeeIso", 50,-5.,100.);
  TH1F *zeeeIso = new TH1F("zeeeIso","zeeeIso", nbins,isobin);
  TH1F *zeehIso = new TH1F("zeehIso","zeehIso", 50,0.,100.);
  TH1F *zeetIso = new TH1F("zeetIso","zeetIso", 50,0.,100.);
  TH1F *zeehOe = new TH1F("zeehOe","zeehOe", 25,0.,1.);
  //TH1F *zeesie = new TH1F("zeesie","zeesie", 100,0.,0.07);
  TH1F *zeesie = new TH1F("zeesie","zeesie", 25,0.,0.03);

  //EWK
  /*TH1F *ewkphoPt = new TH1F("ewkphoPt","ewkphoPt", 18,20.,200.);
  TH1F *ewkelePt = new TH1F("ewkelePt","ewkelePt", 18,20.,200.);
  TH1F *ewkMaxMass = new TH1F("ewkMaxMass","ewkMaxMass", 48,20.,500.);
  TH1F *ewkMee = new TH1F("ewkMee","ewkMee", 48,20.,500.);

  ///isolation, ID
  TH1F *ewkeIso = new TH1F("ewkeIso","ewkeIso", 50,0.,100.);
  TH1F *ewkhIso = new TH1F("ewkhIso","ewkhIso", 50,0.,100.);
  TH1F *ewktIso = new TH1F("ewktIso","ewktIso", 50,0.,100.);
  TH1F *ewkhOe = new TH1F("ewkhOe","ewkhOe", 25,0.,1.);
  //TH1F *ewksie = new TH1F("ewksie","ewksie", 100,0.,0.07);
  TH1F *ewksie = new TH1F("ewksie","ewksie", 25,0.,0.03);
  */


  //ZZ
  TH1F *zzphoPt = new TH1F("zzphoPt","zzphoPt", 18,20.,200.);
  TH1F *zzelePt = new TH1F("zzelePt","zzelePt", 18,20.,200.);
  TH1F *zzMaxMass = new TH1F("zzMaxMass","zzMaxMass", 48,20.,500.);
  TH1F *zzMee = new TH1F("zzMee","zzMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *zzeIso = new TH1F("zzeIso","zzeIso", 50,0.,100.);
  //TH1F *zzeIso = new TH1F("zzeIso","zzeIso", 50,-5.,100.);
  TH1F *zzeIso = new TH1F("zzeIso","zzeIso", nbins,isobin);
  TH1F *zzhIso = new TH1F("zzhIso","zzhIso", 50,0.,100.);
  TH1F *zztIso = new TH1F("zztIso","zztIso", 50,0.,100.);
  TH1F *zzhOe = new TH1F("zzhOe","zzhOe", 25,0.,1.);
  //TH1F *zzsie = new TH1F("zzsie","zzsie", 100,0.,0.07);
  TH1F *zzsie = new TH1F("zzsie","zzsie", 25,0.,0.03);

  //WZ
  TH1F *wzphoPt = new TH1F("wzphoPt","wzphoPt", 18,20.,200.);
  TH1F *wzelePt = new TH1F("wzelePt","wzelePt", 18,20.,200.);
  TH1F *wzMaxMass = new TH1F("wzMaxMass","wzMaxMass", 48,20.,500.);
  TH1F *wzMee = new TH1F("wzMee","wzMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *wzeIso = new TH1F("wzeIso","wzeIso", 50,0.,100.);
  //TH1F *wzeIso = new TH1F("wzeIso","wzeIso", 50,-5.,100.);
  TH1F *wzeIso = new TH1F("wzeIso","wzeIso", nbins,isobin);
  TH1F *wzhIso = new TH1F("wzhIso","wzhIso", 50,0.,100.);
  TH1F *wztIso = new TH1F("wztIso","wztIso", 50,0.,100.);
  TH1F *wzhOe = new TH1F("wzhOe","wzhOe", 25,0.,1.);
  //TH1F *wzsie = new TH1F("wzsie","wzsie", 100,0.,0.07);
  TH1F *wzsie = new TH1F("wzsie","wzsie", 25,0.,0.03);

  //WW
  TH1F *wwphoPt = new TH1F("wwphoPt","wwphoPt", 18,20.,200.);
  TH1F *wwelePt = new TH1F("wwelePt","wwelePt", 18,20.,200.);
  TH1F *wwMaxMass = new TH1F("wwMaxMass","wwMaxMass", 48,20.,500.);
  TH1F *wwMee = new TH1F("wwMee","wwMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *wweIso = new TH1F("wweIso","wweIso", 50,0.,100.);
  //TH1F *wweIso = new TH1F("wweIso","wweIso", 50,-5.,100.);
  TH1F *wweIso = new TH1F("wweIso","wweIso", nbins,isobin);
  TH1F *wwhIso = new TH1F("wwhIso","wwhIso", 50,0.,100.);
  TH1F *wwtIso = new TH1F("wwtIso","wwtIso", 50,0.,100.);
  TH1F *wwhOe = new TH1F("wwhOe","wwhOe", 25,0.,1.);
  //TH1F *wwsie = new TH1F("wwsie","wwsie", 100,0.,0.07);
  TH1F *wwsie = new TH1F("wwsie","wwsie", 25,0.,0.03);
  
  //WENU
  TH1F *wenuphoPt = new TH1F("wenuphoPt","wenuphoPt", 18,20.,200.);
  TH1F *wenuelePt = new TH1F("wenuelePt","wenuelePt", 18,20.,200.);
  TH1F *wenuMaxMass = new TH1F("wenuMaxMass","wenuMaxMass", 48,20.,500.);
  TH1F *wenuMee = new TH1F("wenuMee","wenuMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *wenueIso = new TH1F("wenueIso","wenueIso", 50,0.,100.);
  //TH1F *wenueIso = new TH1F("wenueIso","wenueIso", 50,-5.,100.);
  TH1F *wenueIso = new TH1F("wenueIso","wenueIso", nbins,isobin);
  TH1F *wenuhIso = new TH1F("wenuhIso","wenuhIso", 50,0.,100.);
  TH1F *wenutIso = new TH1F("wenutIso","wenutIso", 50,0.,100.);
  TH1F *wenuhOe = new TH1F("wenuhOe","wenuhOe", 25,0.,1.);
  //TH1F *wenusie = new TH1F("wenusie","wenusie", 100,0.,0.07);
  TH1F *wenusie = new TH1F("wenusie","wenusie", 25,0.,0.03);

  
  //DIPHO
  TH1F *diphophoPt = new TH1F("diphophoPt","diphophoPt", 18,20.,200.);
  TH1F *diphoelePt = new TH1F("diphoelePt","diphoelePt", 18,20.,200.);
  TH1F *diphoMaxMass = new TH1F("diphoMaxMass","diphoMaxMass", 48,20.,500.);
  TH1F *diphoMee = new TH1F("diphoMee","diphoMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *diphoeIso = new TH1F("diphoeIso","diphoeIso", 50,0.,100.);
  //TH1F *diphoeIso = new TH1F("diphoeIso","diphoeIso", 50,-5.,100.);
  TH1F *diphoeIso = new TH1F("diphoeIso","diphoeIso", nbins,isobin);
  TH1F *diphohIso = new TH1F("diphohIso","diphohIso", 50,0.,100.);
  TH1F *diphotIso = new TH1F("diphotIso","diphotIso", 50,0.,100.);
  TH1F *diphohOe = new TH1F("diphohOe","diphohOe", 25,0.,1.);
  //TH1F *diphosie = new TH1F("diphosie","diphosie", 100,0.,0.07);
  TH1F *diphosie = new TH1F("diphosie","diphosie", 25,0.,0.03);




  //Ztt
  TH1F *zttphoPt = new TH1F("zttphoPt","zttphoPt", 18,20.,200.);
  TH1F *zttelePt = new TH1F("zttelePt","zttelePt", 18,20.,200.);
  TH1F *zttMaxMass = new TH1F("zttMaxMass","zttMaxMass", 48,20.,500.);
  TH1F *zttMee = new TH1F("zttMee","zttMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *ztteIso = new TH1F("ztteIso","ztteIso", 50,0.,100.);
  //  TH1F *ztteIso = new TH1F("ztteIso","ztteIso", 50,-5.,100.);
  TH1F *ztteIso = new TH1F("ztteIso","ztteIso", nbins,isobin);
  TH1F *ztthIso = new TH1F("ztthIso","ztthIso", 50,0.,100.);
  TH1F *ztttIso = new TH1F("ztttIso","ztttIso", 50,0.,100.);
  TH1F *ztthOe = new TH1F("ztthOe","ztthOe", 25,0.,1.);
  //TH1F *zttsie = new TH1F("zttsie","zttsie", 100,0.,0.07);
TH1F *zttsie = new TH1F("zttsie","zttsie", 25,0.,0.03);


  //ttbar
  TH1F *ttphoPt = new TH1F("ttphoPt","ttphoPt", 18,20.,200.);
  TH1F *ttelePt = new TH1F("ttelePt","ttelePt", 18,20.,200.);
  TH1F *ttMaxMass = new TH1F("ttMaxMass","ttMaxMass", 48,20.,500.);
  TH1F *ttMee = new TH1F("ttMee","ttMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *tteIso = new TH1F("tteIso","tteIso", 50,0.,100.);
  //  TH1F *tteIso = new TH1F("tteIso","tteIso", 50,-5.,100.);
  TH1F *tteIso = new TH1F("tteIso","tteIso", nbins,isobin);
  TH1F *tthIso = new TH1F("tthIso","tthIso", 50,0.,100.);
  TH1F *tttIso = new TH1F("tttIso","tttIso", 50,0.,100.);
  TH1F *tthOe = new TH1F("tthOe","tthOe", 25,0.,1.);
  //TH1F *ttsie = new TH1F("ttsie","ttsie", 100,0.,0.07);
TH1F *ttsie = new TH1F("ttsie","ttsie", 25,0.,0.03);

  //QCD
  TH1F *qcdphoPt = new TH1F("qcdphoPt","qcdphoPt", 18,20.,200.);
  TH1F *qcdelePt = new TH1F("qcdelePt","qcdelePt", 18,20.,200.);
  TH1F *qcdMaxMass = new TH1F("qcdMaxMass","qcdMaxMass", 48,20.,500.);
  TH1F *qcdMee = new TH1F("qcdMee","qcdMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *qcdeIso = new TH1F("qcdeIso","qcdeIso", 50,0.,100.);
  //TH1F *qcdeIso = new TH1F("qcdeIso","qcdeIso", 50,-5.,100.);
  TH1F *qcdeIso = new TH1F("qcdeIso","qcdeIso", nbins,isobin);
  TH1F *qcdhIso = new TH1F("qcdhIso","qcdhIso", 50,0.,100.);
  TH1F *qcdtIso = new TH1F("qcdtIso","qcdtIso", 50,0.,100.);
  TH1F *qcdhOe = new TH1F("qcdhOe","qcdhOe", 25,0.,1.);
  //TH1F *qcdsie = new TH1F("qcdsie","qcdsie", 100,0.,0.07);
  TH1F *qcdsie = new TH1F("qcdsie","qcdsie", 25,0.,0.03);


  //Zjets
  TH1F *zjetsphoPt = new TH1F("zjetsphoPt","zjetsphoPt", 18,20.,200.);
  TH1F *zjetselePt = new TH1F("zjetselePt","zjetselePt", 18,20.,200.);
  TH1F *zjetsMaxMass = new TH1F("zjetsMaxMass","zjetsMaxMass", 48,20.,500.);
  TH1F *zjetsMee = new TH1F("zjetsMee","zjetsMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *zjetseIso = new TH1F("zjetseIso","zjetseIso", 50,0.,100.);
  //TH1F *zjetseIso = new TH1F("zjetseIso","zjetseIso", 50,-5.,100.);
  TH1F *zjetseIso = new TH1F("zjetseIso","zjetseIso", nbins,isobin);
  TH1F *zjetshIso = new TH1F("zjetshIso","zjetshIso", 50,0.,100.);
  TH1F *zjetstIso = new TH1F("zjetstIso","zjetstIso", 50,0.,100.);
  TH1F *zjetshOe = new TH1F("zjetshOe","zjetshOe", 25,0.,1.);
  //TH1F *zjetssie = new TH1F("zjetssie","zjetssie", 100,0.,0.07);
  TH1F *zjetssie = new TH1F("zjetssie","zjetssie", 25,0.,0.03);

  //Wjets
  TH1F *wjetsphoPt = new TH1F("wjetsphoPt","wjetsphoPt", 18,20.,200.);
  TH1F *wjetselePt = new TH1F("wjetselePt","wjetselePt", 18,20.,200.);
  TH1F *wjetsMaxMass = new TH1F("wjetsMaxMass","wjetsMaxMass", 48,20.,500.);
  TH1F *wjetsMee = new TH1F("wjetsMee","wjetsMee", 48,20.,500.);

  ///isolation, ID
  //TH1F *wjetseIso = new TH1F("wjetseIso","wjetseIso", 50,0.,100.);
  //TH1F *wjetseIso = new TH1F("wjetseIso","wjetseIso", 50,-5.,100.);
  TH1F *wjetseIso = new TH1F("wjetseIso","wjetseIso", nbins,isobin);
  TH1F *wjetshIso = new TH1F("wjetshIso","wjetshIso", 50,0.,100.);
  TH1F *wjetstIso = new TH1F("wjetstIso","wjetstIso", 50,0.,100.);
  TH1F *wjetshOe = new TH1F("wjetshOe","wjetshOe", 25,0.,1.);
  //TH1F *wjetssie = new TH1F("wjetssie","wjetssie", 100,0.,0.07);
  TH1F *wjetssie = new TH1F("wjetssie","wjetssie", 25,0.,0.03);



//   In a ROOT session, you can do:
//      Root > .L TwoTrees.C
//      Root > TwoTrees t
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
  //if (fChain == 0) return;

  
  /////DATA - MC
  TFile * df = new TFile("data_mc.root");
  TTree * dt = (TTree*)df->Get("ntuple");
  Init(dt);
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //cout<<"jentry "<<jentry<<endl;
    
     if(region=="EB_EB_EB")
      flag_region = EB_EB_EB();

    if(region=="EB_EE_EB")
      flag_region = EB_EE_EB();
      
    if(region=="all")
      flag_region = 1;

    /////DATA
    
    if(flag_region){
      if(itype==0)
	{
	  //cout<<"DATA : w : "<<w<<endl;
	  dataphoPt->Fill(finalPt); 
	  dataelePt->Fill(ele1Pt); 
	  dataelePt->Fill(ele2Pt); 
	  dataMaxMass->Fill(finalMaxMass); 
	  dataMee->Fill(finalMee);
	  ////ISO, ID
	  dataeIso->Fill(phoecalIso);
	  datahIso->Fill(phohcalIso);
	  datatIso->Fill(photrkIso);
	  datahOe->Fill(phoHoE);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  datasie->Fill(phosieie);
	  //////////////////////////////////
	}
      

      double nw = w;

      if(itype==-23)  // from HEEP note
	 nw=w*(1656.14/1655.91)*(49293./153001.);
	 

      if(itype==-16)  // Ffomr HEEP note
	 nw=w*(8.65536/8.83704);

      if(itype==-17) // from HEep note
	 nw=w*(1.16685/1.21376);
       
      if(itype==-18) // from HEEP note
	 nw=w*(0.0296832/0.0303007);

      if(itype==-19) //from HEEP note
	 nw=w*(0.0041088/0.00399327);
       
       if(itype==-14)  //from ZGamma pre-approval
	 nw=w*(3048./2475.);

       //wjets
       if(itype==-22)  //from ZGamma pre-approval 
	 nw=w*(31313./27770.);

       //if(itype==-20)
       //nw=w*(33.2/13.79);
       

      ////ZEE - Zee
       //if(itype==-23 ||itype==-16 || itype==-17 || itype==-18 || itype==-19)
       if(itype==-20)
	{
	  
	  //cout<<"gluonpt "<<gluonPt<<endl;
	  //if(gluonPt<gluonPtcut)
	    {
	      //cout<<"ZEE : w :"<<w<<endl;
	      
	      zeephoPt->Fill(finalPt,nw); 
	      zeeelePt->Fill(ele1Pt,nw); 
	      zeeelePt->Fill(ele2Pt,nw); 
	      zeeMaxMass->Fill(finalMaxMass,nw); 
	      zeeMee->Fill(finalMee,nw); 
	      ////ISO, ID
	      zeeeIso->Fill(phoecalIso,nw);
	      zeehIso->Fill(phohcalIso,nw);
	      zeetIso->Fill(photrkIso,nw);
	      zeehOe->Fill(phoHoE,nw);
	      
	      if(phoecalIso<0.0)
		cout<<"zee: ecal Iso is -ve:"<<phoecalIso<<endl;

	      ////////June 10 - sieie checks - after tight cuts///////////////////
	      double et = finalPt;
	      //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	      zeesie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////////
	    }//if(gluonPt<gluonPtcut)
	}
      
      
       
      //-1:ZZ
      if(itype==-1)
	{
	  //cout<<"ZZ : w :"<<w<<endl;
	  zzphoPt->Fill(finalPt,nw); 
	  zzelePt->Fill(ele1Pt,nw); 
	  zzelePt->Fill(ele2Pt,nw); 
	  zzMaxMass->Fill(finalMaxMass,nw); 
	  zzMee->Fill(finalMee,nw);
	  ////ISO, ID
	  zzeIso->Fill(phoecalIso,nw);
	  zzhIso->Fill(phohcalIso,nw);
	  zztIso->Fill(photrkIso,nw);
	  zzhOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	    zzsie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	
	}
      
      //-2:WW
      if(itype==-2)
	{
	  wwphoPt->Fill(finalPt,nw); 
	  wwelePt->Fill(ele1Pt,nw); 
	  wwelePt->Fill(ele2Pt,nw); 
	  wwMaxMass->Fill(finalMaxMass,nw); 
	  wwMee->Fill(finalMee,nw);
	  ////ISO, ID
	  wweIso->Fill(phoecalIso,nw);
	  wwhIso->Fill(phohcalIso,nw);
	  wwtIso->Fill(photrkIso,nw);
	  wwhOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  wwsie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	    
	}

      //-3:WZ
      if(itype==-3)
	{
	  wzphoPt->Fill(finalPt,nw); 
	  wzelePt->Fill(ele1Pt,nw); 
	  wzelePt->Fill(ele2Pt,nw); 
	  wzMaxMass->Fill(finalMaxMass,nw); 
	  wzMee->Fill(finalMee,nw);
	  ////ISO, ID
	  wzeIso->Fill(phoecalIso,nw);
	  wzhIso->Fill(phohcalIso,nw);
	  wztIso->Fill(photrkIso,nw);
	  wzhOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  wzsie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	    
	}

      //-4:Wenu
      //-21: wgamma
      //if(itype==-4)
      if(itype==-21)
	{
	  wenuphoPt->Fill(finalPt,nw); 
	  wenuelePt->Fill(ele1Pt,nw); 
	  wenuelePt->Fill(ele2Pt,nw); 
	  wenuMaxMass->Fill(finalMaxMass,nw); 
	  wenuMee->Fill(finalMee,nw);
	  ////ISO, ID
	  wenueIso->Fill(phoecalIso,nw);
	  wenuhIso->Fill(phohcalIso,nw);
	  wenutIso->Fill(photrkIso,nw);
	  wenuhOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  wenusie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	}

      //-7:born25to250, -8:born250toinf, -9:box25to250, -10:bpx250toinf
      if(itype==-7 || itype==-8 || itype==-9 || itype==-10)
	{
	  diphophoPt->Fill(finalPt,nw); 
	  diphoelePt->Fill(ele1Pt,nw); 
	  diphoelePt->Fill(ele2Pt,nw); 
	  diphoMaxMass->Fill(finalMaxMass,nw); 
	  diphoMee->Fill(finalMee,nw);
	  ////ISO, ID
	  diphoeIso->Fill(phoecalIso,nw);
	  diphohIso->Fill(phohcalIso,nw);
	  diphotIso->Fill(photrkIso,nw);
	  diphohOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  diphosie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	}
      

      //Z->tautau
      //-5: ztotautau
      if(itype==-5)
	{
	  //cout<<"ZTT : w :"<<w<<endl;
	  zttphoPt->Fill(finalPt,nw); 
	  zttelePt->Fill(ele1Pt,nw); 
	  zttelePt->Fill(ele2Pt,nw); 
	  zttMaxMass->Fill(finalMaxMass,nw); 
	  zttMee->Fill(finalMee,nw); 
	  ////ISO, ID
	  ztteIso->Fill(phoecalIso,nw);
	  ztthIso->Fill(phohcalIso,nw);
	  ztttIso->Fill(photrkIso,nw);
	  ztthOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	    zttsie->Fill(phosieie,nw);
	  /////////////////////////////////////////////////////////////////
	}
      
      //ttbar
      //-6: ttbar
      if(itype==-6)
	{
	  //cout<<"TT : w :"<<nw<<endl;
	  ttphoPt->Fill(finalPt,nw); 
	  ttelePt->Fill(ele1Pt,nw); 
	  ttelePt->Fill(ele2Pt,nw); 
	  ttMaxMass->Fill(finalMaxMass,nw); 
	  ttMee->Fill(finalMee,nw); 
	  ////ISO, ID
	  tteIso->Fill(phoecalIso,nw);
	  tthIso->Fill(phohcalIso,nw);
	  tttIso->Fill(photrkIso,nw);
	  tthOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	    ttsie->Fill(phosieie,nw);
	  ////////////////////////////////////////////////////////
	}
      
      //QCD
      //-11:qcdbc20to30, -12:qcdbc30to80, -13:qcd80to170
      if(itype==-11 || itype==-12 || itype==-13)
	{
	  //cout<<"QCD : w :"<<w<<endl;
	  qcdphoPt->Fill(finalPt,nw); 
	  qcdelePt->Fill(ele1Pt,nw); 
	  qcdelePt->Fill(ele2Pt,nw); 
	  qcdMaxMass->Fill(finalMaxMass,nw); 
	  qcdMee->Fill(finalMee,nw); 
	  ////ISO, ID
	  qcdeIso->Fill(phoecalIso,nw);
	  qcdhIso->Fill(phohcalIso,nw);
	  qcdtIso->Fill(photrkIso,nw);
	  qcdhOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	    qcdsie->Fill(phosieie,nw);
	  //////////////////////////////////////////////////////////
	}
      
      
      ////ZJETS - Zjets
      if(itype==-14)
      //if(itype==-20)
	{
	  
	  if( !(TMath::Abs(phoMother)>=1&&TMath::Abs(phoMother)<=5) )
	    {
	      //cout<<"ZJETS : w :"<<w<<endl;
	      zjetsphoPt->Fill(finalPt,nw); 
	      zjetselePt->Fill(ele1Pt,nw); 
	      zjetselePt->Fill(ele2Pt,nw); 
	      zjetsMaxMass->Fill(finalMaxMass,nw); 
	      zjetsMee->Fill(finalMee,nw); 
	      ////ISO, ID
	      zjetseIso->Fill(phoecalIso,nw);
	      zjetshIso->Fill(phohcalIso,nw);
	      zjetstIso->Fill(photrkIso,nw);
	      zjetshOe->Fill(phoHoE,nw);
	      
	      ////////June 10 - sieie checks - after tight cuts///////////////////
	      double et = finalPt;
	      //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	      zjetssie->Fill(phosieie,nw);
	      
	  //////////////////////////////////////////////////////////////
	    }//if( abs(phoMother)>=1&&abs(phoMother)<=5 )
	  
	}



      ////WJETS - Wjets
      if(itype==-22)
	{
	  
	  if( !(TMath::Abs(phoMother)>=1&&TMath::Abs(phoMother)<=5) )
	    {
	      //cout<<"WJETS : w :"<<w<<endl;
	      wjetsphoPt->Fill(finalPt,nw); 
	      wjetselePt->Fill(ele1Pt,nw); 
	      wjetselePt->Fill(ele2Pt,nw); 
	      wjetsMaxMass->Fill(finalMaxMass,nw); 
	      wjetsMee->Fill(finalMee,nw); 
	      ////ISO, ID
	      wjetseIso->Fill(phoecalIso,nw);
	      wjetshIso->Fill(phohcalIso,nw);
	      wjetstIso->Fill(photrkIso,nw);
	      wjetshOe->Fill(phoHoE,nw);
	      
	      ////////June 10 - sieie checks - after tight cuts///////////////////
	      double et = finalPt;
	      //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	      wjetssie->Fill(phosieie,nw);
	      
	  //////////////////////////////////////////////////////////////
	    }//if( abs(phoMother)>=1&&abs(phoMother)<=5 )
	  
	}

      
    }//if(flag_region)
      
  }//for (Long64_t jentry=0; jentry<nentries;jentry++)

  //////////////////////OPEN SECOND FILE/////////////////////////////
  /////second file which contains information from 
  ////box25to250, box250Toinf, qcd20To30,qcd30to80 and qcd80to170. First file doesnt have that info

  TFile * df1 = new TFile("box_qcd.root");
  TTree * dt1 = (TTree*)df1->Get("ntuple");
  Init(dt1);
  
  nentries = fChain->GetEntriesFast();
  
  nbytes = 0;
  nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    double nw=w;
    
    if(itype==-9 || itype==-10)
	{
	  diphophoPt->Fill(finalPt,nw); 
	  diphoelePt->Fill(ele1Pt,nw); 
	  diphoelePt->Fill(ele2Pt,nw); 
	  diphoMaxMass->Fill(finalMaxMass,nw); 
	  diphoMee->Fill(finalMee,nw);
	  ////ISO, ID
	  diphoeIso->Fill(phoecalIso,nw);
	  diphohIso->Fill(phohcalIso,nw);
	  diphotIso->Fill(photrkIso,nw);
	  diphohOe->Fill(phoHoE,nw);
	  
	  ////////June 10 - sieie checks - after tight cuts///////////////////
	  double et = finalPt;
	  //if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	  diphosie->Fill(phosieie,nw);
	  ///////////////////////////////////////////////////////////////
	}
    
    //QCD
    //-11:qcdbc20to30, -12:qcdbc30to80, -13:qcd80to170
    if(itype==-11 || itype==-12 || itype==-13)
      {
	//cout<<"QCD : w :"<<w<<endl;
	qcdphoPt->Fill(finalPt,nw); 
	qcdelePt->Fill(ele1Pt,nw); 
	qcdelePt->Fill(ele2Pt,nw); 
	qcdMaxMass->Fill(finalMaxMass,nw); 
	qcdMee->Fill(finalMee,nw); 
	////ISO, ID
	qcdeIso->Fill(phoecalIso,nw);
	qcdhIso->Fill(phohcalIso,nw);
	qcdtIso->Fill(photrkIso,nw);
	qcdhOe->Fill(phoHoE,nw);
	
	////////June 10 - sieie checks - after tight cuts///////////////////
	double et = finalPt;
	//if(phoecalIso<4.2+0.006*et && phohcalIso<2.2+0.0025*et && phoHoE<0.05 && photrkIso<2.0+0.001*et)
	qcdsie->Fill(phosieie,nw);
	//////////////////////////////////////////////////////////
      }
    

    //cout<<"jentry "<<jentry<<endl;
  }//for (Long64_t jentry=0; jentry<nentries;jentry++)



  /*dataMee->Scale(1./dataMee->Integral());
   zeeMee->Scale(1./zeeMee->Integral());
   zzMee->Scale(1./zzMee->Integral());
   wzMee->Scale(1./wzMee->Integral());
   wwMee->Scale(1./wwMee->Integral());
   wenuMee->Scale(1./wenuMee->Integral());
   diphoMee->Scale(1./diphoMee->Integral());
   zttMee->Scale(1./zttMee->Integral());
   ttMee->Scale(1./ttMee->Integral());
   qcdMee->Scale(1./qcdMee->Integral());
   zjetsMee->Scale(1./zjetsMee->Integral());
  */

  ////some checks
  float datain  = dataMee->Integral();
   float zeein   = zeeMee->Integral();
   float zzin   = zzMee->Integral();
   float wzin   = wzMee->Integral();
   float wwin   = wwMee->Integral();
   float wenuin = wenuMee->Integral();
   float diphoin   = diphoMee->Integral();
   
   float dyttin  = zttMee->Integral();
   float ttin    = ttMee->Integral();
   float qcdin   = qcdMee->Integral();
   float zjetsin   = zjetsMee->Integral();
   float wjetsin   = wjetsMee->Integral();
   //zjetsin = 0;
   float totexp = zeein+zzin+wzin+wwin+wenuin+diphoin+dyttin+ttin+qcdin+zjetsin+wjetsin;
   cout<<"data: "<<datain<<":"<<"Zee: "<<zeein<<":"<<"ZZ: "<<zzin<<":"<<"WZ: "<<wzin<<":"<<"WW: "<<wwin<<":"<<"Wenu: "<<wenuin<<":"<<"Dipho: "<<diphoin<<":"<<"dytt: "<<dyttin<<":"<<"ttbar:"<<ttin<<":"<<"QCD:"<<qcdin<<":"<<"zjets: "<<zjetsin<<":"<<"wjets: "<<wjetsin<<":"<<"TOT EXP: "<<totexp<<endl;
  


   
   

  /*int bin = dataMee->FindBin(0.);
  cout<<"bin: "<<bin<<endl;
  
  float datain  = dataMee->Integral(bin,48);
  float zeein   = zeeMee->Integral(bin,48);
  float zzin   = zzMee->Integral(bin,48);
  float wzin   = wzMee->Integral(bin,48);
  float wwin   = wwMee->Integral(bin,48);
  float wenuin = wenuMee->Integral(bin,48);
  float diphoin   = diphoMee->Integral(bin,48);
  
  float dyttin  = zttMee->Integral(bin,48);
  float ttin    = ttMee->Integral(bin,48);
  float qcdin   = qcdMee->Integral(bin,48);
  float zjetsin   = zjetsMee->Integral(bin,48);
  float totexp = zeein+zzin+wzin+wwin+wenuin+diphoin+dyttin+ttin+qcdin+zjetsin;
  cout<<"data: "<<datain<<":"<<"Zee: "<<zeein<<":"<<"ZZ: "<<zzin<<":"<<"WZ: "<<wzin<<":"<<"WW: "<<wwin<<":"<<"Wenu: "<<wenuin<<":"<<"Dipho: "<<diphoin<<":"<<"dytt: "<<dyttin<<":"<<"ttbar:"<<ttin<<":"<<"QCD:"<<qcdin<<":"<<"zjets: "<<zjetsin<<":"<<"TOT EXP: "<<totexp<<endl;
  
  */


  /*float datain  = dataMaxMass->Integral();
   float zeein   = zeeMaxMass->Integral();
   float zzin   = zzMaxMass->Integral();
   float wzin   = wzMaxMass->Integral();
   float wwin   = wwMaxMass->Integral();
   float wenuin = wenuMaxMass->Integral();
   float diphoin   = diphoMaxMass->Integral();
   
   float dyttin  = zttMaxMass->Integral();
   float ttin    = ttMaxMass->Integral();
   float qcdin   = qcdMaxMass->Integral();
   float zjetsin   = zjetsMaxMass->Integral();
   float totexp = zeein+zzin+wzin+wwin+wenuin+diphoin+dyttin+ttin+qcdin+zjetsin;
   cout<<"data: "<<datain<<":"<<"Zee: "<<zeein<<":"<<"ZZ: "<<zzin<<":"<<"WZ: "<<wzin<<":"<<"WW: "<<wwin<<":"<<"Wenu: "<<wenuin<<":"<<"Dipho: "<<diphoin<<":"<<"dytt: "<<dyttin<<":"<<"ttbar:"<<ttin<<":"<<"QCD:"<<qcdin<<":"<<"zjets: "<<zjetsin<<":"<<"TOT EXP: "<<totexp<<endl;
  */

   // USAGE: drawCanvas(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5, char *xtitle, char *dirName,char *gifName)
  //////QCD, ttbar, ztotautau, ZZ, Zee, data
  
   //PhoPt          
   /*char *outputName[1] = {(char*)("phoPt"+region).c_str()};      
   char *xtitle[1] = {"Pt_{#gamma} (GeV/c)"};  

   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdphoPt,ttphoPt,zttphoPt,zzphoPt, zeephoPt, zjetsphoPt,dataphoPt, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdphoPt,ttphoPt,zttphoPt,zzphoPt,wwphoPt,wzphoPt,wenuphoPt,diphophoPt, zeephoPt,zjetsphoPt,dataphoPt, xtitle[0], "plots", outputName[0]);
   */

   //ElePt          
   /*char *outputName[1] = {("elePt"+region).c_str()};      
   char *xtitle[1] = {"Pt_{e} (GeV/c)"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdelePt,ttelePt,zttelePt,zzelePt, zeeelePt, zjetselePt,dataelePt, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdelePt,ttelePt,zttelePt,zzelePt,wwelePt,wzelePt,wenuelePt,diphoelePt, zeeelePt,zjetselePt,dataelePt, xtitle[0], "plots", outputName[0]);
   

   //Massegma
   char *outputName[1] = {("Megam"+region).c_str()};      
   char *xtitle[1] = {"M_{e#gamma} (GeV/c^{2})"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdMaxMass,ttMaxMass,zttMaxMass,zzMaxMass, zeeMaxMass, zjetsMaxMass,dataMaxMass, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdMaxMass,ttMaxMass,zttMaxMass,zzMaxMass,wwMaxMass,wzMaxMass,wenuMaxMass,diphoMaxMass, zeeMaxMass,zjetsMaxMass,dataMaxMass, xtitle[0], "plots", outputName[0]);
   */

   //Mee
   char *outputName[1] = {(char*)("Mee"+region).c_str()};      
   char *xtitle[1] = {"M_{ee} (GeV/c^{2})"};  
   cout<<"outputname = "<<*outputName[0]<<endl;
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdMee,ttMee,zttMee,zzMee, zeeMee, zjetsMee, dataMee, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdMee,ttMee,zttMee,zzMee,wwMee,wzMee,wenuMee,diphoMee, zeeMee,zjetsMee,dataMee, xtitle[0], "plots", outputName[0]);
   drawCanvasForUL12(qcdMee,ttMee,zttMee,zzMee,wwMee,wzMee,wenuMee,diphoMee, zeeMee,zjetsMee,wjetsMee,dataMee, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL10(qcdMee,ttMee,zttMee,zzMee,wwMee,wzMee,wenuMee,diphoMee, zeeMee,dataMee, xtitle[0], "plots", outputName[0]);
   

   
   //ISO, ID
   //ecalIso        
   /* char *outputName[1] = {(char*)("ecalIso"+region).c_str()};      
   char *xtitle[1] = {"Ecal Isolation (GeV/c)"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdeIso,tteIso,ztteIso,zzeIso, zeeeIso,zjetseIso, dataeIso, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdeIso,tteIso,ztteIso,zzeIso,wweIso,wzeIso,wenueIso,diphoeIso, zeeeIso,zjetseIso,dataeIso, xtitle[0], "plots", outputName[0]);
   drawCanvasForUL12(qcdeIso,tteIso,ztteIso,zzeIso,wweIso,wzeIso,wenueIso,diphoeIso, zeeeIso,zjetseIso,wjetseIso,dataeIso, xtitle[0], "plots", outputName[0]);
   */
   
   //hcalIso        
   /*char *outputName[1] = {("hcalIso"+region).c_str()};      
   char *xtitle[1] = {"Hcal Isolation (GeV/c)"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdhIso,tthIso,ztthIso,zzhIso, zeehIso,zjetshIso,datahIso, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdhIso,tthIso,ztthIso,zzhIso,wwhIso,wzhIso,wenuhIso,diphohIso, zeehIso,zjetshIso,datahIso, xtitle[0], "plots", outputName[0]);

   //trkIso        
   char *outputName[1] = {("trkIso"+region).c_str()};      
   char *xtitle[1] = {"Track Isolation (GeV/c)"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdtIso,tttIso,ztttIso,zztIso, zeetIso,zjetstIso,datatIso, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdtIso,tttIso,ztttIso,zztIso,wwtIso,wztIso,wenutIso,diphotIso, zeetIso,zjetstIso,datatIso, xtitle[0], "plots", outputName[0]);

   //H/E       
   char *outputName[1] = {("HoverE"+region).c_str()};      
   char *xtitle[1] = {"H/E"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdhOe,tthOe,ztthOe,zzhOe, zeehOe,zjetshOe,datahOe, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdhOe,tthOe,ztthOe,zzhOe,wwhOe,wzhOe,wenuhOe,diphohOe, zeehOe,zjetshOe,datahOe, xtitle[0], "plots", outputName[0]);
   */
   //sie        
   /*char *outputName[1] = {(char*)("sie"+region).c_str()};      
   char *xtitle[1] = {"#sigma_{i#etai#eta}"};  
   //cout<<"drawcanvas called"<<endl;   
   //drawCanvasForUL7(qcdsie,ttsie,zttsie,zzsie, zeesie, zjetssie,datasie, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL11(qcdsie,ttsie,zttsie,zzsie,wwsie,wzsie,wenusie,diphosie, zeesie,zjetssie,datasie, xtitle[0], "plots", outputName[0]);
   //drawCanvasForUL10(qcdsie,ttsie,zttsie,zzsie,wwsie,wzsie,wenusie,diphosie, zeesie,datasie, xtitle[0], "plots", outputName[0]);
   drawCanvasForUL12(qcdsie,ttsie,zttsie,zzsie,wwsie,wzsie,wenusie,diphosie, zeesie,zjetssie,wjetssie,datasie, xtitle[0], "plots", outputName[0]);
*/
   

   //////CUMULATIVE HISTO - drawCanvas2 is called from there
   //phopt 
   /*TH1F *clonehisto = (TH1F*)zeephoPt->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdphoPt);
   addhist.push_back(ttphoPt);
   addhist.push_back(zttphoPt);
   addhist.push_back(zzphoPt);
   addhist.push_back(wzphoPt);
   addhist.push_back(wwphoPt);
   addhist.push_back(wenuphoPt);
   addhist.push_back(diphophoPt);
   addhist.push_back(zeephoPt);
   addhist.push_back(zjetsphoPt);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("phoPt_cumul"+region).c_str()};      
   char *xtitle[1] = {"Pt_{#gamma} (GeV/c)"};  
   //cout<<" called"<<endl;   
   //cumulativeHistoLR(tothist,dataphoPt, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoRL(tothist,dataphoPt, xtitle[0], "plots", outputName[0]);
   


   
   //Mee
   TH1F *clonehisto = (TH1F*)zeeMee->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdMee);
   addhist.push_back(ttMee);
   addhist.push_back(zttMee);
   addhist.push_back(zzMee);
   addhist.push_back(wzMee);
   addhist.push_back(wwMee);
   addhist.push_back(wenuMee);
   addhist.push_back(diphoMee);

   addhist.push_back(zeeMee);
   addhist.push_back(zjetsMee);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("Mee_cumul"+region).c_str()};      
   char *xtitle[1] = {"M_{ee} (GeV/c)"};  
   //cout<<" called"<<endl;   
   //cumulativeHistoLR(tothist,dataMee, xtitle[0], "plots", outputName[0]);
   double integ = tothist->Integral();
   //tothist->Scale(1./integ);
   //integ = dataMee->Integral();
   //dataMee->Scale(1./integ);
   //cumulativeHistoRL(tothist,dataMee, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoLR(tothist,dataMee, xtitle[0], "plots", outputName[0]);

   //MaxMass
   TH1F *clonehisto = (TH1F*)zeeMaxMass->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdMaxMass);
   addhist.push_back(ttMaxMass);
   addhist.push_back(zttMaxMass);
   addhist.push_back(zzMaxMass);
   addhist.push_back(wzMaxMass);
   addhist.push_back(wwMaxMass);
   addhist.push_back(wenuMaxMass);
   addhist.push_back(diphoMaxMass);

   addhist.push_back(zeeMaxMass);
   addhist.push_back(zjetsMaxMass);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("MaxMass_cumul"+region).c_str()};      
   char *xtitle[1] = {"M_{e#gamma} (GeV/c)"};  
   //cout<<" called"<<endl;   
   //cumulativeHistoLR(tothist,dataMaxMass, xtitle[0], "plots", outputName[0]);
   //double integ = tothist->Integral();
   //tothist->Scale(1./integ);
   //integ = dataMaxMass->Integral();
   //dataMaxMass->Scale(1./integ);
   //cumulativeHistoRL(tothist,dataMaxMass, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoLR(tothist,dataMaxMass, xtitle[0], "plots", outputName[0]);
   
   
   */
   
   
   ///sieie
   /*TH1F *clonehisto = (TH1F*)zeesie->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdsie);
   addhist.push_back(ttsie);
   addhist.push_back(zttsie);
   addhist.push_back(zzsie);
   addhist.push_back(wzsie);
   addhist.push_back(wwsie);
   addhist.push_back(wenusie);
   addhist.push_back(diphosie);
   
   addhist.push_back(zeesie);
   addhist.push_back(zjetssie);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("sie_cumul"+region).c_str()};      
   char *xtitle[1] = {"Pt_{#gamma} (GeV/c)"};  
   //cout<<" called"<<endl;   
   cumulativeHistoLR(tothist,datasie, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoRL(tothist,datasie, xtitle[0], "plots", outputName[0]);
   */



   //ISO, ID
   //ecal Iso
   /*TH1F *clonehisto = (TH1F*)zeeeIso->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdeIso);
   addhist.push_back(tteIso);
   addhist.push_back(ztteIso);
   addhist.push_back(zzeIso);
   addhist.push_back(wzeIso);
   addhist.push_back(wweIso);
   addhist.push_back(wenueIso);
   addhist.push_back(diphoeIso);
   addhist.push_back(zeeeIso);
   addhist.push_back(zjetseIso);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("eIso_cumul"+region).c_str()};      
   char *xtitle[1] = {"Ecal Isolation(#gamma) (GeV/c)"};  
   //cout<<" called"<<endl;   
   //cumulativeHistoLR(tothist,dataeIso, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoRL(tothist,dataeIso, xtitle[0], "plots", outputName[0]);

   TH1F *clonehisto = (TH1F*)zeeeIso->Clone();
   
   ////H/E
   vector<TH1F*> addhist;
   addhist.push_back(qcdhOe);
   addhist.push_back(tthOe);
   addhist.push_back(ztthOe);
   addhist.push_back(zzhOe);
   addhist.push_back(wzhOe);
   addhist.push_back(wwhOe);
   addhist.push_back(wenuhOe);
   addhist.push_back(diphohOe);
   addhist.push_back(zeehOe);
   addhist.push_back(zjetshOe);

   TH1F *tothist = addHisto(addhist);
   
   char *outputName[1] = {("hOe_cumul"+region).c_str()};      
   char *xtitle[1] = {"H/E(#gamma) (GeV/c)"};  
   //cout<<" called"<<endl;   
   //cumulativeHistoLR(tothist,datahOe, xtitle[0], "plots", outputName[0]);
   //cumulativeHistoRL(tothist,datahOe, xtitle[0], "plots", outputName[0]);
   

   */


   //do ratio plots also
   //photon pt
   /*char *outputName[1] = {("phoPt_ratio"+region).c_str()};      
   char *xtitle[1] = {"Pt_{#gamma} (GeV/c)"};  
   //cout<<" called"<<endl;   
   drawCanvasForUL6Ratio(qcdphoPt,ttphoPt,zttphoPt,zzphoPt, zeephoPt, dataphoPt, xtitle[0], "plots", outputName[0]);
   */
   

   //////STAT TEST
   ///ISO_ID vars
   
   //1. sieie
   /*TH1F *clonehisto = (TH1F*)zeesie->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdsie);
   addhist.push_back(ttsie);
   addhist.push_back(zttsie);
   addhist.push_back(zzsie);
   addhist.push_back(wzsie);
   addhist.push_back(wwsie);
   addhist.push_back(wenusie);
   addhist.push_back(diphosie);
   
   addhist.push_back(zeesie);

   TH1F *tothist = addHisto(addhist);
   
   ///USAGE: tothistexp,hdata
   kstest(tothist,datasie);;
   */


   //2. tIso
   /*TH1F *clonehisto = (TH1F*)zeetIso->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdtIso);
   addhist.push_back(tttIso);
   addhist.push_back(ztttIso);
   addhist.push_back(zztIso);
   addhist.push_back(wztIso);
   addhist.push_back(wwtIso);
   addhist.push_back(wenutIso);
   addhist.push_back(diphotIso);
   
   addhist.push_back(zeetIso);

   TH1F *tothist = addHisto(addhist);
   
   ///USAGE: tothistexp,hdata
   kstest(tothist,datatIso);;
   */

   //3. eIso
   /*TH1F *clonehisto = (TH1F*)zeeeIso->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdeIso);
   addhist.push_back(tteIso);
   addhist.push_back(ztteIso);
   addhist.push_back(zzeIso);
   addhist.push_back(wzeIso);
   addhist.push_back(wweIso);
   addhist.push_back(wenueIso);
   addhist.push_back(diphoeIso);
   
   addhist.push_back(zeeeIso);

   TH1F *tothist = addHisto(addhist);
   
   ///USAGE: tothistexp,hdata
   kstest(tothist,dataeIso);;
   */
   
   //4. hIso
   /*TH1F *clonehisto = (TH1F*)zeehIso->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdhIso);
   addhist.push_back(tthIso);
   addhist.push_back(ztthIso);
   addhist.push_back(zzhIso);
   addhist.push_back(wzhIso);
   addhist.push_back(wwhIso);
   addhist.push_back(wenuhIso);
   addhist.push_back(diphohIso);

   addhist.push_back(zeehIso);

   TH1F *tothist = addHisto(addhist);
   
   ///USAGE: tothistexp,hdata
   kstest(tothist,datahIso);;
   */
   
   
   //5. H/E
   TH1F *clonehisto = (TH1F*)zeehOe->Clone();
   
   vector<TH1F*> addhist;
   addhist.push_back(qcdhOe);
   addhist.push_back(tthOe);
   addhist.push_back(ztthOe);
   addhist.push_back(zzhOe);
   addhist.push_back(wzhOe);
   addhist.push_back(wwhOe);
   addhist.push_back(wenuhOe);
   addhist.push_back(diphohOe);

   addhist.push_back(zeehOe);

   TH1F *tothist = addHisto(addhist);
   
   ///USAGE: tothistexp,hdata
   kstest(tothist,datahOe);
}
