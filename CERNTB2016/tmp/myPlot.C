#define myPlot_cxx
#include "myPlot.h"

#include "headers.h"

bool myPlot::spikeclean(int scind)
{
  if( sce2e9[scind]<0.95 )
    return true;
  else 
    return false;
}




void myPlot::Loop(string fileToOpen, string fileToWrite)
{
  using namespace std;
  TStopwatch t1 ;
  t1.Start() ;


  //TFile *file   = new TFile("histo.root", "RECREATE");
  TFile *file   = new TFile(fileToWrite.c_str(), "RECREATE");
  TTree *tree   = new TTree("ntuple","a tree with histograms");
  setbranch(tree);

  TFile* f;
  TTree* t;

  int fileindex;
  fileindex = 0;

  ///////before running, verify the values/////////////
  double weight = 1.0;
  int ftype = 1;
  int reqNele = 2;
  int reqNpho = 1;
  double energycorr1 = 0.0;
  double energycorr2 = 0.0;
  double ptcut = 0.0;
  ////////////////////////////////////////////////////////

  //start opening the files                    
  ifstream infile;
  infile.open("inputFilenames.list", ifstream::in );
  char infilename[200];
  while(!infile.eof()){
    infile >> infilename;
    cout<<"infilename = "<<infilename<<endl;
    if(strncmp(infilename,"#",1)==0)
      {
	//fileindex++;                                                                                                                                             
	continue;
      }
    TChain *mychain;
    string str(infilename);
    string key(".");
    size_t found;
    found           = str.rfind(key);
    string filename = str.substr(0,found);
    cout<<"filename = "<<filename<<endl;
    char file_num[100];
    sprintf(file_num, "%d", fileindex);
    std::string i_file = file_num;
    cout<<"i_file = "<<i_file<<endl;

    //f = new TFile(infilename,"READ");
    //t = (TTree*)f->Get("myEvent");
    //Init(t);

    if( filename == "data" )
      {
	ifstream datafile;
	//datafile.open("dataFilenames.list", ifstream::in );
	datafile.open(fileToOpen.c_str(), ifstream::in );
	char datafilename[200];
	mychain = (TChain *)new TChain("myEvent");
	while(!datafile.eof()){
	  datafile >> datafilename;
	  cout<<"datafilename = "<<datafilename<<endl;
	  if(strncmp(datafilename,"#",1)==0)
	    {
	      //fileindex++; 
	      continue;
	    }
	  mychain->Add(datafilename);
	}//while(!datafile.eof())    
	Init(mychain);
      }//if( filename == "data.root" )


    if( filename == "estar" )
      {
	ifstream estarfile;
	//estarfile.open("estarFilenames.list", ifstream::in );
	estarfile.open(fileToOpen.c_str(), ifstream::in );
	char estarfilename[200];
	mychain = (TChain *)new TChain("myEvent");
	while(!estarfile.eof()){
	  estarfile >> estarfilename;
	  cout<<"estarfilename = "<<estarfilename<<endl;
	  if(strncmp(estarfilename,"#",1)==0)
	    {
	      //fileindex++; 
	      continue;
	    }
	  mychain->Add(estarfilename);
	}//while(!estarfile.eof())    
	Init(mychain);
      }//if( filename == "estar.root" )

    if( filename == "z" )
      {
	ifstream zfile;
	//zfile.open("zFilenames.list", ifstream::in );
	zfile.open(fileToOpen.c_str(), ifstream::in );
	char zfilename[200];
	mychain = (TChain *)new TChain("myEvent");
	while(!zfile.eof()){
	  zfile >> zfilename;
	  cout<<"zfilename = "<<zfilename<<endl;
	  if(strncmp(zfilename,"#",1)==0)
	    {
	      //fileindex++; 
	      continue;
	    }
	  mychain->Add(zfilename);
	}//while(!zfile.eof())    
	Init(mychain);
      }//if( filename == "z.root" )

    if( filename == "bkg" )
      {
	ifstream bkgfile;
	//bkgfile.open("bkgFilenames.list", ifstream::in );
	bkgfile.open(fileToOpen.c_str(), ifstream::in );
	char bkgfilename[200];
	mychain = (TChain *)new TChain("myEvent");
	while(!bkgfile.eof()){
	  bkgfile >> bkgfilename;
	  cout<<"bkgfilename = "<<bkgfilename<<endl;
	  if(strncmp(bkgfilename,"#",1)==0)
	    {
	      //fileindex++; 
	      continue;
	    }
	  mychain->Add(bkgfilename);
	}//while(!bkgfile.eof())    
	Init(mychain);
      }//if( filename == "bkg.root" )

    
    if (fChain == 0) return;
    
    
    
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t mynentries = fChain->GetEntries();
    cout<<"mynentries = "<<mynentries<<endl;
    //nentries = 60881;
    Long64_t nbytes = 0, nb = 0;
    
    
    int ntot        = 0;
    int ev_passjson = 0;
    int ev_goodvtx = 0;
    int heepev = 0;
    double heepevwt = 0;
    int phoidev = 0;
    double phoidevwt = 0;
    int invev   = 0;
    double invevwt   = 0;
    int hltev   = 0;
    double hltevwt   = 0;
    int m150    = 0;
    double m150wt    = 0;
    int m350    = 0;
    double m350wt    = 0;
    int m550    = 0;
    double m550wt    = 0;

    TLorentzVector *sc_p4;
    //nentries = 94801;

    ////for conversion --- using conversion tools                               
    int nele1ConvTool=0, nele2ConvTool=0;
    int nanyeleConvTool=0, nbotheleConvTool=0;
    int nConvel=0;
    
    ///////////////fChain->SetBranchStatus("hlt*",0);
    fChain->SetBranchStatus("*",0);
    /*fChain->SetBranchStatus("run",1);
    fChain->SetBranchStatus("event",1);
    fChain->SetBranchStatus("lumis",1);
    fChain->SetBranchStatus("bx",1);
    fChain->SetBranchStatus("gentrk*",1);
    fChain->SetBranchStatus("patel*",1);
    fChain->SetBranchStatus("samspateleheepid*",1);
    fChain->SetBranchStatus("patpho*",1);
    fChain->SetBranchStatus("sc*",1);
    fChain->SetBranchStatus("hlt*",1);
    fChain->SetBranchStatus("vtx*",1);
    */

    ///for setting status and addbranchtoCache
    string brnames[] = {"run","event","lumis","bx*","gentrk*","patel*","samspateleheepid*","patpho*","sc*","hlt*","vtx*","rho","pvi"};
    std::vector<std::string> brname(brnames, brnames + sizeof(brnames) / sizeof(string));
    
    branchStatus1(brname);
    //std::cout<<"branchstatus of scsize "<<fChain->GetBranchStatus("scsize")<<endl;

    ///for fast processing
    addbranchtoCache(brname);


    if(filename!="data")
      {
	fChain->SetBranchStatus("pileup_*",1);
	fChain->SetBranchStatus("gen*",1);
	//cache --- since space is already allocated in                 
	fChain->AddBranchToCache("pileup_*");
        fChain->AddBranchToCache("gen*");
	fChain->AddBranchToCache("pvi*");

      }

    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      ////////////////////////////FOR FAST PROCESSING/////////////////////////////////////////////////
      ///1. require here to have more than 1 elec else continue

      b_patelesize->GetEntry(ientry);
      if(patelesize<reqNele) 
	{
	  continue;
	}

      b_patphosize->GetEntry(ientry);
      if(patphosize<reqNpho)
        {
          continue;
        }
      


      ///2. now require it to have passed json
      b_run->GetEntry(ientry);
      b_lumis->GetEntry(ientry);
      b_event->GetEntry(ientry);
      
      if( filename == "data" ){
	//b_run->GetEntry(ientry);
	//b_lumis->GetEntry(ientry);
	if(!lumiRunSelection(jentry) )
	  {
	    continue;
	  }
      }//if( filename == "data" )

      if(filename == "bkg" || filename == "estar")
	{
	  b_pileup_nvtx->GetEntry(ientry);
	  b_pileup_bunchXing->GetEntry(ientry);
	  b_pvi->GetEntry(ientry);
	}
      //////////////////////////////////////NOW MOVE TO THE NEXT STEPS IF PASSED/////////////////////////

      b_rho->GetEntry(ientry);
      
      //cout<<"rho:"<<rho<<endl;
      //cout<<"====NEW EVENT====="<<endl;
      //cout<<"entered into the loop"<<endl;
      //nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<"bytes:"<<nbytes<<endl;
      
      
      if (Cut(ientry) < 0) continue;
      //cout<<"jentry = "<<jentry<<endl;
      //cout<<"RUN = "<<run<<endl;
      //cout<<"event = "<<event<<endl;
      
      //cout<<"RUN: event :  "<<run<<" : "<<event<<endl;
      //weight
      if(filename == "data")
	{
	  weight = 1.0;
	  ftype  = 0;
	}

      string fname(fChain->GetCurrentFile()->GetName());
      //cout<<"fname = "<<fname<<endl;
      
      //double lumi = 643.654659608;
      //double lumi = 406.123521242;
      //double lumi = 521.005103417;
      double lumi=getLumi();
      double scale = 1; 
      //double scale = 0.966; //normalizing to Z peak in data
      
      //cout<<"pile up:"<<pileup_nvtx<<endl;


      if(filename!="data")
	{
	  ///call reading of gen br here
	  readGenbr(ientry);
	  
	  setweight(fname, lumi, scale, weight, ftype);
	}//if(filename!="data")
      
      //HERE STARTS MY SELECTION 
      //1. for data: JSON -----> to be put after what Elizabeth says
       bool goodlumi = true;
       if(filename == "data")
	 goodlumi = lumiRunSelection(jentry);
       
       if( goodlumi )
	 {
	   ev_passjson++ ; 
	   
	   int noscrape = 1;
	   
	   //if(filename == "data")
	   readTrkVtxbr(ientry);
	   
	   noscrape = noscraping();
	   
	   
	   int foundzgam = 1;  // true for every file except for zee for which we call a func and make a chk
	   
	   //if( filename == "z" )
	   //foundzgam = selectzgamma();
	   
	   if( foundzgam )
	     ntot++;
	   
	   ///cutshat

	   bool shatcut = true;

	   string fname(fChain->GetCurrentFile()->GetName());
	   int file;

	   file = fname.find("/mydyee20/",0);
	   if( file!=string::npos)
	     {
	       //cout<<"file name = "<<fname<<"so shat upper cut = 120"<<endl;
	       shatcut = cutshat(120);
	     }



	   file = fname.find("/dyee120/",0);
	   if( file!=string::npos)
	     {
	       //cout<<"file name = "<<fname<<"so shat upper cut = 200"<<endl;                                                                                      
	       shatcut = cutshat(200);
	     }

	   file = fname.find("/mydyee200/",0);
	   if( file!=string::npos)
	     {
	       //cout<<"file name = "<<fname<<"so shat upper cut = 500"<<endl;                                                                                      
	       shatcut = cutshat(500);
	     }

	   file = fname.find("/dyee500/",0);
	   if( file!=string::npos)
	     {
	       //cout<<"file name = "<<fname<<"so shat upper cut = 800"<<endl;                                                                                      
	       shatcut = cutshat(800);
	     }

	   ////call gen here.
	   
	  if(goodvertex() && noscrape && foundzgam && shatcut)
	    {
	      ev_goodvtx++;
	      //2. 2 HEEP electrons
	      int ele_select = 0;
	      int ele_index[100]; //array to store index of those electrons which pass the HEEPID cut
	      //cout<<"elesize = "<<patelesize<<endl;
	      // cout<<"sc size = "<<scsize<<endl;
	      double eleeffFac[100];


	      //read ele and sc br here
	      readElebr(ientry);
	      //readPhobr(ientry);
	      //readHLTbr(ientry);
	      readScbr(ientry);
	      
	      /////Apply scale energy correction   to MC only                 
	      /*if(filename !="data" )
		{
		  
		  
		  for(int iele=0; iele<patelesize; iele++)
		    {
		      
		      TLorentzVector *ele = (TLorentzVector*)patelp4->At(iele);
		      double energy = ele->Energy();
		      double px = ele->Px();
		      
		      //cout<<"=======BEfore correction======="<<endl;      
		      //cout<<"energy = "<<energy<<endl;                          
		      //cout<<"px     = "<<px<<endl;                                 
		      
		      
		      scaleE(ele,energycorr1,energycorr2);
		    }
		  
		  
		  for(int ipho=0; ipho<patphosize; ipho++)
		    {
		      TLorentzVector *pho = (TLorentzVector*)patphop4->At(ipho);
		      scaleE(pho,energycorr1,energycorr2);
		    }
		  
		}//if(filename !=data)
	      */
	      //cout<<""<<endl;

	      /*for(int isc=0;isc<scsize;isc++)
		{
		  TLorentzVector *sp4 = (TLorentzVector*)scp4->At(isc);
		  double sceta = sp4->Eta();
		  cout<<"isc : sceta = "<<isc<<" : "<<sceta<<endl;
		}
	      */
	      

	      for(int iele=0; iele<patelesize; iele++)
		{
		  //cout<<"iele : pt "<<iele<<endl;
		  int scind = patelescind[iele];
		  //cout<<"scind    = "<<scind<<endl;
		  //cout<<"scsize     = "<<scsize<<endl;
				  
		  int flag_scind = 1;
		  double eleeta;
		  
		  
		  //if(scind<sc_n && !(scind<0))
		  //{
		  //  TLorentzVector *scp4 = (TLorentzVector*)sc_p4->At(scind);
		  //  eleeta = scp4->Eta();
		  //}

		  //if(filename == "data")
		  //{
		  flag_scind = scind<scsize && !(scind<0);
		  if(flag_scind)
		    {
		      TLorentzVector *sp4 = (TLorentzVector*)scp4->At(scind);
		      eleeta = sp4->Eta();
		      
		      //TLorentzVector *elsc = (TLorentzVector*)patelsc->At(iele);
		      //eleeta = elsc->Eta();
		      //cout<<"eleeta : patelesceta : "<<eleeta<<" " << patelesceta<<endl;
		    }
		  
		  
		  //////changed defn of elp4
		  //TLorentzVector *elp4 = (TLorentzVector*)patelp4->At(iele);
		  TLorentzVector *elp4 = new TLorentzVector();
		  elp4 = setPtEtaPhiE(iele);

		  //cout<<"iele : pt : flag_scind : "<<iele<<" : "<<elp4->Pt()<<" : "<<flag_scind<<endl;
		  
		  if( flag_scind )
		    {
		      int nospike = 1;
		      
		      if( fabs(eleeta)  < 1.442 )
			nospike = spikeclean(scind);
			
		      
		      bool isConv = 0;
		      
		      ConvTool(iele, isConv);
		      
		      //cout<<"nospike : isConv: "<<nospike<<" : "<<!isConv<<endl;
		      if( nospike && !isConv ) 
			{
                          //if( samspateleheepid[iele] && elp4->Pt() >= ptcut )
			  if( eleID(iele) && elp4->Pt() >= ptcut )
			    {
			      ele_index[ele_select] = iele;
			      ele_select++;
			    }
			}//if( nospike )
		    }//if(reco_ind<sc_n)
		}//for(int iele=0; iele<pat_el_n; iele++)
	      //cout<<"out of electron loop"<<endl;
	    
	      bool elemcmatch = true;
              if(filename!="data" && filename!="estar" && ele_select>=2)
                {
                  int ele1 = ele_index[0];
                  int ele2 = ele_index[1];
		  //chnaged defn - 28th sept
                  //TLorentzVector *el1 = (TLorentzVector*)patelp4->At(ele1);
                  //TLorentzVector *el2 = (TLorentzVector*)patelp4->At(ele2);

		  TLorentzVector *el1 = new TLorentzVector();
		  el1 = setPtEtaPhiE(ele_index[0]);
		  
		  TLorentzVector *el2 = new TLorentzVector();
		  el2 = setPtEtaPhiE(ele_index[1]);
		  
                  elemcmatch = selectee(el1, el2);
		  //cout<<"elematch = "<<elemcmatch<<endl;
                }
	      
	      if(elemcmatch){
		if(ele_select>=2 ) {
		  heepev++; 
		  
		  //cout<<"found 2 HEEPs"<<endl;
		  int ele1 = ele_index[0];
		  int ele2 = ele_index[1];
		  int scind1 = patelescind[ele1];
		  int scind2 = patelescind[ele2];
		  
		  TLorentzVector *elscp4_1 = (TLorentzVector*)scp4->At(scind1);
		  TLorentzVector *elscp4_2 = (TLorentzVector*)scp4->At(scind2);
		  
		  double eleeta1 = elscp4_1->Eta();
		  double eleeta2 = elscp4_2->Eta();
		  
		  if( fabs(eleeta1)<1.442 )
		    eleeffFac[0] = 0.978;
		  
		  if( fabs(eleeta1)<2.5 && fabs(eleeta1)>1.56 )
		    eleeffFac[0] = 0.994;
		  
		  
		  if( fabs(eleeta2)<1.442 )
		    eleeffFac[1] = 0.978;
		  
		  if( fabs(eleeta2)<2.5 && fabs(eleeta2)>1.56 )
		    eleeffFac[1] = 0.994;
		  
		  heepevwt = eleeffFac[0]*eleeffFac[1]+heepevwt;

		  TLorentzVector *el1 = (TLorentzVector*)patelp4->At(ele1);
                  TLorentzVector *el2 = (TLorentzVector*)patelp4->At(ele2);
		  //		  cout<<"pt of el1 : el2 : "<<el1->Pt()<<" : "<<el2->Pt()<<endl;
		}
		
		
		/////Conversion rejection tools - to chk how many times the selected electrons are from conversion - as a xcheck to above (MC match)
                
		bool ele1ConvTool = false, ele2ConvTool = false;
                if(ele_select>=2)
                  {
                    int ele1 = ele_index[0];
                    int ele2 = ele_index[1];
                    TLorentzVector *el1 = (TLorentzVector*)patelp4->At(ele1);  
                    TLorentzVector *el2 = (TLorentzVector*)patelp4->At(ele2);       
                    ConvTool(ele1, ele1ConvTool);
                    ConvTool(ele2, ele2ConvTool);

                    if(ele1ConvTool)
                      {
			nele1ConvTool++;
			//cout<<"CONVERSION ELECTRON FOUND: run:lumis:event: ele1Pt"<<run<<":"<<lumis<<":"<<event<<":"<<el1->Pt();
			nConvel++;
		      }
                    //second electron             
                    if(ele2ConvTool)
		      {
			nele2ConvTool++;
			//cout<<"CONVERSION ELECTRON FOUND: run:lumis:event: ele2Pt"<<run<<":"<<lumis<<":"<<event<<":"<<el2->Pt();
			nConvel++;
		      }

                    //any electron              
                    if(ele1ConvTool||ele2ConvTool)
		      {
			nanyeleConvTool++;
		      }
                    //both   
                    if(ele1ConvTool && ele2ConvTool)
		      nbotheleConvTool++;
		    
                  }//if(ele_select>=2)

		//tight photon selection
		int pho_select=0;
		int phoindex[100]; //to store indices of photons passing tight photonID 
		if(ele_select>=2)
		  {
		    ///read photon br here:
		    readPhobr(ientry);

		    //cout<<""<<endl;
		    //cout<<"ievent:"<<event<<endl;

		    //cout<<"patphosize = "<<patphosize<<endl;
		    for(int ipho=0; ipho<patphosize; ipho++)
		      {
			//cout<<"ipho = "<<ipho<<endl;
			TVector3 * phop3;
			phop3 = (TVector3 *) patphocalopos->At(ipho);
			//phop3 = (TVector3 *) pho_calopos->At(ipho);
			double phoeta=phop3->Eta();
			//cout<<"phoeta = "<<phoeta<<endl;
			//should not match with electron SC
			//cout<<"recoind1 = "<<recoind1<<endl;
			
			////defn. now changed - sept 28
			//TLorentzVector *ele1p4 = (TLorentzVector*)patelp4->At(ele_index[0]);
			//TLorentzVector *ele2p4 = (TLorentzVector*)patelp4->At(ele_index[1]);
			
			TLorentzVector *ele1p4 = new TLorentzVector();
			ele1p4 = setPtEtaPhiE(ele_index[0]);
			
			TLorentzVector *ele2p4 = new TLorentzVector();
			ele2p4 = setPtEtaPhiE(ele_index[1]);

			TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(ipho);
			
			//cout<<"myPlot.C before scaling:ipho : pt : "<<ipho<<" : "<<pho_p4->Pt()<<endl;
			////////scale photons energy here if it is MC
			/*if(filename=="bkg")
			  {
			    pho_p4 = scalephoE(ipho);
			  }
			*/

			//cout<<"myPlot.C after scaling:ipho : pt : "<<ipho<<" : "<<pho_p4->Pt()<<endl;

			double dR1 = pho_p4->DeltaR(*ele1p4);
			double dR2 = pho_p4->DeltaR(*ele2p4);
			
			int nomatch;
			nomatch = patelescind[ele_index[0]]!=patphoscind[ipho] && dR1>0.5 && patelescind[ele_index[1]]!=patphoscind[ipho] && dR2>0.5;
			//nomatch = dR1>0.5 &&dR2>0.5;
			
			//cout<<"patelescind[ele_index[0]] = "<<patelescind[ele_index[0]]<<endl;
			//cout<<"patelescind[ele_index[1]] = "<<patelescind[ele_index[1]]<<endl;
			//cout<<"patphoscind[ipho]  = "<<patphoscind[ipho] <<endl;
			
			//if( (fabs(phoeta)<1.4442) || (fabs(phoeta)>1.566 && fabs(phoeta)<2.4) )
			
			//cout<<"nomatch = "<<nomatch<<endl;

			if(nomatch && pho_p4->Pt()> 20.)
			  {
			    //cout<<"PASSED NOMATCH"<<endl;
			    if( fabs(phoeta)<1.4442 )  // Elizabeth is not using EE photons for the time being
			      {
				bool nospike=true;
				
				int scind_pho = patphoscind[ipho];
				
				if( (fabs(phoeta)<1.4442) )
				  nospike = spikeclean(scind_pho);
				  
				//cout<<"nospike : "<<nospike<<endl;
				if( nospike )
				  {
				    //cout<<"passed nospike"<<endl;
				    //cout<<"phoID : "<<phoID(ipho)<<endl;
				    if( phoID(ipho) )
				      {
					//cout<<"passed PHOID"<<endl;
					phoindex[pho_select] = ipho;
					pho_select++;
				      }
				  }//if( spikeclean(phorecoind) )
			      }//if( (fabs(phoeta)<1.4442) || (fabs(phoeta)>1.566 && fabs(phoeta)<2.4) )
			  }//if(nomatch)
		      }//for(int ipho=0; ipho<pat_pho_n; ipho++)
		  }//if(ele_select>=2)
		
		bool phomcmatch = true;
		if(filename!="data" && filename!="estar" && pho_select>=1)
		  {
		    int phorecoind = phoindex[0];
		    TLorentzVector *phop4  = (TLorentzVector*)patphop4->At(phorecoind);
		    phomcmatch = selectgam(phop4);
		    //phomcmatch = 1;                                   
		    //cout<<"phomcmatch : "<<phomcmatch<<endl;
		  }
		
		if(phomcmatch)
		  {
		    if(pho_select>=1)
		      {
			phoidev++;
			phoidevwt = eleeffFac[0]*eleeffFac[1]*0.98 + phoidevwt;
		      }
		    
		    int invflag = 0;
		    double invm_e1e2;
		    double invm_e1pho;
		    double invm_e2pho;
		    double invm_eegam;
		    
		    if(pho_select>=1 && ele_select>=2)
		      {
			//Z window rejection
			TVector3 *phop3 = (TVector3 *) patphocalopos->At(phoindex[0]);
			
			
			/////////now use the Et = SC energy*sin(Gsftrk theta) - 28th sept
			TLorentzVector *ele1p4 = new TLorentzVector();
			ele1p4 = setPtEtaPhiE(ele_index[0]);
			
			TLorentzVector *ele2p4 = new TLorentzVector();
			ele2p4 = setPtEtaPhiE(ele_index[1]);

                        TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(phoindex[0]);
			
			
			invm_e1e2     = ( (*ele1p4)+(*ele2p4) ).M();
			invm_e1pho    = ( (*ele1p4)+(*pho_p4) ).M();
			invm_e2pho    = ( (*ele2p4)+(*pho_p4) ).M();
			invm_eegam    = ( (*ele1p4)+(*ele2p4)+(*pho_p4) ).M();

			

			//invflag = !(invm_e1e2>81.19 && invm_e1e2<101.19) && !(invm_e1pho>81.19 && invm_e1pho<101.19) && !(invm_e2pho>81.19 && invm_e2pho<101.19) && !(invm_eegam>81.19 && invm_eegam<101.19);
			invflag = invm_e1e2>60.;
			
			//invflag = invm_e1e2>60. && !(invm_e1e2>(91.2-15)&&invm_e1e2<(91.2+15));
			//cout<<"invflag : "<<invflag<<endl;
		      }// if(pho_select>=1 && ele_select>=2)
		    
		    if(invflag)
		      {
			invev++;
			invevwt = eleeffFac[0]*eleeffFac[1]*0.98 + invevwt;
		      }
		    

		    ///read HLT br here
		    readHLTbr(ientry);
		    
		    //call HLT
		    bool passhlt = false;
		    TLorentzVector *ele1p4 = new TLorentzVector();
		    TLorentzVector *ele2p4 = new TLorentzVector();
		    TLorentzVector *pho_p4;
		    
		    if(pho_select>=1 && ele_select>=2 && invflag)
		      {

			
			ele1p4 = (TLorentzVector*)patelp4->At(ele_index[0]);
			ele2p4 = (TLorentzVector*)patelp4->At(ele_index[1]);

			//chnaged the defn of Etso as to be consistent wit hteh HEEP note
			ele1p4 = setPtEtaPhiE(ele_index[0]);
			ele2p4 = setPtEtaPhiE(ele_index[1]);
			
			pho_p4  = (TLorentzVector*)patphop4->At(phoindex[0]);
			//cout<<"myPlot.C after scaling: inside some other function: pt : "<<pho_p4->Pt()<<endl;
			//cout<<"event = "<<event<<endl;
			//cout<<"calling HLT"<<endl;
			
			bool passhlt_e1 = findHLTs(ele1p4, run, filename); 
			//cout<<"passhlt_e1 = "<<passhlt_e1<<endl;
			bool passhlt_e2 = findHLTs(ele2p4, run, filename); 
			//cout<<"passhlt_e2 = "<<passhlt_e2<<endl;
			bool passhlt_pho = findHLTs(pho_p4, run, filename);
			//cout<<"passhlt_pho = "<<passhlt_pho<<endl;
			//if(passhlt_e1 || passhlt_e2 || passhlt_pho)

			
			bool evpasshlt = findHLTPaths(run, filename);
			
			if( (passhlt_e1 || passhlt_e2) && evpasshlt)
			  passhlt = true;
			else
			  passhlt= false;
			
			
			//cout<<"passHLT = "<<passhlt<<endl;
			//cout<<"returned from data HLT"<<endl;
		      }
		    
		    if(filename=="bkg")
		      {
			passhlt = true;
		      }
		    
		    if(passhlt)
		      {
			hltev++;
			hltevwt =  eleeffFac[0]*eleeffFac[1]*0.98 + hltevwt;
			//cout<<"run no    = "<<run<<endl;
			//cout<<"event no. = "<<event<<endl;
		      }
		    
		    if(passhlt && pho_select>=1 && ele_select>=2 && invflag)
		      {
			
			double maxmass = TMath::Max(invm_e1pho,invm_e2pho);
			double minmass = TMath::Min(invm_e1pho,invm_e2pho);
			
			if(maxmass>180.)
			  {
			    m150 ++;
			    m150wt = eleeffFac[0]*eleeffFac[1]*0.98 + m150wt;
			  }
			
			if(maxmass>350.)
			  {
			    m350 ++;
			    m350wt = eleeffFac[0]*eleeffFac[1]*0.98 + m350wt; 
			  }
			
			if(maxmass>500.)
			  {
			    m550 ++;
			    m550wt = eleeffFac[0]*eleeffFac[1]*0.98 + m550wt;
			  }
			
			//filling of tree
			finalPt          = pho_p4->Pt();
			finalMaxMass     = maxmass; 
			finalMee         = invm_e1e2;
			finalMeegam      = invm_eegam;
			finalMinMass     = minmass;
			w                = weight;
			itype            = ftype;
			ele1Pt = ele1p4->Pt();  //// so now here corrected Pt consistent with the HEEP defn of Pt is going
			ele2Pt = ele2p4->Pt();  ////// same here
			
			ele1Px = ele1p4->Px();
			ele1Py = ele1p4->Py();
			ele1Pz = ele1p4->Pz();
			ele1E = ele1p4->E();
			
			ele2Px = ele2p4->Px();
			ele2Py = ele2p4->Py();
			ele2Pz = ele2p4->Pz();
			ele2E = ele2p4->E();
			
			phoE = pho_p4->E();
			phoPx = pho_p4->Px();
			phoPy = pho_p4->Py();
			phoPz = pho_p4->Pz();
			
			runno = run;
			eventno = event;
			lumino = lumis;
			TVector3 *phop3 = (TVector3 *) patphocalopos->At(phoindex[0]);
			finalEta = phop3->Eta();
			
			int scind = patelescind[ele_index[0]];
			TLorentzVector *sp4 = (TLorentzVector*)scp4->At(scind);
			double ele1eta = sp4->Eta();
			ele1Eta = ele1eta;

			int scind2 = patelescind[ele_index[1]];
			TLorentzVector *sp24 = (TLorentzVector*)scp4->At(scind2);
			double ele2eta = sp24->Eta();
			ele2Eta = ele2eta;


			
			int scind1   = patelescind[ele_index[0]];
			scind2   = patelescind[ele_index[1]];
			
			TLorentzVector *elscp4_1 = (TLorentzVector*)scp4->At(scind1);
			TLorentzVector *elscp4_2 = (TLorentzVector*)scp4->At(scind2);
			double eleta1 = elscp4_1->Eta();
			double eleta2 = elscp4_2->Eta();
			
			
			if( fabs(eleta1)<1.442 )
			  ele1effFac = 0.978;
			
			if( fabs(eleta1)<2.5 && fabs(eleta1)>1.56 )
			  ele1effFac =0.994;
			
			if( fabs(eleta2)<1.442 )
			  ele2effFac = 0.978;
			
			if( fabs(eleta2)<2.5 && fabs(eleta2)>1.56 )
			  ele2effFac =0.994;
			
			
			////store information from PU in case of MC
			if(filename=="bkg" || filename=="estar")
			  {
			    pileupvtx = pileup_nvtx[1];

			    ///for 3-D PU reweighting
			    if(pileup_bunchXing[0]==-1)
			      npm1 = pileup_nvtx[0];

			    if(pileup_bunchXing[1]==0)
			      np0 = pileup_nvtx[1];

			    if(pileup_bunchXing[2]==1)
			      npp1 = pileup_nvtx[2];

			    
			  }


			ngoodvtx = nvtx();
			
			tree->Fill();
		      }
		  }//if(phomcmatch)
	      }//if(elemcmatch)
	      
	    }//if(goodvertex() && noscraping())
	}//if( lumiRunSelection(jentry) )
       
       
    }//for (Long64_t jentry=0; jentry<nentries;jentry++) 
    cout<<"ntot = "<<ntot<<endl;
    cout<<"ev_passjson = "<<ev_passjson<<endl;
    cout<<"ev_goodvtx = "<<ev_goodvtx<<endl;
    cout<<"heepev = "<<heepev<<endl;
    cout<<"phoidev= "<<phoidev<<endl;
    cout<<"invev = "<<invev<<endl;
    cout<<"hltev = "<<hltev<<endl;
    cout<<"m150 = "<<m150<<endl;
    cout<<"m350 = "<<m350<<endl;
    cout<<"m550 = "<<m550<<endl;
    cout<<"weight = "<<w<<endl;
    cout<<"=====NOW eff wt====="<<endl;
    cout<<"heepevwt = "<<heepevwt<<endl;
    cout<<"phoidevwt= "<<phoidevwt<<endl;
    cout<<"invevwt = "<<invevwt<<endl;
    cout<<"hltevwt = "<<hltevwt<<endl;
    cout<<"m150wt = "<<m150wt<<endl;
    cout<<"m350wt = "<<m350wt<<endl;
    cout<<"m550wt = "<<m550wt<<endl;
    
    cout<<"======using conversion tools====="<<endl;
    cout<<"nele1ConvTool = "<<nele1ConvTool<<endl;
    cout<<"nele2ConvTool = "<<nele2ConvTool<<endl;
    cout<<"nanyeleConvTool = "<<nanyeleConvTool<<endl;
    cout<<"nbotheleConvTool = "<<nbotheleConvTool<<endl;
    cout<<"nConvel = "<<nConvel<<endl;
    cout<<"Thus using conversion tools fraction of times we find atleast one elctron converted in e* sample = "<<(double)nanyeleConvTool/(double)heepev<<endl;

    fileindex++;
    

  }//while(!infile.eof())
  //infile.close();
  file->cd();
  tree->Write();
  file->Write();
  file->Close();
  t1.Print() ;

}
