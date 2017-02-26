float myPlot::getLumi(){
 
  //double lumi = 902.978045209; //- till 12th sept - wrong lumi by lumiCalc.py
  //double lumi = 948.735; //-  righ lumi by lumiCalc2.py
  //double lumi = 1916.;
  //double lumi = 2141.;
  //double lumi = 2152.; //from corrected Lumicalc2.py
  //double lumi = 4604.; ///as on 5th Nov
  //double lumi = 4679.; //as of 24th nov
  double lumi = 4680.; //as of 15th dec
  
 return lumi;
}

bool myPlot::goodvertex()
{

  bool dec = false;
  //do scraping and vtx filters                                                                                                                                        
  //   vertex selection                                                                                                                                                
  //turn off when doing pho fakerate                                                                                                                                   
  int goodVertex = -1;
  for(int iVtx=0; iVtx<vtx_std_n; ++iVtx){
    TVector3 * pos_vtx = (TVector3 *) vtx_std_xyz->At(iVtx);

    ////isFake: see in vertex.h - comments above
    //if vtx_std_x2dof==0 && vtx_std_ndof == 0
    if(vtx_std_ndof[iVtx]>=4 && fabs(pos_vtx->Z())<=24. && fabs(pos_vtx->Perp())<=2. && !(vtx_std_x2dof==0 && vtx_std_ndof == 0) )
      {
	//cout<<"found a good vertex"<<endl;
	goodVertex = iVtx;
	dec = true;
	break;
      }//if(vtx_std_ndof[iVtx]>=4 && fabs(pos_vtx->Z())<=24. && fabs(pos_vtx->Perp())<=2.)                                                                           
  }//for(int iVtx=0; iVtx<vtx_std_n; ++iVtx)                                                                                                                         

  if(goodVertex == -1)
    {
      //cout<<"didnt find a good vertex"<<endl;                                                                                                                      
      dec = false;
      return dec;
    }

  return dec;

}//bool myPlot::goodvertex()               

bool myPlot::noscraping()
{
  // scraping events                                                                                                                                                   
  bool passscrapefilter =true;

  if (gentrksize > 10) {
    Int_t goodTracks = 0;
    for (Int_t tks=0; tks<gentrksize; tks++) {
      //      cout<<"quality = "<<tk_quality[tks]<<endl;    
      if( (gentrkqualityMask[tks]&2)==2){// || gentrkquality[tks] == 111){           
	goodTracks++;
      }
    }

    //cout<<"ntrks:frac of tracks:"<<gentrksize<<":"<<(float)goodTracks/(float)gentrksize<<endl;
    
    if ((float)goodTracks/(float)gentrksize > .25)
      passscrapefilter = true;
    else
      passscrapefilter = false;
  }

  else
    passscrapefilter = true; // if(ntracks<10) --> pass anyway       

  if(passscrapefilter == false)
    {
      //cout<<"i'm a scraper"<<endl;                                                                                                                                   
      return passscrapefilter;
    }

  return passscrapefilter;
}//bool myPlot::noscraping()


bool myPlot::findHLTs(TLorentzVector *p4, int run, string filename) {
  bool dbg=false;
  bool mydbg=false;
  //  load the labels corresponding to the photon10 and photon15 triggers
  if(dbg)cout<<"-----------------------------------------------  entering findHLTS  ------------------------------------------------"<<endl;
  std::vector<std::string> hltLabels;
  std::string s="abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz";
  //for data choose trigger by run number
  if(filename == "data"){
    //cout<<"data HLT called"<<endl;
    string hlt = "aaa";

    
    //if(run>=160404 && run<=163869)
    if(run>=160404 && run<=165000)
      {
	hlt="hltDoublePhoton33EgammaLHEDoubleFilter";
	hltLabels.push_back(hlt);
      }
    
    //if(run>=165088 && run<=178380)
    if(run>165000)
      {
	hlt="hltDiEle33CaloIdLPixelMatchDoubleFilter";
	hltLabels.push_back(hlt);
	hlt = "hltDiEle45CaloIdLPixelMatchDoubleFilter";
	hltLabels.push_back(hlt);
	hlt = "hltDiEle33CaloIdTPixelMatchDoubleFilter";
	hltLabels.push_back(hlt);
      }
    
      //if(run>178380)
      //hlt="hltDiEle45CaloIdLPixelMatchDoubleFilter";

    
    //hltLabels.push_back(hlt);

  }
  //for MC don't use it for now - because HEEP paper is also not using
  else{
    //hltLabels.push_back("hltDoublePhoton33EgammaLHEDoubleFilter");
    //hltLabels.push_back("hltDiEle33CaloIdLPixelMatchDoubleFilter");
    return true;
    
  }

  //cout<<"hltLabels[0] = "<<hltLabels[0]<<endl;
  ///////////////////////DEBUG///////////////////////////////
  /*for(int ii=0; ii<hlt_label_names_1->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_1)[ii])
        cout<<"FOUND MY HLT 1"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_2->size();ii++)
      {
      if(hltLabels[0]==(*hlt_label_names_2)[ii])
        cout<<"FOUND MY HLT 2"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_3->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_3)[ii])
        cout<<"FOUND MY HLT 3"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_4->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_4)[ii])
        cout<<"FOUND MY HLT 4"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_5->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_5)[ii])
        cout<<"FOUND MY HLT 5"<<endl;
	    }

  for(int ii=0; ii<hlt_label_names_6->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_6)[ii])
        cout<<"FOUND MY HLT 6"<<endl;
    }
  */

  //////////////////////END OF DEBUG///////////////////////

  unsigned int nHLTs = hltLabels.size();
  if(dbg) {for(int l=0; l<nHLTs;  l++) cout<<l<<" "<<nHLTs<<"  label="<<hltLabels[l];  cout<<endl;}

  //  Unfortunately, the way things are stored makes this extremely complicated, so I have to loop over all the candidates and look for triggers.
  //  In the futrue, I think we should just have the triggers and leave the candidate matching to the lepton list or to the ntuple analysis.
  //  The current cure seems worse than the disease.
  if(dbg) cout<<"hlt_n = "<<hlt_n<<endl;
  for(unsigned int iCandidate=0; iCandidate<hlt_n; ++iCandidate) {
    if(dbg) cout<<"hlt_n = "<<hlt_n<<endl;
    TLorentzVector * p4_hlt = (TLorentzVector *) hlt_p4->At(iCandidate);
    if(dbg) cout<<"ele eta = "<<p4->Eta()<<" ele phi = "<<p4->Phi()<<" trig eta = "<<p4_hlt->Eta()<<" trig phi = "<<p4_hlt->Phi()<<" hlt_dr = "<< p4->DeltaR(*p4_hlt) <<endl;
    if(p4->DeltaR(*p4_hlt)<0.2) ///matching
      {
        // only look at triggers that match in eta phi
        //  loop over trigger bits: b {loop over candidates {  loop over the triggers we accept: l  (and we have to do this for the four trigger words, hence 4 lines of code)
      //  dbg =1;
	if(dbg) cout<<"close cand "<<endl;
	ULong64_t a = 1;
	for (int b=0; b<64;  b++) { if(hlt_candpath_1[iCandidate] >> b & a) { s=(*hlt_label_names_1)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"1  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
	for (int b=0; b<64;  b++) { if(hlt_candpath_2[iCandidate] >> b & a) { s=(*hlt_label_names_2)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"2  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
	for (int b=0; b<64;  b++) { if(hlt_candpath_3[iCandidate] >> b & a) { s=(*hlt_label_names_3)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"3  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
	for (int b=0; b<64;  b++) { if(hlt_candpath_4[iCandidate] >> b & a) { s=(*hlt_label_names_4)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"4  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
	for (int b=0; b<64;  b++) { if(hlt_candpath_5[iCandidate] >> b & a) { s=(*hlt_label_names_5)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"5  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
	for (int b=0; b<64;  b++) { if(hlt_candpath_6[iCandidate] >> b & a) { s=(*hlt_label_names_6)[b];  if(dbg)cout<<"my i cand = "<<iCandidate<<b<<"6  "<<s<<endl; \
				      for(int l=0; l<nHLTs;  l++) {if(s == hltLabels[l]) { if(mydbg) cout<<"FOUND..."<<endl; return true;} } } }
      }//if(p4->DeltaR(*p4_hlt)<0.2)
  }//for(unsigned int iCandidate=0; iCandidate<hlt_n; ++iCandidate)
  return false;
}//bool myPlot::findHLTs(TLorentzVector *p4, int run, int itype)


////////////6th Dec - only for telling if event has passed HLT
bool myPlot::findHLTPaths(int run, string filename) {
  bool dbg=false;
  bool mydbg=false;
  //  load the labels corresponding to the photon10 and photon15 triggers
  if(dbg)cout<<"-----------------------------------------------  entering findHLTS  ------------------------------------------------"<<endl;
  std::vector<std::string> hltLabels;
  std::string s="abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz";
  //for data choose trigger by run number
  if(filename == "data"){
    //cout<<"data HLT called"<<endl;
    string hlt = "aaa";

    ///here give the names of the HLT path
    //if(run>=160404 && run<=163869)
    if(run>=160404 && run<=165000)
      {
	hlt="DoublePhoton33_";
	hltLabels.push_back(hlt);
      }
    
    //if(run>=165088 && run<=178380)
    if(run>165000)
      {
	hlt="DoubleEle33_CaloIdL_v";
	hltLabels.push_back(hlt);
	hlt = "DoubleEle33_CaloIdT_v";
	hltLabels.push_back(hlt);
	hlt = "DoubleEle45_CaloIdL_v";
	hltLabels.push_back(hlt);
      }
    
      //if(run>178380)
      //hlt="hltDiEle45CaloIdLPixelMatchDoubleFilter";

    
    //hltLabels.push_back(hlt);

  }
  //for MC don't use it for now - because HEEP paper is also not using
  else{
    //hltLabels.push_back("hltDoublePhoton33EgammaLHEDoubleFilter");
    //hltLabels.push_back("hltDiEle33CaloIdLPixelMatchDoubleFilter");
    return true;
    
  }

  //cout<<"hltLabels[0] = "<<hltLabels[0]<<endl;
  ///////////////////////DEBUG///////////////////////////////
  /*for(int ii=0; ii<hlt_label_names_1->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_1)[ii])
        cout<<"FOUND MY HLT 1"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_2->size();ii++)
      {
      if(hltLabels[0]==(*hlt_label_names_2)[ii])
        cout<<"FOUND MY HLT 2"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_3->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_3)[ii])
        cout<<"FOUND MY HLT 3"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_4->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_4)[ii])
        cout<<"FOUND MY HLT 4"<<endl;
    }

  for(int ii=0; ii<hlt_label_names_5->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_5)[ii])
        cout<<"FOUND MY HLT 5"<<endl;
	    }

  for(int ii=0; ii<hlt_label_names_6->size();ii++)
    {
      if(hltLabels[0]==(*hlt_label_names_6)[ii])
        cout<<"FOUND MY HLT 6"<<endl;
    }
  */

  //////////////////////END OF DEBUG///////////////////////

  if(dbg){
    cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
  }

  unsigned int nHLTs = hltLabels.size();
  if(dbg) {for(int l=0; l<nHLTs;  l++) cout<<l<<" "<<nHLTs<<"  paths "<<hltLabels[l];  cout<<endl;}

  ULong64_t a = 1;

  if(dbg)
    {
      cout<<"HLT1 size:"<<hlt_path_names_HLT1_1->size()<<endl;
      cout<<"HLT2 size:"<<hlt_path_names_HLT1_2->size()<<endl;
      cout<<"HLT3 size:"<<hlt_path_names_HLT1_3->size()<<endl;
      cout<<"HLT4 size:"<<hlt_path_names_HLT1_4->size()<<endl;
      cout<<"HLT5 size:"<<hlt_path_names_HLT1_5->size()<<endl;
      cout<<"HLT6 size:"<<hlt_path_names_HLT1_6->size()<<endl;
    }

  for(int ipath=0; ipath<64; ipath++)
    {
      int fstr;
      string hltname = "aaa";

      
      
      for(int ii=0; ii<nHLTs; ii++)
	{
	  
	  //////HLT1
	  if(hlt_path_names_HLT1_1->size()>=(ipath+1) && hlt_path_names_HLT1_1->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_1)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 1:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_1 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		    
		    return true;
		  }//if(hlt1_bit_1 >> ipath & a)
	      }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_1->size()>=ipath)
	  
	  //////HLT2
	  if(hlt_path_names_HLT1_2->size()>=(ipath+1) && hlt_path_names_HLT1_2->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_2)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 2:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_2 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		  
		    return true;
		  }//if(hlt1_bit_2 >> ipath & a)
	      }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_2->size()>=(ipath+1))
	  
	  //////HLT3
	  if(hlt_path_names_HLT1_3->size()>=(ipath+1) && hlt_path_names_HLT1_3->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_3)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 3:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_3 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		    
		    return true;
		  }//if(hlt1_bit_3 >> ipath & a)
	      }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_3->size()>=(ipath+1))

	  //////HLT4
	  if(hlt_path_names_HLT1_4->size()>=(ipath+1) && hlt_path_names_HLT1_4->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_4)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 4:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_4 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		    
		    return true;
		  }//if(hlt1_bit_4 >> ipath & a)
	      }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_4->size()>=(ipath+1))
	  
	  //////HLT5
	  if(hlt_path_names_HLT1_5->size()>=(ipath+1) && hlt_path_names_HLT1_5->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_5)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 5:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_5 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		    
		    return true;
		  }//if(hlt1_bit_5 >> ipath & a)
	      }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_5->size()>=(ipath+1))
	  
	  //////HLT6
	  if(hlt_path_names_HLT1_6->size()>=(ipath+1) && hlt_path_names_HLT1_6->size()!=0){
	    hltname =  (*hlt_path_names_HLT1_6)[ipath] ;
	    if(dbg)
	      {
		cout<<"HLT name in 6:hlt name in vector:"<<hltname<<":"<<hltLabels[ii]<<endl;
	      }
	    
	    fstr = hltname.find(hltLabels[ii],0);
	    if(fstr!=string::npos)
	      {
		if(dbg)
		  {
		    cout<<"found this HLT:"<<hltLabels[ii]<<" in "<<hltname<<endl;
		  }
		
		if(hlt1_bit_6 >> ipath & a)
		  {
		    if(dbg)
		      {
			cout<<"hurray!passed HLT (:)"<<hltLabels[ii]<<endl;
			cout<<"run:lumis:event:"<<run<<":"<<lumis<<":"<<event<<endl;
		      }
		    
		    return true;
		  }//if(hlt1_bit_6 >> ipath & a)
	    }// if(fstr!=string::npos) 
	  }//if(hlt_path_names_HLT1_6->size()>=(ipath+1))

	}//for(int ii=0; ii<nHLTs; ii++)
    }//for(int ipath=0; ipath<64; ipath++)

  return false;
  
}//bool myPlot::findHLTPaths(int run, string filename)
       
  


TLorentzVector* myPlot::scalephoE(int ipho)
{
  double perc = 0;
  
  TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(ipho);
  TVector3 * phop3;
  phop3 = (TVector3 *) patphocalopos->At(ipho);
  double phoeta=phop3->Eta();
  
  int scind = patphoscind[ipho];
  double e3x3 = sce3x3[scind];
  double rawenergy = scrawEnergy[scind];
  double r9 = e3x3/rawenergy;
  
  if(r9 > 0.94)
    {
      if( fabs(phoeta)<1.4442 )
	perc = 0.02;
      
      if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 )
	perc = 0.0;
      
    }
  
  if(r9 < 0.94)
    {
      if( fabs(phoeta)<1.4442 )
	perc = 0.03;
      
      if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 )
	perc = 0.02;
      
    }

  //perc = 1;

  //cout<<"commonFunctions before scaling:ipho : pt : "<<ipho<<" : "<<pho_p4->Pt()<<endl;
  double num = (1.0 + perc/100.);
  //cout<<"num = "<<num<<endl;
  pho_p4->SetE(pho_p4->E()*num);
  pho_p4->SetPx(pho_p4->Px()*num);
  pho_p4->SetPy(pho_p4->Py()*num);
  pho_p4->SetPz(pho_p4->Pz()*num);
  //cout<<"commonFunctions after scaling:ipho : pt : "<<ipho<<" : "<<pho_p4->Pt()<<endl;
  
  return pho_p4;
  
}

TLorentzVector* myPlot::scaleeleE(int iele)
{
  TLorentzVector *elp4 = new TLorentzVector();
  elp4 = setPtEtaPhiE(iele);
  
  int scind = patelescind[iele];
  int flag_scind = 1;
  double eleeta;
  flag_scind = scind<scsize && !(scind<0);
  if(flag_scind)
    {
      TLorentzVector *sp4 = (TLorentzVector*)scp4->At(scind);
      eleeta = sp4->Eta();
    }

  bool elebarrel = fabs(eleeta) < 1.442;
  bool eleendcap = fabs(eleeta) > 1.56 && fabs(eleeta) < 2.5;
  
  double corr = 1;
  double perc = 0;
  
  if( elebarrel)
    {
      perc = 0;
    }

  if( eleendcap)
    {
      perc = 0;
    }
  
  corr = 1+perc/100.;

  TLorentzVector *elscp4 = (TLorentzVector*)patelsc->At(iele);
  TLorentzVector *elep4 = (TLorentzVector*)patelp4->At(iele);
  double eleet = elscp4->Energy()*sin(elep4->Theta());
  //cout<<"inside setptetaphi:elet:"<<eleet<<endl;                                                                                                                     
  elp4->SetPtEtaPhiE( eleet*corr, elep4->Eta(), elep4->Phi(),eleet*corr*cosh(elep4->Eta()) );
  
  return elp4;

}

bool myPlot::cutshat(double shat)
{
  bool passev = false;

  int ele1=-1, ele2=-1;
  for (int genpart = 0; genpart < gensize; genpart++){
    if (genpdgid[genpart]== 11 && genstatus[genpart]==3 && ele1<0)
      ele1=genpart;
    if (genpdgid[genpart]==-11 && genstatus[genpart]==3 && ele2<0)
      ele2=genpart;
  }
  if (ele1 >= 0 && ele2 >= 0) {
    TLorentzVector * genp4_1 = (TLorentzVector *) genp4->At(ele1);
    TLorentzVector * genp4_2 = (TLorentzVector *) genp4->At(ele2);
    double mass = ( (*genp4_1)+(*genp4_2) ).M();
    if(mass<shat)
      passev = true;
  }
  return passev;
  
}//int myPlot::selectee()


bool myPlot::selectee(TLorentzVector *ele1, TLorentzVector *ele2)
{
  bool ele1MCMatch = false;
  bool ele2MCMatch = false;
  bool passMCeleCheck = false;
  int match1 = 999;
  int match2 = 999;

  //  cout<<"gp_n = "<<gp_n<<endl;
  for (int genpart = 0; genpart < gensize; genpart++){

    //cout<<"genpart1 = "<<genpart<<endl; 
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    //cout<<"got gen_p41"<<endl;

    //match to any electron or any photon ONLY coming from a gluon or quark or electron
    //if( gp_status[genpart] ==1 && ( gen_p4->DeltaR(*ele1) < 0.1 ) && ( (fabs(gp_pdgid[genpart])==11 && (fabs(gp_pdgid[gp_mother[genpart]])==23 || fabs(gp_pdgid[gp_mother[genpart]])==24|| fabs(gp_pdgid[gp_mother[genpart]])==11)) || (fabs(gp_pdgid[genpart])==22 && ( (abs(gp_pdgid[gp_mother[genpart]])>= 1 && abs(gp_pdgid[gp_mother[genpart]])<= 5 ) || abs(gp_pdgid[gp_mother[genpart]])==21 || abs(gp_pdgid[gp_mother[genpart]])==11 ))) ){
      
    
    if(( (fabs(genpdgid[genpart])==11 &&  (fabs(genpdgid[genmother[genpart]])==23 || fabs(genpdgid[genmother[genpart]])==24|| fabs(genpdgid[genmother[genpart]])==11)) || (genpdgid[genpart]==22 && ( (abs(genpdgid[genmother[genpart]])>= 1 &&abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 ||abs(genpdgid[genmother[genpart]])==11 ))) && genstatus[genpart]==1  && gen_p4->DeltaR(*(ele1)) < 0.1 ){
    ele1MCMatch = true;
      match1 = genpart;
    }//if( fabs(genpdgid[genpart])==11 && ...)
    

  }//for (int genpart = 0; genpart < genn; genpart++)

  //cout<<"genn2 = "<<genn<<endl;
  for (int genpart = 0; genpart < gensize; genpart++){
    //cout<<"genpart1 = "<<genpart<<endl; 
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    //cout<<"got gen_p42"<<endl;

    //      if( genstatus[genpart] ==1 && ( gen_p4->DeltaR(*ele2) < 0.1 ) && ( (fabs(genpdgid[genpart])==11 && (fabs(genpdgid[genmother[genpart]])==23 || fabs(genpdgid[genmother[genpart]])==24|| fabs(genpdgid[genmother[genpart]])==11)) || (fabs(genpdgid[genpart])==22 && ( (abs(genpdgid[genmother[genpart]])>= 1 && abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 || abs(genpdgid[genmother[genpart]])==11 ))) ){

        if(( (fabs(genpdgid[genpart])==11 &&  (fabs(genpdgid[genmother[genpart]])==23 || fabs(genpdgid[genmother[genpart]])==24|| fabs(genpdgid[genmother[genpart]])==11)) || (genpdgid[genpart]==22 && ( (abs(genpdgid[genmother[genpart]])>= 1 &&abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 ||abs(genpdgid[genmother[genpart]])==11 ))) && genstatus[genpart]==1  && gen_p4->DeltaR(*(ele2)) < 0.1 ){
      ele2MCMatch = true;
      match2 = genpart;
    }//if( fabs(genpdgid[genpart])==11 && ...)

  }//for (int genpart = 0; genpart < genn; genpart++)                                    
  
  //if(ele1MCMatch && ele2MCMatch && (match1!=match2)) passMCeleCheck = true;
  if(ele1MCMatch && ele2MCMatch) passMCeleCheck = true;

  return passMCeleCheck;
  
}//int myPlot::selectee()



bool myPlot::selecte(TLorentzVector *ele1)
{
  bool ele1MCMatch = false;
  bool passMCeleCheck = false;
  int match1 = 999;
  
  //  cout<<"gp_n = "<<gp_n<<endl;
  for (int genpart = 0; genpart < gensize; genpart++){

    //cout<<"genpart1 = "<<genpart<<endl; 
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    //cout<<"got gen_p41"<<endl;

    //match to any electron or any photon ONLY coming from a gluon or quark or electron
    //if( gp_status[genpart] ==1 && ( gen_p4->DeltaR(*ele1) < 0.1 ) && ( (fabs(gp_pdgid[genpart])==11 && (fabs(gp_pdgid[gp_mother[genpart]])==23 || fabs(gp_pdgid[gp_mother[genpart]])==24|| fabs(gp_pdgid[gp_mother[genpart]])==11)) || (fabs(gp_pdgid[genpart])==22 && ( (abs(gp_pdgid[gp_mother[genpart]])>= 1 && abs(gp_pdgid[gp_mother[genpart]])<= 5 ) || abs(gp_pdgid[gp_mother[genpart]])==21 || abs(gp_pdgid[gp_mother[genpart]])==11 ))) ){
      
    
    if(( (fabs(genpdgid[genpart])==11 &&  (fabs(genpdgid[genmother[genpart]])==23 || fabs(genpdgid[genmother[genpart]])==24|| fabs(genpdgid[genmother[genpart]])==11)) || (genpdgid[genpart]==22 && ( (abs(genpdgid[genmother[genpart]])>= 1 &&abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 ||abs(genpdgid[genmother[genpart]])==11 ))) && genstatus[genpart]==1  && gen_p4->DeltaR(*(ele1)) < 0.1 ){
    ele1MCMatch = true;
      match1 = genpart;
    }//if( fabs(genpdgid[genpart])==11 && ...)
    

  }//for (int genpart = 0; genpart < genn; genpart++)

   if(ele1MCMatch) passMCeleCheck = true;

  return passMCeleCheck;
  
}//int myPlot::selecte()


bool myPlot::selectgam(TLorentzVector *pho)
{
  bool passMCphoCheck = false;
  
  for (int genpart = 0; genpart < gensize; genpart++){
      TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
        
      //
      //if( (fabs(genpdgid[genpart])==22 || fabs(genpdgid[genpart])==11)  && genstatus[genpart]==1 && (gen_p4->DeltaR(*pho) < 0.1) )



      //if( genstatus[genpart]==1 && (gen_p4->DeltaR(*pho) < 0.1) && genpdgid[genmother[genpart]] != 111 && ((genpdgid[genpart]==22 && ( (abs(genpdgid[genmother[genpart]])>= 1&&abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 || abs(genpdgid[genmother[genpart]])==11 || abs(genpdgid[genmother[genpart]])==23))   || (fabs(genpdgid[genpart])==11&&  (fabs(genpdgid[genmother[genpart]])==23 ||fabs(genpdgid[genmother[genpart]])==24|| fabs(genpdgid[genmother[genpart]])==11) )) )
      if( ( (genpdgid[genpart]==22 && ( (abs(genpdgid[genmother[genpart]])>= 1 &&abs(genpdgid[genmother[genpart]])<= 5 ) || abs(genpdgid[genmother[genpart]])==21 ||abs(genpdgid[genmother[genpart]])==11 || abs(genpdgid[genmother[genpart]])==23 ) ) || (fabs(genpdgid[genpart])==11&&(fabs(genpdgid[genmother[genpart]])==23 || fabs(genpdgid[genmother[genpart]])==24||fabs(genpdgid[genmother[genpart]])==11) )) && genstatus[genpart]>=1 && genpdgid[genmother[genpart]] != 111 && gen_p4->DeltaR(*(pho)) < 0.1 ) 
	passMCphoCheck = true;
      
      
      
      
    }

  return passMCphoCheck;
  
}//int myPlot::selectzgamma()

TLorentzVector* myPlot::setPtEtaPhiE(int iele)
{
  TLorentzVector *ele = new TLorentzVector();
  TLorentzVector *elscp4 = (TLorentzVector*)patelsc->At(iele);
  TLorentzVector *elep4 = (TLorentzVector*)patelp4->At(iele);
  double eleet = elscp4->Energy()*sin(elep4->Theta());
  //cout<<"inside setptetaphi:elet:"<<eleet<<endl;
  ele->SetPtEtaPhiE( eleet, elep4->Eta(), elep4->Phi(),eleet*cosh(elep4->Eta()) );
  // cout<<"inside set pt eta : ele->Pt():Energy:"<<ele->Pt()<<":"<<ele->E()<<endl;
  
  return ele;

}


void myPlot::branchStatus1( vector<string> &brname){
  
  for(int ibr=0; ibr<(int)brname.size(); ibr++)
    fChain->SetBranchStatus(brname[ibr].c_str(),1);
}



void myPlot::ConvTool(int iele, bool &eleConv)
{
  eleConv = false;
  
  if( pateleExpectednumberOfHits[iele]>0 && (fabs(pateleconvDist[iele])<0.02 && fabs(pateleconvDcot[iele])<0.02) )
    {
      //cout<<"pateleExpectednumberOfHits:convDist:convDcot "<<pateleExpectednumberOfHits[iele]<<":"<<pateleconvDist[iele]<<":"<<pateleconvDcot[iele]<<endl;
      eleConv = true;
    }
  
  return;
  
}//int myPlot::ConvTool()




bool myPlot::phoID(int ipho)
{
  double HoE          = patphohadronicOverEm[ipho];
  double ecalIso      = patphoecalRecHitSumEtConeDR04[ipho];
  double hcalIso      = patphohcalTowerSumEtConeDR04[ipho];
  double trkIso       = patphotrkSumPtHollowConeDR04[ipho];
  double sieie        = patphosigmaIetaIeta[ipho];
  //int haspixel        = pho_haspixseed[ipho];

  double cutHoE       = 0.05;

  TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(ipho);
  double et              = pho_p4->Pt();

  double cuttrkIso    = 2.0+0.001*et;
  double cutecalIso   = 4.2+0.006*et;
  double cuthcalIso      = 2.2+0.0025*et;

  double cutsieie;

  TVector3 * phop3;
  phop3 = (TVector3 *) patphocalopos->At(ipho);
  double phoeta=phop3->Eta();
  if( (fabs(phoeta)<1.4442) ) cutsieie = 0.013;
  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 ) cutsieie = 0.03;

  bool flag = false;
  //flag = HoE<cutHoE && trkIso<cuttrkIso && ecalIso<cutecalIso && hcalIso<cuthcalIso && sieie<cutsieie  && patphohasPixelSeed[ipho]==0;
  flag = HoE<cutHoE && trkIso<cuttrkIso && ecalIso<cutecalIso && hcalIso<cuthcalIso && sieie<cutsieie;
  //flag = HoE<cutHoE && hcalIso<cuthcalIso && sieie<cutsieie;

  return flag;
}


bool myPlot::loosephoID(int ipho)
{
  
  bool passpho = false;
  
  int pho = ipho;   //passed ipho is that of reco

  TVector3 *phop3 = (TVector3 *) patphocalopos->At(pho);
  double phoeta=phop3->Eta();
  TLorentzVector *phop4  = (TLorentzVector*)patphop4->At(pho);

  double cutecal     = TMath::Min(5*(4.2 + 0.006*phop4->Pt()),0.2*phop4->Pt());
  double cuthcal     = TMath::Min(5*(2.2 + 0.0025*phop4->Pt()),0.2*phop4->Pt());  
  double cuttrksio   = TMath::Min(5*(3.5 + 0.001*phop4->Pt()),0.2*phop4->Pt());
  double cuthoe      = 0.05;
  
  //bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe && patphohasPixelSeed[ipho]==0) ;
  bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe) ;

  
  double fcuttrksio     = (3.5 + 0.001*phop4->Pt());
  double fcutecal       = (4.2 + 0.006*phop4->Pt());
  double fcuthcal       = (2.2 + 0.0025*phop4->Pt());
  double fcutsieie;
  
  
  if( fabs(phoeta)<1.4442 )
    fcutsieie = 0.013;

  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 )
    fcutsieie =0.03;

  bool flipiso =  ( patphotrkSumPtHollowConeDR04[pho] > fcuttrksio ) || ( patphoecalRecHitSumEtConeDR04[pho] > fcutecal ) || ( patphohcalTowerSumEtConeDR04[pho] > fcuthcal ) || (patphosigmaIetaIeta[pho] > fcutsieie);  
  
  if( looseid && flipiso )
    passpho = true;
  
  return passpho;


}



bool myPlot::phoIDforloosesel(int ipho)
{
  double HoE          = patphohadronicOverEm[ipho];
  double ecalIso      = patphoecalRecHitSumEtConeDR04[ipho];
  double hcalIso      = patphohcalTowerSumEtConeDR04[ipho];
  double trkIso       = patphotrkSumPtHollowConeDR04[ipho];
  double sieie        = patphosigmaIetaIeta[ipho];
  //int haspixel        = pho_haspixseed[ipho];

  double cutHoE       = 0.05;

  TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(ipho);
  double et              = pho_p4->Pt();

  double cuttrkIso    = 2.0+0.001*et;
  double cutecalIso   = 4.2+0.006*et;
  double cuthcalIso      = 2.2+0.0025*et;

  double cutsieie;

  TVector3 * phop3;
  phop3 = (TVector3 *) patphocalopos->At(ipho);
  double phoeta=phop3->Eta();
  if( (fabs(phoeta)<1.4442) ) cutsieie = 0.0105;
  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 ) cutsieie = 0.03;

  bool flag = false;
  //flag = HoE<cutHoE && trkIso<cuttrkIso && ecalIso<cutecalIso && hcalIso<cuthcalIso && sieie<cutsieie  && patphohasPixelSeed[ipho]==0;
  //flag = HoE<cutHoE && hcalIso<cuthcalIso && sieie<cutsieie && patphohasPixelSeed[ipho]==0;
  flag = HoE<cutHoE && hcalIso<cuthcalIso && sieie<cutsieie ;

  return flag;
}


bool myPlot::loosephoIDforcontrolR(int ipho)
{
  
  bool passpho = false;
  
  int pho = ipho;   //passed ipho is that of reco

  TVector3 *phop3 = (TVector3 *) patphocalopos->At(pho);
  double phoeta=phop3->Eta();
  TLorentzVector *phop4  = (TLorentzVector*)patphop4->At(pho);

  double cutecal     = TMath::Min(5*(4.2 + 0.006*phop4->Pt()),0.2*phop4->Pt());
  double cuthcal     = TMath::Min(5*(2.2 + 0.0025*phop4->Pt()),0.2*phop4->Pt());  
  double cuttrksio   = TMath::Min(5*(3.5 + 0.001*phop4->Pt()),0.2*phop4->Pt());
  double cuthoe      = 0.05;
  
  //bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe && patphohasPixelSeed[ipho]==0) ;
  //bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe) ;
  bool looseid = ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphohadronicOverEm[pho]< cuthoe) ;
  
  double fcuttrksio     = (3.5 + 0.001*phop4->Pt());
  double fcutecal       = (4.2 + 0.006*phop4->Pt());
  double fcuthcal       = (2.2 + 0.0025*phop4->Pt());
  double fcutsieie;
  
  
  if( fabs(phoeta)<1.4442 )
    fcutsieie = 0.0105;

  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 )
    fcutsieie =0.03;

  //bool flipiso =  ( patphotrkSumPtHollowConeDR04[pho] > fcuttrksio ) || ( patphoecalRecHitSumEtConeDR04[pho] > fcutecal ) || ( patphohcalTowerSumEtConeDR04[pho] > fcuthcal ) || (patphosigmaIetaIeta[pho] > fcutsieie);  
  bool flipiso =  (patphosigmaIetaIeta[pho] > fcutsieie);  
  
  if( looseid && flipiso )
    passpho = true;
  
  return passpho;


}


/////VgammaID
bool myPlot::vgphoID(int ipho)
{
  double HoE          = patphohadronicOverEm[ipho];
  double ecalIso      = patphoecalRecHitSumEtConeDR04[ipho];
  double hcalIso      = patphohcalTowerSumEtConeDR04[ipho];
  double trkIso       = patphotrkSumPtHollowConeDR04[ipho];
  double sieie        = patphosigmaIetaIeta[ipho];

  double cutHoE       = 0.05;
  double rho25        = rho;

  TLorentzVector *pho_p4  = (TLorentzVector*)patphop4->At(ipho);
  double et              = pho_p4->Pt();

  double cuttrkIso    = 2.0+0.001*et + 0.0167*rho25;
  double cutecalIso   = 4.2+0.006*et + 0.183*rho25;
  double cuthcalIso   = 2.2+0.0025*et + 0.062*rho25;

  double cutsieie;

  TVector3 * phop3;
  phop3 = (TVector3 *) patphocalopos->At(ipho);
  double phoeta=phop3->Eta();
  if( (fabs(phoeta)<1.4442) ) cutsieie = 0.011;
  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 ) cutsieie = 0.03;

  bool flag = false;
  //flag = HoE<cutHoE && trkIso<cuttrkIso && ecalIso<cutecalIso && hcalIso<cuthcalIso && sieie<cutsieie  && patphohasPixelSeed[ipho]==0;
  flag = HoE<cutHoE && trkIso<cuttrkIso && ecalIso<cutecalIso && hcalIso<cuthcalIso && sieie<cutsieie;
  //flag = HoE<cutHoE && hcalIso<cuthcalIso && sieie<cutsieie;

  return flag;
}


bool myPlot::vgloosephoID(int ipho)
{
  
  bool passpho = false;
  
  double rho25 = rho;

  int pho = ipho;   //passed ipho is that of reco

  TVector3 *phop3 = (TVector3 *) patphocalopos->At(pho);
  double phoeta=phop3->Eta();
  TLorentzVector *phop4  = (TLorentzVector*)patphop4->At(pho);

  double cutecal     = TMath::Min(5*(4.2 + 0.006*phop4->Pt() + 0.183*rho25),0.2*phop4->Pt());
  double cuthcal     = TMath::Min(5*(2.2 + 0.0025*phop4->Pt() + 0.062*rho25),0.2*phop4->Pt());  
  double cuttrksio   = TMath::Min(5*(3.5 + 0.001*phop4->Pt() + 0.0167*rho25),0.2*phop4->Pt());
  double cuthoe      = 0.05;
  
  double cutsieie = 0.014;
  //bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe && patphohasPixelSeed[ipho]==0) ;
  bool looseid = ( patphoecalRecHitSumEtConeDR04[pho] < cutecal ) && ( patphohcalTowerSumEtConeDR04[pho] < cuthcal ) && ( patphotrkSumPtHollowConeDR04[pho] < cuttrksio ) && ( patphohadronicOverEm[pho]< cuthoe) && (patphosigmaIetaIeta[pho] < cutsieie);

  
  double fcuttrksio     = (3.5 + 0.001*phop4->Pt() + 0.0167*rho25);
  double fcutecal       = (4.2 + 0.006*phop4->Pt() + 0.183*rho25);
  double fcuthcal       = (2.2 + 0.0025*phop4->Pt()+ 0.062*rho25);
  double fcutsieie;
  
  
  if( fabs(phoeta)<1.4442 )
    fcutsieie = 0.011;

  if( fabs(phoeta)>1.566 && fabs(phoeta)<2.5 )
    fcutsieie =0.03;

  bool flipiso =  ( patphotrkSumPtHollowConeDR04[pho] > fcuttrksio ) || ( patphoecalRecHitSumEtConeDR04[pho] > fcutecal ) || ( patphohcalTowerSumEtConeDR04[pho] > fcuthcal ) || (patphosigmaIetaIeta[pho] > fcutsieie);  
  
  if( looseid && flipiso )
    passpho = true;
  
  return passpho;


}

///////////VgammaID

bool myPlot::eleID(int iele)
{
  double HoE           = patelehadronicOverEm[iele];
  double ecalIso       = pateledr03EcalRecHitSumEt[iele];
  double hcal1Iso      = pateledr03HcalDepth1TowerSumEt[iele];
  double hcal2Iso      = pateledr03HcalDepth2TowerSumEt[iele];
  double trkIso        = pateledr03TkSumPt[iele];
  double sieie         = patelesigmaIetaIeta[iele];
  double deta          = pateledeltaEtaSuperClusterTrackAtVtx[iele];
  double dphi          = pateledeltaPhiSuperClusterTrackAtVtx[iele];
  double e1bye5        = patelee1x5[iele]/patelee5x5[iele];
  double e2bye5        = patelee2x5Max[iele]/patelee5x5[iele];
  
  double cutHoE       = 0.05;
  
  int mishits = pateleExpectednumberOfHits[iele];
  
  TLorentzVector *elep4  = (TLorentzVector*)patelp4->At(iele);
  double eleet              = elep4->Pt();

  int scind = patelescind[iele];
  TLorentzVector *sc4 = (TLorentzVector *) scp4->At(scind);
  double eleeta=sc4->Eta();
  //double et              = sc4->Pt();
  double et              = sc4->Energy()*sin(elep4->Theta());

  bool flag = false;

  //EB
  if( (fabs(eleeta)<1.442) )
    {
      flag = et>35. && pateleecalDrivenSeed[iele] && fabs(deta)<0.005 && fabs(dphi)<0.06 && HoE<cutHoE && trkIso<5.0 && (ecalIso+hcal1Iso)<(2+0.03*et) && (e1bye5>0.83|| e2bye5>0.94) && mishits==0;
    }

  //EE
  if( fabs(eleeta)>1.56 && fabs(eleeta)<2.5 )
    {
      double ehcalcut = 9999.;
      if(et<50.) ehcalcut = 2.5;
      if(et>50.) ehcalcut = (2.5+0.03*(et-50.));

      flag = et>40. && pateleecalDrivenSeed[iele] && fabs(deta)<0.007 && fabs(dphi)<0.06 && HoE<cutHoE && trkIso<5. && (ecalIso+hcal1Iso)<ehcalcut /*&& hcal2Iso<0.5*/ &&sieie<0.03 && mishits==0;
    }
  return flag;
}


bool myPlot::eleIDforHeepfake(int iele)
{
  double HoE           = patelehadronicOverEm[iele];
  double ecalIso       = pateledr03EcalRecHitSumEt[iele];
  double hcal1Iso      = pateledr03HcalDepth1TowerSumEt[iele];
  double hcal2Iso      = pateledr03HcalDepth2TowerSumEt[iele];
  double trkIso        = pateledr03TkSumPt[iele];
  double sieie         = patelesigmaIetaIeta[iele];
  double deta          = pateledeltaEtaSuperClusterTrackAtVtx[iele];
  double dphi          = pateledeltaPhiSuperClusterTrackAtVtx[iele];
  double e1bye5        = patelee1x5[iele]/patelee5x5[iele];
  double e2bye5        = patelee2x5Max[iele]/patelee5x5[iele];
  
  double cutHoE       = 0.05;
  
  int mishits = pateleExpectednumberOfHits[iele];
  
  TLorentzVector *elep4  = (TLorentzVector*)patelp4->At(iele);
  //double et              = elep4->Pt();


  int scind = patelescind[iele];
  TLorentzVector *sc4 = (TLorentzVector *) scp4->At(scind);
  double et              = sc4->Energy()*sin(elep4->Theta());
  double eleeta=sc4->Eta();

  bool flag = false;

  //EB
  if( (fabs(eleeta)<1.442) )
    {
      //flag = et>=35. && pateleecalDrivenSeed[iele] && fabs(deta)<0.005 && fabs(dphi)<0.09 && HoE<cutHoE && trkIso<7.5 && (ecalIso+hcal1Iso)<(2+0.03*et) && (e1bye5>0.83|| e2bye5>0.94) && mishits==0;
      flag = pateleecalDrivenSeed[iele] && fabs(deta)<0.005 && fabs(dphi)<0.06 && HoE<cutHoE && trkIso<5. && (ecalIso+hcal1Iso)<(2+0.03*et) && (e1bye5>0.83|| e2bye5>0.94) && mishits==0;
    }

  //EE
  if( fabs(eleeta)>1.56 && fabs(eleeta)<2.5 )
    {
      double ehcalcut = 9999.;
      if(et<50.) ehcalcut = 2.5;
      if(et>50.) ehcalcut = (2.5+0.03*(et-50.));

      //flag = et>=40. && pateleecalDrivenSeed[iele] && fabs(deta)<0.007 && fabs(dphi)<0.09 && HoE<cutHoE && trkIso<15. && (ecalIso+hcal1Iso)<ehcalcut && hcal2Iso<0.5 &&sieie<0.03 && mishits==0;
      flag = pateleecalDrivenSeed[iele] && fabs(deta)<0.007 && fabs(dphi)<0.06 && HoE<cutHoE && trkIso<5. && (ecalIso+hcal1Iso)<ehcalcut && /*hcal2Iso<0.5 &&*/ sieie<0.03 && mishits==0;
    }
  return flag;
}


bool myPlot::selectpi0(TLorentzVector *pho)
{
  bool passMCphoCheck = false;

  for (int genpart = 0; genpart < gensize; genpart++){
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    
    if( genpdgid[genpart]==22 && genstatus[genpart]>=1 && abs(genpdgid[genmother[genpart]]) > 100 && gen_p4->DeltaR(*(pho)) < 0.1 )
      passMCphoCheck = true;

}

return passMCphoCheck;

}//int myPlot::selectpi0()


bool myPlot::foundHardGluon(double &mygluonPt)
{
  using namespace std;

  bool passMCgluCheck = false;

  for (int genpart = 0; genpart < gensize; genpart++){
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);

    int grandmom = genmother[genmother[genpart]];

    if( genpdgid[genpart]==21 && (abs(genpdgid[genmother[genpart]])>= 1 &&abs(genpdgid[genmother[genpart]])<= 5) && abs(genpdgid[grandmom])== 2212 && gen_p4->Pt()>=0.)
      {
	passMCgluCheck = true;
	mygluonPt = gen_p4->Pt();
	//cout<<"mygluonPt "<<mygluonPt <<endl;
      }

  }

  return passMCgluCheck;

}//int myPlot::selectzgamma()


int myPlot::phomother(TLorentzVector *pho)
{
  int phoMother = -999;
  int nmoth = 0;

  for (int genpart = 0; genpart < gensize; genpart++){
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    
    if( genpdgid[genpart]==22 && genstatus[genpart]>=1 && gen_p4->DeltaR(*(pho)) < 0.1 )
      {
	phoMother=genpdgid[genmother[genpart]];
	nmoth++;
      }
    
  }//int myPlot::selectzgamma()

  /*if(nmoth>1)
    {
      //cout<<"WARNING!!!mother of photon is > 1"<<endl;
    }
  */

  return phoMother;
}



void myPlot::addbranchtoCache( vector<string> &brname){
  
  fChain->SetCacheSize(10000000);
  for(int ibr=0; ibr<(int)brname.size(); ibr++)
    fChain->AddBranchToCache(brname[ibr].c_str());
}

void myPlot::setweight(string fname, double lumi, double scale, double &weight, int &ftype)
{
  int fstr;

  fstr = fname.find("m200",0);
  if(fstr!=string::npos)
    {
      weight = (0.19719*lumi)/22000.0;
      ftype  = 2;
    }// if(fstr!=string::npos)

  fstr = fname.find("m400",0);
  if(fstr!=string::npos)
    {
      weight = (0.08089*lumi)/22000.0;
      ftype  = 3;
    }// if(fstr!=string::npos)

  fstr = fname.find("m600",0);
  if(fstr!=string::npos)
    {
      weight = (0.03694*lumi)/21753.0;
      ftype  = 4;
    }// if(fstr!=string::npos)

  fstr = fname.find("m800",0);
  if(fstr!=string::npos)
    {
      weight = (0.016651*lumi)/20200.0;
      ftype  = 5;
    }// if(fstr!=string::npos)

  fstr = fname.find("m1000",0);
  if(fstr!=string::npos)
    {
      weight = (0.007404*lumi)/22000.0;
      ftype  = 6;
    }// if(fstr!=string::npos)

  fstr = fname.find("m1200",0);
  if(fstr!=string::npos)
    {
      weight = (0.00252*lumi)/22000.0;
      ftype  = 7;
    }// if(fstr!=string::npos)
  
  fstr = fname.find("m1500",0);
  if(fstr!=string::npos)
    {
      weight = (0.000944*lumi)/19540.0;
      ftype  = 8;
    }// if(fstr!=string::npos)

  fstr = fname.find("m2000",0);
  if(fstr!=string::npos)
    {
      weight = (0.00252*lumi)/20592.0;
      ftype  = 9;
    }// if(fstr!=string::npos)


  fstr = fname.find("/zz/",0);
  if(fstr!=string::npos)
    {
      //cout<<"found ZZ"<<endl;
      //weight = (4.297*scale*lumi)/2108608.0;
      weight = (5.9*scale*lumi)/4187885.;
      ftype  = -1;
      //cout<<"weigth = "<<weight<<endl;
    }// if(fstr!=string::npos)

  fstr = fname.find("/ww/",0);
  if(fstr!=string::npos)
    {
      //weight = (27.79*lumi)/110000.0;
      weight = (43.*lumi)/4225916.0;
      ftype  = -2;
    }// if(fstr!=string::npos)

  
  fstr = fname.find("/wz/",0);
  if(fstr!=string::npos)
    {
      //weight = (10.4*scale*lumi)/110000.0;
      //weight = (18.2*scale*lumi)/4265243.0;
      weight = (18.2*scale*lumi)/4009331.0;
      ftype  = -3;
    }// if(fstr!=string::npos)

  fstr = fname.find("/wenu/",0);
  if(fstr!=string::npos)
    {
      //weight = (7899.0*lumi)/4336757.0;
      weight = ( (31314.0/3.)*lumi)/5204108.0;
      ftype  = -4;
    }// if(fstr!=string::npos)

  fstr = fname.find("/dytt/",0);
  if(fstr!=string::npos)
    {
      //weight = (1666.0*scale*lumi)/15078879.0; /// for PYTHIA sameple, 1666 is the NNLO xsec
      weight = (1627.0*scale*lumi)/19937479.0; /// for MADGRAPH sameple, 1627 is the NNLO xsec
      ftype  = -5;
    }// if(fstr!=string::npos)

  fstr = fname.find("/mytt/",0);
  if(fstr!=string::npos)
    {
      weight = (157.5*lumi)/3701947.;
      ftype  = -6;
    }// if(fstr!=string::npos)

  fstr = fname.find("/born25to250/",0);
  if(fstr!=string::npos)
    {
      weight = (22.37*lumi)/532864.0;
      ftype  = -7;
    }// if(fstr!=string::npos)
  
  fstr = fname.find("/born250toinf/",0);
  if(fstr!=string::npos)
    {
      //weight = (8.072e-03*lumi)/526240.0;
      weight = (8.072e-03*lumi)/526240.0;
      ftype  = -8;
    }// if(fstr!=string::npos)

  fstr = fname.find("/box25to250/",0);
  if(fstr!=string::npos)
    {
      weight = (12.37*lumi)/518288.0;
      ftype  = -9;
    }// if(fstr!=string::npos)

  fstr = fname.find("/box250toinf/",0);
  if(fstr!=string::npos)
    {
      weight = (2.08e-04*lumi)/514514.0;
      ftype  = -10;
    }// if(fstr!=string::npos)
  
  fstr = fname.find("/qcdbctoe20to30/",0);
  if(fstr!=string::npos)
    {
      //weight = (132160.0*scale*lumi)/1734634.0;
      weight = ((2.361e+8*5.9e-4)*scale*lumi)/1734634.0;
      ftype  = -11;
    }// if(fstr!=string::npos)                                                                                                                           \


    fstr = fname.find("/qcdbctoe30to80/",0);
    if(fstr!=string::npos)
      {
	//weight = (54980000.0*0.00230*scale*lumi)/790577.0;
	weight = ((5.944e+7*0.00242)*0.00230*scale*lumi)/790577.0;
	ftype  = -12;
      }// if(fstr!=string::npos)                                                                                                                           \



      fstr = fname.find("/qcdbctoe80to170/",0);
      if(fstr!=string::npos)
	{
	  weight = ((898200*0.0105)*0.0140*scale*lumi)/1082691.0;
	  ftype  = -13;
	}// if(fstr!=string::npos)                                                                                                                           \



	fstr = fname.find("/myzjets/",0);
	if(fstr!=string::npos)
	  {
	    //weight = (2475.*scale*lumi)/115341.;
	    weight = (2475.*scale*lumi)/34101287.;
	    ftype  = -14;
	  }// if(fstr!=string::npos)

	//////DYtoEE
	fstr = fname.find("/mydyee20/",0);
	if(fstr!=string::npos)
	  {
	    //cout<<"fname = "<<fname<<endl;
	    //weight = (1628.*scale*lumi)/2249449.0;  //POWHEG
	    weight = (1656.14*scale*lumi)/2249449.0;  //PYTHIA
	    ftype  = -15;
	  }// if(fstr!=string::npos)
	
	fstr = fname.find("/dyee120/",0);
	if(fstr!=string::npos)
	  {
	    weight = (8.65536*scale*lumi)/54550.0;
	    ftype  = -16;
	  }// if(fstr!=string::npos)      

	  fstr = fname.find("/mydyee200/",0);
	  if(fstr!=string::npos)
	    {
	      weight = (1.16685*scale*lumi)/55000.0;
	      ftype  = -17;
	    }// if(fstr!=string::npos)               
	  
	  fstr = fname.find("/dyee500/",0);
	    if(fstr!=string::npos)
	    {
	      weight = (0.0296832*scale*lumi)/54698.0;
	      ftype  = -18;
	    }// if(fstr!=string::npos)            

	    fstr = fname.find("/mydyee800/",0);
	    if(fstr!=string::npos)
	      {
		weight = (0.0041088*scale*lumi)/55000.0;
		ftype  = -19;
	      }// if(fstr!=string::npos)


	    fstr = fname.find("/vgamma/myzgamma",0);
	    if(fstr!=string::npos)
	      {
		//weight = (13.79*scale*lumi)/(62317.+62399.+62506.);
		weight = (13.79*scale*lumi)/(62317.+62399.+62506.+61801.+62130.+62682.);
		ftype  = -20;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/vgamma/mywgamma",0);
	    if(fstr!=string::npos)
	      {
		//weight = (21.41*scale*lumi)/(63639.+52079.+64807.+64687.+65017.);
		weight = (21.41*scale*lumi)/(63639.+52079.+64807.+64687.+65017.+64697.);
		ftype  = -21;
	      }// if(fstr!=string::npos)

	    fstr = fname.find("/wjets/",0);
	    if(fstr!=string::npos)
	      {
		weight = (27770.*scale*lumi)/81352581.0;
		ftype  = -22;
	      }// if(fstr!=string::npos)

	    
	    fstr = fname.find("/GJet30To50/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (16690.*scale*lumi)/2187260.0;  //PYTHIA
		ftype  = -23;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/GJet50To80/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (2722.*scale*lumi)/2036704.0;  //PYTHIA
		ftype  = -24;
	      }// if(fstr!=string::npos)


	    fstr = fname.find("/GJet80To120/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (447.2*scale*lumi)/2046637.0;  //PYTHIA
		ftype  = -25;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/GJet120To170/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (84.17*scale*lumi)/2022712.0;  //PYTHIA
		ftype  = -26;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/GJet170To300/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (22.64*scale*lumi)/2069161.0;  //PYTHIA
		ftype  = -27;
	      }// if(fstr!=string::npos)

	    fstr = fname.find("/GJet300To470/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (1.493*scale*lumi)/2076880.0;  //PYTHIA
		ftype  = -28;
	      }// if(fstr!=string::npos)

	    fstr = fname.find("/GJet470To800/",0);
	    if(fstr!=string::npos)
	      {
		//weight = (1655.91*scale*lumi)/49293.0;  //PYTHIA
		weight = (0.1323*scale*lumi)/1345504.0;  //PYTHIA
		ftype  = -29;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/myoffzgam_wgam/zgamma",0);
	    if(fstr!=string::npos)
	      {
		double kfactor = zkFactor();
		//cout<<"fname : kfactor : "<<fname<<" : "<<kfactor<<endl;
		//weight = (27.*lumi*kfactor)/323881.;
		weight = (34.16*lumi)/323881.;
		ftype  = -30;
	      }// if(fstr!=string::npos)

	    fstr = fname.find("/myoffzgam_wgam/wgamma",0);
	    if(fstr!=string::npos)
	      {
		//double kfactor = wkFactor();
		//cout<<"fname : kfactor : "<<fname<<" : "<<kfactor<<endl;
		//weight = (100.*lumi*kfactor)/524503.;
		weight = (114.7*lumi)/524503.;
		ftype  = -31;
	      }// if(fstr!=string::npos)

	    //////////POWHEG samples - added 13th jan
	    //////DYtoEE
	    fstr = fname.find("/mypowheg20_jan13/",0);
	    if(fstr!=string::npos)
	      {
		weight = (1666*scale*lumi)/16913054.0; 
		ftype  = -32;
	      }// if(fstr!=string::npos)
	    
	    fstr = fname.find("/powhdyee120/",0);
	    if(fstr!=string::npos)
	      {
		weight = (10.3*scale*lumi)/45798.0;
		ftype  = -33;
	      }// if(fstr!=string::npos)      
	    
	    fstr = fname.find("/powhdyee200/",0);
	    if(fstr!=string::npos)
	      {
		weight = (1.28*scale*lumi)/47911.;
		ftype  = -34;
	      }// if(fstr!=string::npos)               
	    
	    fstr = fname.find("/powhdyee500/",0);
	    if(fstr!=string::npos)
	      {
		weight = (0.0284*scale*lumi)/37090.;
		ftype  = -35;
	      }// if(fstr!=string::npos)            
	    
	    fstr = fname.find("/powhdyee800/",0);
	    if(fstr!=string::npos)
	      {
		weight = (0.00415*scale*lumi)/42520.;
		ftype  = -36;
	      }// if(fstr!=string::npos)
	    
	    
	    
	    


}//void myPlot::setweight(string fname, double &weight, int &ftype)


double myPlot::zkFactor()
{
  using namespace std;
  double pt = 0;
  int pidmomgammagen = 22;
  double kfactor = 1.2;
  for (int genpart = 0; genpart < gensize; genpart++){
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    if( abs(genpdgid[genpart])==22 )
      {
        pt = gen_p4->Pt();
        pidmomgammagen = genpdgid[genmother[genpart]];
        //cout<<"myPlot::kFactor: pidmomgammagen : status : "<<pidmomgammagen<<" : "<<genstatus[genpart]<<endl;
        break;
      }//if( abs(genpdgid[genpart])==22 )      
  }//for (int genpart = 0; genpart < genn; genpart++)                      
  
  if( pidmomgammagen != 23){
    if( pt > 150) {
      pt = 150;
    }
    kfactor = 1.3756 +   1.32178e-03 * pt  + 3.72292e-05 * pow(pt,2) - 2.09795e-07 * pow(pt,3);
  }//if( pidmomgammagen != 23) 
  
  //cout<<"zgamma:myPlot::kFactor: "<<kfactor<<endl;
  return kfactor;
}//double myPlot::kFactor()          


double myPlot::wkFactor()
{
  using namespace std;
  double pt = 0;
  int pidmomgammagen = 22;
  double kfactor = 1.2;
  for (int genpart = 0; genpart < gensize; genpart++){
    TLorentzVector * gen_p4 = (TLorentzVector *) genp4->At(genpart);
    if( abs(genpdgid[genpart])==22 )
      {
        pt = gen_p4->Pt();
        pidmomgammagen = genpdgid[genmother[genpart]];
        //cout<<"myPlot::kFactor: pidmomgammagen : status : "<<pidmomgammagen<<" : "<<genstatus[genpart]<<endl;
        break;
      }//if( abs(genpdgid[genpart])==22 )      
  }//for (int genpart = 0; genpart < genn; genpart++)                      
  

  if( abs(pidmomgammagen) != 24){
    if( pt > 200) {
      pt = 200;
    }
    kfactor = 1.419+0.0113*pt-0.00001011*pt*pt;
  }//if( abs(pidmomgammagen) != 24)
       
  //cout<<"wgamma:myPlot::kFactor: "<<kfactor<<endl;
  return kfactor;
}//double myPlot::kFactor()          


int myPlot::nvtx()
{

  int nvertex = 0;
  for(int iVtx=0; iVtx<vtx_std_n; ++iVtx){
    TVector3 * pos_vtx = (TVector3 *) vtx_std_xyz->At(iVtx);

    ////isFake: see in vertex.h - comments above
    //if vtx_std_x2dof==0 && vtx_std_ndof == 0
    if(vtx_std_ndof[iVtx]>=4 && fabs(pos_vtx->Z())<=24. && fabs(pos_vtx->Perp())<=2. && !(vtx_std_x2dof==0 && vtx_std_ndof == 0) )
      {
	//cout<<"found a good vertex"<<endl;                                                                                                                           
        nvertex++;
      }//if(vtx_std_ndof[iVtx]>=4 && fabs(pos_vtx->Z())<=24. && fabs(pos_vtx->Perp())<=2.)
  }//for(int iVtx=0; iVtx<vtx_std_n; ++iVtx)  
  return nvertex;

}//bool myPlot::goodvertex()                  
//////////////////////////////////////////////////////////////////////////  

void myPlot::readElebr(int ientry){
  //cout<<"inside readElebr.C"<<endl;
  b_patelesize->GetEntry(ientry);
  b_patelsc->GetEntry(ientry);  
  b_patelp4->GetEntry(ientry);  
  b_patelescind->GetEntry(ientry); 
  b_pateletrackerDrivenSeed->GetEntry(ientry);
  b_pateleecalDrivenSeed->GetEntry(ientry); 
  b_pateledr03EcalRecHitSumEt->GetEntry(ientry); 
  b_pateledr03HcalDepth1TowerSumEt->GetEntry(ientry);
  b_pateledr03HcalDepth2TowerSumEt->GetEntry(ientry);
  b_pateledr03HcalTowerSumEt->GetEntry(ientry); 
    b_patelee1x5->GetEntry(ientry); 
  b_patelee2x5Max->GetEntry(ientry); 
  b_patelee5x5->GetEntry(ientry);
  b_patelehadronicOverEm->GetEntry(ientry);
  b_patelesigmaIetaIeta->GetEntry(ientry);
  b_pateledeltaEtaSuperClusterTrackAtVtx->GetEntry(ientry);
  b_pateledeltaPhiSuperClusterTrackAtVtx->GetEntry(ientry);
  b_pateledr03TkSumPt->GetEntry(ientry);
  b_patelee2e9->GetEntry(ientry);
    b_pateleisEB->GetEntry(ientry);
  b_pateleisEBEEGap->GetEntry(ientry);
  b_pateleisEBEtaGap->GetEntry(ientry);
  b_pateleisEBGap->GetEntry(ientry);
  b_pateleisEBPhiGap->GetEntry(ientry); 
  b_pateleisEE->GetEntry(ientry); 
  b_pateleisEEDeeGap->GetEntry(ientry); 
  b_pateleisEEGap->GetEntry(ientry);
  b_pateleisEERingGap->GetEntry(ientry); 
   b_pateleExpectednumberOfHits->GetEntry(ientry); 
 b_pateleconvDist->GetEntry(ientry); 
 b_pateleconvDcot->GetEntry(ientry); 

 /*  b_patelmomvtx->GetEntry(ientry);
  b_patelmomvtxconst->GetEntry(ientry);
  b_patelmomcalo->GetEntry(ientry); 
  b_patelmomout->GetEntry(ientry); 
  b_patelposvtx->GetEntry(ientry); 
  b_patelposcalo->GetEntry(ientry);
  b_patelefbrem->GetEntry(ientry); 
  b_patelecharge->GetEntry(ientry); 
   b_pateledr04EcalRecHitSumEt->GetEntry(ientry);
  b_pateledr04HcalDepth1TowerSumEt->GetEntry(ientry);
  b_pateledr04HcalDepth2TowerSumEt->GetEntry(ientry);
  b_pateledr04HcalTowerSumEt->GetEntry(ientry); 
  b_pateleeEleClusterOverPout->GetEntry(ientry);
  b_pateleeSeedClusterOverP->GetEntry(ientry); 
  b_pateleeSeedClusterOverPout->GetEntry(ientry);
  b_pateleeSuperClusterOverP->GetEntry(ientry);
    b_pateledeltaEtaEleClusterTrackAtCalo->GetEntry(ientry); 
  b_pateledeltaEtaSeedClusterTrackAtCalo->GetEntry(ientry);
  b_pateledeltaPhiEleClusterTrackAtCalo->GetEntry(ientry); 
  b_pateledeltaPhiSeedClusterTrackAtCalo->GetEntry(ientry);
  b_pateledeltaPhiSuperClusterTrackAtVtx->GetEntry(ientry);
  b_patelemva->GetEntry(ientry); 
  b_patelenumberOfTracks->GetEntry(ientry); 
  
  b_pateledr04TkSumPt->GetEntry(ientry);
  b_patelepin->GetEntry(ientry); 
  b_patelepout->GetEntry(ientry);
  b_patelepfeta->GetEntry(ientry);
  b_patelepfphi->GetEntry(ientry);
  b_patelepfE->GetEntry(ientry);
  b_patelemaxEnergyXtal->GetEntry(ientry);
  b_pateleswissCross->GetEntry(ientry); 
  b_pateleswissBasedspikevar->GetEntry(ientry);
  b_pateleseedtime->GetEntry(ientry); 
  b_patelerecoFlag->GetEntry(ientry); 
  //b_pateleseverityLevel->GetEntry(ientry);
  
  b_pateletrkpt->GetEntry(ientry); 
  b_pateletrkcharge->GetEntry(ientry); 
  b_pateletrkchi2->GetEntry(ientry);
  b_pateletrketa->GetEntry(ientry); 
  b_pateletrknumberOfLostHits->GetEntry(ientry); 
  b_pateletrknumberOfValidHits->GetEntry(ientry);
  b_pateletrklost->GetEntry(ientry); 
  b_pateletrkd0->GetEntry(ientry);
  b_pateletrkdxy->GetEntry(ientry);
  b_pateletrkdz->GetEntry(ientry); 
  b_pateletrkptin->GetEntry(ientry); 
  b_pateletrkptout->GetEntry(ientry);
b_pateletrkfbrem->GetEntry(ientry);
 b_pateletrkqoverp->GetEntry(ientry); 
 b_pateletrkvx->GetEntry(ientry); 
 b_pateletrkvy->GetEntry(ientry); 
 b_pateletrkvz->GetEntry(ientry); 
b_pateletrkphi->GetEntry(ientry); 
 b_pateletrkndof->GetEntry(ientry); 
 b_pateletrkrecHitsSize->GetEntry(ientry);
 b_pateletrktheta->GetEntry(ientry);
 b_pateletrkqualityMask->GetEntry(ientry); 
b_pateletrkouterX->GetEntry(ientry); 
 b_pateletrkouterY->GetEntry(ientry);
 b_pateletrkouterZ->GetEntry(ientry);
 b_pateletrkouterRadius->GetEntry(ientry);
 b_pateletrkinnerX->GetEntry(ientry); 
 b_pateletrkinnerY->GetEntry(ientry); 
 b_pateletrkinnerZ->GetEntry(ientry); 
b_pateleheepcutword->GetEntry(ientry);
 b_pateleheepid->GetEntry(ientry); 
 */

 //cout<<"patelesize:"<<patelesize<<endl;
return;
}


void myPlot::readPhobr(int ientry){
  b_patphosize->GetEntry(ientry);
  b_patphop4->GetEntry(ientry);  
  b_patphocalopos->GetEntry(ientry);
  b_patphoscp4->GetEntry(ientry); 
  b_patphoscind->GetEntry(ientry);
  b_patphoecalRecHitSumEtConeDR04->GetEntry(ientry);
  b_patphohcalDepth1TowerSumEtConeDR04->GetEntry(ientry);
  b_patphohcalDepth2TowerSumEtConeDR04->GetEntry(ientry);
  b_patphohcalTowerSumEtConeDR04->GetEntry(ientry);
  b_patphotrkSumPtHollowConeDR04->GetEntry(ientry);
  b_patphotrkSumPtSolidConeDR04->GetEntry(ientry); 
  b_patphonTrkHollowConeDR04->GetEntry(ientry); 
  b_patphonTrkSolidConeDR04->GetEntry(ientry); 
  b_patphohadronicOverEm->GetEntry(ientry); 
  b_patphosigmaIetaIeta->GetEntry(ientry);  
  b_patphohasPixelSeed->GetEntry(ientry);
  b_patphoeta->GetEntry(ientry);  

  /*b_patphoe1x5->GetEntry(ientry); 
  b_patphoe2x5->GetEntry(ientry); 
  b_patphoe5x5->GetEntry(ientry); 
  b_patphonumberOfTracks->GetEntry(ientry); 
  b_patphoecalRecHitSumEtConeDR03->GetEntry(ientry);  
  b_patphohcalDepth1TowerSumEtConeDR03->GetEntry(ientry);
  b_patphohcalDepth2TowerSumEtConeDR03->GetEntry(ientry);
  b_patphohcalTowerSumEtConeDR03->GetEntry(ientry);
  b_patphotrkSumPtHollowConeDR03->GetEntry(ientry);
  b_patphotrkSumPtSolidConeDR03->GetEntry(ientry); 
  b_patphonTrkHollowConeDR03->GetEntry(ientry);
  b_patphonTrkSolidConeDR03->GetEntry(ientry); 
  b_patphoisConvertedPhoton->GetEntry(ientry); 
  b_patphomaxEnergyXtal->GetEntry(ientry); 
  b_patphoTightIDcutword->GetEntry(ientry);
  b_patphotightid->GetEntry(ientry);
  b_patphoisEB->GetEntry(ientry);
  b_patphoisEBEEGap->GetEntry(ientry);
  b_patphoisEBEtaGap->GetEntry(ientry);
  b_patphoisEBGap->GetEntry(ientry); 
  b_patphoisEBPhiGap->GetEntry(ientry); 
  b_patphoisEE->GetEntry(ientry);
  b_patphoisEEDeeGap->GetEntry(ientry);
  b_patphoisEEGap->GetEntry(ientry); 
  b_patphoisEERingGap->GetEntry(ientry);
  b_patphoswissCross->GetEntry(ientry); 
  b_patphoswissBasedspikevar->GetEntry(ientry); 
  b_patphoseedtime->GetEntry(ientry); 
  b_patphorecoFlag->GetEntry(ientry); 
  //b_patphoseverityLevel->GetEntry(ientry);
  b_patphoe2e9->GetEntry(ientry);  
  b_patphoconvsize->GetEntry(ientry); 
  b_patphohasConversionTracks->GetEntry(ientry);
  b_patphoconvtxX->GetEntry(ientry); 
  b_patphoconvtxY->GetEntry(ientry); 
  b_patphoconvtxZ->GetEntry(ientry); 
  b_patphoconvtxR->GetEntry(ientry); 
  b_patphotrkpt->GetEntry(ientry);  
  b_patphotrkcharge->GetEntry(ientry); 
  b_patphotrkchi2->GetEntry(ientry); 
  b_patphotrketa->GetEntry(ientry);  
  b_patphotrknumberOfLostHits->GetEntry(ientry); 
  b_patphotrknumberOfValidHits->GetEntry(ientry);
  b_patphotrklost->GetEntry(ientry); 
  b_patphotrkd0->GetEntry(ientry);   
  b_patphotrkdxy->GetEntry(ientry);  
  b_patphotrkdz->GetEntry(ientry);   
  b_patphotrkptin->GetEntry(ientry); 
  b_patphotrkptout->GetEntry(ientry);
  b_patphotrkfbrem->GetEntry(ientry);
  b_patphotrkqoverp->GetEntry(ientry);
  b_patphotrkvx->GetEntry(ientry); 
  b_patphotrkvy->GetEntry(ientry); 
  b_patphotrkvz->GetEntry(ientry); 
  b_patphotrkphi->GetEntry(ientry);
  b_patphotrkndof->GetEntry(ientry);
  b_patphotrkrecHitsSize->GetEntry(ientry);
  b_patphotrktheta->GetEntry(ientry);
  b_patphotrkqualityMask->GetEntry(ientry);
  b_patphotrkouterX->GetEntry(ientry); 
  b_patphotrkouterY->GetEntry(ientry); 
  b_patphotrkouterZ->GetEntry(ientry); 
  b_patphotrkouterRadius->GetEntry(ientry);
  b_patphotrkinnerX->GetEntry(ientry); 
  b_patphotrkinnerY->GetEntry(ientry); 
  b_patphotrkinnerZ->GetEntry(ientry); 
  b_patphomGenisJet->GetEntry(ientry); 
  b_patphomGenisPhoton->GetEntry(ientry); 
  b_patphomGenisElectron->GetEntry(ientry); 
  b_patphomGenpdgId->GetEntry(ientry); 
  b_patphomGenstatus->GetEntry(ientry);
  b_patphonummoth->GetEntry(ientry); 
  b_patphomGenmompdgId->GetEntry(ientry);
  b_patphomGengranpdgId->GetEntry(ientry);
  b_patphomGentheta->GetEntry(ientry); 
  b_patphomGeneta->GetEntry(ientry); 
  b_patphomGenphi->GetEntry(ientry); 
  b_patphomGenpt->GetEntry(ientry);  
  b_patphomGenpx->GetEntry(ientry);  
  b_patphomGenpy->GetEntry(ientry);  
  b_patphomGenpz->GetEntry(ientry);
  b_patphomGenenergy->GetEntry(ientry);
  */  
  return;
}


void myPlot::readHLTbr(int ientry){
b_hlt1_bit_1->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit_2->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit_3->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit_4->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit_5->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit_6->GetEntry(ientry);   //!                                                                                                                                 
b_hlt1_bit->GetEntry(ientry);   //!                                                                                                                                   
b_hlt_path_names_HLT1_1->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1_2->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1_3->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1_4->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1_5->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1_6->GetEntry(ientry);   //!                                                                                                                      
b_hlt_path_names_HLT1->GetEntry(ientry);   //!                                                                                                                        
b_hlt_n->GetEntry(ientry);   //!                                                                                                                                      
b_hlt_p4->GetEntry(ientry);   //!                     
b_hlt_candpath_1->GetEntry(ientry);   //!                                                                                                                             
b_hlt_candpath_2->GetEntry(ientry);   //!                                                                                                                             
b_hlt_candpath_3->GetEntry(ientry);   //!                                                                                                                             
b_hlt_candpath_4->GetEntry(ientry);   //!                                                                                                                             
b_hlt_candpath_5->GetEntry(ientry);   //!                                                                                                                             
b_hlt_candpath_6->GetEntry(ientry);   //!                                                                                                                             
//b_hlt_candpath->GetEntry(ientry);   //!                                                                                                                             
b_hlt_label_names_1->GetEntry(ientry);   //!                                                                                                                          
b_hlt_label_names_2->GetEntry(ientry);   //!                                                                                                                          
b_hlt_label_names_3->GetEntry(ientry);   //!  
b_hlt_label_names_4->GetEntry(ientry);   //!                                                                                                                          
b_hlt_label_names_5->GetEntry(ientry);   //!                                                                                                                          
b_hlt_label_names_6->GetEntry(ientry);   //!   

 return;
}


void myPlot::readScbr(int ientry){
  b_scsize->GetEntry(ientry);
  b_scp4->GetEntry(ientry);
  b_sce2e9->GetEntry(ientry);
  b_sce3x3->GetEntry(ientry);
  b_scrawEnergy->GetEntry(ientry);
  
  return;
}


void myPlot::readTrkVtxbr(int ientry){
  
  b_vtx_std_n->GetEntry(ientry);   //!  
  b_vtx_std_xyz->GetEntry(ientry);   //!
  b_vtx_std_x2dof->GetEntry(ientry);   //!                                                                                                                              
  b_vtx_std_ndof->GetEntry(ientry);   //! 
  b_vtx_std_ntks->GetEntry(ientry);
  b_gentrksize->GetEntry(ientry);
  b_gentrkqualityMask->GetEntry(ientry);
  
  
  return;
}

void myPlot::readGenbr(int ientry){
  b_gensize->GetEntry(ientry);   //!                                             
  b_genp4->GetEntry(ientry);   //!                                               
  b_genvtx->GetEntry(ientry);   //!                                              
  b_genstatus->GetEntry(ientry);   //!                                           
  b_gencharge->GetEntry(ientry);   //!                                           
  b_genpdgid->GetEntry(ientry);   //!                                            
  b_genmother->GetEntry(ientry);   //!                                           
  b_genndau->GetEntry(ientry);   //!                                             
  b_gennmoth->GetEntry(ientry);   //!  

}
