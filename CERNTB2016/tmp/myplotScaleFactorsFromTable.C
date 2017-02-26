#include "TMath.h"
#include <iostream>

double err( double n, double d, double err_n, double err_d){

  double err = (n/d)*sqrt( pow(err_n/n,2) + pow(err_d/d, 2) );

  return err;

}


void setTCanvasNicev1(TCanvas *can0){
  can0->SetFillColor(0);
  can0->SetBorderMode(0);
  can0->SetBorderSize(2);
  can0->SetTickx(1);
  can0->SetTicky(1);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
  return;
}



void savePlot(TH2F *h, string name, bool plotErr){

  TCanvas *c = new TCanvas("c","c",600,600);
  h->SetMinimum(-0.001);
  h->SetMaximum(1.001);
  //h->SetEntries(16);
  //  h->SetDirectory(0);
  h->SetStats(0);
  //h->SetOption("colztexte");
  //TCanvas *c = new TCanvas("c","c",600,600);
  if(plotErr) h->Draw("colztexte");
  else h->Draw("colztext");

  h->SetContour(20);
  h->SetContourLevel(0,-0.001);
  h->SetContourLevel(1,0.0491);
  h->SetContourLevel(2,0.0992);
  h->SetContourLevel(3,0.1493);
  h->SetContourLevel(4,0.1994);
  h->SetContourLevel(5,0.2495);
  h->SetContourLevel(6,0.2996);
  h->SetContourLevel(7,0.3497);
  h->SetContourLevel(8,0.3998);
  h->SetContourLevel(9,0.4499);
  h->SetContourLevel(10,0.5);
  h->SetContourLevel(11,0.5501);
  h->SetContourLevel(12,0.6002);
  h->SetContourLevel(13,0.6503);
  h->SetContourLevel(14,0.7004);
  h->SetContourLevel(15,0.7505);
  h->SetContourLevel(16,0.8006);
  h->SetContourLevel(17,0.8507);
  h->SetContourLevel(18,0.9008);
  h->SetContourLevel(19,0.9509);

  

  //h->Draw("colztexte");
  char filename[100];
  sprintf(filename,"%s.gif",name.c_str());
  c->Print(filename);

  sprintf(filename,"%s.C",name.c_str());
  c->Print(filename);




}


////save canvas of many histograms 
void savePlot(map<string, TH1F*> h, int startbin, int endbin, string name, string legname[], string plotname){


  int W = 800;
  int H = 600;

  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R                                                                                                                                                     
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TCanvas* c = new TCanvas("c","c",50,50,W,H);

  c->SetLeftMargin( L/W );
  c->SetRightMargin( R/W );
  c->SetTopMargin( T/H );
  c->SetBottomMargin( B/H );
  
  gStyle->SetOptStat(0);
  

  setTCanvasNicev1(c);

  //TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
  //TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.198,22);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0,21);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.3,22);
  pad1->SetFillColor(0);
  pad2->SetFillColor(0);


  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetGrid();
  if(name=="etabin") pad1->SetLogx();

  gStyle->SetOptStat(0);
  //h->SetMinimum(-0.001);
  //h->SetMaximum(1.001);

  float maxf = -999;
  float minf = 999;
  for(int i=startbin; i<=endbin; i++){
  
    maxf = std::max(maxf,(float)h[Form("%s%d",name.c_str(),i)]->GetMaximum());
    minf = std::min(minf,(float)h[Form("%s%d",name.c_str(),i)]->GetMinimum());
  }

  
  maxf = 1.2;
  minf = 0.4;
  
  TLegend *leg = new TLegend(0.6080402,0.7125436,0.8994975,0.8954704,NULL,"brNDC");
  //TLegend *leg = new TLegend(0.5752508,0.3704348,0.8662207,0.5530435,NULL,"brNDC");
  
  for(int i=startbin; i<=endbin; i++){
    
    if((i==1 || i==6) && name=="etabin") continue; //gap region

    h[Form("dataeff_%s%d",name.c_str(),i)]->SetName("");
    h[Form("dataeff_%s%d",name.c_str(),i)]->SetTitle("");

    h[Form("dataeff_%s%d",name.c_str(),i)]->SetMaximum(maxf);
    h[Form("dataeff_%s%d",name.c_str(),i)]->SetMinimum(minf);
    if(i==0) h[Form("dataeff_%s%d",name.c_str(),i)]->Draw();
    else h[Form("dataeff_%s%d",name.c_str(),i)]->Draw("same");
 
    h[Form("mceff_%s%d",name.c_str(),i)]->SetName("");
    h[Form("mceff_%s%d",name.c_str(),i)]->SetTitle("");

    h[Form("mceff_%s%d",name.c_str(),i)]->SetMaximum(maxf);
    h[Form("mceff_%s%d",name.c_str(),i)]->SetMinimum(minf);

    h[Form("mceff_%s%d",name.c_str(),i)]->Draw("same");
   
    leg->AddEntry(h[Form("dataeff_%s%d",name.c_str(),i)],(legname[i]+" (Data)").c_str(),"lp");
    leg->AddEntry(h[Form("mceff_%s%d",name.c_str(),i)],(legname[i]+" (MC)").c_str(),"lp");

    //cout<<"leg name is "<<legname[i]<<endl;
  }


  

  leg->Draw();
  pad1->Modified();
  pad1->Update();

  pad2->cd();
  if(name=="etabin") pad2->SetLogx();
  pad2->SetGrid();
  for(int i=startbin; i<=endbin; i++){
    
    if((i==1 || i==6) && name=="etabin") continue; //gap region

    h[Form("%s%d",name.c_str(),i)]->SetName("");
    h[Form("%s%d",name.c_str(),i)]->SetTitle("");

    h[Form("%s%d",name.c_str(),i)]->SetMaximum(1.2);

    h[Form("%s%d",name.c_str(),i)]->SetMinimum(0.8);
    h[Form("%s%d",name.c_str(),i)]->GetYaxis()->SetTitle("Scale Factor");
    h[Form("%s%d",name.c_str(),i)]->GetYaxis()->SetTitleSize(0.08);
    h[Form("%s%d",name.c_str(),i)]->GetXaxis()->SetTitleSize(0.08);
    h[Form("%s%d",name.c_str(),i)]->GetYaxis()->SetTitleOffset(.5);
    h[Form("%s%d",name.c_str(),i)]->GetXaxis()->SetTitleOffset(.5);

    h[Form("%s%d",name.c_str(),i)]->GetYaxis()->SetLabelSize(0.08);
    h[Form("%s%d",name.c_str(),i)]->GetXaxis()->SetLabelSize(0.06);

    if(i==0) h[Form("%s%d",name.c_str(),i)]->Draw();
    else h[Form("%s%d",name.c_str(),i)]->Draw("same");

  }


  //h->Draw("colztexte");
  char filename[100];
  sprintf(filename,"plots/%s.gif",plotname.c_str());
  c->Print(filename);

  sprintf(filename,"plots/%s.C",plotname.c_str());
  c->Print(filename);




}



void myplotScaleFactorsFromTable(int type = 0){

  string file;
  
  if(type==0) file="Loose";
  if(type==1) file="Medium";
  if(type==2) file="Tight";
  if(type==3) file="MVA";
  if(type==4) file="Presel";

  ifstream datafile;
  //datafile.open( "data_"+file+"_GJ.txt");
  //datafile.open( "data_"+file+".txt");
  datafile.open( "passing"+file+".txt");
  
  ifstream mcfile;
  //mcfile.open("MCtruth_"+file+".txt");
  mcfile.open("passing"+file+"_MC.txt");
  //mcfile.open("data_"+file+"_SJ.txt");

  cout<<"MC file is "<<("passing"+file+"_MC.txt")<<endl;
  TFile *fout = new TFile( (file+"_sf.root").c_str(), "RECREATE");

  //double sf[100], sferr[100], ptmin[100], ptmax[100], etamin[100], etamax[100];


  /*int nptbins=0;
  int netabins = 0;
  while(!mcfile.eof()){
    
    double mcetamin, mcetamax, mcptmin, mcptmax, mceff, mcerr;
    mcfile>>mcetamin>>mcetamax>>mcptmin>>mcptmax>>mceff>>mcerr;
    
    if(nptbins==0){
      ptmin[nptbins] = mcptmin;
      ptmax[nptbins] = mcptmax;
      nptbins++;
    }
    
    if(netabins==0){
      etamin[netabins] = mcetamin;
      etamax[netabins] = mcetamax;
      netabins++;
    }

    if(nptbins!=0){
      if(mcptmin!=ptmin[nptbins-1]){
	ptmin[nptbins] = mcptmin;
	ptmax[nptbins] = mcptmax;
	nptbins++;
      }
    }

    if(netabins!=0){
      if(mcetamin!=etamin[netabins-1]){
      etamin[netabins] = mcetamin;
      etamax[netabins] = mcetamax;
      netabins++;
      }
    }
    
    
  }//while(!mcfile.eof())
    
  ptmin[nptbins+1] = ptmax[nptbins];
  etamin[netabins+1] = etamax[netabins];
  
  mcfile.close();

  for(int i=0; i<nptbins+1; i++){
    cout<<"pt["<<i<<"] "<<ptmin[i]<<endl;
  }

  for(int i=0; i<netabins+1; i++){
    cout<<"eta["<<i<<"] "<<etamin[i]<<endl;
  }


  

  TH2F *h = new TH2F("scalefactors","",nptbins+1,ptmin,netabins+1,etamin);
  mcfile.open("MCtruth_"+file+".txt");
  */

  /*
  const int nptbins = 3;
  const int netabins = 4;
  double ptmin[nptbins+1] = {25.,40,50,200};
  double etamin[netabins+1] = {0,1.,1.4,1.6,2.5};
  */

  ///if change here, then change in the function below: printSystematics
  const int nptbins = 4;
  const int netabins = 8;
  double ptmin[nptbins+1] = {20.,30.,40,50,200};
  double etamin[netabins+1] = {-2.5,-1.566,-1.4442,-1,0,1.,1.4442,1.566,2.5};

  
  //TH2F *h = new TH2F("scalefactors","",nptbins,ptmin,netabins,etamin);
  TH2F *h = new TH2F( "sf","",netabins,etamin,nptbins,ptmin);
  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("#p_{T} [GeV]");


  TH2F *hdeff = new TH2F( "hdeff","",netabins,etamin,nptbins,ptmin);
  hdeff->GetXaxis()->SetTitle("#eta");
  hdeff->GetYaxis()->SetTitle("#p_{T} [GeV]");

  TH2F *hmeff = new TH2F( "hmeff","",netabins,etamin,nptbins,ptmin);
  hmeff->GetXaxis()->SetTitle("#eta");
  hmeff->GetYaxis()->SetTitle("#p_{T} [GeV]");

  while(!datafile.eof()){

    double detamin, detamax, dptmin, dptmax, deff, derr;
    double mcetamin, mcetamax, mcptmin, mcptmax, mceff, mcerr;

    //datafile>>detamin>>detamax>>dptmin>>dptmax>>deff>>derr;
    datafile>>dptmin>>dptmax>>detamin>>detamax>>deff>>derr;
    cout<<"DATA === "<<detamin<<"  "<<detamax<<"  "<<dptmin<<"  "<<dptmax<<"  "<<deff<<"  "<<derr<<endl;

    //mcfile>>mcetamin>>mcetamax>>mcptmin>>mcptmax>>mceff>>mcerr;
    mcfile>>mcptmin>>mcptmax>>mcetamin>>mcetamax>>mceff>>mcerr;
    cout<<"MC===="<<mcetamin<<"  "<<mcetamax<<"  "<<mcptmin<<"  "<<mcptmax<<"  "<<mceff<<"  "<<mcerr<<endl;
    
    double scaleFactor = deff/mceff;
    double scaleFactorerr = err(deff, mceff, derr, mcerr);

    cout<<"detamin : detamax "<<fabs(detamin)<<" : "<<fabs(detamax)<<endl;
    if( (detamin>=1.4442 && detamax<=1.566) || (detamin>=-1.566 && detamax<=-1.4442) ) {
      cout<<"==============INSIDE=============="<<endl;
      cout<<detamin<<"  "<<detamax<<"  "<<dptmin<<"  "<<dptmax<<"  "<<deff<<"  "<<derr<<endl;
      scaleFactor = 0;
      deff = 0;
      mceff = 0;
      scaleFactorerr = 0;
      derr = 0;
      mcerr = 0;
    }
    

    int iptbin = h->GetYaxis()->FindBin(dptmin);
    int ietabin = h->GetXaxis()->FindBin(detamin);
    
    cout<<""<<endl;
    cout<<iptbin<<" " <<ietabin<<endl;
    hdeff->SetBinContent(ietabin,iptbin,deff);
    hdeff->SetBinError(ietabin,iptbin,derr);

    hmeff->SetBinContent(ietabin,iptbin,mceff);
    hmeff->SetBinError(ietabin,iptbin,mcerr);

    fout->cd();
    h->SetBinContent(ietabin,iptbin,scaleFactor);
    h->SetBinError(ietabin,iptbin,scaleFactorerr);

    cout<<"eta bin ; etbin : "<<mcetamin<<" "<<mcetamax<<" "<<mcptmin<<" "<<mcptmax<<endl;

    cout<<"Dataeff : MCeff : "<<deff<<" : "<<mceff<<endl;
    cout<<"SF is "<<scaleFactor<<"+/-"<<scaleFactorerr<<endl;
  }


  savePlot(hdeff, ("dataeff_"+file).c_str(),true );
  savePlot(hmeff, ("mceff_"+file).c_str(),true );
  savePlot(h,file, true);
  
  
  cout<<"Now writing the latex table"<<endl;
  
  int nbinsx = h->GetNbinsX();
  cout<<"nbinsx "<<nbinsx<<endl;
  
  int nbinsy = h->GetNbinsY();
  cout<<"nbinsy "<<nbinsy<<endl;
  
  
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\label{tab:EvtYields2010}"<<endl;
  cout<<"\\end{table}"<<endl;                     
  
  cout<<""<<endl;
  cout<<""<<endl;
  cout<<""<<endl;
  //////////////////////try to print in a table format///////////////
  //1. eff and scale factors//////////////////////////
  cout<<"\\begin{table}[htp]"<<endl;
  cout<<"\\begin{center}"<<endl;
  cout<<"\\begin{tabular}{|c|c|c|c|c|c|c|}"<<endl;
  cout<<"\\hline"<<endl;
  //cout<<"$\\eta$ & $\\p_T$ & baseline SF & MINOS err. & sys: BKG fit   & sys: SIG fit & sys: fit range & total error\\\\"<<endl;
  cout<<"$\\eta$ & p$_T$ & Eff in Data & Eff in MC & SF $\\pm$ toterrr\\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){

      double xlow = h->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = h->GetXaxis()->GetBinUpEdge(ix);
      
      double ylow = h->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = h->GetYaxis()->GetBinUpEdge(iy);

      double binc = h->GetBinContent(ix,iy);
      double binerr = h->GetBinError(ix,iy);
      
      double dataeff_base = hdeff->GetBinContent(ix,iy);
      double mceff_base = hmeff->GetBinContent(ix,iy);

      if(xlow>=1.4 && xhigh<=1.6) {
	continue;
      }
      
      cout<<xlow<<"-"<<xhigh<<" & "<<ylow<<"-"<<yhigh
	//<<binc<<" & "<<binerr<<" & "<<sys_bkg<<" & "<<sys_sig<<" & "<<sys_fit<<" & "<<toterr<<"\\\\"<<endl;
	  <<" & "<<dataeff_base<<" & "<<mceff_base<<" & "<<binc<<" $\\pm$ "<<binerr<<"\\\\"<<endl;
      cout<<"\\hline"<<endl;
      
    }
  }

  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\label{tab:sf_presel}"<<endl;
  cout<<"\\end{table}"<<endl;                     

  fout->Write();
  fout->Close();
  

}



//void printSystematics(string baseline_dir="sf_baseline", string bkgsys_dir="sf_sysbkgPOLY", string sigsys_dir="sf_sigsys", string fitsys_dir = "sf_FITrange", 
void printSystematics(string baseline_dir="sf_baseline", string bkgsys_dir="sf_sysbkgPOLY", string sigsys_dir="sf_sigsys", 
		      string tagpt_dir = "sf_tagpT35", string tagmedid_dir = "sf_tagMedID", string gen_dir = "sf_mcATnlo", string pud_dir = "sf_puDown", 
		      string puu_dir = "sf_puUp",
		      int type=0){
  
  string file;
  if(type==0) file="Loose";
  if(type==1) file="Medium";
  if(type==2) file="Tight";
  if(type==3) file="MVA";
  if(type==4) file="Presel";

  TFile *fin_base   = TFile::Open((baseline_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_bkgsys = TFile::Open((bkgsys_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_sigsys = TFile::Open((sigsys_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_fitsys = TFile::Open((fitsys_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_tagptsys = TFile::Open((tagpt_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_tagidsys = TFile::Open((tagmedid_dir+"/"+file+"_sf.root").c_str());

  TFile *fin_gensys = TFile::Open((gen_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_pudsys = TFile::Open((pud_dir+"/"+file+"_sf.root").c_str());
  TFile *fin_puusys = TFile::Open((puu_dir+"/"+file+"_sf.root").c_str());


  TH2F *hdeff_base = (TH2F*)fin_base->Get("hdeff");
  TH2F *hdeff_bkgsys = (TH2F*)fin_bkgsys->Get("hdeff");
  TH2F *hdeff_sigsys = (TH2F*)fin_sigsys->Get("hdeff");
  TH2F *hdeff_fitsys = (TH2F*)fin_fitsys->Get("hdeff");
  TH2F *hdeff_ptsys = (TH2F*)fin_tagptsys->Get("hdeff");
  TH2F *hdeff_idsys = (TH2F*)fin_tagidsys->Get("hdeff");

  TH2F *hdeff_gensys = (TH2F*)fin_gensys->Get("hdeff");
  TH2F *hdeff_pudsys = (TH2F*)fin_pudsys->Get("hdeff");
  TH2F *hdeff_puusys = (TH2F*)fin_puusys->Get("hdeff");

  TH2F *hmeff_base = (TH2F*)fin_base->Get("hmeff");
  TH2F *hmeff_bkgsys = (TH2F*)fin_bkgsys->Get("hmeff");
  TH2F *hmeff_sigsys = (TH2F*)fin_sigsys->Get("hmeff");
  TH2F *hmeff_fitsys = (TH2F*)fin_fitsys->Get("hmeff");
  TH2F *hmeff_ptsys = (TH2F*)fin_tagptsys->Get("hmeff");
  TH2F *hmeff_idsys = (TH2F*)fin_tagidsys->Get("hmeff");

  TH2F *hmeff_gensys = (TH2F*)fin_gensys->Get("hmeff");
  TH2F *hmeff_pudsys = (TH2F*)fin_pudsys->Get("hmeff");
  TH2F *hmeff_puusys = (TH2F*)fin_puusys->Get("hmeff");


  TH2F *hsf_base = (TH2F*)fin_base->Get("sf");
  TH2F *hsf_bkgsys = (TH2F*)fin_bkgsys->Get("sf");
  TH2F *hsf_sigsys = (TH2F*)fin_sigsys->Get("sf");
  TH2F *hsf_fitsys = (TH2F*)fin_fitsys->Get("sf");
  TH2F *hsf_ptsys = (TH2F*)fin_tagptsys->Get("sf");
  TH2F *hsf_idsys = (TH2F*)fin_tagidsys->Get("sf");

  TH2F *hsf_gensys = (TH2F*)fin_gensys->Get("sf");
  TH2F *hsf_pudsys = (TH2F*)fin_pudsys->Get("sf");
  TH2F *hsf_puusys = (TH2F*)fin_puusys->Get("sf");

  
  int nbinsx = hsf_base->GetNbinsX();
  int nbinsy = hsf_base->GetNbinsY();

  //TFile *fout = new TFile( (file+"_withsys_sf.root").c_str(), "RECREATE");
  TH2F *hsf_fullsys = (TH2F*)hsf_base->Clone();
  
  TH2F *h_stat = (TH2F*)hsf_base->Clone();
  TH2F *h_bkgsys = (TH2F*)hsf_base->Clone();
  TH2F *h_sigsys = (TH2F*)hsf_base->Clone();
  TH2F *h_fitsys = (TH2F*)hsf_base->Clone();
  TH2F *h_ptsys = (TH2F*)hsf_base->Clone();
  TH2F *h_idsys = (TH2F*)hsf_base->Clone();

  TH2F *h_gensys = (TH2F*)hsf_base->Clone();
  TH2F *h_pudsys = (TH2F*)hsf_base->Clone();
  TH2F *h_puusys = (TH2F*)hsf_base->Clone();
  TH2F *h_pusys = (TH2F*)hsf_base->Clone();

  TFile *ffinal = new TFile("ffinal.root","RECREATE");
  TH2F *hsf_final = (TH2F*)hsf_base->Clone();
  hsf_final->SetName("hsf_final");

  map<string,TH1F*>hsf_func_eta;
  map<string,TH1F*>hsf_func_pt;

  const int nptbins = 4; //==nbinsy
  const int netabins = 8; //==nbinsx
  double ptmin[nptbins+1] = {20.,30.,40,50,200};
  double etamin[netabins+1] = {-2.5,-1.566,-1.4442,-1,0,1.,1.4442,1.566,2.5};

  int col[9] = {1,2,4,6,8,9,40,16,28};
  int sty[9] = {20,21,22,23,24,25,26,27,28};

  

  for(int i=0; i<nptbins; i++){

    hsf_func_eta[Form("ptbin%d",i)] = new TH1F(Form("ptbin%d",i),Form("ptbin%d",i),netabins,etamin);
    hsf_func_eta[Form("ptbin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_eta[Form("ptbin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_eta[Form("ptbin%d",i)]->SetLineColor(col[i]);
    hsf_func_eta[Form("ptbin%d",i)]->SetLineWidth(3);
    hsf_func_eta[Form("ptbin%d",i)]->GetXaxis()->SetTitle("#eta");
    hsf_func_eta[Form("ptbin%d",i)]->GetYaxis()->SetTitle("Scale Factor");

    //hsf_func_eta[Form("ptbin%d",i)]->SetMaximum(2);
    //hsf_func_eta[Form("ptbin%d",i)]->SetMinimum(1);

    ///data Eff
    hsf_func_eta[Form("dataeff_ptbin%d",i)] = new TH1F(Form("dataeff_ptbin%d",i),Form("dataeff_ptbin%d",i),netabins,etamin);
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->SetLineColor(col[i]);
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->SetLineWidth(3);
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->GetXaxis()->SetTitle("#eta");
    hsf_func_eta[Form("dataeff_ptbin%d",i)]->GetYaxis()->SetTitle("Efficiency");

    ///mc Eff
    hsf_func_eta[Form("mceff_ptbin%d",i)] = new TH1F(Form("mceff_ptbin%d",i),Form("mceff_ptbin%d",i),netabins,etamin);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->SetLineColor(col[i]);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->SetLineStyle(2);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->SetLineWidth(3);
    hsf_func_eta[Form("mceff_ptbin%d",i)]->GetXaxis()->SetTitle("#eta");
    hsf_func_eta[Form("mceff_ptbin%d",i)]->GetYaxis()->SetTitle("Efficiency");
    
  }


  for(int i=0; i<netabins; i++){
    hsf_func_pt[Form("etabin%d",i)] = new TH1F(Form("etabin%d",i),Form("etabin%d",i),nptbins,ptmin);

    hsf_func_pt[Form("etabin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_pt[Form("etabin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_pt[Form("etabin%d",i)]->SetLineColor(col[i]);
    hsf_func_pt[Form("etabin%d",i)]->SetLineWidth(3);
    hsf_func_pt[Form("etabin%d",i)]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hsf_func_pt[Form("etabin%d",i)]->GetYaxis()->SetTitle("Scale Factor");

    ///data Eff
    hsf_func_pt[Form("dataeff_etabin%d",i)] = new TH1F(Form("dataeff_etabin%d",i),Form("dataeff_etabin%d",i),nptbins,ptmin);
    hsf_func_pt[Form("dataeff_etabin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_pt[Form("dataeff_etabin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_pt[Form("dataeff_etabin%d",i)]->SetLineColor(col[i]);
    hsf_func_pt[Form("dataeff_etabin%d",i)]->SetLineWidth(3);
    hsf_func_pt[Form("dataeff_etabin%d",i)]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hsf_func_pt[Form("dataeff_etabin%d",i)]->GetYaxis()->SetTitle("Efficiency");


    ///mc Eff
    hsf_func_pt[Form("mceff_etabin%d",i)] = new TH1F(Form("mceff_etabin%d",i),Form("mceff_etabin%d",i),nptbins,ptmin);
    hsf_func_pt[Form("mceff_etabin%d",i)]->SetMarkerStyle(sty[i]);
    hsf_func_pt[Form("mceff_etabin%d",i)]->SetMarkerColor(col[i]);
    hsf_func_pt[Form("mceff_etabin%d",i)]->SetLineColor(col[i]);
    hsf_func_pt[Form("mceff_etabin%d",i)]->SetLineStyle(2);
    hsf_func_pt[Form("mceff_etabin%d",i)]->SetLineWidth(3);
    hsf_func_pt[Form("mceff_etabin%d",i)]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hsf_func_pt[Form("mceff_etabin%d",i)]->GetYaxis()->SetTitle("Efficiency");

  }


  
  cout<<"nbins X : Y : "<<nbinsx<<" : "<<nbinsy<<endl;

  
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){

      double xlow = hsf_base->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = hsf_base->GetXaxis()->GetBinUpEdge(ix);

      double ylow = hsf_base->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = hsf_base->GetYaxis()->GetBinUpEdge(iy);

      cout<<""<<endl;
      cout<<"xlow : xhigh : ylow : yhigh : "<<xlow<<" : "<<xhigh<<" : "<<ylow<<" : "<<yhigh<<endl;
      
      double binc = hsf_base->GetBinContent(ix,iy);
      double binerr = hsf_base->GetBinError(ix,iy);
      
      double binc_bkgsys = hsf_bkgsys->GetBinContent(ix,iy); 
      double binc_sigsys = hsf_sigsys->GetBinContent(ix,iy); 
      double binc_fitsys = hsf_fitsys->GetBinContent(ix,iy); 
      double binc_ptsys = hsf_ptsys->GetBinContent(ix,iy); 
      double binc_idsys = hsf_idsys->GetBinContent(ix,iy); 

      double binc_gensys = hsf_gensys->GetBinContent(ix,iy); 
      double binc_pudsys = hsf_pudsys->GetBinContent(ix,iy); 
      double binc_puusys = hsf_puusys->GetBinContent(ix,iy); 
      

      double sys_bkg = fabs(binc - binc_bkgsys);
      double sys_sig = fabs(binc - binc_sigsys);
      double sys_fit = fabs(binc - binc_fitsys);
      double sys_pt = fabs(binc - binc_ptsys);
      double sys_id = fabs(binc - binc_idsys);

      double sys_gen = fabs(binc - binc_gensys);
      double sys_pud = fabs(binc - binc_pudsys);
      double sys_puu = fabs(binc - binc_puusys);
      double sys_pu = (sys_pud + sys_puu)/2.;
      
      cout<<"baseline_sf : bkg_sf : sig_sf : fit_sf : "<<binc <<" : "<<binc_bkgsys<< " : "<<binc_sigsys<< " : "<<binc_fitsys<<endl;
      cout<<"binc : binerr : sys_bkg : sys_sig : fit_sys : "<<binc<<" : "<<binerr<<" : "<<sys_bkg <<" : " << sys_sig<<" : "<<sys_fit<<endl;

      double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) + pow(sys_gen,2) + pow(sys_pu,2) );
      
      /*if(xlow>=1.4 && xhigh<=1.6) {
	binc = 0;
	toterr = 0;
      }
      */

      if( (xlow>=1.4442 && xhigh<=1.566) || (xlow>=-1.566 && xhigh<=-1.4442) ) {

	binc = 0;
	toterr = 0;

	
      }
  
      

      hsf_fullsys->SetBinContent(ix,iy,binc);
      hsf_fullsys->SetBinError(ix,iy,toterr);

      h_stat->SetBinContent(ix,iy,binerr);
      h_bkgsys->SetBinContent(ix,iy,sys_bkg);
      h_sigsys->SetBinContent(ix,iy,sys_sig);
      h_fitsys->SetBinContent(ix,iy,sys_fit);

      h_ptsys->SetBinContent(ix,iy,sys_pt);
      h_idsys->SetBinContent(ix,iy,sys_id);

      h_gensys->SetBinContent(ix,iy,sys_gen);

      h_pudsys->SetBinContent(ix,iy,sys_pud);
      h_puusys->SetBinContent(ix,iy,sys_puu);
      h_pusys->SetBinContent(ix,iy,sys_pu);


      //cout<<"sys_sig and bincontent from sigsys is "<<sys_sig<<":"<<h_sigsys->GetBinContent(ix,iy)<<endl;

      cout<<"New histogram : binc : binerr : "<<hsf_fullsys->GetBinContent(ix,iy)<<" : "<<hsf_fullsys->GetBinError(ix,iy)<<endl;

    }//for(int iy=1; iy<=nbinsy; iy++)
  }//for(int ix=1; ix<=nbinsx; ix++)
  
  //fout->Write();




  savePlot(hsf_fullsys, (file+"_fullsys").c_str(),true );
  savePlot(h_stat, (file+"_stat").c_str(),false );
  savePlot(h_bkgsys, (file+"_bkgsys").c_str(),false );
  savePlot(h_sigsys, (file+"_sigsys").c_str(),false );
  savePlot(h_fitsys, (file+"_fitsys").c_str(),false );
  savePlot(h_ptsys, (file+"_ptsys").c_str(),false );
  savePlot(h_idsys, (file+"_idsys").c_str(),false );
  savePlot(h_gensys, (file+"_gensys").c_str(),false );  
  savePlot(h_pusys, (file+"_pusys").c_str(),false );  


  ////now save the 1D plots in eta and pt/////
  
  ///as a function of eta
  for(int ipt=0; ipt<nptbins; ipt++){
    for(int ieta=1; ieta<=netabins; ieta++){
      
      double binc = hsf_fullsys->GetBinContent(ieta,ipt+1);
      double binerr = hsf_fullsys->GetBinError(ieta,ipt+1);

      hsf_func_eta[Form("ptbin%d",ipt)]->SetBinContent(ieta,binc);
      hsf_func_eta[Form("ptbin%d",ipt)]->SetBinError(ieta,binerr);

      //eff in data
      double deff     = hdeff_base->GetBinContent(ieta,ipt+1);
      double deff_binerr = hdeff_base->GetBinError(ieta,ipt+1);
      double binc_bkgsys = hdeff_bkgsys->GetBinContent(ieta,ipt+1); 
      double binc_sigsys = hdeff_sigsys->GetBinContent(ieta,ipt+1); 
      double binc_fitsys = hdeff_fitsys->GetBinContent(ieta,ipt+1); 

      double binc_ptsys = hdeff_ptsys->GetBinContent(ieta,ipt+1); 
      double binc_idsys = hdeff_idsys->GetBinContent(ieta,ipt+1); 

      double binc_gensys = hdeff_gensys->GetBinContent(ieta,ipt+1); 
      double binc_pudsys = hdeff_pudsys->GetBinContent(ieta,ipt+1); 
      double binc_puusys = hdeff_puusys->GetBinContent(ieta,ipt+1); 
      
      double sys_bkg = fabs(deff - binc_bkgsys);
      double sys_sig = fabs(deff - binc_sigsys);
      double sys_fit = fabs(deff - binc_fitsys);

      double sys_pt = fabs(deff - binc_ptsys);
      double sys_id = fabs(deff - binc_idsys);
      
      double sys_gen = fabs(deff - binc_gensys);
      double sys_pud = fabs(deff - binc_pudsys);
      double sys_puu = fabs(deff - binc_puusys);
      
      double sys_pu = (sys_pud + sys_puu)/2.;
      
      double toterr = sqrt( pow(deff_binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) + pow(sys_gen,2) + pow(sys_pu,2) );
      
      hsf_func_eta[Form("dataeff_ptbin%d",ipt)]->SetBinContent(ieta,deff);
      hsf_func_eta[Form("dataeff_ptbin%d",ipt)]->SetBinError(ieta,toterr);


      //eff in MC
      double mceff     = hmeff_base->GetBinContent(ieta,ipt+1);
      double mceff_binerr = hmeff_base->GetBinError(ieta,ipt+1);
      hsf_func_eta[Form("mceff_ptbin%d",ipt)]->SetBinContent(ieta,mceff);
      hsf_func_eta[Form("mceff_ptbin%d",ipt)]->SetBinError(ieta,mceff_binerr);
      
      

    }
  }

 
  ///as a function of pt
  for(int ieta=0; ieta<netabins; ieta++){
    for(int ipt=1; ipt<=nptbins; ipt++){
      
      double binc = hsf_fullsys->GetBinContent(ieta+1,ipt);
      double binerr = hsf_fullsys->GetBinError(ieta+1,ipt);

      hsf_func_pt[Form("etabin%d",ieta)]->SetBinContent(ipt,binc);
      hsf_func_pt[Form("etabin%d",ieta)]->SetBinError(ipt,binerr);
      

      //eff in data
      double deff     = hdeff_base->GetBinContent(ieta+1,ipt);
      double deff_binerr = hdeff_base->GetBinError(ieta+1,ipt);
      double binc_bkgsys = hdeff_bkgsys->GetBinContent(ieta+1,ipt); 
      double binc_sigsys = hdeff_sigsys->GetBinContent(ieta+1,ipt); 
      double binc_fitsys = hdeff_fitsys->GetBinContent(ieta+1,ipt); 

      double binc_ptsys = hdeff_ptsys->GetBinContent(ieta+1,ipt); 
      double binc_idsys = hdeff_idsys->GetBinContent(ieta+1,ipt); 
      
      double binc_gensys = hdeff_gensys->GetBinContent(ieta+1,ipt); 
      double binc_pudsys = hdeff_pudsys->GetBinContent(ieta+1,ipt); 
      double binc_puusys = hdeff_puusys->GetBinContent(ieta+1,ipt); 
            
      double sys_bkg = fabs(deff - binc_bkgsys);
      double sys_sig = fabs(deff - binc_sigsys);
      double sys_fit = fabs(deff - binc_fitsys);

      double sys_pt = fabs(deff - binc_ptsys);
      double sys_id = fabs(deff - binc_idsys);

      double sys_gen = fabs(binc - binc_gensys);
      double sys_pud = fabs(binc - binc_pudsys);
      double sys_puu = fabs(binc - binc_puusys);
      double sys_pu = (sys_pud + sys_puu)/2.;
      

      //double toterr = sqrt( pow(deff_binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) );
      double toterr = sqrt( pow(deff_binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) );
      
      hsf_func_pt[Form("dataeff_etabin%d",ieta)]->SetBinContent(ipt,deff);
      hsf_func_pt[Form("dataeff_etabin%d",ieta)]->SetBinError(ipt,toterr);


      //eff in MC
      double mceff     = hmeff_base->GetBinContent(ieta+1,ipt);
      double mceff_binerr = hmeff_base->GetBinError(ieta+1,ipt);
      hsf_func_pt[Form("mceff_etabin%d",ieta)]->SetBinContent(ipt,mceff);
      hsf_func_pt[Form("mceff_etabin%d",ieta)]->SetBinError(ipt,mceff_binerr);
      
      
    }
  }


  string name_ptbins[nptbins];
  string name_etabins[netabins];

  for(int i=0; i<nptbins; i++){
    char tmp[100];
    sprintf(tmp,"%2.0f < p_{T} <%2.0f",ptmin[i],ptmin[i+1]);
    name_ptbins[i] = tmp;
    //cout<<"pt bin name "<<name_ptbins[i]<<endl;      
  }


  for(int i=0; i<netabins; i++){
    char tmp[100];
    sprintf(tmp, "%2.1f < #eta < %2.1f",etamin[i],etamin[i+1]);
    name_etabins[i] = tmp;
  }
  
  for(int i=0; i<nptbins; i++){
    cout<<"AFTER pt bin name "<<name_ptbins[i]<<endl;      
  }

  /*
  savePlot(hsf_func_eta,nptbins,"ptbin",name_ptbins,file+"_sf_funcEta");
  savePlot(hsf_func_eta,nptbins,"dataeff_ptbin",name_ptbins,file+"_dataeff_funcEta");
  savePlot(hsf_func_eta,nptbins,"mceff_ptbin",name_ptbins,file+"_mceff_funcEta");

  ///divide this
  savePlot(hsf_func_pt,netabins,"etabin",name_etabins,file+"_sf_funcPt");
  savePlot(hsf_func_pt,netabins,"dataeff_etabin",name_etabins,file+"_dataeff_funcPt");
  savePlot(hsf_func_pt,netabins,"mceff_etabin",name_etabins,file+"_mceff_funcPt");
  */


  savePlot(hsf_func_eta,0,nptbins-1,"ptbin",name_ptbins,file+"_sf_funcEta");
  savePlot(hsf_func_pt,0,3,"etabin",name_etabins,file+"_sf_funcPt_NegativeEtaBins");
  savePlot(hsf_func_pt,4,7,"etabin",name_etabins,file+"_sf_funcPt_PositiveEtaBins");

  ///divide this
  /*savePlot(hsf_func_pt,netabins,"etabin",name_etabins,file+"_sf_funcPt");
  savePlot(hsf_func_pt,netabins,"dataeff_etabin",name_etabins,file+"_dataeff_funcPt");
  savePlot(hsf_func_pt,netabins,"mceff_etabin",name_etabins,file+"_mceff_funcPt");
  */


  /////Write all the numbers in a txt file for Dec jamboree///////////////
  std::ofstream num ( (file+"numbers.txt").c_str(), std::ofstream::out);
  
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){

      double xlow = hsf_base->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = hsf_base->GetXaxis()->GetBinUpEdge(ix);
      
      double ylow = hsf_base->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = hsf_base->GetYaxis()->GetBinUpEdge(iy);

      double binc = hsf_base->GetBinContent(ix,iy);
      double binerr = hsf_base->GetBinError(ix,iy);
      
      double sys_bkg = h_bkgsys->GetBinContent(ix,iy);
      double sys_sig = h_sigsys->GetBinContent(ix,iy);
      double sys_fit = h_fitsys->GetBinContent(ix,iy);
      
      double sys_pt = h_ptsys->GetBinContent(ix,iy);
      double sys_id = h_idsys->GetBinContent(ix,iy);

      double sys_gen = h_gensys->GetBinContent(ix,iy);
      double sys_pu = h_pusys->GetBinContent(ix,iy);
      
      double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) + pow(sys_gen,2) + pow(sys_pu,2) );

      double dataeff_base = hdeff_base->GetBinContent(ix,iy);
      double dataeff_err = hdeff_base->GetBinError(ix,iy);

      double mceff_base = hmeff_base->GetBinContent(ix,iy);
      double mceff_err = hmeff_base->GetBinError(ix,iy);

      //////////////Pure numbers
      double binc_bkgsys = hsf_bkgsys->GetBinContent(ix,iy); 
      double binc_sigsys = hsf_sigsys->GetBinContent(ix,iy); 
      double binc_fitsys = hsf_fitsys->GetBinContent(ix,iy); 
      double binc_ptsys = hsf_ptsys->GetBinContent(ix,iy); 
      double binc_idsys = hsf_idsys->GetBinContent(ix,iy); 

      double binc_gensys = hsf_gensys->GetBinContent(ix,iy); 
      double binc_pudsys = hsf_pudsys->GetBinContent(ix,iy); 
      double binc_puusys = hsf_puusys->GetBinContent(ix,iy); 
      
      double binc_pu = (binc_pudsys+binc_puusys)/2.;


      ///
      //eff in data
      double deff     = hdeff_base->GetBinContent(ix,iy);
      double deff_binerr = hdeff_base->GetBinError(ix,iy);
      binc_bkgsys = hdeff_bkgsys->GetBinContent(ix,iy); 
      binc_sigsys = hdeff_sigsys->GetBinContent(ix,iy); 
      binc_fitsys = hdeff_fitsys->GetBinContent(ix,iy); 
      
      binc_ptsys = hdeff_ptsys->GetBinContent(ix,iy); 
      binc_idsys = hdeff_idsys->GetBinContent(ix,iy); 
      
      ///in MC
      double meff        = hmeff_base->GetBinContent(ix,iy);
      double meff_binerr = hmeff_base->GetBinError(ix,iy);
      double mbinc_ptsys = hmeff_ptsys->GetBinContent(ix,iy); 
      double mbinc_idsys = hmeff_idsys->GetBinContent(ix,iy); 
      double mbinc_gensys = hmeff_gensys->GetBinContent(ix,iy);
      double mbinc_pudsys = hmeff_pudsys->GetBinContent(ix,iy);
      double mbinc_puusys = hmeff_puusys->GetBinContent(ix,iy);



    if( (xlow>=1.4442 && xhigh<=1.566) || (xlow>=-1.566 && xhigh<=-1.4442) ) {
	continue;
      }

    //minEta   maxEta   minPt   maxPt   effData    statError effMC   statError   systBkgShape    systSigShape   systFitRange   systNLOvsLO   systPU
    //num << xlow << "\t" << xhigh << "\t" << ylow << "\t" << yhigh <<"\t" << dataeff_base <<"\t" << dataeff_err << "\t" << mceff_base <<"\t" << mceff_err << "\t" << sys_bkg <<"\t" << sys_sig << "\t" << sys_fit <<"\t" << "-1" << "\t" << "-1" <<endl;

    //minEta   maxEta   minPt   maxPt   effData    statError effMC   statError   systBkgShape    systSigShape   systFitRange   systNLOvsLO   systPU  systID  systpT  
    //num << xlow << "\t" << xhigh << "\t" << ylow << "\t" << yhigh <<"\t" << dataeff_base <<"\t" << dataeff_err << "\t" << mceff_base <<"\t" << mceff_err << "\t" << sys_bkg <<"\t" << sys_sig << "\t" << sys_fit <<"\t" << sys_gen << "\t" << sys_pu << "\t" << sys_id <<"\t" << sys_pt << endl;

    //num << xlow << "\t" << xhigh << "\t" << ylow << "\t" << yhigh <<"\t" << dataeff_base <<"\t" << dataeff_err << "\t" << mceff_base <<"\t" << mceff_err << "\t" << binc_bkgsys <<"\t" << binc_sigsys << "\t" << binc_fitsys <<"\t" << binc_gensys << "\t" << binc_pu << "\t" << binc_idsys <<"\t" << binc_ptsys << endl;
   


    num << xlow << "\t" << xhigh << "\t" << ylow << "\t" << yhigh <<"\t" << deff <<"\t" << deff_binerr << "\t" << meff <<"\t" << meff_binerr << "\t" << binc_bkgsys <<"\t" << binc_sigsys << "\t" << binc_fitsys <<"\t" << binc_ptsys <<"\t" << mbinc_ptsys <<"\t"<< binc_idsys<<"\t"<< mbinc_idsys <<"\t" << mbinc_gensys << "\t" << mbinc_pudsys << "\t" << mbinc_puusys << endl;

      //cout<<"Eff_basline : Eff_sigsys: sys : "<<dataeff_base<<":"<<
          
    }
  }
 
  num.close();
  

  //////////////////////try to print in a table format///////////////
  
  //1.systematics//////////////////////////
  //cout<<"\\begin{table}[htp]"<<endl;
  cout<<"\\begin{sidewaystable}[htp]"<<endl;
  cout<<"\\begin{center}"<<endl;
  cout<<"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}"<<endl;
  cout<<"\\hline"<<endl;
  //cout<<"$\\eta$ & p$_T$ & MINOS err. & sys: BKG fit   & sys: SIG fit & sys: fit range & total error\\\\"<<endl;
  cout<<"$\\eta$ & p$_T$ & MINOS err. & sys: BKG fit   & sys: SIG fit & sys: fit range & sys:tagPt & sys:tagMedID & sys:MCGEN & sys:PU & total error\\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){

      double xlow = hsf_base->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = hsf_base->GetXaxis()->GetBinUpEdge(ix);
      
      double ylow = hsf_base->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = hsf_base->GetYaxis()->GetBinUpEdge(iy);

      double binc = hsf_base->GetBinContent(ix,iy);
      double binerr = hsf_base->GetBinError(ix,iy);
      
      double sys_bkg = h_bkgsys->GetBinContent(ix,iy);
      double sys_sig = h_sigsys->GetBinContent(ix,iy);
      double sys_fit = h_fitsys->GetBinContent(ix,iy);

      double sys_pt = h_ptsys->GetBinContent(ix,iy);
      double sys_id = h_idsys->GetBinContent(ix,iy);
      
      double sys_gen = h_gensys->GetBinContent(ix,iy);
      double sys_pu = h_pusys->GetBinContent(ix,iy);
      
      //double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) );

      //double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) );

      double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) + pow(sys_gen,2) + pow(sys_pu,2) );
      
      //double dataeff_base = hdeff_base->GetBinContent(ix,iy);

      //if(xlow>=1.4 && xhigh<=1.6) {
      if( (xlow>=1.4442 && xhigh<=1.566) || (xlow>=-1.566 && xhigh<=-1.4442) ) {
	continue;
      }
      
  
      
      cout<<xlow<<"-"<<xhigh<<" & "<<ylow<<"-"<<yhigh
	//<<" & "<<binerr<<" & "<<sys_bkg<<" & "<<sys_sig<<" & "<<sys_fit<<" & "<<toterr<<"\\\\"<<endl;
	  <<" & "<<binerr<<" & "<<sys_bkg<<" & "<<sys_sig<<" & "<<sys_fit<<" & " <<
	sys_pt <<" & " << sys_id <<" & "
	  <<sys_gen << " & " <<sys_pu << " & " <<toterr<<"\\\\"<<endl;

      cout<<"\\hline"<<endl;
      
    }
  }

  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\label{tab:EvtYields2010}"<<endl;
  cout<<"\\end{sidewaystable}"<<endl;                     
  
  cout<<""<<endl;
  cout<<""<<endl;
  cout<<""<<endl;
  //////////////////////try to print in a table format///////////////
  //2. eff and scale factors//////////////////////////
  cout<<"\\begin{table}[htp]"<<endl;
  cout<<"\\begin{center}"<<endl;
  cout<<"\\begin{tabular}{|c|c|c|c|c|c|c|}"<<endl;
  cout<<"\\hline"<<endl;
  //cout<<"$\\eta$ & $\\p_T$ & baseline SF & MINOS err. & sys: BKG fit   & sys: SIG fit & sys: fit range & total error\\\\"<<endl;
  cout<<"$\\eta$ & p$_T$ & Eff in Data & Eff in MC & SF $\\pm$ toterrr\\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){

      double xlow = hsf_base->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = hsf_base->GetXaxis()->GetBinUpEdge(ix);
      
      double ylow = hsf_base->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = hsf_base->GetYaxis()->GetBinUpEdge(iy);

      double binc = hsf_base->GetBinContent(ix,iy);
      double binerr = hsf_base->GetBinError(ix,iy);
      
      double sys_bkg = h_bkgsys->GetBinContent(ix,iy);
      double sys_sig = h_sigsys->GetBinContent(ix,iy);
      double sys_fit = h_fitsys->GetBinContent(ix,iy);

      double sys_pt = h_ptsys->GetBinContent(ix,iy);
      double sys_id = h_idsys->GetBinContent(ix,iy);

      double sys_gen = h_gensys->GetBinContent(ix,iy);
      double sys_pu = h_pusys->GetBinContent(ix,iy);

      //double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) );
      //double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) );

      double toterr = sqrt( pow(binerr,2) + pow(sys_bkg,2) + pow(sys_sig,2) + pow(sys_fit,2) + pow(sys_pt,2) + pow(sys_id,2) + pow(sys_gen,2) + pow(sys_pu,2) );

      double dataeff_base = hdeff_base->GetBinContent(ix,iy);
      double mceff_base = hmeff_base->GetBinContent(ix,iy);

      //if(xlow>=1.4 && xhigh<=1.6) {
      if( (xlow>=1.4442 && xhigh<=1.566) || (xlow>=-1.566 && xhigh<=-1.4442) ) {
	continue;
      }
      
    
      hsf_final->SetBinContent(ix,iy,binc);
      hsf_final->SetBinError(ix,iy,toterr);


      cout<<xlow<<"-"<<xhigh<<" & "<<ylow<<"-"<<yhigh
	//<<binc<<" & "<<binerr<<" & "<<sys_bkg<<" & "<<sys_sig<<" & "<<sys_fit<<" & "<<toterr<<"\\\\"<<endl;
	  <<" & "<<dataeff_base<<" & "<<mceff_base<<" & "<<binc<<" $\\pm$ "<<toterr<<"\\\\"<<endl;
      cout<<"\\hline"<<endl;
      
    }
  }

  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\label{tab:EvtYields2010}"<<endl;
  cout<<"\\end{table}"<<endl;                     


  ffinal->cd();
  hsf_final->Write();
  ffinal->Write();


  ///just chk once
  
  for(int ix=1; ix<=nbinsx; ix++){
    for(int iy=1; iy<=nbinsy; iy++){
      
      double xlow = hsf_final->GetXaxis()->GetBinLowEdge(ix);
      double xhigh = hsf_final->GetXaxis()->GetBinUpEdge(ix);
      
      double ylow = hsf_final->GetYaxis()->GetBinLowEdge(iy);
      double yhigh = hsf_final->GetYaxis()->GetBinUpEdge(iy);

      double binc = hsf_final->GetBinContent(ix,iy);
      double binerr = hsf_final->GetBinError(ix,iy);

      cout<<"xlow : xhigh : ylow : yhigh : SF : err "<<xlow<<" : "<<xhigh<<" : "<<ylow<<" : "<<yhigh<<" : "<<binc<<" : "<<binerr<<endl;

    }//for(int iy=1; iy<=nbinsy; iy++)
  }//for(int ix=1; ix<=nbinsx; ix++)
  ///

  ffinal->Close();



}

