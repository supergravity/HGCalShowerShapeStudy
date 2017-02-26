#include <map>
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TH1D.h"
#include "TFile.h"
#include <iostream>
#include "TROOT.h"
#include "math.h"
#include <fstream>
#include "TStyle.h"
#include "TPaveStats.h"
#include "TColor.h"
#include "TFractionFitter.h"
#include<sys/stat.h>
#include<sys/types.h>
#include "THStack.h"
#include "TMath.h"

//#include <sys/stat.h>

int main(int argc, char *argv[])
{
   using namespace std;
  //=========Macro generated from canvas: c1/transparent pad
  //=========  (Wed Feb  4 16:17:30 2009) by ROOT version5.18/00

   cout<<"argc:"<<argc<<endl;
   //cout<<"Check the code related to signal templates in"<<endl;

   if(argc < 2)
     {
       //cout<<"req 2 more arguments: filename followed by operation on the file (append or newfile)"<<endl;
       cout<<"req 2 more arguments: filename followed the HLT needed"<<endl;
     }
   
   //bool changepar = false;
   bool changepar = true;

   ////Photon triggered dataset
   //string hlt = "0"; //30
   //string hlt = "1"; //50
   //string hlt = "2"; //75
   //string hlt = "3"; //90
   //string hlt = "4"; //135
   string hlt = "5"; //150

   /// efficiency of signal for sieie cut of 0.011
   ///here enter the percentage of pure photon events which will fail the cut of sieie > 0.011
   //double sig_perc_gt = 1-0.99; /// make it pass thru the file if pt dependent
   double sig_perc_gt = 0; /// make it pass thru the file if pt dependent
   
   ofstream myfile;
   //myfile.open("tempfrac.txt");
   string outfile(argv[1]);
   
   ///in which mode to open the file
   //string op(argv[2]);
   
   ////which hlt
   hlt=argv[2];
   cout<<"hlt:"<<hlt<<endl;

   //      Check if a file exists
   //@param[in] filename - the name of the file to check
   //@return    true if the file exists, else false

   struct stat buf;
   if (stat( (outfile+hlt+".txt").c_str(), &buf) != -1)
     {
       //cout<<"G8!!!file exists"<<endl;
       cout<<"WARNING!!!!FILE EXISTS - WILL APPEN IN THE FILE"<<endl;
       myfile.open( (outfile+hlt+".txt").c_str(), ios::app );
     }

   else
     {
       cout<<"FILE DOESNT EXIST - WILL OPEN A NEW FILE"<<endl;
       myfile.open( (outfile+hlt+".txt").c_str() );
     }
   
   /*cout<<"op:"<<op<<endl;
   
   if( op == "append" )
     {
       myfile.open( (outfile+".txt").c_str(), ios::app );
     }

   if( op == "newfile" )
     {
       myfile.open( (outfile+".txt").c_str() );
     }
   */
   
   
   //create a directory here
   if(mkdir(("plotshlt"+hlt).c_str(),0777)==-1)//creating a directory
     {
       //cerr<<"Error :  "<<strerror(errno)<<endl;
       cerr<<"Error in creating directory" <<endl;
       //exit(1);
     }

   gStyle->SetErrorX(1);
   
   TFile *f1;
   TFile *f2; 
   TH1D  *hdata;
   TH1D  *hqcd;
   TH1D  *hphojet;
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   TFractionFitter* fit;
   
   std::map<string,TCanvas*> c1;
   std::map<string,TLegend*> leg;
   
   
   TString filename1;
   TString filename2;
   //filename1 = "file_dR0.300000.root";
   //filename1 = "file_dR0.300000_nosafeveto.root";
   filename1 = "file_dR0.300000.root";
   //filename1 = "file_dR0.300000_chargediso50.root";
   filename2 = "nopixelphojet.root";  //phojet file name

   string shistname;
   TString histname;
   TString datahistname;
   TString qcdhistname; //trkflip - > histos will go here
   TString fullhistname; //phojet
   TString sighistname;
  

  char cxvariable[200];
  ifstream infile_xvar;
  infile_xvar.open("xvariables.list", ifstream::in );
  while(!infile_xvar.eof()){
    infile_xvar >>cxvariable;
    string xvariable(cxvariable);
    if(strncmp(cxvariable,"#",1)==0)
      {
	continue;
      }
    histname  = xvariable;
    shistname = histname; 
    datahistname = "datasigietahlt"+hlt+"_ptbin"+shistname; //data histos without trk flip
    fullhistname = "sigmaieta"+shistname; //phojet histos without trk flip
    qcdhistname  = "bkgsigietahlt"+hlt+"_ptbin"+shistname; //data histos with trkflip

    
    sighistname = "sigietaEBptbin"+shistname;

    /*if( shistname=="300_400" || shistname=="400_1000" || shistname=="1000_2000" )
      sighistname = "sigietaEBptbin200_300";
    */

    
    cout<<"datahistname = "<<datahistname<<endl;
    cout<<"fullhistname = "<<fullhistname<<endl;
    cout<<"sigistname = "<<sighistname<<endl; 

    c1[xvariable] = new TCanvas(("c"+xvariable).c_str(), "transparent pad",179,30,698,498);
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(3);
    gStyle->SetPadGridY(3);
    gStyle->SetGridStyle(3);
    c1[xvariable]->Range(0,0,1,1);
    c1[xvariable]->SetFillColor(0);
    c1[xvariable]->SetLogy();
    c1[xvariable]->SetBorderMode(0);
    c1[xvariable]->SetBorderSize(2);
    c1[xvariable]->SetFrameFillStyle(0);
    c1[xvariable]->SetFrameBorderMode(0);
    c1[xvariable]->Divide(2);
    
    //legend
    //leg[xvariable] = new TLegend(0.2,0.6,0.4,0.8);
    //leg[xvariable] = new TLegend(0.4755043,0.5212766,0.8227666,0.8574468);
    leg[xvariable] = new TLegend(0.5360231,0.7723404,0.9971182,0.9957447,NULL,"brNDC");
    leg[xvariable]->SetFillColor(0);
    //gPad->SetLogy(1);
    f1 = TFile::Open(filename1) ;
    f1->cd();

    hdata = (TH1D*)f1->Get(datahistname)->Clone();
    hqcd  = (TH1D*)f1->Get(qcdhistname)->Clone();
    //hqcd  = (TH1D*)f1->Get("tsigmaietaEE_eb100to120")->Clone();
    
    double integ_data = hdata->Integral();
    double integ_qcd = hqcd->Integral();
    cout<<"integ_data = "<<integ_data<<endl;
    cout<<"integ_qcd = "<<integ_qcd<<endl;

    //hqcd->Scale(integ_data/integ_qcd); 

    f2 = TFile::Open(filename2) ;
    f2->cd();
  
    /*TH1D *htemp = (TH1D*)f2->Get("sigmaietaEB_eb30to40")->Clone();
    int mynbins = htemp->GetNbinsX();
    double xlow = htemp->GetXaxis()->GetXmin();
    double xhigh = htemp->GetXaxis()->GetXmax();
    */
      
    
    TH1D *htemp = (TH1D*)f2->Get("sigietaEBptbin30_40")->Clone();
    int mynbins = htemp->GetNbinsX();
    double xlow = htemp->GetXaxis()->GetXmin();
    double xhigh = htemp->GetXaxis()->GetXmax();

    TH1D *hphojet = new TH1D("hphojet","signal",mynbins,xlow,xhigh);

    if( shistname=="30_50" )
      {
	TH1D *h1 = (TH1D*)f2->Get("sigietaEBptbin30_40")->Clone();
	TH1D *h2 = (TH1D*)f2->Get("sigietaEBptbin40_50")->Clone();
	hphojet->Add(h1);
	hphojet->Add(h2);
      }


    if( shistname=="50_75" )
      {
	TH1D *h1 = (TH1D*)f2->Get("sigietaEBptbin50_60")->Clone();
	TH1D *h2 = (TH1D*)f2->Get("sigietaEBptbin60_75")->Clone();
	hphojet->Add(h1);
	hphojet->Add(h2);
      }

    if( shistname=="75_90" )
      {
	hphojet = (TH1D*)f2->Get("sigietaEBptbin75_90")->Clone();
      }

    if( shistname=="90_135" )
      {
	TH1D *h1 = (TH1D*)f2->Get("sigietaEBptbin90_110")->Clone();
	TH1D *h2 = (TH1D*)f2->Get("sigietaEBptbin110_135")->Clone();
	hphojet->Add(h1);
	hphojet->Add(h2);
      }

    if( shistname=="135_150" )
      {
	TH1D *hphojet = (TH1D*)f2->Get("sigietaEBptbin135_150")->Clone();
      }

    if( shistname=="150_400" )
      {
	TH1D *h1 = (TH1D*)f2->Get("sigietaEBptbin150_200")->Clone();
	TH1D *h2 = (TH1D*)f2->Get("sigietaEBptbin200_300")->Clone();
	TH1D *h3 = (TH1D*)f2->Get("sigietaEBptbin300_400")->Clone();
	hphojet->Add(h1);
	hphojet->Add(h2);
	hphojet->Add(h3);
      }

    if( shistname=="400_2000" )
      {
	TH1D *h1 = (TH1D*)f2->Get("sigietaEBptbin400_1000")->Clone();
	TH1D *h2 = (TH1D*)f2->Get("sigietaEBptbin1000_2000")->Clone();
	hphojet->Add(h1);
	hphojet->Add(h2);

      }


  //hphojet  = (TH1D*)f2->Get(fullhistname)->Clone();
    
    double integ_phojet = hphojet->Integral();
    cout<<"integ_phojet = "<<integ_phojet<<endl;

    //hphojet->Scale(integ_data/integ_phojet);
    int nbins = hphojet->GetNbinsX();
    cout<<"nbins = "<<nbins<<endl;
    
    ////try smoothening here
    //hphojet->Smooth(10);
    //hqcd->Smooth(10);
    //hphojet->Smooth(10);
    //hqcd->Smooth(9);
    //////////////////////

    //till this point I have defined all the three histos

    /*hdata->Rebin(4);
    hqcd->Rebin(4);
    hphojet->Rebin(4);
    */

    int binsTorebin = 1;
    //int binsTorebin = 4;

    TH1D *hndata = (TH1D*)hdata->Rebin(binsTorebin,"hndata");
    TH1D *hnqcd = (TH1D*)hqcd->Rebin(binsTorebin,"hnqcd");
    TH1D *hnphojet = (TH1D*)hphojet->Rebin(binsTorebin,"hnphojet");

    /////for the TFractionFitter MC prediction
    //TH1D *hfqcd = (TH1D*)hqcd->Rebin(binsTorebin,"hfqcd");
    //TH1D *hfphojet = (TH1D*)hphojet->Rebin(binsTorebin,"hfphojet");
    ///////////////////////////////////////////////////SAMS way/////////////////////////////////////////////
    cout<<"bin corresponding ot 0.011 : "<<hphojet->FindBin(0.011)-1<<endl;
    int myreqbin = hphojet->FindBin(0.011)-1;  ////should be 22 for 0.011; 26 if 0.013

    cout<<"nbins:hdata:hndata:"<<hdata->GetNbinsX()<<":"<<hndata->GetNbinsX()<<endl;
    nbins = hndata->GetNbinsX();
    
    
    //Now get the fraction of photons in 0.011 range

    double data_lt     = hdata->Integral(0,myreqbin);
    double data_gt     = hdata->Integral(myreqbin+1, nbins);
    double data_tot = hdata->Integral(1, nbins);
    
    double uqcd_gt     = hqcd->Integral(myreqbin+1, nbins);
    double uqcd_lt     = hqcd->Integral(1,myreqbin);
    double uqcd_tot = hqcd->Integral(1,nbins);
    
   double corr_scfac = ( data_gt - sig_perc_gt*data_tot )/( uqcd_gt - sig_perc_gt*uqcd_tot ); 
   
   ///scaled QCD
   double sqcd_lt = corr_scfac*uqcd_lt;
   
   cout<<"====scale factor: "<<corr_scfac<<endl;
    
    ///this is true if sig contribution for > region = 0
    //double sqcd_lt = (data_gt/uqcd_gt)*uqcd_lt;
    //double sigperc_lt = (data_lt-sqcd_lt)/data_lt;
    //double qcdperc_lt = (sqcd_lt)/data_lt;
    
    ///for non-contribution fro signal in hte region of > 
    ///(chk the derivation in the notes and code in hte laptop here:
    ////home/shilpi/excitedLepton/myAna2011/phoFR/study_2012/afterremovingbeamhalo_badchannelsMC_5oct/templates_samsway
    

    double sigperc_lt = (data_lt-sqcd_lt)/data_lt;
    double qcdperc_lt = (sqcd_lt)/data_lt;
    
    cout<<"======SAMSWAY======"<<endl;
    cout<<"sigperc_lt:"<<sigperc_lt
	<<" :qcdperc_lt:"<<qcdperc_lt<<endl;

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    mc->Add(hnqcd);
    mc->Add(hnphojet);
    fit = new TFractionFitter(hndata, mc);    // initialise
    
    ///try to change the fit parameters
    if(changepar){
      TVirtualFitter* vfit = fit->GetFitter();
      vfit->SetParameter(0, "qcd", 0.5, 0.0001, 0., 1.);
      vfit->SetParameter(1, "realpho", 0.5, 0.0001, 0., 1.);
    }
    
    //fit = new TFractionFitter(hndata, mc, "V");    // initialise
    fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
    fit->SetRangeX(1,nbins);
    
    //hnqcd->Smooth(1);

    ////so that the fit happens
    /*if(integ_qcd>=10 && integ_data>=10){
      Int_t status = fit->Fit();               // perform the fit
      cout << "fit status: " << status << endl;
      if (status == 0) {                       // check on fit status
	TH1F* result = (TH1F*) fit->GetPlot();
	
	  } //Fit Status
    }//if(integ_qcd>30 && integdata>30)
    
    TH1F* result;
    
	cout<<"RESULTS : "<<endl;
	double value1 = 0, error1 = 0;
	double value2 = 0, error2 = 0;
	
	fit-> GetResult(0,value1,error1);
	cout<<"qcd       = "<<value1<<endl;
	cout<<"qcd error = "<<error1<<endl;
	
	fit-> GetResult(1,value2,error2);
	cout<<"phojet       = "<<value2<<endl;
	cout<<"phojet error = "<<error2<<endl;
	
	double frac_qcd     = (value1*hndata->Integral());
	double frac_realpho = (value2*hndata->Integral());
	
	//scale MC template histos ot the frac.
	//qcd
	double ninteg_qcd = hqcd->Integral();
	hqcd->Scale(frac_qcd/ninteg_qcd);
	hnqcd->Scale(frac_qcd/ninteg_qcd);
	
	//phojet
	double ninteg_es = hphojet->Integral();
	hphojet->Scale(frac_realpho/ninteg_es);
	hnphojet->Scale(frac_realpho/ninteg_es);

	cout<<"bin corresponding ot 0.011 : "<<hphojet->FindBin(0.011)-1<<endl;
	int reqbin = hphojet->FindBin(0.011)-1;  ////should be 22 for 0.011; 26 if 0.013
	
	//Now get the fraction of photons in 0.011 range
	double frac_realpho_range011 = hphojet->Integral(0,reqbin);
	double frac_qcd_range011     = hqcd->Integral(0,reqbin);
	double data_range011         = hdata->Integral(0,reqbin);
    
	
	cout<<"========== USE THIS================="<<endl;
	cout<<"frac of real pho = "<<frac_realpho_range011/(frac_realpho_range011+frac_qcd_range011)<<endl;
        cout<<"frac of qcd      = "<<frac_qcd_range011/(frac_realpho_range011+frac_qcd_range011)<<endl;

	cout<<"fraction sum in whole range : "<<value1+value2<<endl;
        cout<<"total = frac_qcd_range011 + frac_realpho_range011 = "<<(frac_qcd_range011+frac_realpho_range011)/(frac_qcd_range011+frac_realpho_range011)<<endl;
        cout<<"========================================="<<endl;
	

	
	hndata->SetMarkerStyle(20);
	
	Int_t ci;   // for color index setting
	ci = 4;
	
	hnqcd->SetLineColor(ci);
	hnqcd->SetFillStyle(3005);
	hnqcd->SetFillColor(ci);
	
	TH1F *hfqcd = (TH1F*)fit->GetMCPrediction(0);
        hfqcd->SetLineColor(ci);
        hfqcd->SetFillStyle(3005);
        hfqcd->SetFillColor(ci);
	double finteg_qcd = hfqcd->Integral();
        hfqcd->Scale(frac_qcd/finteg_qcd);

	
	
	ci = 1;
	hnphojet->SetLineColor(ci);
	hnphojet->SetFillStyle(6);
	hnphojet->SetFillColor(3);
	TH1F *hfphojet = (TH1F*)fit->GetMCPrediction(1);
        hfphojet->SetLineColor(ci);
        hfphojet->SetFillStyle(6);
        hfphojet->SetFillColor(3);
	double finteg_phojet = hfphojet->Integral();
        hfphojet->Scale(frac_realpho/finteg_phojet);



	result->SetLineColor(2);
	result->SetLineWidth(2);
	c1[xvariable]->cd();
	
	THStack *myStack;
	myStack = new THStack("myStack","myStack");
	//myStack->Add(hnphojet);
	//myStack->Add(hnqcd);
	
	./fn	
	myStack->Add(hfphojet);
	myStack->Add(hfqcd);

	
	float max = TMath::Max(myStack->GetMaximum(),hndata->GetMaximum());
	myStack->SetMaximum(max*1.2);
	//myStack->SetMinimum(0.01);
	myStack->SetMinimum(0.5);
	hndata->SetMaximum(max*1.2);
	//hndata->SetMinimum(0.01);
	hndata->SetMinimum(0.5);
	
	///to suppress the errors along X direction
	gStyle->SetErrorX(0.);
	
	myStack->Draw();
	
	//hnphojet->Draw();
	//hnqcd->Draw("same");
	//hndata->Draw("PEsame");
	hndata->Draw("samePE1");
	result->Draw("same");
	
	leg[xvariable]->AddEntry(hnphojet,"fit to signal photons : #gamma+jet MC","f");
	leg[xvariable]->AddEntry(hnqcd,"fit to fake photons : Data","f");
	leg[xvariable]->AddEntry(hndata,"Data","PE");
	leg[xvariable]->AddEntry(result,"fit to Data","l");
	
	leg[xvariable]->Draw();
	c1[xvariable]->Modified();
	c1[xvariable]->cd();
	c1[xvariable]->SetSelected(c1[xvariable]);
	c1[xvariable]->SaveAs(("plotshlt"+hlt+"/"+shistname+".gif").c_str());
	c1[xvariable]->SaveAs(("plotshlt"+hlt+"/"+shistname+".C").c_str());
	
	/////from hte getMCpredicted histos  
	cout<<"bin corresponding ot 0.011 : "<<hphojet->FindBin(0.011)-1<<endl;
	//Now get the fraction of photons in 0.011 range
	double ffrac_realpho_range011 = hfphojet->Integral(0,reqbin);
	double ffrac_qcd_range011     = hfqcd->Integral(0,reqbin);
	double fdata_range011         = hdata->Integral(0,reqbin);

        cout<<"=======================FROM MC PREDICTED HISTSO======================"<<endl;
        //cout<<"frac of real pho = "<<ffrac_realpho_range011/fdata_range011<<endl;
	//cout<<"frac of qcd      = "<<ffrac_qcd_range011/fdata_range011<<endl;

	cout<<"frac of real pho = "<<ffrac_realpho_range011/(ffrac_realpho_range011+ffrac_qcd_range011)<<endl;
	cout<<"frac of qcd      = "<<ffrac_qcd_range011/(ffrac_realpho_range011+ffrac_qcd_range011)<<endl;

	double netdiffQCD = ffrac_qcd_range011/(ffrac_realpho_range011+ffrac_qcd_range011) - frac_qcd_range011/(frac_realpho_range011+frac_qcd_range011);
	
	cout<<"fraction sum in whole range : "<<value1+value2<<endl;
	cout<<"total = frac_qcd_range011 + frac_realpho_range011 = "<<(ffrac_qcd_range011+ffrac_realpho_range011)/fdata_range011<<endl;
	
	cout<<" netdiffQCD : "<< netdiffQCD <<endl;
	
	*/
	
	double samserr = 0.01; //asssuming 0.02% of the signal fails sieie > 0.01
	
	///draw the templates here
	hnqcd->SetMinimum(0.1);
	hnqcd->SetLineColor(40);
	hnqcd->SetFillColor(40);
	hnqcd->SetFillStyle(3001);
	hnqcd->SetLineWidth(3);
	
	hdata->SetMinimum(0.1);
	hdata->SetLineColor(1);
	hdata->SetFillColor(1);
	hdata->SetMarkerColor(1);
	hdata->SetMarkerStyle(20);
	
	double max[3] = {0.};
	max[0] = hdata->GetMaximum();
	max[1] = hnqcd->GetMaximum();



	
	c1[xvariable]->cd(1);
	gPad->SetLogy();
	double maxf = TMath::MaxElement(2,max);
	
	

	//hnqcd->SetMaximum(maxf*1.2);
	//hdata->SetMaximum(maxf*1.2);
	
	hnqcd->SetMinimum(0.1);
	hnqcd->SetTitle("Original template");
	//hdata->Draw("P");
	hnqcd->Draw("");
	
	leg[xvariable] = new TLegend(0.5255163,0.6721398,0.9878122,0.8972458,NULL,"brNDC");
	leg[xvariable]->SetFillColor(0);
	//leg[xvariable]->AddEntry(hdata,"Data","lp");
	leg[xvariable]->AddEntry(hnqcd,"Unscaled background","f");
	leg[xvariable]->Draw();
  

	c1[xvariable]->cd(2);
	gPad->SetLogy();
	hnqcd->Scale(corr_scfac);
	hnphojet->Add(hndata,hnqcd,1,-1);
	
	////no signal beyong 0.011
	for( int ibin=myreqbin+1; ibin<nbins; ibin++){
	  
	  hnphojet->SetBinContent(ibin,0);
	}
	
	hnphojet->SetMinimum(0.1);
	hnphojet->SetLineColor(42);
	hnphojet->SetFillColor(42);
	hnphojet->SetFillStyle(3305);
	hnphojet->SetLineWidth(3);

	hnqcd->SetMinimum(0.1);
	hdata->SetMinimum(0.1);

	THStack *myStack;
	myStack = new THStack("myStack","Template Fits");
	myStack->Add(hnphojet);
	myStack->Add(hnqcd);

	maxf = TMath::Max(myStack->GetMaximum(),hdata->GetMaximum());
	myStack->SetMaximum(maxf*1.2);
	hdata->SetMaximum(maxf*1.2);

	myStack->Draw();
	//hdata->Draw("PEsame");
	hdata->Draw("Psame");
	
	//hnqcd->Draw("same");
	//hnphojet->Draw("same");
	

	leg[(xvariable+"1")] = new TLegend(0.5255163,0.6721398,0.9878122,0.8972458,NULL,"brNDC");
	leg[(xvariable+"1")]->SetFillColor(0);
	leg[(xvariable+"1")]->AddEntry(hdata,"Data","lp");
	leg[(xvariable+"1")]->AddEntry(hnqcd,"Scaled background","f");
	leg[(xvariable+"1")]->AddEntry(hnphojet,"Estimated Signal","f");
	leg[(xvariable+"1")]->Draw();

	
	myfile<< qcdperc_lt <<"\n";
	myfile<< samserr <<"\n";
	
	c1[xvariable]->Modified();
	c1[xvariable]->cd();
	c1[xvariable]->SetSelected(c1[xvariable]);
	c1[xvariable]->SaveAs(("plotshlt"+hlt+"/"+xvariable+"_all4.gif").c_str());
	c1[xvariable]->SaveAs(("plotshlt"+hlt+"/"+xvariable+"_all4.C").c_str());
	
	//  } //Fit Status
	//}//if(integ_qcd>30 && integdata>30)
    
	/*  else 
      {
	cout<<"!!!!!!!CANNOT FIT!!! integral < 10"<<endl;
	cout<<"take frac of real pho = 1 and QCD = 0"<<endl;
      }
	*/

	//mc->Clear();
  
  } //Loop over xvariable.list lines

  //cout<<"WARNING!!!!fix sighistname sigietaEBptbin200_300 for pt > 300"<<endl;
  
}//main() ends here


#ifdef __MAKECINT__                           
#pragma link C++ class map<string,TCanvas*>;
#pragma link C++ class map<string,TPad*>;
#pragma link C++ class map<string,TLegend*>;
#pragma link C++ class map<string,TGraphErrors*>;
#pragma link C++ class map<string,TLine*>;
#endif
