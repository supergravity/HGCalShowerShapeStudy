#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
//#include "utils.h"



//======================================================================//
// Global Functions:
// SigmaE/E for calorimeter with stoch and const term only
double resolutionf1(double *x,double *par) {
  // x[0]   - Beam Energy (just 6 points for now) 
  // par[0] - Stochastic Term
  // par[1] - Constant Term


  //cout << x[0] << "," << par[0] << ", " << par[1] << endl;
  double func_value=sqrt(par[0]*par[0]/x[0]+par[1]*par[1]);
  return func_value;
}
//======================================================================//


void profilePlots(){

  gROOT->Reset();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();



  //TFile* fresol = new TFile("outputMC100GeV.root");
  TFile* fresol = new TFile("output.root");
  //TCanvas *c = new TCanvas("c","c",600,400);
  //*********************************************************
  //* 200 GeV
  //*********************************************************

  char histoDraw[100];
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.95);
  gStyle->SetPalette(1);

  //atlasStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(6.0);
  //atlasStyle->SetPadBottomMargin(0.16);
  //atlasStyle->SetPadLeftMargin(0.);

  bool doScale = true;


  double detaCut[31] = {
    0.025 ,0.030, 0.030,
    0.050 ,0.050, 0.050,,0.050, 0.070, 0.070, 0.070, 0.070, 
    0.070, 0.070, 0.070, 0.070, 0.070, 0.070, 0.090, 0.090, 0.090, 0.080,
    0.080, 0.080, 0.070, 0.070, 0.070, 0.050, 0.050, 0.025, 0.020, 0.020
  };
  double dphiCut[31] = {
    0.025 ,0.030, 0.030,
    0.050 ,0.050, 0.050,,0.050, 0.070, 0.070, 0.070, 0.070, 
    0.070, 0.070, 0.070, 0.070, 0.070, 0.070, 0.090, 0.090, 0.090, 0.080,
    0.080, 0.080, 0.070, 0.070, 0.070, 0.050, 0.050, 0.025, 0.020, 0.020
  };


  //
  //* Draw eta projections in one graph
  //
  TCanvas p;
  p.SetLogz();
  gPad->SetMargin(0.1,0.18,0.15,0.1);

  //TH1D* detaVlay = detaVlayer->Clone();
  //detaVlay->Draw("colz");
  //detaVlay->SetXTitle("#eta_{hit}-#eta_{hitmax}");
  //detaVlay->SetYTitle("Layer");
  //detaVlay->SetZTitle("E(MIPs)/event");


  //p.Update();
  //getchar();

  TH1D* dRVlay = dRVlayer->Clone();

  dRVlay->Draw("colz");
  dRVlay->SetXTitle("dR(cm)");
  dRVlay->SetYTitle("Layer");
  dRVlay->SetZTitle("E(MIPs)/pad/event");
  p.Print("dRxyVlayer.png");
  p.Update();
  getchar();



  p.SetLogy();


  //First fit the deta shape:
  gPad->SetMargin(0.15,0.1,0.15,0.1);
  TH1D* detaShower = dRprofAll->Clone();
  detaShower->Sumw2();
  detaShower->Scale(1./detaShower->GetMaximum());
  detaShower->Draw();
  //TF1 *func = new TF1("func","[0]*exp([1]*x*x)/(1.+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x)",0.0,0.15);
  TF1 *func = new TF1("func","[0]*exp([1]*x)*(1.+[2]*x+[3]*x*x)+([4]*x+[5]*x*x+[6]*x*x*x+[7]*x*x*x*x)",-0.5,14.0);


  //func->FixParameter(0,1);
  //func->FixParameter(1,-1.20682);
  //func->FixParameter(2,-0.609984);
  //func->FixParameter(3,0.151971);
  //func->FixParameter(4,5.65297e-05);
  //func->FixParameter(5,-1.4194e-05);
  //func->FixParameter(6,1.16217e-06);
  //func->FixParameter(7,-3.12166e-08);
  
  //100GeV MC fit:
  func->FixParameter(0,1);
  func->FixParameter(1,-1.20682);
  func->FixParameter(2,-0.609984);
  func->FixParameter(3,0.151971);
  func->FixParameter(4,5.65297e-05);
  func->FixParameter(5,-1.4194e-05);
  func->FixParameter(6,1.16217e-06);
  func->FixParameter(7,-3.12166e-08);


  detaShower->Fit("func","","",1.5,14.);
  //detaShower->Fit("func","","",1.5,24);//for 3rd compartment
  detaShower->SetXTitle("dR (cm)");
  detaShower->SetYTitle("#frac{1}{E} #frac{dE}{dR}");

  //func->Draw();
  for(int i=0; i<8; i++)
    cout <<"func->SetParameter("<<i<<","<<func->GetParameter(i)<< ");"<< endl;

  func->Draw("same");
  TLegend leg(0.15,0.65,0.65,0.94);
  leg.SetFillStyle(0);  
  leg.SetBorderSize(0);
  p.Print("FitCompt1.png");
  p.Update();
  getchar();
  

  //First Compartment
  gStyle->SetOptFit(0);
  TLegend legen(0.60,0.64,0.94,0.89);
  legen.SetFillStyle(0);  
  legen.SetBorderSize(0);
  int firstLayer=0;//meaning firstlayer is first (start from 0)
  int lastLayer=8;//meaning lastlayer-1 is last
  for (int i=firstLayer; i<lastLayer; i++) {

    int icolor=i-firstLayer+1;
    if(i>=9) icolor+=1;
    if(doScale) {
      sprintf(histoDraw,"dRprof_%d->Scale(1./dRprof_%d->GetMaximum())",i,i,i);
      //sprintf(histoDraw,"dRprof_%d->Scale(1./dRprof_%d->GetEntries())",i,i,i);
      gROOT->ProcessLine(histoDraw);
    }
    if(i==firstLayer) {
      sprintf(histoDraw,"dRprof_%d->Draw(\"\")",i,i);
      gROOT->ProcessLine(histoDraw);
      sprintf(histoDraw,"dRprof_%d->SetXTitle(\"dR (cm)\")",i,i);
      gROOT->ProcessLine(histoDraw);
      sprintf(histoDraw,"dRprof_%d->SetYTitle(\"#frac{1}{E} #frac{dE}{dR}\")",i,i);
      gROOT->ProcessLine(histoDraw);
      //legen.SetHeader("something");
    }
    else {
      sprintf(histoDraw,"dRprof_%d->Draw(\"same\")",i,i);
      gROOT->ProcessLine(histoDraw);
    }
    atlasStyle->SetOptStat(0);
    sprintf(histoDraw,"dRprof_%d->SetLineColor(%d)",i,icolor);
    gROOT->ProcessLine(histoDraw);
    sprintf(histoDraw,"dRprof_%d->SetLineWidth(2.)",i);
    gROOT->ProcessLine(histoDraw);
    //sprintf(histoDraw,"legen.AddEntry(dRprof_%d,\"layer %d\",\"L\")",i,i);
    //gROOT->ProcessLine(histoDraw);

    if(i==lastLayer-1) func->Draw("same");
    cout << "Layer: " << i << endl;
    p.Update();
    getchar();
    
  }
  legen.AddEntry(dRprof_0,"layer 1","L");
  legen.AddEntry(dRprof_1,"layer 2","L");
  legen.AddEntry(dRprof_2,"layer 3","L");
  legen.AddEntry(dRprof_3,"layer 4","L");
  legen.AddEntry(dRprof_4,"layer 5","L");
  legen.AddEntry(dRprof_5,"layer 6","L");
  legen.AddEntry(dRprof_6,"layer 7","L");
  legen.AddEntry(dRprof_7,"layer 8","L");
  legen.Draw();
  p.Print("Compt1Layers.png");
  getchar();



#ifdef doDETA

  //First Compartment
  int firstLayer=3;
  //int lastLayer=11;
  int lastLayer=26;
  for (int i=firstLayer; i<lastLayer; i++) {
    sprintf(histoDraw,"proj%d = hetaphiprof_%d->ProjectionX()",i,i);
    gROOT->ProcessLine(histoDraw);

    if(doScale) {
      sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetMaximum())",i,i,i);
      //sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetEntries())",i,i,i);
      gROOT->ProcessLine(histoDraw);
    }
    if(i==firstLayer) {
      sprintf(histoDraw,"proj%d->Draw(\"\")",i,i);
      gROOT->ProcessLine(histoDraw);
    }
    else {
      sprintf(histoDraw,"proj%d->Draw(\"same\")",i,i);
      gROOT->ProcessLine(histoDraw);
    }

    if(i==lastLayer-1) func->Draw("same");
    cout << "Layer: " << i << endl;
    p.Update();
    getchar();
    
  }

  //Second Compartment
  firstLayer=11;
  lastLayer=20;
  for (int i=firstLayer; i<lastLayer; i++) {
    sprintf(histoDraw,"proj%d = hetaphiprof_%d->ProjectionX()",i,i);
    gROOT->ProcessLine(histoDraw);

    if(doScale) {
      sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetMaximum())",i,i,i);
      //sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetEntries())",i,i,i);
      gROOT->ProcessLine(histoDraw);
    }
    if(i==firstLayer) {
      sprintf(histoDraw,"proj%d->Draw(\"\")",i,i);
      gROOT->ProcessLine(histoDraw);
    }
    sprintf(histoDraw,"proj%d->Draw(\"same\")",i,i);
    gROOT->ProcessLine(histoDraw);


    cout << "Layer: " << i << endl;
    if(i==lastLayer-1) func->Draw("same");
    p.Update();
    getchar();
    
  }

  //Third Compartment
  firstLayer=20;
  lastLayer=27; //ignore last layers
  for (int i=firstLayer; i<lastLayer; i++) {
    sprintf(histoDraw,"proj%d = hetaphiprof_%d->ProjectionX()",i,i);
    gROOT->ProcessLine(histoDraw);

    if(doScale) {
      sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetMaximum())",i,i,i);
      //sprintf(histoDraw,"proj%d->Scale(1./proj%d->GetEntries())",i,i,i);
      gROOT->ProcessLine(histoDraw);
    }
    if(i==firstLayer) {
      sprintf(histoDraw,"proj%d->Draw(\"\")",i,i);
      gROOT->ProcessLine(histoDraw);
    }
    sprintf(histoDraw,"proj%d->Draw(\"same\")",i,i);
    gROOT->ProcessLine(histoDraw);

    cout << "Layer: " << i << endl;
    if(i==lastLayer-1) func->Draw("same");
    p.Update();
    getchar();
    
  }

#endif

  cout << "Finished with eta projections... press return to continue." << endl;
  getchar();

  //
  //* Draw the 2D projection of the shower
  //
  TCanvas c;
  c.SetLogz();

  for (int i=0; i<8; i++) {

    TLegend legS(0.15,0.82,0.65,0.94);
    legS.SetFillStyle(0);  
    legS.SetBorderSize(0);

    sprintf(histoDraw,"hetaphiprof_%d->Draw(\"colz\")",i);
    gROOT->ProcessLine(histoDraw);
    sprintf(histoDraw,"hetaphiprof_%d->SetMinimum(0.1)",i);
    gROOT->ProcessLine(histoDraw);

    sprintf(histoDraw,"hetaphiprof_%d->SetXTitle(\"x_{hit}-x_{hitmax}(cm)\")",i);
    gROOT->ProcessLine(histoDraw);
    sprintf(histoDraw,"hetaphiprof_%d->SetYTitle(\"y_{hit}-y_{hitmax}(cm)\")",i);
    gROOT->ProcessLine(histoDraw);

    sprintf(histoDraw,"E=100GeV, MIPs/pad/event, Layer:%d",i+1);
    legS->SetHeader(histoDraw);
    legS->Draw();

    ostringstream buffer ;
    buffer <<"XYEProfileE100.Layer."<<i+1<<".png" ;
    string binM_string=buffer.str(); 
    buffer.str("");
    c.Update();
    c.Print(binM_string.c_str());

    getchar();
  }



#endif

}
