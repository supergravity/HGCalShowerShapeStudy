

void plot()
{

  double m[8] = {200,400,600,800,1000,1200,1500,2000};
  

  //no pixe; veto
  double phoisele_nopix[8] = {0.302,0.247,0.233,0.213,0.204,0.199,0.187,0.181};
  
  TGraph *g = new TGraph(8,m,phoisele_nopix);
  g->SetLineColor(2);
  g->Draw("APL");

  //no pixe; veto
  double phoisele_pix[8] = {0.3072,0.250,0.233,0.212,0.202,0.196,0.181,0.173};
  
  TGraph *g1 = new TGraph(8,m,phoisele_pix);
  g1->SetLineColor(2);
  g1->Draw("PL");
  
  

}
