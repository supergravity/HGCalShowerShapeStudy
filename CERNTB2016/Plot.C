void Plot(){

  gROOT->Reset();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  


  TFile* f = new TFile("output.root");
  TH1F *hLayerMC = hLayerEnergySum->Clone();

  //
  //* Layer totat energy plot:
  //

  //hLayerData->Draw("errors");
  hLayerMC->Draw("histosame");

  ostringstream buffer ;
  float binM=hLayerMC->GetBinWidth(1);
  buffer <<"Events/("<<(binM)<<" GeV)" ;
  string binM_string=buffer.str(); 
  buffer.str("");
  hLayerMC->SetXTitle("m_{#gamma#gamma} (GeV)");
  hLayerMC->SetYTitle(binM_string.c_str());


  c1->Update();
  //c1->Print("hLayerEnergySum.png");
  getchar();


}
