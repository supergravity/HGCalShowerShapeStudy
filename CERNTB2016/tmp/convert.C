void convert(Int_t json=0, Int_t xs=69) {

  // json == 0 (golden)
  // json == 1 (silver)

  Char_t fname[200];
  if (json == 0) sprintf(fname, "PU_13TeV_GoldenJSON_%dmb.root", xs);
  if (json == 1) sprintf(fname, "PU_13TeV_SilverJSON_%dmb.root", xs);
  TFile *fdata = new TFile(fname);
  //TFile *fmc   = new TFile("/data6/ggNtuples/V07_04_14_00/job_spring15_WJetsToLNu_aMCatNLO_miniAOD.root");
  //TFile *fmc   = new TFile("/data6/ggNtuples/V07_04_14_01/job_spring15_DYJetsToLL_m50_miniAOD.root");
  TFile *fmc   = new TFile("/data6/ggNtuples/V07_04_14_01/job_spring15_Zg_aMCatNLO_miniAOD.root");

  //TH1D *hnew = new TH1D("mcwei_run000001", "pileup", 1000, 0, 200);
  TH1D *hnew = new TH1D("mcwei_run000001", "pileup", 200, 0, 200);

  TH1D *hdata = (TH1D*) fdata->Get("pileup");
  TH1D *hmc   = (TH1D*) fmc->Get("ggNtuplizer/hPUTrue");

  hdata->Rebin(5);
  hmc->Rebin(5);

  hdata->Scale(1./hdata->Integral(-1, -1));
  hmc->Scale(1./hmc->Integral(-1, -1));

  hnew->Divide(hdata, hmc);
  cout<<hdata->Integral(-1, -1)<<endl;
  cout<<hmc->Integral(-1, -1)<<endl;
  for (Int_t i=1; i<=200; ++i) {
    cout<<hdata->GetBinContent(i)<<" "<<hmc->GetBinContent(i)<<" "<<hdata->GetBinContent(i)/hmc->GetBinContent(i)<<endl;
  }

  if (json == 0) sprintf(fname, "PU_histo_13TeV_GoldenJSON_%dmb.root", xs);
  if (json == 1) sprintf(fname, "PU_histo_13TeV_SilverJSON_%dmb.root", xs);
  TFile *fout = new TFile(fname, "recreate");
  fout->cd();
  hnew->Write();
  //hdata->Write();
  //hmc->Write();
  fout->Close();

}
