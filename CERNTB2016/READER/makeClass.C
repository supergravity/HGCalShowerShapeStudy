//**********************************************************************
// A Macro to automatically create a skeleton for reading a tree
// and parsing the various leaves.
// When executed, this macro produces a macro "analysis.C" ready to use
//**********************************************************************
{
  gROOT->Reset();

  gROOT->SetStyle("Plain");

  //TChain* tt = new TChain("h1");

  //tt->AddFile("../files/geant3.21/jets/dc1.002001.lumi10.recon.008.00001.hlt.pythia_jet_25.eg8.602.root",-1);
  //tt->AddFile("../files/geant3.21/jets/dc1.002001.lumi10.recon.008.00002.hlt.pythia_jet_25.eg8.602.root",-1);


  //
  // Create a file containing the root object
  //
  //TFile f("/home/paganis/Analysis/ATLAS_Analysis/src/NtupleFilter/D3PD_0.0125.0.root");
  //TFile f("/home/paganis/Data/2011/mc10_7TeV.107652.AlpgenJimmyZeeNp2_pt20.merge.NTUP_HSG2.e737_s933_s946_r2215_r2260_p537_tid312078_00/NTUP_HSG2.312078._000013.root.1");
  //TFile f("../data/Single22_CMSSW_6_2_0_SLHC25_patch6/RECO-PU0/paganis_SimHits_0.root");
  //TFile f("../data/tmp/TBGenSim_142.root");
  TFile f("../data/DATA/Output.root");

  //
  // List the contents of the file (usually it contains a TDirectory)
  //
  f.ls();

  //  
  // List the pwd
  //
  f.pwd();
  //
  // Read Tree h3333 in memory. Notice that h3333 tree is under CBNT dir
  //
  //TTree *tt = (TTree*)f.Get("CollectionTree"); //gets the run info block
  //TTree *tt = (TTree*)f.Get("physics"); //gets the variables
  //TTree *tt = (TTree*)f.Get("h4lTree"); //gets the variables
  //TTree *tt = (TTree*)f.Get("HGCalTBAnalyzer/HGCTB"); //gets the variables
  TTree *tt = (TTree*)f.Get("hgcaltbntuple/HGC_Events");

  tt->Print();



  //
  // Make the automatic parsing code if you like
  //
  //gROOT->LoadMacro("MakeCode.C");
  //MakeCode(tt,"analysis.C");
  tt->MakeClass("analysis");
}
