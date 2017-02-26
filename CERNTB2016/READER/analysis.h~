//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 29 23:09:34 2016 by ROOT version 5.34/21
// from TTree HGC_Events/HGCAL TB variables tree
// found on file: ../data/DATA/Output.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "/build/bellenot/source/root_v5.34.21/root/cint/cint/lib/prec_stl/vector"
#include "/build/bellenot/source/root_v5.34.21/root/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evtID;
   Int_t           nhit;
   vector<int>     *cellID;
   vector<double>  *x;
   vector<double>  *y;
   vector<double>  *z;
   vector<double>  *energy;
   Float_t         thrustX0;
   Float_t         thrustX;
   Float_t         thrustY0;
   Float_t         thrustY;
   vector<double>  *cluster_x;
   vector<double>  *cluster_y;
   vector<double>  *cluster_z;
   vector<double>  *cluster_energy;
   vector<int>     *cluster_size;

   // List of branches
   TBranch        *b_evtID;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_cellID;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_thrustX0;   //!
   TBranch        *b_thrustX;   //!
   TBranch        *b_thrustY0;   //!
   TBranch        *b_thrustY;   //!
   TBranch        *b_cluster_x;   //!
   TBranch        *b_cluster_y;   //!
   TBranch        *b_cluster_z;   //!
   TBranch        *b_cluster_energy;   //!
   TBranch        *b_cluster_size;   //!

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/DATA/Output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/DATA/Output.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../data/DATA/Output.root:/hgcaltbntuple");
      dir->GetObject("HGC_Events",tree);

   }
   Init(tree);
}

analysis::~analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   cellID = 0;
   x = 0;
   y = 0;
   z = 0;
   energy = 0;
   cluster_x = 0;
   cluster_y = 0;
   cluster_z = 0;
   cluster_energy = 0;
   cluster_size = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evtID", &evtID, &b_evtID);
   fChain->SetBranchAddress("nhit", &nhit, &b_nhit);
   fChain->SetBranchAddress("cellID", &cellID, &b_cellID);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("thrustX0", &thrustX0, &b_thrustX0);
   fChain->SetBranchAddress("thrustX", &thrustX, &b_thrustX);
   fChain->SetBranchAddress("thrustY0", &thrustY0, &b_thrustY0);
   fChain->SetBranchAddress("thrustY", &thrustY, &b_thrustY);
   fChain->SetBranchAddress("cluster_x", &cluster_x, &b_cluster_x);
   fChain->SetBranchAddress("cluster_y", &cluster_y, &b_cluster_y);
   fChain->SetBranchAddress("cluster_z", &cluster_z, &b_cluster_z);
   fChain->SetBranchAddress("cluster_energy", &cluster_energy, &b_cluster_energy);
   fChain->SetBranchAddress("cluster_size", &cluster_size, &b_cluster_size);
   Notify();
}

Bool_t analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
