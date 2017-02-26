#define myanalysis_cxx
#include "myanalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "commonFunctions.C"

void myanalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L myanalysis.C
//      Root > myanalysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      int nele = 0; ///SJ
      
      for(int iele=0; iele<patelesize; iele++){

	bool flag = eleID2012(iele);

	if(flag) {
	  //cout<<"ievent: iele : Good electron"<<jentry<<":"<<iele<<endl;
	  nele = nele+1; ///SJ
	}
	
      }//for(int iele=0....)
      cout<<"Number of electrons for event number "<<jentry<<" are "<<nele<<endl; ///SJ
   }
}
