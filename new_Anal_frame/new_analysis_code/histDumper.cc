// vim: set sts=4 sw=4 fdm=marker et:
#include <iostream>
#include <array>
#include "math.h"
#include "cstdio.h"
#include "gROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "CompartmentObject.h"
#include "waferGeom.h"
#include "TH1D.h"
#include "TProfile.h"

const int G_nLayers = 8;
const int G_inputEnergy = 125;


//int nHits;

void regTree()
{//{{{
    // Copy the settings from makePlots
    std::string treeName = "HGCLayer";
    TChain *ch = new TChain(treeName);
    ch -> Add("./HGCLayerTree.root");
    
}//}}}

void histDumper(){
    // Follow the main function in histDumper.py
    regTree();
    TFile *fout = TFile("dumpedHist.root",'RECREATE');
    // Import simple plots to be drawn
    
    TFile *handle = open("./histDef.py");
    exec(handle);
    handle -> Close();


    
    //For complicate plots 
    TH1D *h1_erawHit = new TH1D();
    std::array<int,>ringRadii = ('f',{10,19,29,34,38,51,58,68,77});
    TProfile *pr_twoPtCor_rings = new ("pr_twoPtCor_rings","",9,ringRadii);
    
    Long64_t nentries = ch -> GetEntries();
    for (int evt = 0; ){
        for (){
            for (){
            }
        }
    }

    fout -> Write();
    fout -> Close();
}

int main(){
    histDumper();
    return 0;
}
