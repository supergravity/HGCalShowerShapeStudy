//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 15 17:10:51 2011 by ROOT version 5.28/00b
// from TTree H4l /H4l Tree
// found on file: ../MakePlots/H4l.root
//////////////////////////////////////////////////////////

#ifndef makePlots_h
#define makePlots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <vector>
#include <algorithm> //std::sort
#include <functional>
#include "TProfile.h"
#include "TGraphErrors.h"
//#include "TFitResultPtr.h"

// Local objects:
#include "CompartmentObject.h" //hold the layer(s)
class CompartmentObject;

#include "waferGeom.h"
class HexGeometry;


//#ifndef  BASE4MOM
//#define  BASE4MOM
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/GenVector/LorentzVector.h"
typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > Base4Mom;
#include "P4Helpers.h"
//#endif

//const double ENEPERMIP    =0.000052;//GeV (52KeV/MIP) (electron)
const double ENEPERMIP    =  56.31e-06;  // 125pion MPV response
const int    MAXENEBINS   =  1.5;
const int    NLAYERS      =  8;
const int    Ring         =  7;
const int    HIT_LIMIT    =  2350;
const int    HIT_MAX_ARRAY=  4096;
const double MMtoCM       =  0.1;
const double MEVTOGEV     =  0.001;
const double GEVTOMEV     =  1000;
//const double raw_enefactor = 1.205;//(for MIPtoMEV == 53.8)(for electron))
const double raw_enefactor = 1.115; //(for MIPtoMEV == 56.31)(for pion)

//* Prototype Configuration:
//                     Layer:   0   1    2     3    4    5    6    7
const double leadX0 [NLAYERS]= {4.5,0.0 ,0.0 , 0.0 ,0.  ,0.  ,1.0 ,0.  };
const double tungX0 [NLAYERS]= {0.0,2.85,2.85, 2.17,2.51,0.  ,1.37,2.17};
const double coppX0 [NLAYERS]= {0.4,0.4 ,0.4 , 0.4 ,0.  ,0.8 ,0.  ,0.8 };
const double cuwX0  [NLAYERS]= {0.15,0.15,0.15,0.25,0.  ,0.75,0.  ,1.0 };
const double dRGeo1 [7] = {0., 1.12455, 2.24909, 3.37364, 4.49818, 5.84331, 6.84034};
const double NPADS  [7] = {1., 6., 12., 18., 24., 30., 36. };
const double MIP_cut_number [5] = {2,3,4,5,6};




using namespace std;



class makePlots
{
public :

    makePlots (TChain* inchain);
    ~makePlots();

    //PFOInfoBranches PFO();


    //Two point Correlation
    TH1D     *hEiEj_Multi_normal     [NLAYERS];
    TH1D     *hdist_Two_hits         [NLAYERS];
    TH1D     *helayerRaw             [NLAYERS];
    TH1D     *helayerRawtotal        ;
    TH1D     *hene_total_layer       ;
    TH1D     *Pad_ring_hist          [NLAYERS][Ring];
    //TH2D     *EiEj_multi_dis_2D      [NLAYERS];
    //TGraph   *EiEj_multi_dis         [NLAYERS];
    TH2D   *EiEj_nor_rel_dis         [NLAYERS];
    TH2D   *EiEj_nor_dis_pad         [NLAYERS];
    TH2D   *EiEj_nor_dis_avg         [NLAYERS];
    TProfile *p1                     [NLAYERS];   //[hit_size]
    TProfile *EiEj_prof              [NLAYERS];
    //TProfile *p2                     [NLAYERS];  
    TGraphErrors *EiEj_nor_dis_graph       [NLAYERS];
    //======== Getting interlayer showershape profile ========
    TGraph2D     *h2_prof_mc;
    TGraph2D     *h2_prof_data;



    //TFitResultPtr r                  [NLAYERS];
    //Declaration of Canvas
    TCanvas *c1                    ;


    bool FillResolHistos();
    void PrintGenInfo();
    double  DiffPhi(double dPhi)
    {
        if (fabs(dPhi) > M_PI) return fabs(2.*M_PI-fabs(dPhi));
        return fabs(dPhi);
    }
    double GetAreaInPadUnits(int iL,double dR)
    {
        double area=1.0;
        //double padsize=3*0.65*0.65*sqrt(3.)/2.; //geom1
        for(int i=0; i<7; i++)
        {
            if(fabs(dR-dRGeo1[i])<0.45)
            {
                area = NPADS[i];
                return area;
            }
        }
        if(dR>3.0 && dR<4.0) return NPADS[3];//are these calib pads?
        if(dR>5.0 && dR<6.0) return NPADS[5];
        if(dR>6.0)           return NPADS[6];
        cout << "no matching dR found ... dR=" << dR << endl;
        return area;
    };



    bool doTruth;// doTruth = 1 MC / = 0 Data


    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain


    //MC:
    // Declaration of leaf types
    vector<float>   *simHitLayEn1E;
    vector<float>   *simHitLayEn2E;
    vector<float>   *simHitLayEn1H;
    vector<float>   *simHitLayEn2H;
    Double_t        xBeam;
    Double_t        yBeam;
    Double_t        zBeam;
    Double_t        pBeam;
    vector<unsigned int> *simHitCellIdE;
    vector<float>   *simHitCellEnE;
    vector<unsigned int> *simHitCellIdH;
    vector<float>   *simHitCellEnH;

    // List of branches
    TBranch        *b_simHitLayEn1E;   //!
    TBranch        *b_simHitLayEn2E;   //!
    TBranch        *b_simHitLayEn1H;   //!
    TBranch        *b_simHitLayEn2H;   //!
    TBranch        *b_xBeam;   //!
    TBranch        *b_yBeam;   //!
    TBranch        *b_zBeam;   //!
    TBranch        *b_pBeam;   //!
    TBranch        *b_simHitCellIdE;   //!
    TBranch        *b_simHitCellEnE;   //!
    TBranch        *b_simHitCellIdH;   //!
    TBranch        *b_simHitCellEnH;   //!

    //Data:
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


    virtual Int_t    Cut(Long64_t entry);
    void     Init();
    void     Loop();
    virtual void     Show(Long64_t entry = -1);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual Bool_t   Notify();

    //TFile* outfile;

private:

    Int_t           hit_size;
    Int_t           hit_layer[HIT_MAX_ARRAY];   //[hit_size]
    Double_t        hit_EMenergy[HIT_MAX_ARRAY];   //[hit_size]
    Double_t        hit_Hadenergy[HIT_MAX_ARRAY];   //[hit_size]
    Double_t        hit_Rawenergy[HIT_MAX_ARRAY];   //[hit_size]
    Double_t        hit_Eta[HIT_MAX_ARRAY];   //[hit_Eta]
    Double_t        hit_Phi[HIT_MAX_ARRAY];   //[hit_Phi]
    Double_t        hit_X  [HIT_MAX_ARRAY];   //
    Double_t        hit_Y  [HIT_MAX_ARRAY];   //
    Double_t        hit_Z  [HIT_MAX_ARRAY];   //
    Int_t           cell_N [HIT_MAX_ARRAY];   //

    //Make an NxMIP cut:
    bool keepHit(int ihit)
    {
        double ene   =  0.;
        //int    celln =  cell_N[ihit];
        //if(celln==2) return false;
        if(doTruth)
        {
            ene = simHitCellEnE->at(ihit);//in GeV
        }
        else
        {
            ene = (raw_enefactor)*energy->at(ihit)*MEVTOGEV;//needs conversion
            //ene = energy->at(ihit)*MEVTOGEV;//needs conversion
        }
        //cout << ene << " " << 2*ENEPERMIP << endl;
        if(ene<MIP_cut_number[0]*ENEPERMIP) return false;
        return true;
    };

    //* Layers:
    std::vector<CompartmentObject*> Layer;

    int buildCompObjects()
    {

        // Uncode x,y,z from input
        //cout << "testing icell size:" << simHitCellIdE->size() << endl;
        HexGeometry geomc(false);//127 cells (133 pads some ganged together)
        if(doTruth)
        {
            for(unsigned int icell=0; icell<simHitCellIdE->size(); icell++)
            {
                int cellno = (simHitCellIdE->at(icell)>>0)&0xFF;
                std::pair<double,double> xy = geomc.position(cellno);
                double lx    =  xy.first;
                double ly    =  xy.second;
                double ene   =  simHitCellEnE->at(icell);
                int    layer = ((simHitCellIdE->at(icell)>>19)&0x7F);
                //cout << "testing icell:" << icell << " x:"<<lx<<" y:"<<ly
                //     <<" layer:"<<layer  << " E=" << ene << endl;

                hit_layer    [icell]=layer; //starting from 1
                hit_Rawenergy[icell]=ene;
                hit_X        [icell]=lx;
                hit_Y        [icell]=ly;
                cell_N       [icell]=cellno;
            }
        }
        else
        {
            for(unsigned int icell=0; icell<cellID->size(); icell++)
            {
                //is the following correct for data?
                int cellno = (cellID->at(icell)>>0)&0xFF;
                std::pair<double,double> xy = geomc.position(cellno);
                double lx    =  xy.first;
                double ly    =  xy.second;
                //double z     =  xy.third;
                double ene   =  raw_enefactor*energy->at(icell);
                //double ene   =  energy->at(icell);
                int    layer = ((cellID->at(icell)>>24)&0x7F);
                //cout << "testing icell:" << icell << " x:"<<lx<<" y:"<<ly
                //     <<" layer:"<<layer  << " E=" << ene << endl;
                //cout <<x->at(icell) << " " << y->at(icell) << endl;
                hit_layer    [icell]=layer; //starting from 1
                hit_Rawenergy[icell]=ene*MEVTOGEV;
                hit_X        [icell]=x->at(icell)*10.;
                hit_Y        [icell]=y->at(icell)*10.;
                cell_N       [icell]=cellno;
            }
        }
        //getchar();

        // Find the seed
        for(int iL=0; iL<NLAYERS; ++iL)
        {
            double maxene=0.;
            double maxeta=-999.;
            double maxphi=-999.;
            double maxX  =-999.;
            double maxY  =-999.;
            double maxZ  =-999.;
            if(doTruth)
            {
                for(unsigned int ih=0; ih<simHitCellIdE->size(); ih++)
                {
                    if(!keepHit(ih)) continue;
                    if(hit_layer[ih]!=(iL+1)) continue;
                    if(hit_Rawenergy[ih]>maxene)
                    {
                        maxene=hit_Rawenergy[ih];
                        maxX  =hit_X[ih];
                        maxY  =hit_Y[ih];
                    }
                }
            }
            else
            {
                for(unsigned int ih=0; ih<cellID->size(); ih++)
                {
                    if(!keepHit(ih)) continue;
                    if(hit_layer[ih]!=(iL+1)) continue;
                    if(hit_Rawenergy[ih]>maxene)
                    {
                        maxene=hit_Rawenergy[ih];
                        maxX  =hit_X[ih];
                        maxY  =hit_Y[ih];
                    }
                }
            }

            CompartmentObject *tmp = new CompartmentObject;
            int nhitsInLayer=0;
            //cout << maxene << " " << maxene/ENEPERMIP << endl;
            //tmp->SetEneMax(maxene/ENEPERMIP);
            tmp->SetEneMax(maxene);
            tmp->SetEtaMax(maxeta);
            tmp->SetPhiMax(maxphi);
            tmp->SetXMax  (maxX  );
            tmp->SetYMax  (maxY  );
            tmp->SetZMax  (maxZ  );
            if(doTruth)
            {
                for(unsigned int ih=0; ih<simHitCellIdE->size(); ih++)
                {
                    if(hit_layer[ih]==(iL+1))
                    {
                        tmp->SetCellNo(cell_N[ih],nhitsInLayer);
                        tmp->SetEtaHit(hit_Eta[ih],nhitsInLayer);
                        tmp->SetPhiHit(hit_Phi[ih],nhitsInLayer);
                        tmp->SetXHit(hit_X[ih],nhitsInLayer);
                        tmp->SetYHit(hit_Y[ih],nhitsInLayer);
                        tmp->SetZHit(hit_Z[ih],nhitsInLayer);
                        //tmp->SetErawHit(hit_Rawenergy[ih]/ENEPERMIP,nhitsInLayer);
                        tmp->SetErawHit(hit_Rawenergy[ih],nhitsInLayer);
                        nhitsInLayer++;
                    }
                }
            }
            else
            {
                for(unsigned int ih=0; ih<cellID->size(); ih++)
                {
                    if(hit_layer[ih]==(iL+1))
                    {
                        tmp->SetCellNo(cell_N[ih],nhitsInLayer);
                        tmp->SetEtaHit(hit_Eta[ih],nhitsInLayer);
                        tmp->SetPhiHit(hit_Phi[ih],nhitsInLayer);
                        tmp->SetXHit(hit_X[ih],nhitsInLayer);
                        tmp->SetYHit(hit_Y[ih],nhitsInLayer);
                        tmp->SetZHit(hit_Z[ih],nhitsInLayer);
                        tmp->SetErawHit(hit_Rawenergy[ih],nhitsInLayer);
                        nhitsInLayer++;
                    }
                }
            }
            tmp->SetLeadX0(leadX0[iL]);
            tmp->SetWX0   (tungX0[iL]);
            tmp->SetCuX0  (coppX0[iL]);
            tmp->SetCuWX0 ( cuwX0[iL]);
            tmp->SetNhits(nhitsInLayer);
            Layer.push_back(tmp);
        }
        return 0;
    };

    void Clear()
    {
        //We need to clear all the pointers inside the vectors
        //for (unsigned int i=0;i<PFOs.size();++i){
        //  delete PFOs.at(i);
        //}
        //PFOs.clear();
        for ( int i=0;i<Layer.size(); ++i )
            delete Layer[i];
        Layer.clear();
    }




    //* Ntuple
    TTree *outTree;
    float m_Mgg;
    float m_Dyjj;
    float m_Mjj;
    float m_Etaj1;
    float m_Etastar;
    float m_dPhi2g2j;
    float m_weight;
    float m_BDTG;

private:


};

#endif

