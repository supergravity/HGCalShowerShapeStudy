// vim: set sw=4 sts=4 fdm=marker et:

#include <fstream>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "./CompartmentObject.h"
#include "./waferGeom.h"

class LayerBuilder
{//{{{
    public:
        LayerBuilder(bool isTruth);
        ~LayerBuilder();
        CompartmentObject* build(int);// dependency on layer id

        //const double ENEPERMIP    =0.000052;//GeV (52KeV/MIP) (electron)
        const double ENEPERMIP    =  56.31e-06;  // 125pion MPV response
        const int    MAXENEBINS   =  int(1.5);
        static const int NLAYERS  =  8;
        const int    Ring         =  7;
        const int    HIT_LIMIT    =  2350;
        static const int    HIT_MAX_ARRAY=  4096;
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
            
        //  Variables used in processing
        bool doTruth;
        TChain *ch;
        CompartmentObject* layers[NLAYERS];

    private:
        bool keepHit(int);

        // Branches in inputTree
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

        int    hit_layer     [HIT_MAX_ARRAY];
        double hit_Rawenergy [HIT_MAX_ARRAY];
        double hit_X         [HIT_MAX_ARRAY];
        double hit_Y         [HIT_MAX_ARRAY];
        double hit_Z         [HIT_MAX_ARRAY];
        double hit_Eta       [HIT_MAX_ARRAY];
        double hit_Phi       [HIT_MAX_ARRAY];
        int    cell_N        [HIT_MAX_ARRAY];

};//}}} 

LayerBuilder::LayerBuilder(bool isTruth)
{//{{{
    doTruth = isTruth;
    for(int iL=0; iL<NLAYERS; iL++){
        layers[iL] = new CompartmentObject;
    }

    //Initialization of input tree
    std::string rootDir;
    doTruth ? rootDir = "HGCalTBAnalyzer/HGCTB" : rootDir = "hgcaltbntuple/HGC_Events";
    ch = new TChain(rootDir.c_str());
    
    std::ifstream inTxtFile;
    inTxtFile.open("MC.txt");
    while(!inTxtFile.eof()){
        std::string fileName;
        std::getline(inTxtFile, fileName);
        if (fileName.size()<2){
            break;
        }
        ch->Add(fileName.c_str());
    }

    simHitLayEn1E  = 0;
    simHitLayEn2E  = 0;
    simHitLayEn1H  = 0;
    simHitLayEn2H  = 0;
    simHitCellIdE  = 0;
    simHitCellEnE  = 0;
    simHitCellIdH  = 0;
    simHitCellEnH  = 0;
    cellID         = 0;
    x              = 0;
    y              = 0;
    z              = 0;
    energy         = 0;
    cluster_x      = 0;
    cluster_y      = 0;
    cluster_z      = 0;
    cluster_energy = 0;
    cluster_size   = 0;

    if(doTruth) {
        ch->SetBranchAddress("simHitLayEn1E"  , &simHitLayEn1E ) ;
        ch->SetBranchAddress("simHitLayEn2E"  , &simHitLayEn2E ) ;
        ch->SetBranchAddress("simHitLayEn1H"  , &simHitLayEn1H ) ;
        ch->SetBranchAddress("simHitLayEn2H"  , &simHitLayEn2H ) ;
        ch->SetBranchAddress("xBeam"          , &xBeam         ) ;
        ch->SetBranchAddress("yBeam"          , &yBeam         ) ;
        ch->SetBranchAddress("zBeam"          , &zBeam         ) ;
        ch->SetBranchAddress("pBeam"          , &pBeam         ) ;
        ch->SetBranchAddress("simHitCellIdE"  , &simHitCellIdE ) ;
        ch->SetBranchAddress("simHitCellEnE"  , &simHitCellEnE ) ;
        ch->SetBranchAddress("simHitCellIdH"  , &simHitCellIdH ) ;
        ch->SetBranchAddress("simHitCellEnH"  , &simHitCellEnH ) ;
    } else {
        ch->SetBranchAddress("evtID"          , &evtID         ) ;
        ch->SetBranchAddress("nhit"           , &nhit          ) ;
        ch->SetBranchAddress("cellID"         , &cellID        ) ;
        //ch->SetBranchAddress("cellID"       , &simHitCellIdE ) ;
        ch->SetBranchAddress("x"              , &x             ) ;
        ch->SetBranchAddress("y"              , &y             ) ;
        ch->SetBranchAddress("z"              , &z             ) ;
        ch->SetBranchAddress("energy"         , &energy        ) ;
        //ch->SetBranchAddress("energy"       , &simHitCellIdE ) ;
        ch->SetBranchAddress("thrustX0"       , &thrustX0      ) ;
        ch->SetBranchAddress("thrustX"        , &thrustX       ) ;
        ch->SetBranchAddress("thrustY0"       , &thrustY0      ) ;
        ch->SetBranchAddress("thrustY"        , &thrustY       ) ;
        ch->SetBranchAddress("cluster_x"      , &cluster_x     ) ;
        ch->SetBranchAddress("cluster_y"      , &cluster_y     ) ;
        ch->SetBranchAddress("cluster_z"      , &cluster_z     ) ;
        ch->SetBranchAddress("cluster_energy" , &cluster_energy) ;
        ch->SetBranchAddress("cluster_size"   , &cluster_size  ) ;
    }
}//}}} 

LayerBuilder::~LayerBuilder()
{//{{{
}//}}}

CompartmentObject* LayerBuilder::build(int iL)
{//{{{
    // Uncode x,y,z from input
    //cout << "testing icell size:" << simHitCellIdE->size() << endl;
    HexGeometry geomc(false);//127 cells (133 pads some ganged together)
    if(doTruth){
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
    }else{
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

    CompartmentObject *tmp = layers[iL];
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
    }else{
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

    return tmp;
}//}}}

bool LayerBuilder::keepHit(int iHit)
{//{{{

        double ene   =  0.;
        //int    celln =  cell_N[ihit];
        //if(celln==2) return false;
        if(doTruth){
            ene = simHitCellEnE->at(iHit);//in GeV
        }else{
            ene = (raw_enefactor)*energy->at(iHit)*MEVTOGEV;//needs conversion
        }
        if(ene<MIP_cut_number[0]*ENEPERMIP){
            return false;
        }
        return true;
}//}}}

void layerBuilder(bool isTruth =    true){
    LayerBuilder builder(isTruth);

    TFile *fout = new TFile("HGCLayerTree.root","RECREATE");
    TTree *tout = new TTree("HGCLayer","");
    for(int iL=0; iL<builder.NLAYERS; iL++){
        tout->Branch(TString::Format("layer%d",iL).Data(),builder.layers[iL]);
    } 

    // Event loop
    for(unsigned int entry=0; entry < builder.ch->GetEntries(); entry++){
        builder.ch->GetEntry(entry);
        for(int iL=0; iL<builder.NLAYERS; iL++){
            builder.build(iL);
        }
        tout->Fill();
    }

    fout->Write();
    fout->Close();
}

int main(){
    layerBuilder();
}
