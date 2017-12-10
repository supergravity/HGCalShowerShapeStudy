#ifndef COMPOBJECTH
#define COMPOBJECTH

#include <vector>

const int MAXCHITS=1000;

//Object to hold a compartment of N-calolayers
class CompartmentObject
{

public:


    CompartmentObject()
    {
        X0depth[0]=5.05;
        X0depth[1]=8.45;
        X0depth[2]=11.85;
        X0depth[3]=14.67;
        X0depth[4]=17.18;
        X0depth[5]=18.73;
        X0depth[6]=21.10;
        X0depth[7]=27.07;
    }

    ~CompartmentObject() {}

    //Setters
    void SetNhits (int    d)
    {
        nHits =d;
    };
    void SetEneTot(double d)
    {
        eneTot=d;
    };
    void SetEtaMax(double d)
    {
        etaMax=d;
    };
    void SetPhiMax(double d)
    {
        phiMax=d;
    };
    void SetXMax  (double d)
    {
        XMax  =d;
    };
    void SetYMax  (double d)
    {
        YMax  =d;
    };
    void SetZMax  (double d)
    {
        ZMax  =d;
    };
    void SetEneMax(double d)
    {
        eneMax=d;
    };
    void SetX0    (double d)
    {
        X0    =d;
    };
    void SetShowerDepth(double d)
    {
        shDepthInLayers = d;
    };
    void SetCellNo (int    d, int i)
    {
        cellNo [i]=d;
    };
    void SetEtaHit (double d, int i)
    {
        etaHit [i]=d;
    };
    void SetPhiHit (double d, int i)
    {
        phiHit [i]=d;
    };
    void SetXHit   (double d, int i)
    {
        XHit   [i]=d;
    };
    void SetYHit   (double d, int i)
    {
        YHit   [i]=d;
    };
    void SetZHit   (double d, int i)
    {
        ZHit   [i]=d;
    };
    void SetErawHit(double d, int i)
    {
        erawHit[i]=d;
    };
    void SetLeadX0 (double d)
    {
        leadX0=d;
    };
    void SetWX0    (double d)
    {
        wX0   =d;
    };
    void SetCuX0   (double d)
    {
        cuX0  =d;
    };
    void SetCuWX0  (double d)
    {
        cuwX0 =d;
    };

    //Getters
    int    GetNhits ()
    {
        return nHits;
    };
    double GetEneTot()
    {
        return eneTot;
    };
    double GetEtaMax()
    {
        return etaMax;
    };
    double GetPhiMax()
    {
        return phiMax;
    };
    double GetEneMax()
    {
        return eneMax;
    };
    double GetXMax  ()
    {
        return XMax  ;
    };
    double GetYMax  ()
    {
        return YMax  ;
    };
    double GetXHit  (int i)
    {
        return XHit[i];
    };
    double GetYHit  (int i)
    {
        return YHit[i];
    };
    double GetX0    ()
    {
        return X0    ;
    };
    double GetX0depth(int iLayer)
    {
        return X0depth[iLayer];
    };
    double GetShowerDepthInX0()
    {
        int nL = (int) shDepthInLayers;
        shDepthInX0 = X0depth[nL];
        return shDepthInX0;
    };
    int    GetCellNo (int i)
    {
        return cellNo[i];
    };
    double GetEtaHit (int i)
    {
        return etaHit[i];
    };
    double GetPhiHit (int i)
    {
        return phiHit[i];
    };
    double GetErawHit(int i)
    {
        return erawHit[i];
    };

    double GetDX  (int i)
    {
        return (XHit[i] - XMax);
    };
    double GetDY  (int i)
    {
        return (YHit[i] - YMax);
    };
    double GetDZ  (int i)
    {
        return (ZHit[i] - ZMax);
    };
    double GetDRxy(int i)
    {
        return sqrt(GetDX(i)*GetDX(i)+
                    GetDY(i)*GetDY(i));
    };
    double GetDeta(int i)
    {
        return (etaHit[i] - etaMax);
    };
    double GetDphi(int i)
    {
        return DiffPhi(phiHit[i] - phiMax);
    };
    double GetDR  (int i)
    {
        return sqrt(GetDeta(i)*GetDeta(i)+
                    GetDeta(i)*GetDeta(i));
    };

    double GetLeadX0 ()
    {
        return leadX0;
    };
    double GetWX0    ()
    {
        return wX0;
    };
    double GetCuX0   ()
    {
        return cuX0;
    };
    double GetCuWX0  ()
    {
        return cuwX0;
    };
    //Utils
    double  DiffPhi(double dPhi)
    {
        if (fabs(dPhi) > M_PI) return fabs(2.*M_PI-fabs(dPhi));
        return fabs(dPhi);
    }


private:

    // Total hits in the compartment
    int    nHits;

    // Info about the max Energy cell:
    double etaMax;
    double phiMax;
    double   XMax;
    double   YMax;
    double   ZMax;
    double eneMax;

    // Total energy in the compartment
    double eneTot;

    // X0:
    double X0;
    double X0depth[31];
    double shDepthInLayers;
    double shDepthInX0;
    double leadX0;
    double wX0;
    double cuX0;
    double cuwX0;

    // Raw Hit info:
    double erawHit[MAXCHITS];
    double etaHit [MAXCHITS];
    double phiHit [MAXCHITS];
    double XHit   [MAXCHITS];
    double YHit   [MAXCHITS];
    double ZHit   [MAXCHITS];
    int    cellNo [MAXCHITS];

};


#endif
