#!/usr/bin/env python

import math

from ROOT import gROOT, TChain, TTree, TFile
from ROOT import TH1D, TProfile
from array import array

gROOT.LoadMacro("./CompartmentObject.h")
gROOT.LoadMacro("./waferGeom.h")
from ROOT import CompartmentObject, HexGeometry

G_nLayers       = 8
G_inputEnergy   = 125 

def main(treeName="HGCLayer"):

    # Chain input files.
    ch = TChain(treeName)
    ch.Add("./HGCLayerTree.root")

    # 
    fout = TFile("dumpedHist.root",'RECREATE')

    # Import simple plots to be drawn
    plots = []
    handle = open("./histDef.py")
    exec(handle)
    handle.close()

    for p in plots:
        for step in p[1:-1]:
            ch.Draw(step[0],step[1],"goff")
            pass

    # For complicate plots.
    h1_erawHit = TH1D("h1_erawHit","",8,1,9)
    ringRadii  = array('f',[10,19,29,34,38,50,51,58,68,77])
    pr_twoPtCor_rings = TProfile("pr_twoPtCor_rings","",9,ringRadii)

    for ev in ch:
        for iL in range(G_nLayers):
            layerBranch = getattr(ev,"layer{0}".format(iL))
            for iHit in range(layerBranch.GetNhits()):
                h1_erawHit.Fill(iL+1,layerBranch.GetErawHit(iHit))
                for jHit in range(iHit+1,layerBranch.GetNhits()):
                    dist = math.sqrt((layerBranch.GetXHit(iHit)-layerBranch.GetXHit(jHit))**2
                                    +(layerBranch.GetYHit(iHit)-layerBranch.GetYHit(jHit))**2)
                    pr_twoPtCor_rings.Fill(dist,layerBranch.GetErawHit(iHit)*layerBranch.GetErawHit(jHit))

    fout.Write()
    fout.Close()

if __name__ == '__main__':
    main()
