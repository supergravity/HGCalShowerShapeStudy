#!/usr/bin/env python

h1_nHits = [ 
    TH1D("h1_nHits","",8,1,9),
    ["1 >>  h1_nHits", "layer0.nHits"],
    ["2 >>+ h1_nHits", "layer1.nHits"],
    ["3 >>+ h1_nHits", "layer2.nHits"],
    ["4 >>+ h1_nHits", "layer3.nHits"],
    ["5 >>+ h1_nHits", "layer4.nHits"],
    ["6 >>+ h1_nHits", "layer5.nHits"],
    ["7 >>+ h1_nHits", "layer6.nHits"],
    ["8 >>+ h1_nHits", "layer7.nHits"],
]
plots.append(h1_nHits)


# h1_erawHit = [ 
#     TH1D("h1_erawHit","",8,1,9),
#     ["1 >>  h1_erawHit", "layer0.erawHit"],
#     ["2 >>+ h1_erawHit", "layer1.erawHit"],
#     ["3 >>+ h1_erawHit", "layer2.erawHit"],
#     ["4 >>+ h1_erawHit", "layer3.erawHit"],
#     ["5 >>+ h1_erawHit", "layer4.erawHit"],
#     ["6 >>+ h1_erawHit", "layer5.erawHit"],
#     ["7 >>+ h1_erawHit", "layer6.erawHit"],
#     ["8 >>+ h1_erawHit", "layer7.erawHit"],
# ]
# plots.append(h1_erawHit)
