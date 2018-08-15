cp  ./MC_electron125GeV.txt MC.txt
make 
./makePlots
#mv output.root output_electrontwocorrMC100GeV.root 
mv output.root output_test_MC_125GeV.root 
