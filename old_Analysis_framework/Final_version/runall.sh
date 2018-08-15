cp  ./MC_electron100GeV.txt MC.txt
make 
./makePlots
#mv output.root output_electrontwocorrMC100GeV.root 
mv output.root output_test_MC_3.root 

cp  ./Data_electron100GeV.txt  MC.txt
make 
./makePlots -t 0 
#mv output.root output_electronDatatwocorr100GeV.root 
mv output.root output_test_data_3.root 

