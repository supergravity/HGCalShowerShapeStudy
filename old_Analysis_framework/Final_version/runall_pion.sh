cp  ./MC_pion125GeV.txt MC.txt
make 
./makePlots
mv output.root output_piontwocorrMC125GeV_2.root 

cp  ./Data_pion125GeV.txt  MC.txt
make 
./makePlots -t 0 
mv output.root output_piontwocorrData125GeV_2.root 


