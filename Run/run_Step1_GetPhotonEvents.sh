###---------------------------------
### JETM4 data
###---------------------------------

datapath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/JETM4/JETM4_Data/JETM4_Data_SEPT12_v1/merged/'

echo "Doing data15-16..."
rm -f logfiles/run_Step1_GetPhotonEvents_data15-16.log
root -l -b -q '../Root/GetPhotonEvents.C+("data15-16_merged_processed","gdata","'${datapath}'",1,"data15-16")' > logfiles/run_Step1_GetPhotonEvents_data15-16.log 2>&1 &

echo "Doing data17..."
rm -f logfiles/run_Step1_GetPhotonEvents_data17.log
root -l -b -q '../Root/GetPhotonEvents.C+("data17_merged_processed","gdata","'${datapath}'",1,"data17")' > logfiles/run_Step1_GGetPhotonEvents_data17.log 2>&1 &

echo "Doing data18..."
rm -f logfiles/run_Step1_GetPhotonEvents_data18.log
root -l -b -q '../Root/GetPhotonEvents.C+("data18_merged_processed","gdata","'${datapath}'",1,"data18")' > logfiles/run_Step1_GGetPhotonEvents_data18.log 2>&1 &

