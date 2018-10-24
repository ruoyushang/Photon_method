#!/bin/bash

#int smearing_method = 0; // No smearing
#int smearing_method = 1; // MC smearing
#int smearing_method = 2; // Data smearing
#int smearing_method = 3; // Truth smearing


#for period in data15-16 data17 data18
#for period in data15-16_Vgamma

	      
#for period in data15-16 data15-16_Vgamma

for period in data17 data17_Vgamma data18 data18_Vgamma
do
    
    echo "Smearing" $period "ee channel NoSmear..."
    rm -f logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_ee.log
    root -l -b -q '../Root/GetPhotonSmearing.C+("'${period}'_merged_processed","ee",1,"'${period}'",0)' > logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_ee.log 2>&1 &
    
    echo "Smearing" $period "mm channel NoSmear..."
    rm -f logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_NoSmear.log
    root -l -b -q '../Root/GetPhotonSmearing.C+("'${period}'_merged_processed","mm",1,"'${period}'",0)' > logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_NoSmear.log 2>&1 &
    
    echo "Smearing" ${period} "mm channel McSmear..."
    rm -f logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_McSmear.log
    root -l -b -q '../Root/GetPhotonSmearing.C+("'${period}'_merged_processed","mm",1,"'${period}'",4)' > logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_McSmear.log 2>&1 &
    
    echo "Smearing" ${period} "mm channel DataSmear..."
    rm -f logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_DataSmear.log
    root -l -b -q '../Root/GetPhotonSmearing.C+("'${period}'_merged_processed","mm",1,"'${period}'",5)' > logfiles/run_Step3_GetPhotonSmearing_JETM4_Data_${period}_mm_DataSmear.log 2>&1 &
done

