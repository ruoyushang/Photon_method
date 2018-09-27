###---------------------------------
### SUSY2 data
###---------------------------------

datapath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/SUSY2/SUSY2_Data/'

for period in data15-16 data17 data18 data15-18
do
    echo "Doing SUSY2 data " $period
    rm -f logfiles/run_Step2_GetBaseLineEvents_SUSY2_Data_${period}.log
    treename='data'
    root -l -b -q '../Root/GetBaseLineEvents.C+("'${period}'_merged_processed","zdata","'${datapath}'",1,"'${treename}'")' > logfiles/run_Step2_GetBaseLineEvents_SUSY2_Data_${period}.log 2>&1 &
done


###---------------------------------
### JETM4 mc16a
###---------------------------------

mc16apath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/SUSY2/SUSY2_Bkgs_mc16a/'

#for sample in diboson
#

#for sample in ttbar_410472_dilep
for sample in Zjets lowMassDY ttbar diboson
do
    echo "Doing JETM4 mc16a " $sample
    rm -f logfiles/run_Step2_GetBaseLineEvents_SUSY2_mc16a_${sample}.log

    # ttbar dilep only
    #root -l -b -q '../Root/GetBaseLineEvents.C+("'${sample}'_processed","zmc16a","'${mc16apath}'",0,"ttbar_dilep_NoSys")' > logfiles/run_Step2_GetBaseLineEvents_SUSY2_mc16a_${sample}.log 2>&1 &

    # nominal
    root -l -b -q '../Root/GetBaseLineEvents.C+("'${sample}'_merged_processed","zmc16a","'${mc16apath}'",0,"'${sample}'_NoSys")' > logfiles/run_Step2_GetBaseLineEvents_SUSY2_mc16a_${sample}.log 2>&1 &
done


###---------------------------------
### JETM4 mc16cd
###---------------------------------

mc16cdpath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/SUSY2/SUSY2_Bkgs_mc16cd/'

for sample in Zjets lowMassDY ttbar ttbar_dilep ttbar_nonallhad diboson
do
    echo "Doing JETM4 mc16cd " $sample
    rm -f logfiles/run_Step2_GetBaseLineEvents_SUSY2_mc16cd_${sample}.log
    root -l -b -q '../Root/GetBaseLineEvents.C+("'${sample}'_merged_processed","zmc16cd","'${mc16cdpath}'",0,"'${sample}'_NoSys")' > logfiles/run_Step2_GetBaseLineEvents_SUSY2_mc16cd_${sample}.log 2>&1 &
done




