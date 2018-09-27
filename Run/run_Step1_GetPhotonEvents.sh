###---------------------------------
### JETM4 data
###---------------------------------

datapath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/JETM4_lowpt/JETM4_Data/JETM4_Data_SEPT12_v2/merged/'

for period in data15-16 data17 data18
do
    echo "Doing JETM4 data " $period
    treename=$period
    #treename='tree_NoSys'
    #treename='data'
    #treename='tree_NoSys'
    rm -f logfiles/run_Step1_GetPhotonEvents_JETM4_Data_${period}.log
    root -l -b -q '../Root/GetPhotonEvents.C+("'${period}'_merged_processed","gdata","'${datapath}'",1,"'${treename}'")' > logfiles/run_Step1_GetPhotonEvents_JETM4_Data_${period}.log 2>&1 &
done


###---------------------------------
### JETM4 mc16a
###---------------------------------

mc16apath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/JETM4_lowpt/JETM4_mc16a/JETM4_mc16a_SEPT12_v2/merged/'


#for sample in SinglePhoton222
for sample in Vgamma SinglePhoton222 SinglePhoton211
do
    echo "Doing JETM4 mc16a " $sample
    rm -f logfiles/run_Step1_GetPhotonEvents_JETM4_mc16a_${sample}.log
    root -l -b -q '../Root/GetPhotonEvents.C+("'${sample}'_merged_processed","gmc16a","'${mc16apath}'",0,"'${sample}'_NoSys")' > logfiles/run_Step1_GetPhotonEvents_JETM4_mc16a_${sample}.log 2>&1 &
done


###---------------------------------
### JETM4 mc16cd
###---------------------------------

mc16cdpath='/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.3/JETM4_lowpt/JETM4_mc16cd/JETM4_mc16cd_SEPT12_v1_COPY/merged/'

for sample in Vgamma SinglePhoton222 SinglePhoton211
do
    echo "Doing JETM4 mc16cd " $sample
    rm -f logfiles/run_Step1_GetPhotonEvents_JETM4_mc16cd_${sample}.log
    root -l -b -q '../Root/GetPhotonEvents.C+("'${sample}'_merged_processed","gmc16cd","'${mc16cdpath}'",0,"'${sample}'_NoSys")' > logfiles/run_Step1_GetPhotonEvents_JETM4_mc16cd_${sample}.log 2>&1 &
done


