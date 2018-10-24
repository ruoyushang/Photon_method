#!/bin/bash

#cd ../Root 

#smearing_mode="DataSmear"
smearing_mode="NoSmear"
#variable="Z_pt"
#variable="METl"
#variable="METt"
echo "Doing smearing mode" $smearing_mode

#for variable in "Z_pt"

#for period in data15-16 data17
#for period in data18
#for period in data15-16
#for period in data15-16 data15-16_Vgamma data17 data17_Vgamma data18 data18_Vgamma

for period in data15-16_Vgamma data17_Vgamma data18_Vgamma	      
do
    #for variable in "nVtx" "mu"
    #for variable in "MET"
    #for variable in "Z_pt" "MET"
    for variable in "bjet_n" "jet_n" "METl" "METt" "MET" 
    #for variable in "MET"
    do
	#for channel in "ee" 
	for channel in "ee" "mm" "em"
	#for channel in "mm"
	do	    
	    #echo "Variable channel" $variable $channel

	    ### no smearing
	    echo ${period} ${channel} ${variable} 
	    rm -f logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}.log
	    root -l -b -q '../Root/quickDraw_Data.C+("'${period}'","'${channel}'","'${variable}'")' > logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}.log 2>&1 &


#	    if [ "${channel}" == 'mm' ]
#	    then
#		### MC smearing
#		echo ${period} ${channel} ${variable} "McSmear"
#		rm -f logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}_McSmear.log
#		root -l -b -q '../Root/quickDraw_Data.C+("'${period}'","'${channel}'","'${variable}'","McSmear")' > logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}_McSmear.log 2>&1 &
#		
#		### data smearing
#		echo ${period} ${channel} ${variable} "DataSmear"
#		rm -f logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}_DataSmear.log
#		root -l -b -q '../Root/quickDraw_Data.C+("'${period}'","'${channel}'","'${variable}'","DataSmear")' > logfiles/run_Step5_quickDraw_Data_${period}_${channel}_${variable}_DataSmear.log 2>&1 &
#	    fi
#


	    #sleep 10
	    
	    #echo ${period} ${channel} ${variable} "no WITH subtraction"
	    #rm -f logfiles/run_gdata_quickDraw_Data_${period}_${channel}_${variable}_VgSubtracted.log
	    #root -l -b -q '../Root/quickDraw_Data.C+("'${period}'","'${channel}'","'${variable}'",true)' > logfiles/run_gdata_quickDraw_Data_${period}_${channel}_${variable}_VgSubtracted.log 2>&1 &
	    #sleep 10
	    
	    #echo "with Vg subtraction"
	    #rm -f logfiles/run_gdata_quickDraw_Data_Data15-16_${channel}_${variable}_VgSubtracted.log
	    #root -l -b -q '../Root/quickDraw_Data.C+("data15-16","'${channel}'","'${variable}'",true)' > logfiles/run_gdata_quickDraw_Data_Data15-16_${channel}_${variable}_VgSubtracted.log 2>&1 &
	    #sleep 15
	    
	done
    done
done


# echo "Doing data15-16 mm channel NoSmear..."
# rm -f logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_NoSmear_${variable}_VgSubtraction.log
# root -l -b -q 'quickDraw_Data.C+("data15-16","mm","'${variable}'","NoSmear")' > logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_NoSmear_${variable}_VgSubtraction.log 2>&1 &
# sleep 15
# 
# echo "Doing data15-16 mm channel McSmear..."
# rm -f logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_McSmear_${variable}_VgSubtraction.log
# root -l -b -q 'quickDraw_Data.C+("data15-16","mm","'${variable}'","McSmear")' > logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_McSmear_${variable}_VgSubtraction.log 2>&1 &
# sleep 15
# 
# echo "Doing data15-16 mm channel DataSmear..."
# rm -f logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_DataSmear_${variable}_VgSubtraction.log
# root -l -b -q 'quickDraw_Data.C+("data15-16","mm","'${variable}'","DataSmear")' > logfiles/runRJR_gdata_quickDraw_Data_Data15-16_mm_DataSmear_${variable}_VgSubtraction.log 2>&1 &
# sleep 15
# 

# echo "Doing data15-16 em channel..."
# rm -f runRJR_gdata_quickDraw_Data_Data15-16_em_${variable}.log
# root -l -b -q 'quickDraw_Data.C+("data15-16","em","'${variable}'")' > runRJR_gdata_quickDraw_Data_Data15-16_em_${variable}.log 2>&1 &
# sleep 30
# 
# echo "Doing data17 ee channel..."
# rm -f runRJR_gdata_quickDraw_Data_Data17_ee_${variable}.log 
# root -l -b -q 'quickDraw_Data.C+("data17","ee","'${variable}'")' > runRJR_gdata_quickDraw_Data_Data17_ee_${variable}.log 2>&1 &
# sleep 30
# 
# echo "Doing data17 mm channel..."
# rm -f runRJR_gdata_quickDraw_Data_Data17_mm_${smearing_mode}_${variable}.log
# root -l -b -q 'quickDraw_Data.C+("data17","mm","'${variable}'","'${smearing_mode}'")' > runRJR_gdata_quickDraw_Data_Data17_mm_${smearing_mode}_${variable}.log 2>&1 &
# sleep 30
# 
# echo "Doing data17 em channel..."
# rm -f runRJR_gdata_quickDraw_Data_Data17_em_${variable}.log
# root -l -b -q 'quickDraw_Data.C+("data17","em","'${variable}'")' > runRJR_gdata_quickDraw_Data_Data17_em_${variable}.log 2>&1 &
# sleep 30

