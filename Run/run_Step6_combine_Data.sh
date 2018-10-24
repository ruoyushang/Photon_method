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
#for variable in "nVtx" "mu"
for variable in "MET" "METl" "METt" "jet_n" "bjet_n" 
#for variable in "Z_pt" "MET" "METl" "METt" "jet_n" "HT"
#for variable in "bjet_n" "jet_n" "MET" "Z_pt" "mll"
#for variable in "MET" "Z_pt" 
do
    for channel in "ee" "mm" "em"
		   #for channel in "ee" "mm" "em"
		   #for channel in "mm"
    do	    
	#echo "Variable channel" $variable $channel
	
	echo ${channel} ${variable} "no Vg subtraction"
	rm -f logfiles/run_Step6_combine_Data_${channel}_${variable}.log
	root -l -b -q '../Root/combine_Data.C+("'${channel}'","'${variable}'")' > logfiles/run_Step6_combine_Data_${channel}_${variable}.log 2>&1 &
	#sleep 10X
	
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

