################################################################################
# script to hadd ntuples from JETM4 data and JETM4 Vgamma MC
################################################################################

path=../OutputNtuples/v1.3_v00/

echo "combining samples at " $path


#hadd -f ${path}/gdata/data15-16_Vgamma_merged_processed.root  ${path}/gdata/data15-16_merged_processed.root ${path}/gmc16a/Vgamma_merged_processed.root >! ${path}/gdata/data15-16_Vgamma_merged_processed.log 2>&1 &

#hadd -f ${path}/gdata/data17_Vgamma_merged_processed.root  ${path}/gdata/data17_merged_processed.root ${path}/gmc16cd/Vgamma_merged_processed.root >! ${path}/gdata/data17_Vgamma_merged_processed.log 2>&1 &

rm -f ${path}/gdata/data18_Vgamma_merged_processed.log
hadd -f ${path}/gdata/data18_Vgamma_merged_processed.root  ${path}/gdata/data18_merged_processed.root ${path}/gmc16e/Vgamma_merged_processed.root > ${path}/gdata/data18_Vgamma_merged_processed.log 2>&1 &

