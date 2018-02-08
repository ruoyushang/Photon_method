hadd -f gmc/gmc_raw.root gmc/3*.root
hadd -f gdata/gdata_raw.root gdata/period*.root Vg/3*.root
hadd -f zjets/zjets_ee.root zjets/3*_ee.root
hadd -f zjets/zjets_mm.root zjets/3*_mm.root
hadd -f zjets_221/zjets_ee.root zjets_221/3*_ee.root
hadd -f zjets_221/zjets_mm.root zjets_221/3*_mm.root
hadd -f zjets_tt/ztt_ee.root zjets_tt/3*_ee.root
hadd -f zjets_tt/ztt_mm.root zjets_tt/3*_mm.root
hadd -f tt/ttee.root tt/4*_ee.root
hadd -f tt/ttmm.root tt/4*_mm.root
hadd -f vv/vvee.root vv/3*_ee.root
hadd -f vv/vvmm.root vv/3*_mm.root
hadd -f data/data_ee.root data/period*_ee.root
hadd -f data/data_mm.root data/period*_mm.root

hadd -f bphys/bphys_raw.root bphys/data.root
