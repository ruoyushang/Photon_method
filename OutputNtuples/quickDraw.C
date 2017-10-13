
void quickDraw() {

	TChain * tree1 = new TChain("BaselineTree");
	tree1->Add("zjets_221/zjets_ee.root");
	//tree1->Add("zjets/zjets_ee.root");
	tree1->Add("tt/ttee.root");
	tree1->Add("vv/vvee.root");
	tree1->Add("zjets_221/zjets_mm.root");
	//tree1->Add("zjets/zjets_mm.root");
	tree1->Add("tt/ttmm.root");
	tree1->Add("vv/vvmm.root");
	TChain * tree2 = new TChain("BaselineTree");
	tree2->Add("susy/392330_ee.root");
	tree2->Add("susy/392330_mm.root");

	TCut mycut = "";
	//mycut += "jet_n==1";
	mycut += "bjet_n==0";

	mycut += "lep_pT[0]>20 && lep_pT[1]>20";
	mycut += "abs(lep_eta[0])<2.5 && abs(lep_eta[1])<2.5";
	mycut += "mll>81 && mll<101";

	mycut += "jet_n>2";
	mycut += "jet_pT[0]>30 && jet_pT[1]>30 && jet_pT[2]>30";
	//mycut += "mWmin>60 && mWmin<100";
	mycut += "mWmin>20 && mWmin<140";
	//mycut += "Z_pt>40";
	//mycut += "Wmin_pt>40";
	//mycut += "MET>120";
	mycut += "Z_pt>20";
	mycut += "Wmin_pt>20";
	mycut += "MET>140";
	//mycut += "abs(DPhi_METNonWminJet)>2.5";
	mycut += "MET/NonWminJet_pT>0.4 && MET/NonWminJet_pT<0.8";

//	mycut += "jet_n==2";
//	mycut += "mj0j1>60 && mj0j1<100";
//	mycut += "Z_pt>50";
//	mycut += "W01_pt>50";
//	mycut += "MET>50";
//	mycut += "min(DPhi_METJetLeading,DPhi_METJetSecond)>0.4";

	TCut weight_z = "totalWeight*36.1*1000";
	TCut weight_g = "totalWeight*36.1*1000";
	//TCut weight_g = "totalWeight*36.1*10000";
	//TCut weight_g = "totalWeight*36.1*100000";

	tree1->Draw("MTWino_Wmin>>h1(20,0,400)",mycut*weight_z);
	tree2->Draw("MTWino_Wmin>>h2(20,0,400)",mycut*weight_g);
	//tree1->Draw("min(MTWino_W01,MTWino_Wmin)>>h1(20,0,400)",mycut*weight_z);
	//tree2->Draw("min(MTWino_W01,MTWino_Wmin)>>h2(20,0,400)",mycut*weight_g);
	//tree1->Draw("DPhi_METLepMin>>h1(10,0,4)",mycut*weight_z);
	//tree2->Draw("DPhi_METLepMin>>h2(10,0,4)",mycut*weight_g);
	//tree1->Draw("mWmin>>h1(10,0,200)",mycut*weight_z);
	//tree2->Draw("mWmin>>h2(10,0,200)",mycut*weight_g);
	//tree1->Draw("MET/NonWminJet_pT>>h1(10,0,2)",mycut*weight_z);
	//tree2->Draw("MET/NonWminJet_pT>>h2(10,0,2)",mycut*weight_g);


	h1->SetLineColor(4);
	h1->SetStats(0);
	h2->SetLineColor(2);
	h2->SetStats(0);
	h1->Draw();
	h2->Draw("same");
	//c1->SetLogy();

}
