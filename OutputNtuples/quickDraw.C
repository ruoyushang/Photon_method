
void quickDraw() {

	TChain * tree1 = new TChain("BaselineTree");
	tree1->Add("zjets/zjets_mm.root");
	TChain * tree2 = new TChain("BaselineTree");
	tree2->Add("gmc/gmc_mm_McSmear.root");

	TCut mycut = "";
	mycut += "bjet_n==0";

	mycut += "lep_pT[0]>20 && lep_pT[1]>20";
	mycut += "abs(lep_eta[0])<2.5 && abs(lep_eta[1])<2.5";

	mycut += "jet_n>=2";
	mycut += "Z_pt>50";

	TCut weight_z = "totalWeight*36.1*1000";
	TCut weight_g = "totalWeight*36.1*1000*ptrw_bveto";

	tree1->Draw("MET>>h1(20,0,400)",mycut*weight_z);
	tree2->Draw("MET>>h2(20,0,400)",mycut*weight_g);
	//tree1->Draw("Z_pt>>h1(20,0,500)",mycut*weight_z);
	//tree2->Draw("gamma_pt>>h2(20,0,500)",mycut*weight_g);


	h1->SetLineColor(4);
	h1->SetStats(0);
	h2->SetLineColor(2);
	h2->SetStats(0);
	h1->Draw();
	h2->Draw("same");
	c1->SetLogy();

}
