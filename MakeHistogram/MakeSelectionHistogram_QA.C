
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip> 

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TSpectrum.h"
#include "TVirtualFFT.h"

TTree* inputTree;
Long64_t EventNumber;
Int_t RunNumber;
double totalWeight = 0.;
double mll;
double MET;
double MET_phi;
double METl = 0.;
double METt = 0.;
double lep0_truthPt= 0.;
double lep0_truthEta = 0.;
double lep0_truthPhi = 0.;
double lep1_truthPt= 0.;
double lep1_truthEta = 0.;
double lep1_truthPhi = 0.;
Int_t jet_n;
Int_t bjet_n;
double HT = 0.;
double truthGamma_pt = 0.;
double truthGamma_phi = 0.;
double gamma_pt = 0.;
double gamma_eta = 0.;
double gamma_phi = 0.;
double Z_pt = 0.;
double Z_eta = 0.;
double Z_phi = 0.;
std::vector<double>* jet_pT = new std::vector<double>(10);
std::vector<double>* jet_eta = new std::vector<double>(10);
std::vector<double>* jet_phi = new std::vector<double>(10);
std::vector<double>* jet_m = new std::vector<double>(10);
std::vector<double>* jet_btag = new std::vector<double>(10);
std::vector<double>* lep_pT = new std::vector<double>(10);
std::vector<double>* lep_eta = new std::vector<double>(10);
std::vector<double>* lep_phi = new std::vector<double>(10);
double ptrw = 0;


void SetBasicBranch() {
	inputTree->SetBranchStatus("EventNumber"		,1);
	inputTree->SetBranchAddress("EventNumber"		,&EventNumber);
	inputTree->SetBranchStatus("RunNumber"			,1);
	inputTree->SetBranchAddress("RunNumber"			,&RunNumber);
	inputTree->SetBranchStatus("totalWeight"		,1);
	inputTree->SetBranchAddress("totalWeight"		,&totalWeight);
	inputTree->SetBranchStatus("MET"			,1);
	inputTree->SetBranchAddress("MET"			,&MET);
	inputTree->SetBranchStatus("MET_phi"			,1);
	inputTree->SetBranchAddress("MET_phi"			,&MET_phi);
	inputTree->SetBranchStatus("METl"			,1);
	inputTree->SetBranchAddress("METl"			,&METl);
	inputTree->SetBranchStatus("METt"			,1);
	inputTree->SetBranchAddress("METt"			,&METt);
	inputTree->SetBranchStatus("mll"			,1);
	inputTree->SetBranchAddress("mll"			,&mll);
	inputTree->SetBranchStatus("HT"				,1);
	inputTree->SetBranchAddress("HT"			,&HT);
	inputTree->SetBranchStatus("jet_n"			,1);
	inputTree->SetBranchAddress("jet_n"			,&jet_n);
	inputTree->SetBranchStatus("bjet_n"			,1);
	inputTree->SetBranchAddress("bjet_n"			,&bjet_n);
	inputTree->SetBranchStatus("gamma_pt"			,1);
	inputTree->SetBranchAddress("gamma_pt"			,&gamma_pt);
	inputTree->SetBranchStatus("gamma_eta"			,1);
	inputTree->SetBranchAddress("gamma_eta"			,&gamma_eta);
	inputTree->SetBranchStatus("gamma_phi"			,1);
	inputTree->SetBranchAddress("gamma_phi"			,&gamma_phi);
	inputTree->SetBranchStatus("Z_pt"			,1);
	inputTree->SetBranchAddress("Z_pt"			,&Z_pt);
	inputTree->SetBranchStatus("Z_eta"			,1);
	inputTree->SetBranchAddress("Z_eta"			,&Z_eta);
	inputTree->SetBranchStatus("Z_phi"			,1);
	inputTree->SetBranchAddress("Z_phi"			,&Z_phi);
	inputTree->SetBranchStatus("jet_pT"			,1);
	inputTree->SetBranchAddress("jet_pT"			,&jet_pT);
	inputTree->SetBranchStatus("jet_eta"			,1);
	inputTree->SetBranchAddress("jet_eta"			,&jet_eta);
	inputTree->SetBranchStatus("jet_phi"			,1);
	inputTree->SetBranchAddress("jet_phi"			,&jet_phi);
	inputTree->SetBranchStatus("jet_m"			,1);
	inputTree->SetBranchAddress("jet_m"			,&jet_m);
	inputTree->SetBranchStatus("lep_pT"			,1);
	inputTree->SetBranchAddress("lep_pT"			,&lep_pT);
	inputTree->SetBranchStatus("lep_eta"			,1);
	inputTree->SetBranchAddress("lep_eta"			,&lep_eta);
	inputTree->SetBranchStatus("lep_phi"			,1);
	inputTree->SetBranchAddress("lep_phi"			,&lep_phi);
	//if (isData!=1) T->SetBranchAddress("truthGamma_pt", &truthGamma_pt);
	//if (isData!=1) T->SetBranchAddress("truthGamma_phi", &truthGamma_phi);
}
bool Strong2LSelection(double unit) {
	if (jet_n<2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (mll<12) return false;
	if (Z_pt<50) return false;
	return true;
}
void Strong2L(TString pathToNtuple, TString pathToOutput, TString TreeName, bool DoReweighting, double scale, double unit) {

	TH1::SetDefaultSumw2();

	std::cout << "open file " << pathToNtuple << std::endl;
	TFile*  inputFile      = TFile::Open(pathToNtuple);
	inputTree              = (TTree*)inputFile->Get(TreeName);
	inputTree->SetBranchStatus("*", 0);
	SetBasicBranch();
	ptrw = 1.;
	if (DoReweighting) {
		inputTree->SetBranchStatus("ptrw_bveto"			,1);
		inputTree->SetBranchAddress("ptrw_bveto"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht200"			,1);
		//inputTree->SetBranchAddress("ptrw_ht200"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht400"			,1);
		//inputTree->SetBranchAddress("ptrw_ht400"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht1200"		,1);
		//inputTree->SetBranchAddress("ptrw_ht1200"		,&ptrw);
	}

	Long64_t nentries = inputTree->GetEntries();

	std::cout << "output file " << pathToOutput << std::endl;
	TFile   outputFile(pathToOutput,"recreate");

	TH1D Hist_MET = TH1D(TString("Hist_MET"),"",15,0,300);
	TH1D Hist_RatioMETJet1st = TH1D(TString("Hist_RatioMETJet1st"),"",16,0,40);
	TH1D Hist_RatioMETJet2nd = TH1D(TString("Hist_RatioMETJet2nd"),"",16,0,40);
	TH1D Hist_RatioMETPhoton = TH1D(TString("Hist_RatioMETPhoton"),"",16,0,40);

	int step = 1;
	for (Long64_t i=0;i<nentries;i+=step) {
		if (fmod(i,1e5)==0) std::cout << pathToNtuple << " " << i << "/" << nentries << " events processed." << std::endl;
		inputTree->GetEntry(i);
		if (!Strong2LSelection(unit)) continue;
		if (MET/unit<50.) continue;

		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector jet1_4vec;
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),0);
		TLorentzVector jet2_4vec;
		jet2_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),0);
		TLorentzVector lept1_4vec;
		lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
		TLorentzVector lept2_4vec;
		lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
		TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
		double DPhiMETJet1st = abs(jet1_4vec.DeltaPhi(met_4vec));
		double RatioMETJet1st = (MET/unit)/pow(abs((jet_pT->at(0)/unit)*TMath::Cos(DPhiMETJet1st)),0.5);
		double DPhiMETJet2nd = abs(jet2_4vec.DeltaPhi(met_4vec));
		double RatioMETJet2nd = (MET/unit)/pow(abs((jet_pT->at(1)/unit)*TMath::Cos(DPhiMETJet2nd)),0.5);
		double DPhiMETPhoton = abs(z_4vec.DeltaPhi(met_4vec));
		double RatioMETPhoton = (MET/unit)/pow(abs((z_4vec.Pt()/unit)*TMath::Cos(DPhiMETPhoton)),0.5);

		
		if (RatioMETJet1st<10) continue;
		//if (RatioMETJet2nd<10) continue;

		Hist_MET.Fill(MET/unit,totalWeight*ptrw*scale);
		Hist_RatioMETJet1st.Fill(RatioMETJet1st,totalWeight*ptrw*scale);
		Hist_RatioMETJet2nd.Fill(RatioMETJet2nd,totalWeight*ptrw*scale);
		Hist_RatioMETPhoton.Fill(RatioMETPhoton,totalWeight*ptrw*scale);
	}

	Hist_MET.Write();
	Hist_RatioMETJet1st.Write();
	Hist_RatioMETJet2nd.Write();
	Hist_RatioMETPhoton.Write();
	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}
void RunStrong2L() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;

	//pathToNtuple = TString("../OutputNtuples/tt/ttee.root");
	//pathToOutput = TString("../OutputHistogram/Top_ee_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	//pathToNtuple = TString("../OutputNtuples/vv/vvee.root");
	//pathToOutput = TString("../OutputHistogram/VV_ee_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	//pathToNtuple = TString("../OutputNtuples/data_2016/data_2016_ee.root");
	//pathToOutput = TString("../OutputHistogram/Data_ee_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);
	//pathToNtuple = TString("../OutputNtuples/gdata_2016/data_2016_ee_mcmetl.root");
	//pathToOutput = TString("../OutputHistogram/GData_ee_Strong2L.root");
	//DoReweighting = true;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);
	//pathToNtuple = TString("../OutputNtuples/tt/ttmm.root");
	//pathToOutput = TString("../OutputHistogram/Top_mm_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	//pathToNtuple = TString("../OutputNtuples/vv/vvmm.root");
	//pathToOutput = TString("../OutputHistogram/VV_mm_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	//pathToNtuple = TString("../OutputNtuples/data_2016/data_2016_mm.root");
	//pathToOutput = TString("../OutputHistogram/Data_mm_Strong2L.root");
	//DoReweighting = false;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);
	//pathToNtuple = TString("../OutputNtuples/gdata_2016/data_2016_mm_mcmetl.root");
	//pathToOutput = TString("../OutputHistogram/GData_mm_Strong2L.root");
	//DoReweighting = true;
	//Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = TString("../OutputNtuples/gmc/gmc_ee_McSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_ee_Strong2L.root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/zjets/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC_ee_Strong2L.root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/gmc/gmc_mm_McSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_mm_Strong2L.root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/zjets/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC_mm_Strong2L.root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	

}
void MakeSelectionHistogram_QA() {

	RunStrong2L();

}
