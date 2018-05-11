//-----------------------------------------------------------------------------------------------
// this script takes the photon ntuples from QuickAna and generates smaller ntuples with information required by photon method
// the parameters of the function GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, bool isData) are:
// 	sampleID: DSID of the MC sample
// 	pathToNtuples: the path to the input ntuple
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetPhotonEvents.C+("361039","gmc","root://eosatlas//eos/atlas/user/r/rshang/MET_template/v21/photon_mc/",false)'
//-----------------------------------------------------------------------------------------------


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
#include "TObject.h"

#include "BasicSetting.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void getPhotonSmearingFunction(TString file, TString histname, TH1D* hist) {
	TFile *fData = TFile::Open( file );
	fData->cd("");
	TH1D* temp = (TH1D*)fData->Get(histname);
	for (int bin=1;bin<=temp->GetNbinsX();bin++) {
		hist->SetBinContent(bin,max(temp->GetBinContent(bin),0.));
	}
	fData->Close();
}
void GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, int isData) {

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	//TH1D*   EventCountHist = (TH1D*) inputFile->Get("EventCountHist");
	//Float_t _nGenEvents    = EventCountHist->GetBinContent(2);
	Float_t _nGenEvents    = 1.;
	TTree*  T              = (TTree*)inputFile->Get("tree_NoSys");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << T->GetEntries() << endl;
	cout << "Total generated events : " << _nGenEvents     << endl;

	//-----------------------------
	// access existing branches
	//-----------------------------
	
	
	Long64_t EventNumber;
	Int_t RunNumber;
	int trigMatch_HLT_g35_loose_L1EM15;
	int trigMatch_HLT_g40_loose_L1EM15;
	int trigMatch_HLT_g45_loose_L1EM15;
	int trigMatch_HLT_g50_loose_L1EM15;
	int trigMatch_HLT_g60_loose;
	int trigMatch_HLT_g70_loose;
	int trigMatch_HLT_g80_loose;
	int trigMatch_HLT_g100_loose;
	int trigMatch_HLT_g120_loose;
	int trigMatch_HLT_g140_loose;
	float trigPrescale_HLT_g35_loose_L1EM15;
	float trigPrescale_HLT_g40_loose_L1EM15;
	float trigPrescale_HLT_g45_loose_L1EM15;
	float trigPrescale_HLT_g50_loose_L1EM15;
	float trigPrescale_HLT_g60_loose;
	float trigPrescale_HLT_g70_loose;
	float trigPrescale_HLT_g80_loose;
	float trigPrescale_HLT_g100_loose;
	float trigPrescale_HLT_g120_loose;
	float trigPrescale_HLT_g140_loose;
	Float_t Mu;
	float MET;
	float MET_phi;
	float MET_softTerm;
	float MET_softPhi;
	Float_t truthMET = 0;
	Float_t truthMET_Phi = 0;
	float DPhi_METJetLeading;
	float DPhi_METJetSecond;
	float HT;
	Int_t jet_n;
	Int_t bjet_n;
	Int_t lep_n;
	std::vector<float>* lep_pT = new std::vector<float>(10);
	std::vector<float>* lep_eta = new std::vector<float>(10);
	std::vector<float>* lep_phi = new std::vector<float>(10);
	//std::vector<float>* jet_btag = new std::vector<float>(10);
	std::vector<float>* jet_m = new std::vector<float>(10);
	std::vector<float>* jet_pT = new std::vector<float>(10);
	std::vector<float>* jet_eta = new std::vector<float>(10);
	std::vector<float>* jet_phi = new std::vector<float>(10);
	std::vector<float>* photon_pT = new std::vector<float>(10);
	std::vector<float>* photon_eta = new std::vector<float>(10);
	std::vector<float>* photon_phi = new std::vector<float>(10);
	//std::vector<long>* photon_passAmbi = new std::vector<long>(10);
	float PhotonPt;
	float PhotonEta;
	float PhotonPhi;
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("EventNumber"              ,1); 
    T->SetBranchStatus("RunNumber"              ,1); 
    T->SetBranchStatus("mu"              ,1); 
    T->SetBranchStatus("met_Et"             ,1); 
    T->SetBranchStatus("met_Phi"         ,1); 
    T->SetBranchStatus("TST_Et"         ,1); 
    T->SetBranchStatus("TST_Phi"         ,1); 
    //T->SetBranchStatus("truthMET"             ,1); 
    //T->SetBranchStatus("truthMET_Phi"         ,1); 
    //T->SetBranchStatus("DPhi_METJetLeading"         ,1); 
    //T->SetBranchStatus("DPhi_METJetSecond"          ,1); 
    T->SetBranchStatus("Ht30"              ,1); 
    T->SetBranchStatus("nBJet30_MV2c10_FixedCutBEff_77"          ,1); 
    T->SetBranchStatus("nLep_signal"           ,1); 
    T->SetBranchStatus("lepPt"          ,1); 
    T->SetBranchStatus("lepEta"         ,1); 
    T->SetBranchStatus("lepPhi"         ,1); 
    T->SetBranchStatus("nJet30"           ,1); 
    //T->SetBranchStatus("jet_btag"          ,1); 
    T->SetBranchStatus("jetM"          ,1); 
    T->SetBranchStatus("jetPt"          ,1); 
    T->SetBranchStatus("jetEta"         ,1); 
    T->SetBranchStatus("jetPhi"         ,1); 
    T->SetBranchStatus("trigMatch_HLT_g35_loose_L1EM15", 1); 
    T->SetBranchStatus("trigMatch_HLT_g40_loose_L1EM15", 1); 
    T->SetBranchStatus("trigMatch_HLT_g45_loose_L1EM15", 1); 
    T->SetBranchStatus("trigMatch_HLT_g50_loose_L1EM15", 1); 
    T->SetBranchStatus("trigMatch_HLT_g60_loose", 1); 
    T->SetBranchStatus("trigMatch_HLT_g70_loose", 1); 
    T->SetBranchStatus("trigMatch_HLT_g80_loose", 1); 
    T->SetBranchStatus("trigMatch_HLT_g100_loose", 1); 
    T->SetBranchStatus("trigMatch_HLT_g120_loose", 1); 
    T->SetBranchStatus("trigMatch_HLT_g140_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g35_loose_L1EM15", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g40_loose_L1EM15", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g45_loose_L1EM15", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g50_loose_L1EM15", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g60_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g70_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g80_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g100_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g120_loose", 1); 
    T->SetBranchStatus("trigPrescale_HLT_g140_loose", 1); 
    T->SetBranchStatus("PhotonPt", 1); 
    T->SetBranchStatus("PhotonEta", 1); 
    T->SetBranchStatus("PhotonPhi", 1); 
    //T->SetBranchStatus("photon_passAmbi", 1); 

	T->SetBranchAddress("EventNumber"              ,&EventNumber               );
	T->SetBranchAddress("RunNumber"              ,&RunNumber               );
	T->SetBranchAddress("mu"              ,&Mu               );
	T->SetBranchAddress("met_Et"             ,&MET              );
	T->SetBranchAddress("met_Phi"         ,&MET_phi          );
	T->SetBranchAddress("TST_Et"         ,&MET_softTerm          );
	T->SetBranchAddress("TST_Phi"         ,&MET_softPhi          );
	//T->SetBranchAddress("truthMET"             ,&truthMET              );
	//T->SetBranchAddress("truthMET_Phi"         ,&truthMET_Phi          );
	//T->SetBranchAddress("DPhi_METJetLeading"         ,&DPhi_METJetLeading          );
	//T->SetBranchAddress("DPhi_METJetSecond"          ,&DPhi_METJetSecond           );
	T->SetBranchAddress("Ht30"              ,&HT               );
	T->SetBranchAddress("nBJet30_MV2c10_FixedCutBEff_77"          ,&bjet_n            );
	T->SetBranchAddress("nLep_signal"           ,&lep_n            );
	T->SetBranchAddress("lepPt"          ,&lep_pT           );
	T->SetBranchAddress("lepEta"         ,&lep_eta          );
	T->SetBranchAddress("lepPhi"         ,&lep_phi          );
	T->SetBranchAddress("nJet30"           ,&jet_n            );
	//T->SetBranchAddress("jet_btag"          ,&jet_btag           );
	T->SetBranchAddress("jetM"          ,&jet_m           );
	T->SetBranchAddress("jetPt"          ,&jet_pT           );
	T->SetBranchAddress("jetEta"         ,&jet_eta          );
	T->SetBranchAddress("jetPhi"         ,&jet_phi          );
	T->SetBranchAddress("trigMatch_HLT_g35_loose_L1EM15", 		&trigMatch_HLT_g35_loose_L1EM15);
	T->SetBranchAddress("trigMatch_HLT_g40_loose_L1EM15", 		&trigMatch_HLT_g40_loose_L1EM15);
	T->SetBranchAddress("trigMatch_HLT_g45_loose_L1EM15", 		&trigMatch_HLT_g45_loose_L1EM15);
	T->SetBranchAddress("trigMatch_HLT_g50_loose_L1EM15", 		&trigMatch_HLT_g50_loose_L1EM15);
	T->SetBranchAddress("trigMatch_HLT_g60_loose", 			&trigMatch_HLT_g60_loose);
	T->SetBranchAddress("trigMatch_HLT_g70_loose", 			&trigMatch_HLT_g70_loose);
	T->SetBranchAddress("trigMatch_HLT_g80_loose", 			&trigMatch_HLT_g80_loose);
	T->SetBranchAddress("trigMatch_HLT_g100_loose", 		&trigMatch_HLT_g100_loose);
	T->SetBranchAddress("trigMatch_HLT_g120_loose", 		&trigMatch_HLT_g120_loose);
	T->SetBranchAddress("trigMatch_HLT_g140_loose", 		&trigMatch_HLT_g140_loose);
	T->SetBranchAddress("trigPrescale_HLT_g35_loose_L1EM15", 	&trigPrescale_HLT_g35_loose_L1EM15);
	T->SetBranchAddress("trigPrescale_HLT_g40_loose_L1EM15", 	&trigPrescale_HLT_g40_loose_L1EM15);
	T->SetBranchAddress("trigPrescale_HLT_g45_loose_L1EM15", 	&trigPrescale_HLT_g45_loose_L1EM15);
	T->SetBranchAddress("trigPrescale_HLT_g50_loose_L1EM15", 	&trigPrescale_HLT_g50_loose_L1EM15);
	T->SetBranchAddress("trigPrescale_HLT_g60_loose", 		&trigPrescale_HLT_g60_loose);
	T->SetBranchAddress("trigPrescale_HLT_g70_loose", 		&trigPrescale_HLT_g70_loose);
	T->SetBranchAddress("trigPrescale_HLT_g80_loose", 		&trigPrescale_HLT_g80_loose);
	T->SetBranchAddress("trigPrescale_HLT_g100_loose", 		&trigPrescale_HLT_g100_loose);
	T->SetBranchAddress("trigPrescale_HLT_g120_loose", 		&trigPrescale_HLT_g120_loose);
	T->SetBranchAddress("trigPrescale_HLT_g140_loose", 		&trigPrescale_HLT_g140_loose);
	T->SetBranchAddress("PhotonPt", &PhotonPt);
	T->SetBranchAddress("PhotonEta", &PhotonEta);
	T->SetBranchAddress("PhotonPhi", &PhotonPhi);
	//T->SetBranchAddress("photon_passAmbi", &photon_passAmbi);

	float sampleWeight;
	Float_t eventWeight;
	//std::vector<float>* truthPhoton_pT = new std::vector<float>(10);
	//std::vector<float>* truthPhoton_eta = new std::vector<float>(10);
	//std::vector<float>* truthPhoton_phi = new std::vector<float>(10);
	//std::vector<float>* photon_truthPt = new std::vector<float>(10);
	//std::vector<float>* photon_truthEta = new std::vector<float>(10);
	//std::vector<float>* photon_truthPhi = new std::vector<float>(10);
	if (isData!=1) {
	        T->SetBranchStatus("sampleWeight",1);
	        T->SetBranchStatus("eventWeight",1);
	        //T->SetBranchStatus("truthPhoton_pT", 1);
	        //T->SetBranchStatus("truthPhoton_eta", 1);
	        //T->SetBranchStatus("truthPhoton_phi", 1);
	        //T->SetBranchStatus("photon_truthPt", 1);
	        //T->SetBranchStatus("photon_truthEta", 1);
	        //T->SetBranchStatus("photon_truthPhi", 1);

		T->SetBranchAddress("sampleWeight",&sampleWeight);
		T->SetBranchAddress("eventWeight",&eventWeight);
		//T->SetBranchAddress("truthPhoton_pT", &truthPhoton_pT);
		//T->SetBranchAddress("truthPhoton_eta", &truthPhoton_eta);
		//T->SetBranchAddress("truthPhoton_phi", &truthPhoton_phi);
		//T->SetBranchAddress("photon_truthPt", &photon_truthPt);
		//T->SetBranchAddress("photon_truthEta", &photon_truthEta);
		//T->SetBranchAddress("photon_truthPhi", &photon_truthPhi);
	}

	//-----------------------------
	// add new branches
	//-----------------------------

	TFile   outputFile(TString(outputPath)+"/"+TString(outputName)+"/"+TString(sampleID.c_str())+".root","recreate");
	TTree BaselineTree("BaselineTree","baseline tree");
	float gamma_pt = 0.;
	float gamma_eta = 0.;
	float gamma_phi = 0.;
	float gamma_dR = 999.;
	//int gamma_passAmbi = 0;
	//float truthGamma_pt = 0.;
	//float truthGamma_eta = 0.;
	//float truthGamma_phi = 0.;
	float METl = 0.;
	float METt = 0.;
	//float truthMETl = 0.;
	//float truthMETt = 0.;
	float DPhi_METPhoton = 0.;
	float MinDR_PhotonJet = 0.;
	float MinDPhi_PhotonJet = 0.;
	int pt = 0;
	int ht = 0;
	int njet = 0;
	int nbjet = 0;
	float MT;
	BaselineTree.Branch("pt",&pt,"pt/I");
	BaselineTree.Branch("ht",&ht,"ht/I");
	BaselineTree.Branch("njet",&njet,"njet/I");
	BaselineTree.Branch("nbjet",&nbjet,"nbjet/I");
	BaselineTree.Branch("Mu",&Mu,"Mu/Float_t");
	BaselineTree.Branch("MET_raw",&MET,"MET_raw/F");
	BaselineTree.Branch("METl_raw",&METl,"METl_raw/F");
	BaselineTree.Branch("METt_raw",&METt,"METt_raw/F");
	BaselineTree.Branch("MET_phi_raw",&MET_phi,"MET_phi_raw/F");
	//BaselineTree.Branch("MET_softTerm",&MET_softTerm,"MET_softTerm/F");
	//BaselineTree.Branch("MET_softPhi",&MET_softPhi,"MET_softPhi/F");
	//BaselineTree.Branch("truthMET",&truthMET,"truthMET/Float_t");
	//BaselineTree.Branch("truthMETl",&truthMETl,"truthMETl/F");
	//BaselineTree.Branch("truthMETt",&truthMETt,"truthMETt/F");
	//BaselineTree.Branch("truthMET_Phi",&MET_phi,"MET_phi/Float_t");
	//BaselineTree.Branch("DPhi_METJetLeading_raw",&DPhi_METJetLeading,"DPhi_METJetLeading_raw/F");
	//BaselineTree.Branch("DPhi_METJetSecond_raw",&DPhi_METJetSecond,"DPhi_METJetSecond_raw/F");
	//BaselineTree.Branch("DPhi_METPhoton_raw",&DPhi_METPhoton,"DPhi_METPhoton_raw/F");
	//BaselineTree.Branch("MinDR_PhotonJet",&MinDR_PhotonJet,"MinDR_PhotonJet/F");
	//BaselineTree.Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/F");
	BaselineTree.Branch("HT",&HT,"HT/F");
	BaselineTree.Branch("MT",&MT,"MT/F");
	BaselineTree.Branch("gamma_pt",&gamma_pt,"gamma_pt/F");
	BaselineTree.Branch("gamma_eta",&gamma_eta,"gamma_eta/F");
	BaselineTree.Branch("gamma_phi",&gamma_phi,"gamma_phi/F");
	//BaselineTree.Branch("gamma_passAmbi",&gamma_passAmbi,"gamma_passAmbi/I");
	BaselineTree.Branch("bjet_n",&bjet_n,"bjet_n/Int_t");
	BaselineTree.Branch("jet_n",&jet_n,"jet_n/Int_t");
	//BaselineTree.Branch("jet_btag","std::vector<float>",&jet_btag);
	BaselineTree.Branch("jet_m","std::vector<float>",&jet_m);
	BaselineTree.Branch("jet_pT","std::vector<float>",&jet_pT);
	BaselineTree.Branch("jet_phi","std::vector<float>",&jet_phi);
	BaselineTree.Branch("jet_eta","std::vector<float>",&jet_eta);
	BaselineTree.Branch("EventNumber",&EventNumber,"EventNumber/Long64_t");
	BaselineTree.Branch("RunNumber",&RunNumber,"RunNumber/Int_t");
	BaselineTree.Branch("lep_n_raw",&lep_n,"lep_n_raw/Int_t");
	BaselineTree.Branch("lep_pT_raw","std::vector<float>",&lep_pT);
	BaselineTree.Branch("lep_phi_raw","std::vector<float>",&lep_phi);
	BaselineTree.Branch("lep_eta_raw","std::vector<float>",&lep_eta);
	float totalWeight = 0.;
	BaselineTree.Branch("totalWeight",&totalWeight,"totalWeight/F");
	if (isData!=1) {
		//BaselineTree.Branch("truthGamma_pt",&truthGamma_pt,"truthGamma_pt/F");
		//BaselineTree.Branch("truthGamma_eta",&truthGamma_eta,"truthGamma_eta/F");
		//BaselineTree.Branch("truthGamma_phi",&truthGamma_phi,"truthGamma_phi/F");
	//	BaselineTree.Branch("gamma_dR",&gamma_dR,"gamma_dR/F");
	}

	TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
	hist_low_njet->SetStats(0);
	TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
	hist_low_nbjet->SetStats(0);
	TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,pt_bin);
	hist_low_pt->SetStats(0);
	TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
	hist_sm_pt->SetStats(0);
	TH1D* hist_low_et = new TH1D("hist_low_et","",bin_size,pt_bin);
	hist_low_et->SetStats(0);
	TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);
	hist_low_ht->SetStats(0);
	TH1D* hist_medium_njet = new TH1D("hist_medium_njet","",bin_size,njet_bin);
	hist_medium_njet->SetStats(0);
	TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet","",bin_size,njet_bin);
	hist_medium_nbjet->SetStats(0);
	TH1D* hist_medium_pt = new TH1D("hist_medium_pt","",bin_size,pt_bin);
	hist_medium_pt->SetStats(0);
	TH1D* hist_medium_et = new TH1D("hist_medium_et","",bin_size,pt_bin);
	hist_medium_et->SetStats(0);
	TH1D* hist_medium_ht = new TH1D("hist_medium_ht","",bin_size,ht_bin);
	hist_medium_ht->SetStats(0);
	TH1D* hist_high_njet = new TH1D("hist_high_njet","",bin_size,njet_bin);
	hist_high_njet->SetStats(0);
	TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet","",bin_size,njet_bin);
	hist_high_nbjet->SetStats(0);
	TH1D* hist_high_pt = new TH1D("hist_high_pt","",bin_size,pt_bin);
	hist_high_pt->SetStats(0);
	TH1D* hist_high_et = new TH1D("hist_high_et","",bin_size,pt_bin);
	hist_high_et->SetStats(0);
	TH1D* hist_high_ht = new TH1D("hist_high_ht","",bin_size,ht_bin);
	hist_high_ht->SetStats(0);
	TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet","",bin_size,njet_bin);
	hist_zmet_njet->SetStats(0);
	TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet","",bin_size,njet_bin);
	hist_zmet_nbjet->SetStats(0);
	TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt","",bin_size,pt_bin);
	hist_zmet_pt->SetStats(0);
	TH1D* hist_zmet_et = new TH1D("hist_zmet_et","",bin_size,pt_bin);
	hist_zmet_et->SetStats(0);
	TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht","",bin_size,ht_bin);
	hist_zmet_ht->SetStats(0);
	TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet","",bin_size,njet_bin);
	hist_bveto_njet->SetStats(0);
	TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet","",bin_size,njet_bin);
	hist_bveto_nbjet->SetStats(0);
	TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt","",bin_size,pt_bin);
	hist_bveto_pt->SetStats(0);
	TH1D* hist_bveto_et = new TH1D("hist_bveto_et","",bin_size,pt_bin);
	hist_bveto_et->SetStats(0);
	TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht","",bin_size,ht_bin);
	hist_bveto_ht->SetStats(0);


	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector z_4vec;
	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = T->GetEntries();
	for (Long64_t i=0;i<nentries;i+=event_interval) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		T->GetEntry(i);

		photon_pT->clear();
		photon_eta->clear();
		photon_phi->clear();
		photon_pT->push_back(PhotonPt);
		photon_eta->push_back(PhotonEta);
		photon_phi->push_back(PhotonPhi);

		if (lep_n>0) continue;
		//std::cout << i << " events processed nlep cut." << std::endl;
		if (jet_n==0) continue;
		//std::cout << i << " events processed njet cut." << std::endl;
		if (photon_pT->size()==0) continue;
		//std::cout << i << " events processed size cut." << std::endl;
		if (photon_pT->at(0)<37.) continue;
		//if (photon_pT->at(0)<50.) continue;
		//if (abs(photon_eta->at(0))>1.5 && abs(photon_eta->at(0))<1.6) continue;

		//std::cout << i << " events processed pt cut." << std::endl;
		// find the trigger and prescale
		float trigWeight = 0;

//		if (!is2016) {
//			if (HLT_g15_loose_L1EM7->at(1)==1 && photon_pT->at(0)>(15+2) && photon_pT->at(0)<(20+2)) trigWeight = HLT_g15_loose_L1EM7->at(0);
//			if (HLT_g20_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(20+2) && photon_pT->at(0)<(35+2)) trigWeight = HLT_g20_loose_L1EM15->at(0);
//			if (HLT_g35_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(35+2) && photon_pT->at(0)<(40+5)) trigWeight = HLT_g35_loose_L1EM15->at(0);
//			if (HLT_g40_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(40+5) && photon_pT->at(0)<(45+5)) trigWeight = HLT_g40_loose_L1EM15->at(0);
//			if (HLT_g45_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(45+5) && photon_pT->at(0)<(50+5)) trigWeight = HLT_g45_loose_L1EM15->at(0);
//			if (HLT_g50_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(50+5) && photon_pT->at(0)<(120+5)) trigWeight = HLT_g50_loose_L1EM15->at(0);
//			if (HLT_g120_loose->at(1)==1 && photon_pT->at(0)>(120+5)) trigWeight = 1.;
//		}
		if (trigMatch_HLT_g35_loose_L1EM15==1 && photon_pT->at(0)>(35+5) && photon_pT->at(0)<(40+5)) trigWeight = trigPrescale_HLT_g35_loose_L1EM15;
		if (trigMatch_HLT_g40_loose_L1EM15==1 && photon_pT->at(0)>(40+5) && photon_pT->at(0)<(45+5)) trigWeight = trigPrescale_HLT_g40_loose_L1EM15;
		if (trigMatch_HLT_g45_loose_L1EM15==1 && photon_pT->at(0)>(45+5) && photon_pT->at(0)<(50+5)) trigWeight = trigPrescale_HLT_g45_loose_L1EM15;
		if (trigMatch_HLT_g50_loose_L1EM15==1 && photon_pT->at(0)>(50+5) && photon_pT->at(0)<(60+5)) trigWeight = trigPrescale_HLT_g50_loose_L1EM15;
		if (trigMatch_HLT_g60_loose==1 && photon_pT->at(0)>(60+5) && photon_pT->at(0)<(70+5)) trigWeight = trigPrescale_HLT_g60_loose;
		if (trigMatch_HLT_g70_loose==1 && photon_pT->at(0)>(70+5) && photon_pT->at(0)<(80+5)) trigWeight = trigPrescale_HLT_g70_loose;
		if (trigMatch_HLT_g80_loose==1 && photon_pT->at(0)>(80+5) && photon_pT->at(0)<(100+5)) trigWeight = trigPrescale_HLT_g80_loose;
		if (trigMatch_HLT_g100_loose==1 && photon_pT->at(0)>(100+5) && photon_pT->at(0)<(140+5)) trigWeight = trigPrescale_HLT_g100_loose;
		if (trigMatch_HLT_g140_loose==1 && photon_pT->at(0)>(140+5)) trigWeight = trigPrescale_HLT_g140_loose;
		if (trigWeight==0) continue;
		//std::cout << i << " events processed trigger cut." << std::endl;

		gamma_pt = photon_pT->at(0);
		gamma_eta = photon_eta->at(0);
		gamma_phi = photon_phi->at(0);
		//gamma_passAmbi = photon_passAmbi->at(0);

		njet = hist_low_njet->FindBin(jet_n)-1;
		if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
		nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
		if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
		pt = hist_low_pt->FindBin(gamma_pt)-1;
		if (gamma_pt>pt_bin[bin_size]) pt = bin_size-1;
		ht = hist_low_ht->FindBin(HT)-1;
		if (HT>ht_bin[bin_size]) ht = bin_size-1;

		totalWeight = trigWeight;

		// here we compute the MET parallel and perpendicular components
		METt = MET*TMath::Sin(MET_phi-gamma_phi);
		METl = MET*TMath::Cos(MET_phi-gamma_phi);

		TLorentzVector gamma_4vec;
		gamma_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,0);
		// hist_dPt_Pt histogram, i.e. the photon truth-reco response function
		if (isData!=1) {

			totalWeight = (sampleWeight*eventWeight)/_nGenEvents;
			if (isData==2) totalWeight = totalWeight*-1.; // V+gamma subtraction
			//float Vgamma_SF = 1.5;  // v02-03
			float Vgamma_SF = 1.0;  // v02-04
			if (outputName=="Vg") totalWeight = totalWeight*lumi*Vgamma_SF;

			//truthGamma_pt =  photon_truthPt->at(0);
			//truthGamma_eta = photon_truthEta->at(0);
			//truthGamma_phi = photon_truthPhi->at(0);

			//TLorentzVector truthGamma_4vec;
			//truthGamma_4vec.SetPtEtaPhiM(truthGamma_pt,truthGamma_eta,truthGamma_phi,0);
			//gamma_dR = gamma_4vec.DeltaR(truthGamma_4vec);

			//int pt_truth = hist_low_pt->FindBin(truthGamma_pt)-1;
			//
			//truthMETt = truthMET*TMath::Sin(truthMET_Phi-truthGamma_phi);
			//truthMETl = truthMET*TMath::Cos(truthMET_Phi-truthGamma_phi);

		}
		totalWeight = totalWeight*event_interval;

		DPhi_METPhoton = fabs(TMath::ATan2(METt,METl));
		MinDR_PhotonJet = 1000.;
		MinDPhi_PhotonJet = 1000.;
		TLorentzVector jet_4vec;
		for (int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			float DR_PhotonJet = jet_4vec.DeltaR(gamma_4vec);
			float DPhi_PhotonJet = jet_4vec.DeltaPhi(gamma_4vec);
			if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
			if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
		}

		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector lep_4vec;
		if (lep_pT->size()>0) lep_4vec.SetPtEtaPhiM(lep_pT->at(0),0,lep_phi->at(0),0);  // only transverse component
		else lep_4vec.SetPtEtaPhiM(0,0,0,0);
		if (lep_pT->size()>0) MT = (met_4vec+lep_4vec).M();
		else MT = 0;



		float resample_size = 1.;
		if (MET>200.) {
			if (totalWeight>1 && isData==1) {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else if (gamma_pt>200. && outputName=="gmc") {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else if (outputName=="Vg") {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else BaselineTree.Fill();
		}
		else BaselineTree.Fill();

	}

	BaselineTree.Write();

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;


}
