//-----------------------------------------------------------------------------------------------
// this script takes the outputs from GetBaseLineEvents.C and GetPhotonEvents.C and GetPhotonSmearing.C, and makes photon reweighting factors.
// the parameters of the function GetPhotonReweighting(string ch, bool isData) are:
// 	ch: which dilepton channel you are smearing the photon events to (ee,mm)
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetPhotonReweighting.C+("ee",true)'
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
#include "TGraphAsymmErrors.h"

#include "BasicSetting.C"
#include "GetReweightingHistogram.C"
#include "GetReweightingHistogram_SF.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void GetPhotonReweighting(string ch, int isData, string year) {

	if (smearing_method == 0) photon_tag = "_NoSmear";
	if (smearing_method == 1) photon_tag = "_McSmear";
	if (smearing_method == 2) photon_tag = "_DataSmear";
	if (smearing_method == 3) photon_tag = "_TruthSmear";

	// retrieve the histograms for pT, Njets, HT-reweighting and on-Z rescaling
	std::cout << "Prepare reweighting histograms..." << std::endl;
	//if (isData) GetReweightingHistogram_SF(ch, isData,lumi, photon_tag);
	//else GetReweightingHistogram(ch, isData,lumi, photon_tag);
	GetReweightingHistogram(ch, isData,lumi, photon_tag);

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename;
	if (isData==1) filename = TString(TString(outputPath)+"gdata/gdata_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==0) filename = TString(TString(outputPath)+"gmc/gmc_"+TString(ch)+TString(photon_tag)+".root"); 
	TFile*  f              = new TFile(filename.c_str(),"update");          
	TTree*  outputTree              = (TTree*)f->Get("BaselineTree");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << outputTree->GetEntries() << endl;

	//if( _nGenEvents < 1.0 ) cout << "ERROR 0 EVENTS -----------------------------------------------------" << endl;

	//-----------------------------
	// access existing branches
	//-----------------------------

	double totalWeight = 0.;
	int pt = 0;
	int ptsm = 0;
	int metsm = 0;
	int ht = 0;
	int njet = 0;
	int nbjet = 0;
	Int_t jet_n = 0;
	double mll = 0.;
	double HT = 0.;
	double gamma_pt = 0.;
	double gamma_pt_smear = 0.;
	double MET_smear = 0.;
	double gamma_dR = 0.;
	std::vector<double>* lep_pT = new std::vector<double>(10);
	outputTree->SetBranchAddress("totalWeight"     ,&totalWeight              );
	outputTree->SetBranchAddress("pt"              ,&pt               );
	outputTree->SetBranchAddress("pt_smear"   ,&ptsm               );
	outputTree->SetBranchAddress("met_smear"   ,&metsm               );
	outputTree->SetBranchAddress("ht"              ,&ht               );
	outputTree->SetBranchAddress("njet"            ,&njet             );
	outputTree->SetBranchAddress("nbjet"            ,&nbjet             );
	outputTree->SetBranchAddress("mll"              ,&mll               );
	outputTree->SetBranchAddress("HT"              ,&HT               );
	outputTree->SetBranchAddress("jet_n"           ,&jet_n            );
	outputTree->SetBranchAddress("gamma_pt", &gamma_pt);
	outputTree->SetBranchAddress("Z_pt", &gamma_pt_smear);
	outputTree->SetBranchAddress("MET", &MET_smear);
	outputTree->SetBranchAddress("lep_pT"          ,&lep_pT           );
	if (isData!=1) outputTree->SetBranchAddress("gamma_dR", &gamma_dR);

	//-----------------------------
	// add new branches
	//-----------------------------

	double pt37_bveto_corr_met = 0;
	double pt37_bveto_corr_lpt = 0;
	double ptrw_bveto = 0;
	double ptsmrw_bveto = 0;
	double etrw_bveto = 0;
	double etsmrw_bveto = 0;
	double htrw_bveto = 0;
	double njetrw_bveto = 0;
	double nbjetrw_bveto = 0;
	double ptrw_err_bveto = 0;
	double ptsmrw_err_bveto = 0;
	double etrw_err_bveto = 0;
	double etsmrw_err_bveto = 0;
	double htrw_err_bveto = 0;
	double njetrw_err_bveto = 0;
	double nbjetrw_err_bveto = 0;
	//TBranch *b_pt37_bveto_corr_met = outputTree->Branch("pt37_bveto_corr_met",&pt37_bveto_corr_met,"pt37_bveto_corr_met/D");
	TBranch *b_pt37_bveto_corr_lpt = outputTree->Branch("pt37_bveto_corr_lpt",&pt37_bveto_corr_lpt,"pt37_bveto_corr_lpt/D");
	TBranch *b_ptrw_bveto = outputTree->Branch("ptrw_bveto",&ptrw_bveto,"ptrw_bveto/D");
	TBranch *b_ptsmrw_bveto = outputTree->Branch("ptsmrw_bveto",&ptsmrw_bveto,"ptsmrw_bveto/D");
	//TBranch *b_etrw_bveto = outputTree->Branch("etrw_bveto",&etrw_bveto,"etrw/D");
	//TBranch *b_etsmrw_bveto = outputTree->Branch("etsmrw_bveto",&etsmrw_bveto,"etsmrw_bveto/D");
	TBranch *b_htrw_bveto = outputTree->Branch("htrw_bveto",&htrw_bveto,"htrw_bveto/D");
	//TBranch *b_njetrw_bveto = outputTree->Branch("njetrw_bveto",&njetrw_bveto,"njetrw_bveto/D");
	//TBranch *b_nbjetrw_bveto = outputTree->Branch("nbjetrw_bveto",&nbjetrw_bveto,"nbjetrw_bveto/D");
	//TBranch *b_ptrw_err_bveto = outputTree->Branch("ptrw_err_bveto",&ptrw_err_bveto,"ptrw_err_bveto/D");
	//TBranch *b_ptsmrw_err_bveto = outputTree->Branch("ptsmrw_err_bveto",&ptsmrw_err_bveto,"ptsmrw_err_bveto/D");
	//TBranch *b_etrw_err_bveto = outputTree->Branch("etrw_err_bveto",&etrw_err_bveto,"etrw_err/D");
	//TBranch *b_etsmrw_err_bveto = outputTree->Branch("etsmrw_err_bveto",&etsmrw_err_bveto,"etsmrw_err_bveto/D");
	//TBranch *b_htrw_err_bveto = outputTree->Branch("htrw_err_bveto",&htrw_err_bveto,"htrw_err_bveto/D");
	//TBranch *b_njetrw_err_bveto = outputTree->Branch("njetrw_err_bveto",&njetrw_err_bveto,"njetrw_err_bveto/D");
	//TBranch *b_nbjetrw_err_bveto = outputTree->Branch("nbjetrw_err_bveto",&nbjetrw_err_bveto,"nbjetrw_err_bveto/D");

	double pt37_2j30_corr_met = 0;
	double pt37_2j30_corr_lpt = 0;
	double ptrw_2j30 = 0;
	double ptsmrw_2j30 = 0;
	double etrw_2j30 = 0;
	double etsmrw_2j30 = 0;
	double htrw_2j30 = 0;
	double njetrw_2j30 = 0;
	double nbjetrw_2j30 = 0;
	double ptrw_err_2j30 = 0;
	double ptsmrw_err_2j30 = 0;
	double etrw_err_2j30 = 0;
	double etsmrw_err_2j30 = 0;
	double htrw_err_2j30 = 0;
	double njetrw_err_2j30 = 0;
	double nbjetrw_err_2j30 = 0;
	//TBranch *b_pt37_2j30_corr_met = outputTree->Branch("pt37_2j30_corr_met",&pt37_2j30_corr_met,"pt37_2j30_corr_met/D");
	TBranch *b_pt37_2j30_corr_lpt = outputTree->Branch("pt37_2j30_corr_lpt",&pt37_2j30_corr_lpt,"pt37_2j30_corr_lpt/D");
	TBranch *b_ptrw_2j30 = outputTree->Branch("ptrw_2j30",&ptrw_2j30,"ptrw_2j30/D");
	TBranch *b_ptsmrw_2j30 = outputTree->Branch("ptsmrw_2j30",&ptsmrw_2j30,"ptsmrw_2j30/D");
	//TBranch *b_etrw_2j30 = outputTree->Branch("etrw_2j30",&etrw_2j30,"etrw/D");
	//TBranch *b_etsmrw_2j30 = outputTree->Branch("etsmrw_2j30",&etsmrw_2j30,"etsmrw_2j30/D");
	TBranch *b_htrw_2j30 = outputTree->Branch("htrw_2j30",&htrw_2j30,"htrw_2j30/D");
	//TBranch *b_njetrw_2j30 = outputTree->Branch("njetrw_2j30",&njetrw_2j30,"njetrw_2j30/D");
	//TBranch *b_nbjetrw_2j30 = outputTree->Branch("nbjetrw_2j30",&nbjetrw_2j30,"nbjetrw_2j30/D");
	//TBranch *b_ptrw_err_2j30 = outputTree->Branch("ptrw_err_2j30",&ptrw_err_2j30,"ptrw_err_2j30/D");
	//TBranch *b_ptsmrw_err_2j30 = outputTree->Branch("ptsmrw_err_2j30",&ptsmrw_err_2j30,"ptsmrw_err_2j30/D");
	//TBranch *b_etrw_err_2j30 = outputTree->Branch("etrw_err_2j30",&etrw_err_2j30,"etrw_err/D");
	//TBranch *b_etsmrw_err_2j30 = outputTree->Branch("etsmrw_err_2j30",&etsmrw_err_2j30,"etsmrw_err_2j30/D");
	//TBranch *b_htrw_err_2j30 = outputTree->Branch("htrw_err_2j30",&htrw_err_2j30,"htrw_err_2j30/D");
	//TBranch *b_njetrw_err_2j30 = outputTree->Branch("njetrw_err_2j30",&njetrw_err_2j30,"njetrw_err_2j30/D");
	//TBranch *b_nbjetrw_err_2j30 = outputTree->Branch("nbjetrw_err_2j30",&nbjetrw_err_2j30,"nbjetrw_err_2j30/D");

	double pt37_ht200_corr_met = 0;
	double pt37_ht200_corr_lpt = 0;
	double ptrw_ht200 = 0;
	double ptsmrw_ht200 = 0;
	double etrw_ht200 = 0;
	double etsmrw_ht200 = 0;
	double htrw_ht200 = 0;
	double njetrw_ht200 = 0;
	double nbjetrw_ht200 = 0;
	double ptrw_err_ht200 = 0;
	double ptsmrw_err_ht200 = 0;
	double etrw_err_ht200 = 0;
	double etsmrw_err_ht200 = 0;
	double htrw_err_ht200 = 0;
	double njetrw_err_ht200 = 0;
	double nbjetrw_err_ht200 = 0;
	//TBranch *b_pt37_ht200_corr_met = outputTree->Branch("pt37_ht200_corr_met",&pt37_ht200_corr_met,"pt37_ht200_corr_met/D");
	TBranch *b_pt37_ht200_corr_lpt = outputTree->Branch("pt37_ht200_corr_lpt",&pt37_ht200_corr_lpt,"pt37_ht200_corr_lpt/D");
	TBranch *b_ptrw_ht200 = outputTree->Branch("ptrw_ht200",&ptrw_ht200,"ptrw_ht200/D");
	TBranch *b_ptsmrw_ht200 = outputTree->Branch("ptsmrw_ht200",&ptsmrw_ht200,"ptsmrw_ht200/D");
	//TBranch *b_etrw_ht200 = outputTree->Branch("etrw_ht200",&etrw_ht200,"etrw/D");
	//TBranch *b_etsmrw_ht200 = outputTree->Branch("etsmrw_ht200",&etsmrw_ht200,"etsmrw_ht200/D");
	TBranch *b_htrw_ht200 = outputTree->Branch("htrw_ht200",&htrw_ht200,"htrw_ht200/D");
	//TBranch *b_njetrw_ht200 = outputTree->Branch("njetrw_ht200",&njetrw_ht200,"njetrw_ht200/D");
	//TBranch *b_nbjetrw_ht200 = outputTree->Branch("nbjetrw_ht200",&nbjetrw_ht200,"nbjetrw_ht200/D");
	//TBranch *b_ptrw_err_ht200 = outputTree->Branch("ptrw_err_ht200",&ptrw_err_ht200,"ptrw_err_ht200/D");
	//TBranch *b_ptsmrw_err_ht200 = outputTree->Branch("ptsmrw_err_ht200",&ptsmrw_err_ht200,"ptsmrw_err_ht200/D");
	//TBranch *b_etrw_err_ht200 = outputTree->Branch("etrw_err_ht200",&etrw_err_ht200,"etrw_err/D");
	//TBranch *b_etsmrw_err_ht200 = outputTree->Branch("etsmrw_err_ht200",&etsmrw_err_ht200,"etsmrw_err_ht200/D");
	//TBranch *b_htrw_err_ht200 = outputTree->Branch("htrw_err_ht200",&htrw_err_ht200,"htrw_err_ht200/D");
	//TBranch *b_njetrw_err_ht200 = outputTree->Branch("njetrw_err_ht200",&njetrw_err_ht200,"njetrw_err_ht200/D");
	//TBranch *b_nbjetrw_err_ht200 = outputTree->Branch("nbjetrw_err_ht200",&nbjetrw_err_ht200,"nbjetrw_err_ht200/D");

	double pt37_ht400_corr_met = 0;
	double pt37_ht400_corr_lpt = 0;
	double ptrw_ht400 = 0;
	double ptsmrw_ht400 = 0;
	double etrw_ht400 = 0;
	double etsmrw_ht400 = 0;
	double htrw_ht400 = 0;
	double njetrw_ht400 = 0;
	double nbjetrw_ht400 = 0;
	double ptrw_err_ht400 = 0;
	double ptsmrw_err_ht400 = 0;
	double etrw_err_ht400 = 0;
	double etsmrw_err_ht400 = 0;
	double htrw_err_ht400 = 0;
	double njetrw_err_ht400 = 0;
	double nbjetrw_err_ht400 = 0;
	//TBranch *b_pt37_ht400_corr_met = outputTree->Branch("pt37_ht400_corr_met",&pt37_ht400_corr_met,"pt37_ht400_corr_met/D");
	TBranch *b_pt37_ht400_corr_lpt = outputTree->Branch("pt37_ht400_corr_lpt",&pt37_ht400_corr_lpt,"pt37_ht400_corr_lpt/D");
	TBranch *b_ptrw_ht400 = outputTree->Branch("ptrw_ht400",&ptrw_ht400,"ptrw_ht400/D");
	TBranch *b_ptsmrw_ht400 = outputTree->Branch("ptsmrw_ht400",&ptsmrw_ht400,"ptsmrw_ht400/D");
	//TBranch *b_etrw_ht400 = outputTree->Branch("etrw_ht400",&etrw_ht400,"etrw/D");
	//TBranch *b_etsmrw_ht400 = outputTree->Branch("etsmrw_ht400",&etsmrw_ht400,"etsmrw_ht400/D");
	TBranch *b_htrw_ht400 = outputTree->Branch("htrw_ht400",&htrw_ht400,"htrw_ht400/D");
	//TBranch *b_njetrw_ht400 = outputTree->Branch("njetrw_ht400",&njetrw_ht400,"njetrw_ht400/D");
	//TBranch *b_nbjetrw_ht400 = outputTree->Branch("nbjetrw_ht400",&nbjetrw_ht400,"nbjetrw_ht400/D");
	//TBranch *b_ptrw_err_ht400 = outputTree->Branch("ptrw_err_ht400",&ptrw_err_ht400,"ptrw_err_ht400/D");
	//TBranch *b_ptsmrw_err_ht400 = outputTree->Branch("ptsmrw_err_ht400",&ptsmrw_err_ht400,"ptsmrw_err_ht400/D");
	//TBranch *b_etrw_err_ht400 = outputTree->Branch("etrw_err_ht400",&etrw_err_ht400,"etrw_err/D");
	//TBranch *b_etsmrw_err_ht400 = outputTree->Branch("etsmrw_err_ht400",&etsmrw_err_ht400,"etsmrw_err_ht400/D");
	//TBranch *b_htrw_err_ht400 = outputTree->Branch("htrw_err_ht400",&htrw_err_ht400,"htrw_err_ht400/D");
	//TBranch *b_njetrw_err_ht400 = outputTree->Branch("njetrw_err_ht400",&njetrw_err_ht400,"njetrw_err_ht400/D");
	//TBranch *b_nbjetrw_err_ht400 = outputTree->Branch("nbjetrw_err_ht400",&nbjetrw_err_ht400,"nbjetrw_err_ht400/D");

	double pt37_ht1200_corr_met = 0;
	double pt37_ht1200_corr_lpt = 0;
	double ptrw_ht1200 = 0;
	double ptsmrw_ht1200 = 0;
	double etrw_ht1200 = 0;
	double etsmrw_ht1200 = 0;
	double htrw_ht1200 = 0;
	double njetrw_ht1200 = 0;
	double nbjetrw_ht1200 = 0;
	double ptrw_err_ht1200 = 0;
	double ptsmrw_err_ht1200 = 0;
	double etrw_err_ht1200 = 0;
	double etsmrw_err_ht1200 = 0;
	double htrw_err_ht1200 = 0;
	double njetrw_err_ht1200 = 0;
	double nbjetrw_err_ht1200 = 0;
	//TBranch *b_pt37_ht1200_corr_met = outputTree->Branch("pt37_ht1200_corr_met",&pt37_ht1200_corr_met,"pt37_ht1200_corr_met/D");
	TBranch *b_pt37_ht1200_corr_lpt = outputTree->Branch("pt37_ht1200_corr_lpt",&pt37_ht1200_corr_lpt,"pt37_ht1200_corr_lpt/D");
	TBranch *b_ptrw_ht1200 = outputTree->Branch("ptrw_ht1200",&ptrw_ht1200,"ptrw_ht1200/D");
	TBranch *b_ptsmrw_ht1200 = outputTree->Branch("ptsmrw_ht1200",&ptsmrw_ht1200,"ptsmrw_ht1200/D");
	//TBranch *b_etrw_ht1200 = outputTree->Branch("etrw_ht1200",&etrw_ht1200,"etrw/D");
	//TBranch *b_etsmrw_ht1200 = outputTree->Branch("etsmrw_ht1200",&etsmrw_ht1200,"etsmrw_ht1200/D");
	TBranch *b_htrw_ht1200 = outputTree->Branch("htrw_ht1200",&htrw_ht1200,"htrw_ht1200/D");
	//TBranch *b_njetrw_ht1200 = outputTree->Branch("njetrw_ht1200",&njetrw_ht1200,"njetrw_ht1200/D");
	//TBranch *b_nbjetrw_ht1200 = outputTree->Branch("nbjetrw_ht1200",&nbjetrw_ht1200,"nbjetrw_ht1200/D");
	//TBranch *b_ptrw_err_ht1200 = outputTree->Branch("ptrw_err_ht1200",&ptrw_err_ht1200,"ptrw_err_ht1200/D");
	//TBranch *b_ptsmrw_err_ht1200 = outputTree->Branch("ptsmrw_err_ht1200",&ptsmrw_err_ht1200,"ptsmrw_err_ht1200/D");
	//TBranch *b_etrw_err_ht1200 = outputTree->Branch("etrw_err_ht1200",&etrw_err_ht1200,"etrw_err/D");
	//TBranch *b_etsmrw_err_ht1200 = outputTree->Branch("etsmrw_err_ht1200",&etsmrw_err_ht1200,"etsmrw_err_ht1200/D");
	//TBranch *b_htrw_err_ht1200 = outputTree->Branch("htrw_err_ht1200",&htrw_err_ht1200,"htrw_err_ht1200/D");
	//TBranch *b_njetrw_err_ht1200 = outputTree->Branch("njetrw_err_ht1200",&njetrw_err_ht1200,"njetrw_err_ht1200/D");
	//TBranch *b_nbjetrw_err_ht1200 = outputTree->Branch("nbjetrw_err_ht1200",&nbjetrw_err_ht1200,"nbjetrw_err_ht1200/D");

	double pt37_srhigh_corr_met = 0;
	double ptrw_srhigh = 0;
	double ptsmrw_srhigh = 0;
	double etrw_srhigh = 0;
	double etsmrw_srhigh = 0;
	double htrw_srhigh = 0;
	double njetrw_srhigh = 0;
	double nbjetrw_srhigh = 0;
	double ptrw_err_srhigh = 0;
	double ptsmrw_err_srhigh = 0;
	double etrw_err_srhigh = 0;
	double etsmrw_err_srhigh = 0;
	double htrw_err_srhigh = 0;
	double njetrw_err_srhigh = 0;
	double nbjetrw_err_srhigh = 0;
	//TBranch *b_pt37_srhigh_corr_met = outputTree->Branch("pt37_srhigh_corr_met",&pt37_srhigh_corr_met,"pt37_srhigh_corr_met/D");
	TBranch *b_ptrw_srhigh = outputTree->Branch("ptrw_srhigh",&ptrw_srhigh,"ptrw_srhigh/D");
	TBranch *b_ptsmrw_srhigh = outputTree->Branch("ptsmrw_srhigh",&ptsmrw_srhigh,"ptsmrw_srhigh/D");
	//TBranch *b_etrw_srhigh = outputTree->Branch("etrw_srhigh",&etrw_srhigh,"etrw/D");
	//TBranch *b_etsmrw_srhigh = outputTree->Branch("etsmrw_srhigh",&etsmrw_srhigh,"etsmrw_srhigh/D");
	TBranch *b_htrw_srhigh = outputTree->Branch("htrw_srhigh",&htrw_srhigh,"htrw_srhigh/D");
	//TBranch *b_njetrw_srhigh = outputTree->Branch("njetrw_srhigh",&njetrw_srhigh,"njetrw_srhigh/D");
	//TBranch *b_nbjetrw_srhigh = outputTree->Branch("nbjetrw_srhigh",&nbjetrw_srhigh,"nbjetrw_srhigh/D");
	//TBranch *b_ptrw_err_srhigh = outputTree->Branch("ptrw_err_srhigh",&ptrw_err_srhigh,"ptrw_err_srhigh/D");
	//TBranch *b_ptsmrw_err_srhigh = outputTree->Branch("ptsmrw_err_srhigh",&ptsmrw_err_srhigh,"ptsmrw_err_srhigh/D");
	//TBranch *b_etrw_err_srhigh = outputTree->Branch("etrw_err_srhigh",&etrw_err_srhigh,"etrw_err/D");
	//TBranch *b_etsmrw_err_srhigh = outputTree->Branch("etsmrw_err_srhigh",&etsmrw_err_srhigh,"etsmrw_err_srhigh/D");
	//TBranch *b_htrw_err_srhigh = outputTree->Branch("htrw_err_srhigh",&htrw_err_srhigh,"htrw_err_srhigh/D");
	//TBranch *b_njetrw_err_srhigh = outputTree->Branch("njetrw_err_srhigh",&njetrw_err_srhigh,"njetrw_err_srhigh/D");
	//TBranch *b_nbjetrw_err_srhigh = outputTree->Branch("nbjetrw_err_srhigh",&nbjetrw_err_srhigh,"nbjetrw_err_srhigh/D");

	//-----------------------------
	// build TGraph reweighting factors
	//-----------------------------

	TGraph * gr_37_bveto_met_correction = new TGraph(hist_37_bveto_met_correction);
	TGraph * gr_37_bveto_lpt_correction = new TGraph(hist_37_bveto_lpt_correction);
	TGraph * gr_bveto_ptrw_correction = new TGraph(hist_bveto_ptrw_correction);
	TGraph * gr_bveto_ptsmrw_correction = new TGraph(hist_bveto_ptsmrw_correction);
	TGraph * gr_bveto_etrw_correction = new TGraph(hist_bveto_etrw_correction);
	TGraph * gr_bveto_etsmrw_correction = new TGraph(hist_bveto_etsmrw_correction);
	TGraph * gr_bveto_htrw_correction = new TGraph(hist_bveto_htrw_correction);

	TGraph * gr_37_2j30_met_correction = new TGraph(hist_37_2j30_met_correction);
	TGraph * gr_37_2j30_lpt_correction = new TGraph(hist_37_2j30_lpt_correction);
	TGraph * gr_2j30_ptrw_correction = new TGraph(hist_2j30_ptrw_correction);
	TGraph * gr_2j30_ptsmrw_correction = new TGraph(hist_2j30_ptsmrw_correction);
	TGraph * gr_2j30_etrw_correction = new TGraph(hist_2j30_etrw_correction);
	TGraph * gr_2j30_etsmrw_correction = new TGraph(hist_2j30_etsmrw_correction);
	TGraph * gr_2j30_htrw_correction = new TGraph(hist_2j30_htrw_correction);

	TGraph * gr_37_ht200_met_correction = new TGraph(hist_37_ht200_met_correction);
	TGraph * gr_37_ht200_lpt_correction = new TGraph(hist_37_ht200_lpt_correction);
	TGraph * gr_ht200_ptrw_correction = new TGraph(hist_ht200_ptrw_correction);
	TGraph * gr_ht200_ptsmrw_correction = new TGraph(hist_ht200_ptsmrw_correction);
	TGraph * gr_ht200_etrw_correction = new TGraph(hist_ht200_etrw_correction);
	TGraph * gr_ht200_etsmrw_correction = new TGraph(hist_ht200_etsmrw_correction);
	TGraph * gr_ht200_htrw_correction = new TGraph(hist_ht200_htrw_correction);

	TGraph * gr_37_ht400_met_correction = new TGraph(hist_37_ht400_met_correction);
	TGraph * gr_37_ht400_lpt_correction = new TGraph(hist_37_ht400_lpt_correction);
	TGraph * gr_ht400_ptrw_correction = new TGraph(hist_ht400_ptrw_correction);
	TGraph * gr_ht400_ptsmrw_correction = new TGraph(hist_ht400_ptsmrw_correction);
	TGraph * gr_ht400_etrw_correction = new TGraph(hist_ht400_etrw_correction);
	TGraph * gr_ht400_etsmrw_correction = new TGraph(hist_ht400_etsmrw_correction);
	TGraph * gr_ht400_htrw_correction = new TGraph(hist_ht400_htrw_correction);

	TGraph * gr_37_ht1200_met_correction = new TGraph(hist_37_ht1200_met_correction);
	TGraph * gr_37_ht1200_lpt_correction = new TGraph(hist_37_ht1200_lpt_correction);
	TGraph * gr_ht1200_ptrw_correction = new TGraph(hist_ht1200_ptrw_correction);
	TGraph * gr_ht1200_ptsmrw_correction = new TGraph(hist_ht1200_ptsmrw_correction);
	TGraph * gr_ht1200_etrw_correction = new TGraph(hist_ht1200_etrw_correction);
	TGraph * gr_ht1200_etsmrw_correction = new TGraph(hist_ht1200_etsmrw_correction);
	TGraph * gr_ht1200_htrw_correction = new TGraph(hist_ht1200_htrw_correction);

	TGraph * gr_37_srhigh_met_correction = new TGraph(hist_37_srhigh_met_correction);
	TGraph * gr_srhigh_ptrw_correction = new TGraph(hist_srhigh_ptrw_correction);
	TGraph * gr_srhigh_ptsmrw_correction = new TGraph(hist_srhigh_ptsmrw_correction);
	TGraph * gr_srhigh_etrw_correction = new TGraph(hist_srhigh_etrw_correction);
	TGraph * gr_srhigh_etsmrw_correction = new TGraph(hist_srhigh_etsmrw_correction);
	TGraph * gr_srhigh_htrw_correction = new TGraph(hist_srhigh_htrw_correction);

	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = outputTree->GetEntries();

	for (Long64_t i=0;i<nentries;i++) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		outputTree->GetEntry(i);

		// assess the reweighting factor according to the pT/smeared pT/Njets/HT value
		double gamma_et_smear = pow(gamma_pt_smear*gamma_pt_smear+mll*mll,0.5);
		//double gamma_et_smear = gamma_pt_smear;
		pt37_bveto_corr_met = gr_37_bveto_met_correction->Eval(MET_smear);
		pt37_bveto_corr_lpt = gr_37_bveto_lpt_correction->Eval(lep_pT->at(0));
		ptrw_bveto = gr_bveto_ptrw_correction->Eval(gamma_pt);
		ptsmrw_bveto = gr_bveto_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_bveto = gr_bveto_etrw_correction->Eval(gamma_pt);
		etsmrw_bveto = gr_bveto_etsmrw_correction->Eval(gamma_et_smear);
		htrw_bveto = gr_bveto_htrw_correction->Eval(HT);
		njetrw_bveto = hist_bveto_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_bveto = hist_bveto_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_bveto = hist_bveto_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_bveto = hist_bveto_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_bveto = hist_bveto_etrw_correction->GetBinError(pt+1);
		etsmrw_err_bveto = hist_bveto_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_bveto = hist_bveto_htrw_correction->GetBinError(ht+1);
		njetrw_err_bveto = hist_bveto_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_bveto = hist_bveto_nbjetrw_correction->GetBinError(nbjet+1);

		pt37_2j30_corr_met = gr_37_2j30_met_correction->Eval(MET_smear);
		pt37_2j30_corr_lpt = gr_37_2j30_lpt_correction->Eval(lep_pT->at(1));
		ptrw_2j30 = gr_2j30_ptrw_correction->Eval(gamma_pt);
		ptsmrw_2j30 = gr_2j30_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_2j30 = gr_2j30_etrw_correction->Eval(gamma_pt);
		etsmrw_2j30 = gr_2j30_etsmrw_correction->Eval(gamma_et_smear);
		htrw_2j30 = gr_2j30_htrw_correction->Eval(HT);
		njetrw_2j30 = hist_2j30_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_2j30 = hist_2j30_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_2j30 = hist_2j30_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_2j30 = hist_2j30_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_2j30 = hist_2j30_etrw_correction->GetBinError(pt+1);
		etsmrw_err_2j30 = hist_2j30_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_2j30 = hist_2j30_htrw_correction->GetBinError(ht+1);
		njetrw_err_2j30 = hist_2j30_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_2j30 = hist_2j30_nbjetrw_correction->GetBinError(nbjet+1);

		pt37_ht200_corr_met = gr_37_ht200_met_correction->Eval(MET_smear);
		pt37_ht200_corr_lpt = gr_37_ht200_lpt_correction->Eval(lep_pT->at(1));
		ptrw_ht200 = gr_ht200_ptrw_correction->Eval(gamma_pt);
		ptsmrw_ht200 = gr_ht200_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_ht200 = gr_ht200_etrw_correction->Eval(gamma_pt);
		etsmrw_ht200 = gr_ht200_etsmrw_correction->Eval(gamma_et_smear);
		htrw_ht200 = gr_ht200_htrw_correction->Eval(HT);
		njetrw_ht200 = hist_ht200_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_ht200 = hist_ht200_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_ht200 = hist_ht200_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_ht200 = hist_ht200_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_ht200 = hist_ht200_etrw_correction->GetBinError(pt+1);
		etsmrw_err_ht200 = hist_ht200_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_ht200 = hist_ht200_htrw_correction->GetBinError(ht+1);
		njetrw_err_ht200 = hist_ht200_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_ht200 = hist_ht200_nbjetrw_correction->GetBinError(nbjet+1);

		pt37_ht400_corr_met = gr_37_ht400_met_correction->Eval(MET_smear);
		pt37_ht400_corr_lpt = gr_37_ht400_lpt_correction->Eval(lep_pT->at(1));
		ptrw_ht400 = gr_ht400_ptrw_correction->Eval(gamma_pt);
		ptsmrw_ht400 = gr_ht400_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_ht400 = gr_ht400_etrw_correction->Eval(gamma_pt);
		etsmrw_ht400 = gr_ht400_etsmrw_correction->Eval(gamma_et_smear);
		htrw_ht400 = gr_ht400_htrw_correction->Eval(HT);
		njetrw_ht400 = hist_ht400_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_ht400 = hist_ht400_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_ht400 = hist_ht400_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_ht400 = hist_ht400_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_ht400 = hist_ht400_etrw_correction->GetBinError(pt+1);
		etsmrw_err_ht400 = hist_ht400_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_ht400 = hist_ht400_htrw_correction->GetBinError(ht+1);
		njetrw_err_ht400 = hist_ht400_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_ht400 = hist_ht400_nbjetrw_correction->GetBinError(nbjet+1);

		pt37_ht1200_corr_met = gr_37_ht1200_met_correction->Eval(MET_smear);
		pt37_ht1200_corr_lpt = gr_37_ht1200_lpt_correction->Eval(lep_pT->at(1));
		ptrw_ht1200 = gr_ht1200_ptrw_correction->Eval(gamma_pt);
		ptsmrw_ht1200 = gr_ht1200_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_ht1200 = gr_ht1200_etrw_correction->Eval(gamma_pt);
		etsmrw_ht1200 = gr_ht1200_etsmrw_correction->Eval(gamma_et_smear);
		htrw_ht1200 = gr_ht1200_htrw_correction->Eval(HT);
		njetrw_ht1200 = hist_ht1200_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_ht1200 = hist_ht1200_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_ht1200 = hist_ht1200_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_ht1200 = hist_ht1200_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_ht1200 = hist_ht1200_etrw_correction->GetBinError(pt+1);
		etsmrw_err_ht1200 = hist_ht1200_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_ht1200 = hist_ht1200_htrw_correction->GetBinError(ht+1);
		njetrw_err_ht1200 = hist_ht1200_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_ht1200 = hist_ht1200_nbjetrw_correction->GetBinError(nbjet+1);

		pt37_srhigh_corr_met = gr_37_srhigh_met_correction->Eval(MET_smear);
		ptrw_srhigh = gr_srhigh_ptrw_correction->Eval(gamma_pt);
		ptsmrw_srhigh = gr_srhigh_ptsmrw_correction->Eval(gamma_pt_smear);
		etrw_srhigh = gr_srhigh_etrw_correction->Eval(gamma_pt);
		etsmrw_srhigh = gr_srhigh_etsmrw_correction->Eval(gamma_et_smear);
		htrw_srhigh = gr_srhigh_htrw_correction->Eval(HT);
		njetrw_srhigh = hist_srhigh_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_srhigh = hist_srhigh_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_srhigh = hist_srhigh_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_srhigh = hist_srhigh_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_srhigh = hist_srhigh_etrw_correction->GetBinError(pt+1);
		etsmrw_err_srhigh = hist_srhigh_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_srhigh = hist_srhigh_htrw_correction->GetBinError(ht+1);
		njetrw_err_srhigh = hist_srhigh_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_srhigh = hist_srhigh_nbjetrw_correction->GetBinError(nbjet+1);

		//b_pt37_bveto_corr_met->Fill();     
		b_pt37_bveto_corr_lpt->Fill();     
		b_ptrw_bveto->Fill();     
		b_ptsmrw_bveto->Fill();     
		//b_etrw_bveto->Fill();     
		//b_etsmrw_bveto->Fill();     
		b_htrw_bveto->Fill();     
		//b_njetrw_bveto->Fill();     
		//b_nbjetrw_bveto->Fill();     
		//b_ptrw_err_bveto->Fill();     
		//b_ptsmrw_err_bveto->Fill();     
		//b_etrw_err_bveto->Fill();     
		//b_etsmrw_err_bveto->Fill();     
		//b_htrw_err_bveto->Fill();     
		//b_njetrw_err_bveto->Fill();     
		//b_nbjetrw_err_bveto->Fill();     

		//b_pt37_2j30_corr_met->Fill();     
		b_pt37_2j30_corr_lpt->Fill();     
		b_ptrw_2j30->Fill();     
		b_ptsmrw_2j30->Fill();     
		//b_etrw_2j30->Fill();     
		//b_etsmrw_2j30->Fill();     
		b_htrw_2j30->Fill();     
		//b_njetrw_2j30->Fill();     
		//b_nbjetrw_2j30->Fill();     
		//b_ptrw_err_2j30->Fill();     
		//b_ptsmrw_err_2j30->Fill();     
		//b_etrw_err_2j30->Fill();     
		//b_etsmrw_err_2j30->Fill();     
		//b_htrw_err_2j30->Fill();     
		//b_njetrw_err_2j30->Fill();     
		//b_nbjetrw_err_2j30->Fill();     

		//b_pt37_ht200_corr_met->Fill();     
		b_pt37_ht200_corr_lpt->Fill();     
		b_ptrw_ht200->Fill();     
		b_ptsmrw_ht200->Fill();     
		//b_etrw_ht200->Fill();     
		//b_etsmrw_ht200->Fill();     
		b_htrw_ht200->Fill();     
		//b_njetrw_ht200->Fill();     
		//b_nbjetrw_ht200->Fill();     
		//b_ptrw_err_ht200->Fill();     
		//b_ptsmrw_err_ht200->Fill();     
		//b_etrw_err_ht200->Fill();     
		//b_etsmrw_err_ht200->Fill();     
		//b_htrw_err_ht200->Fill();     
		//b_njetrw_err_ht200->Fill();     
		//b_nbjetrw_err_ht200->Fill();     

		//b_pt37_ht400_corr_met->Fill();     
		b_pt37_ht400_corr_lpt->Fill();     
		b_ptrw_ht400->Fill();     
		b_ptsmrw_ht400->Fill();     
		//b_etrw_ht400->Fill();     
		//b_etsmrw_ht400->Fill();     
		b_htrw_ht400->Fill();     
		//b_njetrw_ht400->Fill();     
		//b_nbjetrw_ht400->Fill();     
		//b_ptrw_err_ht400->Fill();     
		//b_ptsmrw_err_ht400->Fill();     
		//b_etrw_err_ht400->Fill();     
		//b_etsmrw_err_ht400->Fill();     
		//b_htrw_err_ht400->Fill();     
		//b_njetrw_err_ht400->Fill();     
		//b_nbjetrw_err_ht400->Fill();     

		//b_pt37_ht1200_corr_met->Fill();     
		b_pt37_ht1200_corr_lpt->Fill();     
		b_ptrw_ht1200->Fill();     
		b_ptsmrw_ht1200->Fill();     
		//b_etrw_ht1200->Fill();     
		//b_etsmrw_ht1200->Fill();     
		b_htrw_ht1200->Fill();     
		//b_njetrw_ht1200->Fill();     
		//b_nbjetrw_ht1200->Fill();     
		//b_ptrw_err_ht1200->Fill();     
		//b_ptsmrw_err_ht1200->Fill();     
		//b_etrw_err_ht1200->Fill();     
		//b_etsmrw_err_ht1200->Fill();     
		//b_htrw_err_ht1200->Fill();     
		//b_njetrw_err_ht1200->Fill();     
		//b_nbjetrw_err_ht1200->Fill();     

		//b_pt37_srhigh_corr_met->Fill();     
		b_ptrw_srhigh->Fill();     
		b_ptsmrw_srhigh->Fill();     
		//b_etrw_srhigh->Fill();     
		//b_etsmrw_srhigh->Fill();     
		b_htrw_srhigh->Fill();     
		//b_njetrw_srhigh->Fill();     
		//b_nbjetrw_srhigh->Fill();     
		//b_ptrw_err_srhigh->Fill();     
		//b_ptsmrw_err_srhigh->Fill();     
		//b_etrw_err_srhigh->Fill();     
		//b_etsmrw_err_srhigh->Fill();     
		//b_htrw_err_srhigh->Fill();     
		//b_njetrw_err_srhigh->Fill();     
		//b_nbjetrw_err_srhigh->Fill();     

	}

	outputTree->Write();
	hist_bveto_ptrw_correction->Write();
	hist_ht200_ptrw_correction->Write();
	hist_ht400_ptrw_correction->Write();
	hist_ht1200_ptrw_correction->Write();
	gr_bveto_ptrw_correction->Write();
	gr_ht200_ptrw_correction->Write();
	gr_ht400_ptrw_correction->Write();
	gr_ht1200_ptrw_correction->Write();

	std::cout << "done." << std::endl;
	delete f;


}
