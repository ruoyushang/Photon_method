
Long64_t EventNumber;
Int_t RunNumber;
double totalWeight = 0.;
int pt = 0; // pT bin for reweighting
int smpt = 0; // pT bin for smearing
int ht = 0;
int njet = 0;
int nbjet = 0;
Int_t jet_n;
Int_t bjet_n;
double HT = 0.;
double truthGamma_pt = 0.;
double truthGamma_phi = 0.;
double gamma_pt = 0.;
double gamma_eta = 0.;
double gamma_phi = 0.;
double gamma_dR = 0.;
double MET;
double METl;
double METt;
double MET_phi;
double DPhi_METJetLeading;
double DPhi_METJetSecond;
double MinDPhi_PhotonJet;
std::vector<double>* jet_pT = new std::vector<double>(10);
std::vector<double>* jet_eta = new std::vector<double>(10);
std::vector<double>* jet_phi = new std::vector<double>(10);
std::vector<double>* jet_m = new std::vector<double>(10);
std::vector<double>* jet_btag = new std::vector<double>(10);

// new variables to be added to ntuple
int lep_n = 0;
int dpt = 0;
int pt_smear = 0;
int met_smear = 0;
int dphi_smear = 0;
double HTincl = 0.;
double mll = 0.;
double smear_shift = 0.;
double gamma_pt_smear = 0.;
double gamma_phi_smear = 0.;
double MET_smear;
double MET_phi_smear;
double DPhi_METJetLeading_smear;
double DPhi_METJetSecond_smear;
double DPhi_METLepLeading_smear;
double DPhi_METLepSecond_smear;
double METl_smear = 0.;
double METt_smear = 0.;
double DPhi_METPhoton_smear = 0.;
double mj1j2 = 0.;
double DR_J1J2 = 0.;
double DR_Wmin2Jet = 0.;
double DR_J0J1 = 0.;
double DR_W80jj = 0.;
double mWmin = 0.;
double Wmin_pt = 0.;
double Wmin_eta = 0.;
double DPhi_METWmin = 0.;
double DPhi_WminZ = 0.;
double mj0j1 = 0.;
double m80jj = 0.;
double W01_pt = 0.;
double DPhi_METW01 = 0.;
double DPhi_W01Z = 0.;
double DPhi_W80Z = 0.;
double DPhi_METNonWminJet = 0.;
double NonWminJet_pT = 0.;
double W80_pt = 0.;
double DPhi_METW80 = 0.;
std::vector<int>* jet_isW80 = new std::vector<int>(10);
std::vector<int>* jet_isW01 = new std::vector<int>(10);
std::vector<int>* jet_isW12 = new std::vector<int>(10);
std::vector<int>* jet_isWmin = new std::vector<int>(10);
std::vector<int>* jet_isWminJ0 = new std::vector<int>(10);
std::vector<int>* jet_isWminJ1 = new std::vector<int>(10);
double DPhi_METNonWJet = 0.;
double NonWJet_pT = 0.;
double DPhi_METNonW12Jet = 0.;
double NonW12Jet_pT = 0.;

std::vector<double>* lep_pT = new std::vector<double>(10);
std::vector<double>* lep_eta = new std::vector<double>(10);
std::vector<double>* lep_phi = new std::vector<double>(10);

float MT2W;
double DR_2Lep = 0.;

