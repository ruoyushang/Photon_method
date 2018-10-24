//TH1D* z_dpt[bin_size];
TH1D* z_dphi[bin_size];
//TH1D* g_dpt[bin_size];
TH1D* z_metl[bin_size];
TH1D* z_metl_2j[bin_size];
TH1D* tt_metl[bin_size];
TH1D* vv_metl[bin_size];
TH1D* g_metl[bin_size];
TH1D* z_jetmetl[bin_size];
TH1D* tt_jetmetl[bin_size];
TH1D* vv_jetmetl[bin_size];

TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,sm_pt_bin);

//double totalWeight = 0.;
//int jet_n = 0;
//int bjet_n = 0;
//float Z_pt = 0.;
//float Z_truthPt = 0.;
//float gamma_pt = 0.;
//float HT = 0.;
//float mll = 0.;
//float METl = 0.;

//Zfloat totalWeightF = 0.;
float METlF = 0.;
float HTF = 0.;

//smearing_method =
// 0 : no smearing
// 1 : R20 MC smearing
// 2 : R20 data smearing
// 3 : truth smearing DO NOT USE
// 4 : R21 MC smearing
// 5 : R21 data smearing

int gchannel;
Float_t gZ_pt;

void GetSmearingHistogram(string ch, float lumi, string photon_tag,string period, int smearing_method) {

  cout << "GetSmearingHistogram : smearing_method " << smearing_method << endl;

  string mcperiod = "";
  string gmcperiod = "";
  string zperiod = "";
  if( TString(period).Contains("data15-16") ) {
    mcperiod  = "zmc16a/";
    gmcperiod = "gmc16a/";
    zperiod   = "data15-16";
  }
  if( TString(period).Contains("data17")    ){
    mcperiod  = "zmc16cd/";
    gmcperiod = "gmc16cd/";
    zperiod   = "data17";
  }
  if( TString(period).Contains("data18")    ){
    mcperiod  = "zmc16cd/";
    gmcperiod = "gmc16cd/";
    zperiod   = "data18";
  }
  
	for (int bin=0;bin<bin_size;bin++) {
	  //z_dpt[bin] = new TH1D(TString("z_dpt_")+TString::Itoa(bin,10),"",40000,-30000,10000);//
	  //g_dpt[bin] = new TH1D(TString("g_dpt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_metl_2j[bin] = new TH1D(TString("z_metl_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_jetmetl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
	}
	if (smearing_method == 5) {  // R21 data-driven smearing function
		std::cout << "Get smearing function from R21 data." << std::endl;

		//--- smearing with R20.7 samples
		//cout << "Opening file           : " << TString(smearingPath)+"data/data_"+TString(ch)+".root"  << endl;
		//TFile fZ( TString(smearingPath)+"data/data_"+TString(ch)+".root" );

		//--- smearing with R21 samples
		string datafilename = smearingPath + "zdata/" + zperiod + "_merged_processed.root";
		cout << "Opening data smearing file   : " << datafilename << endl;
		TFile fZ( datafilename.c_str() );

		TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		//tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&gZ_pt);
		//tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mll);
		tZ->SetBranchAddress("METl" ,&METl);
		tZ->SetBranchAddress("channel" ,&gchannel);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
			if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
			if (gZ_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gZ_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			z_metl[pt]->Fill(METlF,totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METlF,totalWeight);
		}

		//tZ->Close()
		fZ.Close();

		//--- smearing with R20.7 samples
		//cout << "Opening file           : " << TString(smearingPath)+"tt/tt"+TString(ch)+".root"        << endl;
		//TFile ftt( TString(smearingPath)+"tt/tt"+TString(ch)+".root" );

		//--- smearing R21 samples

		string ttfilename = smearingPath + mcperiod + "ttbar_410472_dilep_processed.root";
		if( TString(mcperiod).Contains("zmc16cd") ) ttfilename = smearingPath + mcperiod + "ttbar_dilep_processed.root";

		cout << "Opening tt smearing file   : " << ttfilename << endl;
		TFile ftt( ttfilename.c_str() );

		TTree*  ttt              = (TTree*)ftt.Get("BaselineTree");
		ttt->SetBranchStatus("*", 0);
		ttt->SetBranchStatus("totalWeight", 1);
		ttt->SetBranchStatus("jet_n", 1);
		ttt->SetBranchStatus("bjet_n", 1);
		ttt->SetBranchStatus("Z_pt", 1);
		//ttt->SetBranchStatus("HT", 1);
		ttt->SetBranchStatus("mll", 1);
		ttt->SetBranchStatus("METl", 1);
		ttt->SetBranchAddress("totalWeight" ,&totalWeight);
		ttt->SetBranchAddress("jet_n" ,&jet_n);
		ttt->SetBranchAddress("bjet_n" ,&bjet_n);
		ttt->SetBranchAddress("Z_pt" ,&gZ_pt);
		//ttt->SetBranchAddress("HT" ,&HT);
		ttt->SetBranchAddress("mll" ,&mll);
		ttt->SetBranchAddress("METl" ,&METlF);
		ttt->SetBranchAddress("channel" ,&gchannel);
		ttt->SetBranchAddress("RandomRunNumber" ,&RandomRunNumber);
		for (int entry=0;entry<ttt->GetEntries();entry++) {
			ttt->GetEntry(entry);
			if( TString(period).Contains("data17") && RandomRunNumber > 348000 ) continue; 
			if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
			if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee			
			if (gZ_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gZ_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			z_metl[pt]->Fill(METlF,-1.*lumi*totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METlF,-1.*lumi*totalWeight);
		}
		//ttt->Close();
		ftt.Close();

		//--- smearing with R20.7 samples
		//cout << "Opening file           : " << TString(smearingPath)+"vv/vv"+TString(ch)+".root"       << endl;
		//TFile fvv( TString(smearingPath)+"vv/vv"+TString(ch)+".root" );

		// smearing with R21 samples
		string vvfilename = smearingPath + mcperiod + "diboson_merged_processed.root";
		cout << "Opening VV smearing file   : " << vvfilename << endl;
		TFile fvv( vvfilename.c_str() );

		TTree*  tvv              = (TTree*)fvv.Get("BaselineTree");
		tvv->SetBranchStatus("*", 0);
		tvv->SetBranchStatus("totalWeight", 1);
		tvv->SetBranchStatus("jet_n", 1);
		tvv->SetBranchStatus("bjet_n", 1);
		tvv->SetBranchStatus("Z_pt", 1);
		//tvv->SetBranchStatus("HT", 1);
		tvv->SetBranchStatus("mll", 1);
		tvv->SetBranchStatus("METl", 1);
		tvv->SetBranchAddress("totalWeight" ,&totalWeight);
		tvv->SetBranchAddress("jet_n" ,&jet_n);
		tvv->SetBranchAddress("bjet_n" ,&bjet_n);
		tvv->SetBranchAddress("Z_pt" ,&gZ_pt);
		//tvv->SetBranchAddress("HT" ,&HT);
		tvv->SetBranchAddress("mll" ,&mll);
		tvv->SetBranchAddress("METl" ,&METlF);
		tvv->SetBranchAddress("channel" ,&gchannel);
		tvv->SetBranchAddress("RandomRunNumber" ,&RandomRunNumber);
		for (int entry=0;entry<tvv->GetEntries();entry++) {
			tvv->GetEntry(entry);
			if( TString(period).EqualTo("data17") && RandomRunNumber > 348000 ) continue; 
			if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
			if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
			if (gZ_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gZ_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			z_metl[pt]->Fill(METlF,-1.*lumi*totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METlF,-1.*lumi*totalWeight);
		}
		//tvv->Close();
		fvv.Close();

		//--- smearing with R20.7 samples
		//TFile fPhoton( TString(smearingPath)+"gdata/gdata_raw.root" );
		
		string gperiod = "";
		if( TString(period).EqualTo("data15-16") ) gperiod = "Data15-16";
		if( TString(period).EqualTo("data17")    ) gperiod = "Data17";

		string gfilename = smearingPath + "gdata/" + period + "_merged_processed.root";
		cout << "Opening photon smearing file   : " << gfilename << endl;
		TFile fPhoton( gfilename.c_str() );

		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");

		cout << "Setting photon branches" << endl;
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		//tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		//tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&MET);
		cout << "Done setting photon branches" << endl;
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			g_metl[pt]->Fill(METl,totalWeight);
		}
		//tPhoton->Close();
		fPhoton.Close();
	}
	else if (smearing_method == 1) { // R20 MC-driven smearing function

	  float Z_ptD;
	  float mllD;
	  
		std::cout << "Get smearing function from R20 MC." << std::endl;
		cout << "Opening Z+jets MC smearing file           : " << TString(oldSmearingPath)+"zjets/zjets_"+TString(ch)+".root"  << endl;
		TFile fZ( TString(oldSmearingPath)+"zjets/zjets_"+TString(ch)+".root" );
		TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		//tZ->SetBranchStatus("Z_truthPt", 1);
		tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchStatus("RunNumber", 1);
		tZ->SetBranchStatus("EventNumber", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&Z_ptD);
		//tZ->SetBranchAddress("Z_truthPt" ,&Z_truthPt);
		tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mllD);
		tZ->SetBranchAddress("METl" ,&METl);
		tZ->SetBranchAddress("RunNumber" ,&RunNumber);
		tZ->SetBranchAddress("EventNumber" ,&EventNumber);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if (Z_ptD<50.) continue;
			int pt = hist_low_pt->FindBin(Z_ptD)-1;
			if (jet_n==0) continue;
			if (jet_n>=2) z_metl_2j[pt]->Fill(METl,totalWeight);
			if (jet_n!=1) continue;
			if (EventNumber==62237 && RunNumber==361405) std::cout << "kick out large weight Ev 62237." << std::endl;
			if (EventNumber==62237 && RunNumber==361405) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==1416 && RunNumber==361414) std::cout << "kick out large weight Ev 1416." << std::endl;
			if (EventNumber==1416 && RunNumber==361414) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==3411 && RunNumber==361415) std::cout << "kick out large weight Ev 3411." << std::endl;
			if (EventNumber==3411 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==1449 && RunNumber==361415) std::cout << "kick out large weight Ev 1449." << std::endl;
			if (EventNumber==1449 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			z_metl[pt]->Fill(METl,totalWeight);
			//z_dpt[pt]->Fill(Z_truthPt-Z_pt,totalWeight);
			if (mllD<90 || mllD>92) continue;
			z_jetmetl[pt]->Fill(METl,totalWeight);
		}
		fZ.Close();

		string gmcfilename = oldSmearingPath + "gmc/gmc_raw.root";
		cout << "Opening photon MC smearing file " << gmcfilename << endl;
		
		TFile fPhoton( gmcfilename.c_str() );

		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		//tPhoton->SetBranchStatus("truthGamma_pt", 1);
		tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		//tPhoton->SetBranchAddress("truthGamma_pt" ,&truthGamma_pt);
		tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&METl);
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n==0) continue;
			if (jet_n!=1) continue;
			g_metl[pt]->Fill(METl,totalWeight);
			//g_dpt[pt]->Fill(truthGamma_pt-gamma_pt,totalWeight);
		}
		fPhoton.Close();
	}
	else if (smearing_method == 4) { // R21 MC-driven smearing function
	  
	        float Z_ptD;
	        float mllD;
	  
		std::cout << "Get smearing function from R21 MC." << std::endl;

		string Zfilename = smearingPath + mcperiod + "Zjets_merged_processed.root";
		
		cout << "Opening Z+jets MC smearing file           : " << Zfilename << endl;
		
		TFile fZ( Zfilename.c_str() );

		TTree*  tZ = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		//tZ->SetBranchStatus("Z_truthPt", 1);
		//tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchStatus("RunNumber", 1);
		//tZ->SetBranchStatus("EventNumber", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&Z_ptD);
		//tZ->SetBranchAddress("Z_truthPt" ,&Z_truthPt);
		//tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mllD);
		tZ->SetBranchAddress("METl" ,&METl);
		tZ->SetBranchAddress("RunNumber" ,&RunNumber);
		//tZ->SetBranchAddress("EventNumber" ,&EventNumber);
		tZ->SetBranchAddress("channel" ,&gchannel);
		tZ->SetBranchAddress("RandomRunNumber" ,&RandomRunNumber);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if( TString(period).EqualTo("data17") && RandomRunNumber > 348000 ) continue; 
			if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
			if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee			

			if (Z_ptD<50.) continue;
			int pt = hist_low_pt->FindBin(Z_ptD)-1;
			if (jet_n==0) continue;
			if (jet_n>=2) z_metl_2j[pt]->Fill(METl,totalWeight);
			if (jet_n!=1) continue;
			// if (EventNumber==62237 && RunNumber==361405) std::cout << "kick out large weight Ev 62237." << std::endl;
			// if (EventNumber==62237 && RunNumber==361405) continue;  // this is a large weight event that screws the smearing function
			// if (EventNumber==1416 && RunNumber==361414) std::cout << "kick out large weight Ev 1416." << std::endl;
			// if (EventNumber==1416 && RunNumber==361414) continue;  // this is a large weight event that screws the smearing function
			// if (EventNumber==3411 && RunNumber==361415) std::cout << "kick out large weight Ev 3411." << std::endl;
			// if (EventNumber==3411 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			// if (EventNumber==1449 && RunNumber==361415) std::cout << "kick out large weight Ev 1449." << std::endl;
			// if (EventNumber==1449 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			z_metl[pt]->Fill(METl,totalWeight);
			//z_dpt[pt]->Fill(Z_truthPt-Z_pt,totalWeight);
			if (mllD<90 || mllD>92) continue;
			z_jetmetl[pt]->Fill(METl,totalWeight);
		}
		fZ.Close();

		string gmcfilename = smearingPath + gmcperiod + "SinglePhoton222_merged_processed.root";

		cout << "Opening photon MC smearing file " << gmcfilename << endl;
		
		TFile fPhoton( gmcfilename.c_str() );

		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		//tPhoton->SetBranchStatus("truthGamma_pt", 1);
		//tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		//tPhoton->SetBranchAddress("truthGamma_pt" ,&truthGamma_pt);
		tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&METl);
		tPhoton->SetBranchAddress("RandomRunNumber" ,&RandomRunNumber);
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if( TString(period).EqualTo("data17") && RandomRunNumber > 348000 ) continue; 
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n==0) continue;
			if (jet_n!=1) continue;
			g_metl[pt]->Fill(METl,totalWeight);
			//g_dpt[pt]->Fill(truthGamma_pt-gamma_pt,totalWeight);
		}
		fPhoton.Close();
	}
	else if (smearing_method == 2) {  // data-driven smearing function
		std::cout << "Get smearing function from data." << std::endl;
		cout << "Opening file           : " << TString(oldSmearingPath)+"data/data_"+TString(ch)+".root"  << endl;
		TFile fZ( TString(oldSmearingPath)+"data/data_"+TString(ch)+".root" );
		TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&gZ_pt);
		tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mll);
		tZ->SetBranchAddress("METl" ,&METl);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if (gZ_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gZ_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			z_metl[pt]->Fill(METl,totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METl,totalWeight);
		}
		fZ.Close();
		//cout << "Opening file           : " << TString(smearingPath)+"tt/tt"+TString(ch)+".root"        << endl;
		//TFile ftt( TString(smearingPath)+"tt/tt"+TString(ch)+".root" );
		//TTree*  ttt              = (TTree*)ftt.Get("BaselineTree");
		//ttt->SetBranchStatus("*", 0);
		//ttt->SetBranchStatus("totalWeight", 1);
		//ttt->SetBranchStatus("jet_n", 1);
		//ttt->SetBranchStatus("bjet_n", 1);
		//ttt->SetBranchStatus("Z_pt", 1);
		//ttt->SetBranchStatus("HT", 1);
		//ttt->SetBranchStatus("mll", 1);
		//ttt->SetBranchStatus("METl", 1);
		//ttt->SetBranchAddress("totalWeight" ,&totalWeight);
		//ttt->SetBranchAddress("jet_n" ,&jet_n);
		//ttt->SetBranchAddress("bjet_n" ,&bjet_n);
		//ttt->SetBranchAddress("Z_pt" ,&Z_pt);
		//ttt->SetBranchAddress("HT" ,&HT);
		//ttt->SetBranchAddress("mll" ,&mll);
		//ttt->SetBranchAddress("METl" ,&METl);
		//for (int entry=0;entry<ttt->GetEntries();entry++) {
		//	ttt->GetEntry(entry);
		//	if (Z_pt<50.) continue;
		//	int pt = hist_low_pt->FindBin(Z_pt)-1;
		//	if (jet_n!=1) continue;
		//	if (bjet_n!=0) continue;
		//	z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//	if (mll<90 || mll>92) continue;
		//	z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//}
		//ftt.Close();
		//cout << "Opening file           : " << TString(smearingPath)+"vv/vv"+TString(ch)+".root"       << endl;
		//TFile fvv( TString(smearingPath)+"vv/vv"+TString(ch)+".root" );
		//TTree*  tvv              = (TTree*)fvv.Get("BaselineTree");
		//tvv->SetBranchStatus("*", 0);
		//tvv->SetBranchStatus("totalWeight", 1);
		//tvv->SetBranchStatus("jet_n", 1);
		//tvv->SetBranchStatus("bjet_n", 1);
		//tvv->SetBranchStatus("Z_pt", 1);
		//tvv->SetBranchStatus("HT", 1);
		//tvv->SetBranchStatus("mll", 1);
		//tvv->SetBranchStatus("METl", 1);
		//tvv->SetBranchAddress("totalWeight" ,&totalWeight);
		//tvv->SetBranchAddress("jet_n" ,&jet_n);
		//tvv->SetBranchAddress("bjet_n" ,&bjet_n);
		//tvv->SetBranchAddress("Z_pt" ,&Z_pt);
		//tvv->SetBranchAddress("HT" ,&HT);
		//tvv->SetBranchAddress("mll" ,&mll);
		//tvv->SetBranchAddress("METl" ,&METl);
		//for (int entry=0;entry<tvv->GetEntries();entry++) {
		//	tvv->GetEntry(entry);
		//	if (Z_pt<50.) continue;
		//	int pt = hist_low_pt->FindBin(Z_pt)-1;
		//	if (jet_n!=1) continue;
		//	if (bjet_n!=0) continue;
		//	z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//	if (mll<90 || mll>92) continue;
		//	z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//}
		//fvv.Close();
		TFile fPhoton( TString(smearingPath)+"gdata/gdata_raw.root" );
		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&METl);
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			g_metl[pt]->Fill(METl,totalWeight);
		}
		fPhoton.Close();
	}

}
