#include <stdio.h>
#include <stdlib.h>
#include <string> 

void drawSmearingHistogram() {

	TString index = "10";

	int bin = index.Atoi();
	TFile *fZ = new TFile("zjets/zjets_mm.root");
	TFile *fG = new TFile("gmc/gmc_mm_McSmear.root");
	TH1 *h1 = fG->Get(TString("z_metl_")+index);
	//TH1 *h1 = fG->Get(TString("z_metl_2j_")+index);
	TH1 *h2 = fG->Get(TString("g_metl_")+index);
	TH1 *h3 = fG->Get(TString("g_metl_smear_")+index);
	//TH1 *h3 = fG->Get(TString("g_metl_smear_2j_")+index);
	TH1 *h4 = fG->Get(TString("smear_final_")+index);
	h1->SetStats(0);
	h2->SetStats(0);
	h3->SetStats(0);
	h4->SetStats(0);

	const int bin_size = 22;
	int pt_bin[bin_size+1] ={50, 75,100,125,150 ,175 , 200,250 ,300, 400 ,500 ,700 ,1000,1200,1400,1600,1e10,1e10,1e10,1e10,1e10,1e10,1e10};

	TCanvas *c_both = new TCanvas("c_both","c both", 200, 10, 600, 600);
	c_both->cd();
	int rebin = 1;
	//if (bin<12) {
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	h3->Rebin(rebin);
	//}
	h1->SetLineColor(4);
	h2->SetLineColor(3);
	h3->SetLineColor(2);
	h4->SetLineColor(6);
	h1->GetXaxis()->SetRangeUser(-8*h1->GetRMS(),8*h1->GetRMS());
	h1->SetMaximum(100*h1->GetMaximum());
	h1->GetXaxis()->SetTitle("MET parallel (GeV)");
	h1->DrawNormalized("hist e");
	h2->DrawNormalized("hist e same");
	h3->DrawNormalized("hist e same");
	h4->DrawNormalized("hist same");
	std::cout << "Z RMS = " << h1->GetRMS() << std::endl;
	std::cout << "G RMS = " << h2->GetRMS() << std::endl;

	TLegend *legend = new TLegend(0.55,0.65,0.90,0.90);
	legend->SetBorderSize(0);
	legend->SetTextFont(42);
	legend->SetFillColor(0);
	legend->SetFillStyle(0);
	legend->SetLineColor(0);
	legend->AddEntry(h1,"Z+jets","f");
	legend->AddEntry(h2,"#gamma+jets","f");
	legend->AddEntry(h3,"#gamma+jets smeared","f");
	legend->AddEntry(h4,"smearing function","f");
	legend->Draw("same");

	char * bufPt1 = new char[32];
	char * bufPt2 = new char[32];
	sprintf(bufPt1,"%d",pt_bin[bin]);
	sprintf(bufPt2,"%d",pt_bin[bin+1]);
	TLatex *lumilab = new TLatex(0.2,0.84,TString(bufPt1)+"< Z pT <"+TString(bufPt2)+" GeV" );
	lumilab->SetNDC();
	lumilab->SetTextSize(0.038);
	lumilab->Draw();

	c_both->SetLogy();
	//c_both->SaveAs(TString("smearing_")+index+TString(".pdf"));
}
