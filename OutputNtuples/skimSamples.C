//--------------------------------------------------------------------------------------------------------------------------------
// This script reads in an ntuple, applies some skimming selection, and outputs a new ntuple with events passing the skim
//--------------------------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

void skim(char* input_path, char* output_path, TCut cut, char* sample, char* treename);

void skimSamples(){

//  TCut mll_cut = "";
//  mll_cut += "Jpsi_mll<5";
//
//  TCut cut      = mll_cut;
//  std::cout << cut.Print() << std::endl;
//
//  char* input_path     = "bphys/";
//  char* output_path    = "Jpsi/";
//  char* treename = "BaselineTree";                                                         // name of tree
//  skim(input_path,output_path,cut,"bdata_ee_NoSmear",treename );
//  skim(input_path,output_path,cut,"bdata_mm_NoSmear",treename );
 
  TCut mll_cut = "";
  mll_cut += "Jpsi_mll>5";

  TCut cut      = mll_cut;
  std::cout << cut.Print() << std::endl;

  char* input_path     = "bphys/";
  char* output_path    = "Upsi/";
  char* treename = "BaselineTree";                                                         // name of tree
  skim(input_path,output_path,cut,"bdata_ee_NoSmear",treename );
  skim(input_path,output_path,cut,"bdata_mm_NoSmear",treename );
 

}



void skim(char* input_path, char* output_path, TCut cut, char* sample, char* treename){

  //--------------------------------------------------
  // path of input and output files
  //--------------------------------------------------

  char* infilename              = Form("%s/%s.root",input_path,sample);
  char* outfilename             = Form("%s/%s.root",output_path,sample);

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << endl << endl;
  cout << "Reading in : " << infilename     << endl;
  cout << "Writing to : " << outfilename    << " " << cut   << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
 
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain(treename);
  chain->Add(infilename);

  //--------------------------------------------------
  // output file and tree
  //--------------------------------------------------

  TFile *file_cut = TFile::Open(outfilename, "RECREATE");
  assert( file_cut != 0 );
  TTree* tree_cut = chain->CopyTree( cut );
  tree_cut->Write();
  file_cut->Close();

  //--------------------------------------------------
  // dummy check
  //--------------------------------------------------

  TChain *chin  = new TChain(treename);
  TChain *chout = new TChain(treename);

  chin->Add(infilename);
  chout->Add(outfilename);

  cout << "Infile  total  entries " << chin->GetEntries()           << endl;
  cout << "Infile  cut    entries " << chin->GetEntries(TCut(cut))  << endl;
  cout << "Outfile total  entries " << chout->GetEntries()          << endl;
  cout << "Outfile cut    entries " << chout->GetEntries(TCut(cut)) << endl;


}
