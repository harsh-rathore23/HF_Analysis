#ifndef ANALYZEHFNTUPLES_H
#define ANALYZEHFNTUPLES_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

class AnalyzeHFNtuples : public NtupleVariables{

 public:
  AnalyzeHFNtuples(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  //AnalyzeHFNtuples(const TString &, const char*, const char*); 
  ~AnalyzeHFNtuples();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *);
  void     BookHistogram(const char *);
  
  TFile *oFile;
  
  //Declared Histograms
  TH1F  hnElectrons;
  TH1F *hpt; 
  TH1F *heta;  
  TH1F *hphi; 
  TH1F *h_nHFEMClust;
  TH1F *h_pt_e1;
  TH1F *h_pt_e2;
  TH1F *h_Z_mass;
  TH1F *h_Z_pt;
  TH1F *h_Z_eta;
  TH1F *h_Z_phi;
  TH1F *h_Z_Energy;
  TH1F *h_Zmass_HFEMClust;
};
#endif

#ifdef AnalyzeHFNtuples_cxx

void AnalyzeHFNtuples::BookHistogram(const char *outFileName) {

  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");

  //Define Histograms
  TH1F hnElectrons("nElectrons"          ,     "No. of electrons"  ,                100  ,   1000   , 0);
  h_nHFEMClust         = new TH1F("h_nHFEMClust"        ,     "No. of HF Clusters",                                100  ,   1000   , 0);
  hpt                  = new TH1F("pt"                  ,     "pt distribution"   ,                                 1000 ,   1000   , 3);
  heta                 = new TH1F("eta"                 ,     "eta distribution"  ,                                 500  ,   -6     , 6);
  hphi                 = new TH1F("phi"                 ,     "phi distribution"  ,                                 500  ,   -4     , 4);
  h_pt_e1              = new TH1F("pt_e1"               ,     "pt distribution of high pt e",                       1000 ,   1000   , 0);
  h_pt_e2              = new TH1F("pt_e2"               ,     "pt distribution of low pt e",                        1000 ,   1000   , 0);
  h_Z_mass             = new TH1F("Z_mass"              ,     "mass of Z boson"   ,                                 1000 ,   0   , 1000);
  h_Z_pt               = new TH1F("Z_pt"                ,     "pt distribuiton of Z boson",                         1000 ,   1000   , 0);
  h_Z_eta              = new TH1F("Z_eta"               ,     "eta distribution of Z boson",                        500  ,   -6     , 6);
  h_Z_phi              = new TH1F("Z_phi"               ,     "phi distribution of Z boson",                        500  ,   -4     , 4);
  h_Z_Energy           = new TH1F("Z_Energy"            ,     "Energy distribution of Z boson",                     1000 ,   1000   , 0);
  h_Zmass_HFEMClust    = new TH1F("Z_mass_HFEMClust"    ,     "mass distribution of Z boson using HFEM Cluster",    1000 ,   0   , 1000);
}


AnalyzeHFNtuples::AnalyzeHFNtuples(const TString &inputFileList, const char *outFileName, const char* dataset) {

  TChain *tree = new TChain("nTuplelize/T");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  NtupleVariables::Init(tree);

  BookHistogram(outFileName);
  
}

Bool_t AnalyzeHFNtuples::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeHFNtuples::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

AnalyzeHFNtuples::~AnalyzeHFNtuples() { 
  
  hnElectrons.Write();
  hpt                 -> Write(); 
  heta                -> Write();  
  hphi                -> Write(); 
  h_nHFEMClust        -> Write();
  h_pt_e1             -> Write();
  h_pt_e2             -> Write();
  h_Z_mass            -> Write();
  h_Z_pt              -> Write();
  h_Z_eta             -> Write();
  h_Z_phi             -> Write();
  h_Z_Energy          -> Write();
  h_Zmass_HFEMClust   -> Write();


  h_Z_mass            -> GetXaxis()->SetTitle("mass of Z boson");
  h_Z_pt              -> GetXaxis()->SetTitle("pt of Z boson");
  h_Z_eta             -> GetXaxis()->SetTitle("eta of Z boson");
  h_Z_phi             -> GetXaxis()->SetTitle("phi of Z boson");
  h_Z_Energy          -> GetXaxis()->SetTitle("energy of Z boson");
  h_Zmass_HFEMClust   -> GetXaxis()->SetTitle("mass of Z boson");

  
  TCanvas *c1 = new TCanvas("c1","distribution of mass of Z for n=2 e",3600,4200);
  TCanvas *c2 = new TCanvas("c2","distribution of pt of Z",3600,4200);
  TCanvas *c3 = new TCanvas("c3","distribution of eta of Z",3600,4200);
  TCanvas *c4 = new TCanvas("c4","distribution of phi of Z",3600,4200);
  TCanvas *c5 = new TCanvas("c5","distribution of energy of Z",3600,4200);
  TCanvas *c6 = new TCanvas("c6","distribution of mass of Z for n=1 e",3600,4200);



  
  c1->cd();
  h_Z_mass ->Draw();
  c2->cd();
  h_Z_pt ->Draw();
  c3->cd();
  h_Z_eta ->Draw();
  c4->cd();
  h_Z_phi->Draw();
  c5->cd();
  h_Z_Energy ->Draw();
  c6->cd();
  h_Zmass_HFEMClust->Draw();

  

  

   c1 -> SaveAs("Zmass.png");
   c2 -> SaveAs("Zpt.png");
   c3 -> SaveAs("Zeta.png");
   c4 -> SaveAs("Zphi.png");
   c5 -> SaveAs("Zenergy.png");
   c6 -> SaveAs("ZmassHFEMclust.png");



  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif
