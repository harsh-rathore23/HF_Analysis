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
#include "TStyle.h"
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
  TH1F *h_Ele_Gen_Pt;
  TH1F *h_Ele_Gen_Eta;
  TH1F *h_Ele_Gen_Phi;
  TH1F *h_Ele_Gen_E;
  TH1F *h_Z_Rapi;
  TH2F *Rapi_vs_eta;
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
  h_nHFEMClust         = new TH1F("h_nHFEMClust"        ,     "No. of HF Clusters",                                      200 ,   1000   , 0     );
  hpt                  = new TH1F("pt"                  ,     "pt distribution"   ,                                      200 ,   1000   , 3     );
  heta                 = new TH1F("eta"                 ,     "eta distribution"  ,                                      200 ,   -6     , 6     );
  hphi                 = new TH1F("phi"                 ,     "phi distribution"  ,                                      200 ,   -4     , 4     );
  h_pt_e1              = new TH1F("pt_e1"               ,     "pt distribution of high pt e",                            200 ,   1000   , 0     );
  h_pt_e2              = new TH1F("pt_e2"               ,     "pt distribution of low pt e",                             200 ,   1000   , 0     );
  h_Z_mass             = new TH1F("Z_mass"              ,     "mass of reco Z boson"   ,                                 200 ,   0      , 400   );
  h_Z_pt               = new TH1F("Z_pt"                ,     "pt distribuiton of reco Z boson",                         200 ,   0      , 200   );
  h_Z_eta              = new TH1F("Z_eta"               ,     "eta distribution of reco Z boson",                        200 ,   -10    , 10    );
  h_Z_phi              = new TH1F("Z_phi"               ,     "phi distribution of reco Z boson",                        200 ,   -4     , 4     );
  h_Z_Energy           = new TH1F("Z_Energy"            ,     "Energy distribution of reco Z boson",                     200 ,   0      , 800   );
  h_Zmass_HFEMClust    = new TH1F("Z_mass_HFEMClust"    ,     "mass distribution of reco Z boson using HFEM Cluster",    200 ,   0      , 700   );
  h_Ele_Gen_Pt         = new TH1F("Ele_Gen_Pt"          ,     "pt distribution of gen level e",                          200 ,   0      , 1000  );                                                                                                                                     
  h_Ele_Gen_Phi        = new TH1F("Ele_Gen_Phi"         ,     "phi distribution of gen level e",                         200 ,   -4     , 4     );               
  h_Ele_Gen_E          = new TH1F("Ele_Gen_E"           ,     "energy distribution of gen level e",                      200 ,   0      , 1000  );        
  h_Ele_Gen_Eta        = new TH1F("Ele_Gen_Eta"         ,     "eta distribution of gen level e",                         200 ,   -10    , 10    );  
  h_Z_Rapi             = new TH1F("Z_Rapi"              ,     "Rapidity distribuiton of reco Z Boson ",                  200 ,   -10    , 10    );   
  Rapi_vs_eta          = new TH2F("Rapi_vs_eta"         ,     "Rapidity vs eta distribuiton of reco Z Boson ",           200 ,   -10    , 10    , 200 ,   -10    , 10);         
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
  h_Ele_Gen_Pt        -> Write();
  h_Ele_Gen_Phi       -> Write();
  h_Ele_Gen_E         -> Write();
  h_Ele_Gen_Eta       -> Write();
  h_Z_Rapi            -> Write();
  Rapi_vs_eta         -> Write();


  h_Z_mass            -> GetXaxis()->SetTitle("mass of Z boson");
  h_Z_pt              -> GetXaxis()->SetTitle("pt of Z boson");
  h_Z_eta             -> GetXaxis()->SetTitle("eta of Z boson");
  h_Z_phi             -> GetXaxis()->SetTitle("phi of Z boson");
  h_Z_Energy          -> GetXaxis()->SetTitle("energy of Z boson");
  h_Zmass_HFEMClust   -> GetXaxis()->SetTitle("mass of Z boson");
  h_Ele_Gen_Pt        -> GetXaxis()->SetTitle("pt of gen level electron");
  h_Ele_Gen_Phi       -> GetXaxis()->SetTitle("phi of gen level electron");
  h_Ele_Gen_E         -> GetXaxis()->SetTitle("energy of gen level electron");
  h_Ele_Gen_Eta       -> GetXaxis()->SetTitle("eta of gen level electron");
  h_Z_Rapi            -> GetXaxis()->SetTitle("rapidity of Z boson");
  Rapi_vs_eta         -> GetXaxis()->SetTitle("rapidity vs eta of Z boson");

  
  TCanvas *c1 = new TCanvas("c1","distribution of mass of Z for n=2 e",3600,4200);
  TCanvas *c2 = new TCanvas("c2","distribution of pt of Z",3600,4200);
  TCanvas *c3 = new TCanvas("c3","distribution of eta of Z",3600,4200);
  TCanvas *c4 = new TCanvas("c4","distribution of phi of Z",3600,4200);
  TCanvas *c5 = new TCanvas("c5","distribution of energy of Z",3600,4200);
  TCanvas *c6 = new TCanvas("c6","distribution of mass of Z for n=1 e",3600,4200);
  TCanvas *c7 = new TCanvas("c7","distribution of pt gen level e",3600,4200);
  TCanvas *c8 = new TCanvas("c8","distribution of phi gen level e",3600,4200);
  TCanvas *c9 = new TCanvas("c9","distribution of energy gen level e",3600,4200);
  TCanvas *c10= new TCanvas("c10","distribution of eta gen level e",3600,4200);
  TCanvas *c11= new TCanvas("c11","distribution of rapidity of Z boson",3600,4200);
  c11->Divide(2,1);

  TCanvas *c12= new TCanvas("c12","distribution of rapidity vs eta of Z boson",3600,4200);


  
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
  c7->cd();
  h_Ele_Gen_Pt->Draw();
  c8->cd();
  h_Ele_Gen_Phi->Draw();
  c9->cd();
  h_Ele_Gen_E->Draw();
  c10->cd();
  h_Ele_Gen_Eta->Draw();  
  c11->cd(1);
  h_Z_Rapi->Draw();
  c11->cd(2);
  h_Z_eta ->Draw();
  c12->cd();
  gStyle->SetPalette(1);
  Rapi_vs_eta->Draw("colz");
  c12->Update();

   c1 -> SaveAs("Zmass.pdf");
   c2 -> SaveAs("Zpt.pdf");
   c3 -> SaveAs("Zeta.pdf");
   c4 -> SaveAs("Zphi.pdf");
   c5 -> SaveAs("Zenergy.pdf");
   c6 -> SaveAs("ZmassHFEMclust.pdf");
   c7 -> SaveAs("Ele_Gen_Pt.pdf");                     
   c8 -> SaveAs("Ele_Gen_Phi.pdf");                     
   c9 -> SaveAs("Ele_Gen_E.pdf");                     
   c10-> SaveAs("Ele_Gen_Eta.pdf");    
   c11-> SaveAs("Z_Rapidity.pdf");     
   c12-> SaveAs("Rapi_vs_eta.pdf");            
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif
