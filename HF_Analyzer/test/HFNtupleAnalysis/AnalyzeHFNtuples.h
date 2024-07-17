#ifndef ANALYZEHFNTUPLES_H
#define ANALYZEHFNTUPLES_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TF1.h"
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"


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
  TH2F *h_test;
  TH1F *h_HFEMClust_eLong3x3; 
  TH1F *h_HFEMClust_eShort3x3;
  TH1F *h_HFEMClust_eLong5x5; 
  TH1F *h_HFEMClust_eShort5x5;
  TH1F *h_test_long; 
  TH1F *h_test_short;
  TH1F *R;
  TH1F *angle_e_Zf;
};
#endif

#ifdef AnalyzeHFNtuples_cxx

void AnalyzeHFNtuples::BookHistogram(const char *outFileName) {

  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");

  //Define Histograms
  TH1F hnElectrons("nElectrons"          ,     "No. of electrons"  ,  100  ,   1000   , 0);
  h_nHFEMClust         = new TH1F("h_nHFEMClust"        ,     "No. of HF Clusters",                                       200 ,    1000  , 0     );
  hpt                  = new TH1F("pt"                  ,     "pt distribution"   ,                                       200 ,    1000  , 3     );
  heta                 = new TH1F("#eta"                ,     "#eta distribution"  ,                                      200 ,   -6     , 6     );
  hphi                 = new TH1F("#phi"                ,     "#phi distribution"  ,                                      200 ,   -4     , 4     );
  h_pt_e1              = new TH1F("pt_e1"               ,     "pt distribution of high pt e",                             200 ,    1000  , 0     );
  h_pt_e2              = new TH1F("pt_e2"               ,     "pt distribution of low pt e",                              200 ,    1000  , 0     );
  h_Z_mass             = new TH1F("Z_mass"              ,     "mass of reco Z boson"   ,                                  200 ,    0     , 400   );
  h_Z_pt               = new TH1F("Z_pt"                ,     "pt distribuiton of reco Z boson",                          200 ,    0     , 200   );
  h_Z_eta              = new TH1F("Z_#eta"              ,     "#eta distribution of reco Z boson",                        200 ,   -10    , 10    );
  h_Z_phi              = new TH1F("Z_#phi"              ,     "#phi distribution of reco Z boson",                        200 ,   -4     , 4     );
  h_Z_Energy           = new TH1F("Z_Energy"            ,     "Energy distribution of reco Z boson",                      200 ,    0     , 800   );
  h_Zmass_HFEMClust    = new TH1F("Z_mass_HFEMClust"    ,     "mass distribution of reco Z boson using HFEM Cluster",     200 ,    0     , 700   );
  h_Ele_Gen_Pt         = new TH1F("Ele_Gen_Pt"          ,     "pt distribution of gen level e",                           200 ,    0     , 1000  );                                                                                                                                     
  h_Ele_Gen_Phi        = new TH1F("Ele_Gen_Phi"         ,     "#phi distribution of gen level e",                         200 ,   -4     , 4     );               
  h_Ele_Gen_E          = new TH1F("Ele_Gen_E"           ,     "energy distribution of gen level e",                       200 ,    0     , 1000  );        
  h_Ele_Gen_Eta        = new TH1F("Ele_Gen_Eta"         ,     "#eta distribution of gen level e",                         200 ,   -10    , 10    );  
  h_Z_Rapi             = new TH1F("Z_Rapi"              ,     "Rapidity distribuiton of reco Z Boson ",                   200 ,   -10    , 10    );   
  Rapi_vs_eta          = new TH2F("Rapi_vs_eta"         ,     "Rapidity vs #eta distribuiton of reco Z Boson ",           300 ,   -10    , 10    , 300 ,   -10      , 10 ); 
  h_test               = new TH2F("P/E_vs_#eta/R"       ,     "P/E vs #eta/Rap. distribution of reco Z boson",            300 ,    0.8   , 1.2   , 300 ,    0.8     , 1.2);   
  h_HFEMClust_eLong3x3 = new TH1F("Long3x3"             ,     "HFEMClust_eLong3x3",                                       200 ,   1000   , 0     );   
  h_HFEMClust_eShort3x3= new TH1F("short3x3"            ,     "HFEMClust_eShort3x3",                                      200 ,   1000   , 0     ); 
  h_HFEMClust_eLong5x5 = new TH1F("long5x5"             ,     "HFEMClust_eLong5x5",                                       200 ,   1000   , 0     );
  h_HFEMClust_eShort5x5= new TH1F("short5x5"            ,     "HFEMClust_eShort5x5",                                      200 ,   1000   , 0     );
  h_test_long          = new TH1F("long"                ,     "long(3x3)/(5x5)",                                          300 ,   0      , 1.2   ); 
  h_test_short         = new TH1F("short"               ,     "short(3x3)/(5x5)",                                         300 ,   0      , 1.2   ); 
  R                    = new TH1F("R"                   ,     "R=(L-S)/(L+S)",                                            300 ,   -0.3   , 1.2   ); 
  angle_e_Zf           = new TH1F("angle_e_Zframe"      ,     "angle b/w two e-s in Z Frame",                             200 ,   0     ,  5     );
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
  hpt                   -> Write(); 
  heta                  -> Write();  
  hphi                  -> Write(); 
  h_nHFEMClust          -> Write();
  h_pt_e1               -> Write();
  h_pt_e2               -> Write();
  h_Z_mass              -> Write();
  h_Z_pt                -> Write();
  h_Z_eta               -> Write();
  h_Z_phi               -> Write();
  h_Z_Energy            -> Write();
  h_Zmass_HFEMClust     -> Write();
  h_Ele_Gen_Pt          -> Write();
  h_Ele_Gen_Phi         -> Write();
  h_Ele_Gen_E           -> Write();
  h_Ele_Gen_Eta         -> Write();
  h_Z_Rapi              -> Write();
  Rapi_vs_eta           -> Write();
  h_test                -> Write();
  h_HFEMClust_eLong3x3  -> Write();
  h_HFEMClust_eShort3x3 -> Write();
  h_HFEMClust_eLong5x5  -> Write();
  h_HFEMClust_eShort5x5 -> Write();
  h_test_long           -> Write();
  h_test_short          -> Write();
  R                     -> Write();
  angle_e_Zf            -> Write();



  h_Z_mass                -> GetXaxis()->SetTitle("mass of Z boson");
  h_Z_pt                  -> GetXaxis()->SetTitle("pt of Z boson");
  h_Z_eta                 -> GetXaxis()->SetTitle("#eta of Z boson");
  h_Z_phi                 -> GetXaxis()->SetTitle("#phi of Z boson");
  h_Z_Energy              -> GetXaxis()->SetTitle("energy of Z boson");
  h_Zmass_HFEMClust       -> GetXaxis()->SetTitle("mass of Z boson");
  h_Ele_Gen_Pt            -> GetXaxis()->SetTitle("pt of gen level electron");
  h_Ele_Gen_Phi           -> GetXaxis()->SetTitle("#phi of gen level electron");
  h_Ele_Gen_E             -> GetXaxis()->SetTitle("energy of gen level electron");
  h_Ele_Gen_Eta           -> GetXaxis()->SetTitle("#eta of gen level electron");
  h_Z_Rapi                -> GetXaxis()->SetTitle("rapidity of Z boson");
  Rapi_vs_eta             -> GetXaxis()->SetTitle("rapidity of Z boson");
  Rapi_vs_eta             -> GetYaxis()->SetTitle("#eta of Z boson");
  h_test                  -> GetXaxis()->SetTitle("#eta/rapidity of Z boson");
  h_test                  -> GetYaxis()->SetTitle("momentum/Energy of Z boson");
  h_HFEMClust_eLong3x3    -> GetXaxis()->SetTitle("Energy"); 
  h_HFEMClust_eShort3x3   -> GetXaxis()->SetTitle("Energy"); 
  h_HFEMClust_eLong5x5    -> GetXaxis()->SetTitle("Energy"); 
  h_HFEMClust_eShort5x5   -> GetXaxis()->SetTitle("Energy");  

  
  
  TCanvas *c1 = new TCanvas("c1 "  ,   "distribution of mass of Z for n=2 e",                        3600,4200);
  TCanvas *c2 = new TCanvas("c2 "  ,   "distribution of pt of Z",                                    3600,4200);
  TCanvas *c3 = new TCanvas("c3 "  ,   "distribution of eta of Z",                                   3600,4200);
  TCanvas *c4 = new TCanvas("c4 "  ,   "distribution of phi of Z",                                   3600,4200);
  TCanvas *c5 = new TCanvas("c5 "  ,   "distribution of energy of Z",                                3600,4200);
  TCanvas *c6 = new TCanvas("c6 "  ,   "distribution of mass of Z for n=1 e",                        3600,4200);
  TCanvas *c7 = new TCanvas("c7 "  ,   "distribution of pt gen level e",                             3600,4200);
  TCanvas *c8 = new TCanvas("c8 "  ,   "distribution of phi gen level e",                            3600,4200);
  TCanvas *c9 = new TCanvas("c9 "  ,   "distribution of energy gen level e",                         3600,4200);
  TCanvas *c10= new TCanvas("c10"  ,   "distribution of eta gen level e",                            3600,4200);
  TCanvas *c11= new TCanvas("c11"  ,   "distribution of rapidity of Z boson",                        3600,4200);
  TCanvas *c12= new TCanvas("c12"  ,   "distribution of rapidity vs eta of Z boson",                 3600,4200);
  TCanvas *c13= new TCanvas("c13"  ,   "distribution of #eta/rapidity vs momentum/Energy of Z boson",3600,4200);
  TCanvas *c14= new TCanvas("c14"  ,   "distribution of HFEMClust_eLong3x3",                         3600,4200);
  TCanvas *c15= new TCanvas("c15"  ,   "distribution of HFEMClust_eShort3x3",                        3600,4200);
  TCanvas *c16= new TCanvas("c16"  ,   "distribution of HFEMClust_eLong5x5",                         3600,4200);
  TCanvas *c17= new TCanvas("c17"  ,   "distribution of HFEMClust_eShort5x5",                        3600,4200);
  TCanvas *c18= new TCanvas("c18"  ,   "distribution of long(3x3)/(5x5)",                            3600,4200);
  TCanvas *c19= new TCanvas("c19"  ,   "distribution of short(3x3)/(5x5)",                           3600,4200);  
  TCanvas *c20= new TCanvas("c20"  ,   "distribution of R",                                          3600,4200);
  TCanvas *c21= new TCanvas("c21"  ,   "distribution of angle between two electrons",                3600,4200);

c12->SetLeftMargin(.15);
c12->SetRightMargin(0.15);
c13->SetRightMargin(0.15);
c13->SetLeftMargin(0.15);

// Breit Wigner Fitting
 
  double massMIN1=73;
  double massMAX1=130;
  //TF1 *funct = new TF1("mybw2",mybw,massMIN, massMAX,3);
  TF1 *funct1 = new TF1("mybw1","[0] / (TMath::Pi() * 2) * [2] / ((x - [1])*(x - [1]) + ([2]*[2]/4))", massMIN1, massMAX1);
  funct1->SetParameter(0,22000.0);        funct1->SetParName(0,"const");
  funct1->SetParameter(2,2.4952);         funct1->SetParName(1,"sigma");
  funct1->SetParameter(1,100.0);         funct1->SetParName(2,"mean");


 double massMIN2=80;
 double massMAX2=400;
TF1 *funct2 = new TF1("mybw2","[0] / (TMath::Pi() * 2) * [2] / ((x - [1])*(x - [1]) + ([2]*[2]/4))", massMIN2, massMAX2);
  funct2->SetParameter(0,50000.0);        funct2->SetParName(0,"const");
  funct2->SetParameter(2,5);         funct2->SetParName(1,"sigma");
  funct2->SetParameter(1,92.0);         funct2->SetParName(2,"mean");

// Fitting of Zmass and Zhfem mass 
  h_Zmass_HFEMClust ->Fit("mybw1","QR");
  TF1 *fit1 = h_Zmass_HFEMClust ->GetFunction("mybw1");
  
  h_Z_mass ->Fit("mybw2","QR");
  TF1 *fit2 = h_Z_mass ->GetFunction("mybw2");

  c1->cd();
  h_Z_mass                ->Draw();
  c2->cd();                    
  h_Z_pt                  ->Draw();
  c3->cd();                          
  h_Z_eta                 ->Draw();
  c4->cd();                              
  h_Z_phi                 ->Draw();
  c5->cd();                                    
  h_Z_Energy              ->Draw();
  c6->cd();                                      
  h_Zmass_HFEMClust       ->Draw();
  c7->cd();                                        
  h_Ele_Gen_Pt            ->Draw();
  c8->cd();                                            
  h_Ele_Gen_Phi           ->Draw();
  c9->cd();                                            
  h_Ele_Gen_E             ->Draw();
  c10->cd();                                           
  h_Ele_Gen_Eta           ->Draw();  
  c11->cd();                                           
  h_Z_Rapi                ->Draw();
  c12->cd();                                     
  gStyle->SetPalette();                                  
  Rapi_vs_eta             ->Draw("colz");
  c13->cd();                              
  gStyle->SetPalette(1);                        
  h_test                  ->Draw("colz");
  c14->cd();                              
  h_HFEMClust_eLong3x3    ->Draw();
  c15->cd();
  h_HFEMClust_eShort3x3   ->Draw();
  c16->cd();
  h_HFEMClust_eLong5x5    ->Draw();
  c17->cd();
  h_HFEMClust_eShort5x5   ->Draw();  
  c18->cd();
  h_test_long             ->Draw();
  c19->cd();
  h_test_short            ->Draw();
  c20->cd();
  R                       ->Draw();   
  c21->cd();
  angle_e_Zf              ->Draw();      

c6->cd();
c6->Update(); 
TLegend *legend2 = new TLegend(0.6, 0.6, 0.9, 0.8); // Set the position (x1, y1, x2, y2)
legend2->SetHeader("Legend", "C");
legend2->AddEntry(h_Zmass_HFEMClust, "Data", "l"); 
legend2->AddEntry(fit1, "Breit-Wigner", "l"); 
double peak2    = fit1->GetParameter(1);
double std_dev2 = fit1->GetParameter(2);
double chi2_2   = fit1->GetChisquare();
int ndf2        = fit1->GetNDF();
legend2->AddEntry(fit1, Form("Peak = %.6f", peak2), "");
legend2->AddEntry((TObject*)0, Form("Standard Deviation = %.6f", std_dev2), "");
legend2->AddEntry((TObject*)0, Form("#chi_2/ndf = %.2f/%d", chi2_2, ndf2), "");
legend2->Draw(); 

TPaveStats *stats5= (TPaveStats*)h_Zmass_HFEMClust->FindObject("stats");
stats5->SetX1NDC(0.6);   // Set X1 position (left edge)
stats5->SetY1NDC(0.8);     // Set Y1 position (bottom edge)
stats5->SetX2NDC(0.9);    // Set X2 position (right edge)
stats5->SetY2NDC(0.9);     // Set Y2 position (top edge)
c6->Modified();

 
c1->cd();
c1->Update(); 
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.8); // Set the position (x1, y1, x2, y2)
legend->SetHeader("Legend", "C");
legend->AddEntry(h_Z_mass, "Data", "l"); 
legend->AddEntry(fit2, "Breit-Wigner", "l"); 
double peak    = fit2->GetParameter(1);
double std_dev = fit2->GetParameter(2);
double chi2    = fit2->GetChisquare();
int ndf        = fit2->GetNDF();
legend->AddEntry(fit2, Form("Peak = %.6f", peak), "");
legend->AddEntry((TObject*)0, Form("Standard Deviation = %.6f", std_dev), "");
legend->AddEntry((TObject*)0, Form("#chi_2/ndf = %.2f/%d", chi2, ndf), "");
legend->Draw(); 

TPaveStats *stats4= (TPaveStats*)h_Z_mass->FindObject("stats");
stats4->SetX1NDC(0.6);   // Set X1 position (left edge)
stats4->SetY1NDC(0.8);     // Set Y1 position (bottom edge)
stats4->SetX2NDC(0.9);    // Set X2 position (right edge)
stats4->SetY2NDC(0.9);     // Set Y2 position (top edge)
c1->Modified();
  

//setting the position and size of stat box
c4->cd();
c4->Update();
TPaveStats *stats= (TPaveStats*)h_Z_phi->FindObject("stats");
stats->SetX1NDC(0.7);   // Set X1 position (left edge)
stats->SetY1NDC(0.83);     // Set Y1 position (bottom edge)
stats->SetX2NDC(0.9);    // Set X2 position (right edge)
stats->SetY2NDC(0.9);     // Set Y2 position (top edge)
c4->Modified();

c12->cd();
c12->Update();
TPaveStats *stats2= (TPaveStats*)Rapi_vs_eta->FindObject("stats");
stats2->SetX1NDC(0.711);   // Set X1 position (left edge)
stats2->SetY1NDC(0.7);     // Set Y1 position (bottom edge)
stats2->SetX2NDC(0.85);    // Set X2 position (right edge)
stats2->SetY2NDC(0.9);     // Set Y2 position (top edge)
c12->Modified();

c13->cd();
c13->Update();
TPaveStats *stats3= (TPaveStats*)h_test->FindObject("stats");
stats3->SetX1NDC(0.711);     // Set X1 position (left edge)
stats3->SetY1NDC(0.7);     // Set Y1 position (bottom edge)
stats3->SetX2NDC(0.85);     // Set X2 position (right edge)
stats3->SetY2NDC(0.9);     // Set Y2 position (top edge)
c12->Modified();

/*
c13->cd();
c13->Update();
TPaveStats *stats2= (TPaveStats*)h_test->FindObject("stats");
stats2->SetX1NDC(0.2);     // Set X1 position (left edge)
stats2->SetY1NDC(0.7);     // Set Y1 position (bottom edge)
stats2->SetX2NDC(0.4);     // Set X2 position (right edge)
stats2->SetY2NDC(0.9);     // Set Y2 position (top edge)
c12->Modified();
*/

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
c13-> SaveAs("test.pdf");   
c14-> SaveAs("HFEMClust_Long3x3.pdf"); 
c15-> SaveAs("HFEMClust_Short3x3.pdf"); 
c16-> SaveAs("HFEMClust_Long5x5.pdf"); 
c17-> SaveAs("HFEMClust_Short5x5.pdf"); 
c18-> SaveAs("Test_Long.pdf"); 
c19-> SaveAs("Test_Short.pdf"); 
c20-> SaveAs("R.pdf");
c21-> SaveAs("angle_e_Zframe.pdf");


c2->cd();
h_Z_pt->GetYaxis()->SetTitle("log(Y)");
h_Z_pt->Draw();
c2->SetLogy();
c2->SaveAs("log_pt.pdf");
  
c5->cd();
h_Z_Energy->GetYaxis()->SetTitle("log(Y)");
h_Z_Energy ->Draw();
c5->SetLogy();
c5->SetTitle("energy distribution in log scale of Z boson");
c5->SaveAs("log_energy.pdf");


  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif
