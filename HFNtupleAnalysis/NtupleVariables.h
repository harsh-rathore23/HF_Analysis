//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr  6 05:44:41 2015 by ROOT version 5.34/25
// from TTree tree/
// found on file: 239754_April06.root
//////////////////////////////////////////////////////////

#ifndef NtupleVariables_h
#define NtupleVariables_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"


// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleVariables : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
   vector<float>   *iEtaEle1;
   vector<float>   *iPhiEle1;
   vector<float>   *Hit_ES_Eta_Ele1;
   vector<float>   *Hit_ES_Phi_Ele1;
   vector<float>   *Hit_ES_X_Ele1;
   vector<float>   *Hit_ES_Y_Ele1;
   vector<float>   *Hit_ES_Z_Ele1;
   vector<float>   *ES_RecHitEnEle1;
   vector<float>   *Hit_Eta_Ele1;
   vector<float>   *Hit_Phi_Ele1;
   vector<float>   *Hit_X_Ele1;
   vector<float>   *Hit_Y_Ele1;
   vector<float>   *Hit_Z_Ele1;
   vector<float>   *RecHitEnEle1;
   vector<float>   *RecHitFracEle1;
   vector<float>   *iEtaEle2;
   vector<float>   *iPhiEle2;
   vector<float>   *Hit_ES_Eta_Ele2;
   vector<float>   *Hit_ES_Phi_Ele2;
   vector<float>   *Hit_ES_X_Ele2;
   vector<float>   *Hit_ES_Y_Ele2;
   vector<float>   *Hit_ES_Z_Ele2;
   vector<float>   *ES_RecHitEnEle2;
   vector<float>   *Hit_Eta_Ele2;
   vector<float>   *Hit_Phi_Ele2;
   vector<float>   *Hit_X_Ele2;
   vector<float>   *Hit_Y_Ele2;
   vector<float>   *Hit_Z_Ele2;
   vector<float>   *RecHitEnEle2;
   vector<float>   *RecHitFracEle2;
   Int_t           nElectrons;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *energy;
   vector<float>   *energy_ecal;
   vector<float>   *energy_ecal_mustache;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<int>     *passMVAMediumId;
   vector<float>   *Ele_R9;
   vector<float>   *Ele_S4;
   vector<float>   *Ele_SigIEIE;
   vector<float>   *Ele_SigIPhiIPhi;
   vector<float>   *Ele_SCEtaW;
   vector<float>   *Ele_SCPhiW;
   vector<float>   *Ele_CovIEtaIEta;
   vector<float>   *Ele_CovIEtaIPhi;
   vector<float>   *Ele_ESSigRR;
   vector<float>   *Ele_SCRawE;
   vector<float>   *Ele_SC_ESEnByRawE;
   vector<float>   *Ele_HadOverEm;
   vector<float>   *HFEMClust_pt;
   vector<float>   *HFEMClust_eta;
   vector<float>   *HFEMClust_phi;
   vector<float>   *HFEMClust_energy;
   vector<float>   *HFEMClust_eLong1x1;
   vector<float>   *HFEMClust_eShort1x1;
   vector<float>   *HFEMClust_eLong3x3;
   vector<float>   *HFEMClust_eShort3x3;
   vector<float>   *HFEMClust_eLong5x5;
   vector<float>   *HFEMClust_eShort5x5;
   vector<float>   *HFEMClust_e1x1;
   vector<float>   *HFEMClust_e3x3;
   vector<float>   *HFEMClust_e5x5;
   vector<float>   *HFEMClust_eSeL;
   vector<float>   *HFEMClust_eCOREe9;
   vector<float>   *HFEMClust_e9e25;
   vector<float>   *HFEMClust_eCore;
   vector<float>   *Ele_Gen_Pt;
   vector<float>   *Ele_Gen_Eta;
   vector<float>   *Ele_Gen_Phi;
   vector<float>   *Ele_Gen_E;
   Float_t         rho;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;

   // List of branches
   TBranch        *b_iEtaEle1;   //!
   TBranch        *b_iPhiEle1;   //!
   TBranch        *b_Hit_ES_Eta_Ele1;   //!
   TBranch        *b_Hit_ES_Phi_Ele1;   //!
   TBranch        *b_Hit_ES_X_Ele1;   //!
   TBranch        *b_Hit_ES_Y_Ele1;   //!
   TBranch        *b_Hit_ES_Z_Ele1;   //!
   TBranch        *b_ES_RecHitEnEle1;   //!
   TBranch        *b_Hit_Eta_Ele1;   //!
   TBranch        *b_Hit_Phi_Ele1;   //!
   TBranch        *b_Hit_X_Ele1;   //!
   TBranch        *b_Hit_Y_Ele1;   //!
   TBranch        *b_Hit_Z_Ele1;   //!
   TBranch        *b_RecHitEnEle1;   //!
   TBranch        *b_RecHitFracEle1;   //!
   TBranch        *b_iEtaEle2;   //!
   TBranch        *b_iPhiEle2;   //!
   TBranch        *b_Hit_ES_Eta_Ele2;   //!
   TBranch        *b_Hit_ES_Phi_Ele2;   //!
   TBranch        *b_Hit_ES_X_Ele2;   //!
   TBranch        *b_Hit_ES_Y_Ele2;   //!
   TBranch        *b_Hit_ES_Z_Ele2;   //!
   TBranch        *b_ES_RecHitEnEle2;   //!
   TBranch        *b_Hit_Eta_Ele2;   //!
   TBranch        *b_Hit_Phi_Ele2;   //!
   TBranch        *b_Hit_X_Ele2;   //!
   TBranch        *b_Hit_Y_Ele2;   //!
   TBranch        *b_Hit_Z_Ele2;   //!
   TBranch        *b_RecHitEnEle2;   //!
   TBranch        *b_RecHitFracEle2;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_energy_ecal;   //!
   TBranch        *b_energy_ecal_mustache;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
   TBranch        *b_passMVAMediumId;   //!
   TBranch        *b_Ele_R9;   //!
   TBranch        *b_Ele_S4;   //!
   TBranch        *b_Ele_SigIEIE;   //!
   TBranch        *b_Ele_SigIPhiIPhi;   //!
   TBranch        *b_Ele_SCEtaW;   //!
   TBranch        *b_Ele_SCPhiW;   //!
   TBranch        *b_Ele_CovIEtaIEta;   //!
   TBranch        *b_Ele_CovIEtaIPhi;   //!
   TBranch        *b_Ele_ESSigRR;   //!
   TBranch        *b_Ele_SCRawE;   //!
   TBranch        *b_Ele_SC_ESEnByRawE;   //!
   TBranch        *b_Ele_HadOverEm;   //!
   TBranch        *b_HFEMClust_pt;   //!
   TBranch        *b_HFEMClust_eta;   //!
   TBranch        *b_HFEMClust_phi;   //!
   TBranch        *b_HFEMClust_energy;   //!
   TBranch        *b_HFEMClust_eLong1x1;   //!
   TBranch        *b_HFEMClust_eShort1x1;   //!
   TBranch        *b_HFEMClust_eLong3x3;   //!
   TBranch        *b_HFEMClust_eShort3x3;   //!
   TBranch        *b_HFEMClust_eLong5x5;   //!
   TBranch        *b_HFEMClust_eShort5x5;   //!
   TBranch        *b_HFEMClust_e1x1;   //!
   TBranch        *b_HFEMClust_e3x3;   //!
   TBranch        *b_HFEMClust_e5x5;   //!
   TBranch        *b_HFEMClust_eSeL;   //!
   TBranch        *b_HFEMClust_eCOREe9;   //!
   TBranch        *b_HFEMClust_e9e25;   //!
   TBranch        *b_HFEMClust_eCore;   //!
   TBranch        *b_Ele_Gen_Pt;   //!
   TBranch        *b_Ele_Gen_Eta;   //!
   TBranch        *b_Ele_Gen_Phi;   //!
   TBranch        *b_Ele_Gen_E;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!

   NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~NtupleVariables() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
  // virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //ClassDef(NtupleVariables,0);
};

#endif

#ifdef NtupleVariables_cxx
void NtupleVariables::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   iEtaEle1 = 0;
   iPhiEle1 = 0;
   Hit_ES_Eta_Ele1 = 0;
   Hit_ES_Phi_Ele1 = 0;
   Hit_ES_X_Ele1 = 0;
   Hit_ES_Y_Ele1 = 0;
   Hit_ES_Z_Ele1 = 0;
   ES_RecHitEnEle1 = 0;
   Hit_Eta_Ele1 = 0;
   Hit_Phi_Ele1 = 0;
   Hit_X_Ele1 = 0;
   Hit_Y_Ele1 = 0;
   Hit_Z_Ele1 = 0;
   RecHitEnEle1 = 0;
   RecHitFracEle1 = 0;
   iEtaEle2 = 0;
   iPhiEle2 = 0;
   Hit_ES_Eta_Ele2 = 0;
   Hit_ES_Phi_Ele2 = 0;
   Hit_ES_X_Ele2 = 0;
   Hit_ES_Y_Ele2 = 0;
   Hit_ES_Z_Ele2 = 0;
   ES_RecHitEnEle2 = 0;
   Hit_Eta_Ele2 = 0;
   Hit_Phi_Ele2 = 0;
   Hit_X_Ele2 = 0;
   Hit_Y_Ele2 = 0;
   Hit_Z_Ele2 = 0;
   RecHitEnEle2 = 0;
   RecHitFracEle2 = 0;
   pt = 0;
   eta = 0;
   phi = 0;
   energy = 0;
   energy_ecal = 0;
   energy_ecal_mustache = 0;
   passMediumId = 0;
   passTightId = 0;
   passMVAMediumId = 0;
   Ele_R9 = 0;
   Ele_S4 = 0;
   Ele_SigIEIE = 0;
   Ele_SigIPhiIPhi = 0;
   Ele_SCEtaW = 0;
   Ele_SCPhiW = 0;
   Ele_CovIEtaIEta = 0;
   Ele_CovIEtaIPhi = 0;
   Ele_ESSigRR = 0;
   Ele_SCRawE = 0;
   Ele_SC_ESEnByRawE = 0;
   Ele_HadOverEm = 0;
   HFEMClust_pt = 0;
   HFEMClust_eta = 0;
   HFEMClust_phi = 0;
   HFEMClust_energy = 0;
   HFEMClust_eLong1x1 = 0;
   HFEMClust_eShort1x1 = 0;
   HFEMClust_eLong3x3 = 0;
   HFEMClust_eShort3x3 = 0;
   HFEMClust_eLong5x5 = 0;
   HFEMClust_eShort5x5 = 0;
   HFEMClust_e1x1 = 0;
   HFEMClust_e3x3 = 0;
   HFEMClust_e5x5 = 0;
   HFEMClust_eSeL = 0;
   HFEMClust_eCOREe9 = 0;
   HFEMClust_e9e25 = 0;
   HFEMClust_eCore = 0;
   Ele_Gen_Pt = 0;
   Ele_Gen_Eta = 0;
   Ele_Gen_Phi = 0;
   Ele_Gen_E = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iEtaEle1", &iEtaEle1, &b_iEtaEle1);
   fChain->SetBranchAddress("iPhiEle1", &iPhiEle1, &b_iPhiEle1);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele1", &Hit_ES_Eta_Ele1, &b_Hit_ES_Eta_Ele1);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele1", &Hit_ES_Phi_Ele1, &b_Hit_ES_Phi_Ele1);
   fChain->SetBranchAddress("Hit_ES_X_Ele1", &Hit_ES_X_Ele1, &b_Hit_ES_X_Ele1);
   fChain->SetBranchAddress("Hit_ES_Y_Ele1", &Hit_ES_Y_Ele1, &b_Hit_ES_Y_Ele1);
   fChain->SetBranchAddress("Hit_ES_Z_Ele1", &Hit_ES_Z_Ele1, &b_Hit_ES_Z_Ele1);
   fChain->SetBranchAddress("ES_RecHitEnEle1", &ES_RecHitEnEle1, &b_ES_RecHitEnEle1);
   fChain->SetBranchAddress("Hit_Eta_Ele1", &Hit_Eta_Ele1, &b_Hit_Eta_Ele1);
   fChain->SetBranchAddress("Hit_Phi_Ele1", &Hit_Phi_Ele1, &b_Hit_Phi_Ele1);
   fChain->SetBranchAddress("Hit_X_Ele1", &Hit_X_Ele1, &b_Hit_X_Ele1);
   fChain->SetBranchAddress("Hit_Y_Ele1", &Hit_Y_Ele1, &b_Hit_Y_Ele1);
   fChain->SetBranchAddress("Hit_Z_Ele1", &Hit_Z_Ele1, &b_Hit_Z_Ele1);
   fChain->SetBranchAddress("RecHitEnEle1", &RecHitEnEle1, &b_RecHitEnEle1);
   fChain->SetBranchAddress("RecHitFracEle1", &RecHitFracEle1, &b_RecHitFracEle1);
   fChain->SetBranchAddress("iEtaEle2", &iEtaEle2, &b_iEtaEle2);
   fChain->SetBranchAddress("iPhiEle2", &iPhiEle2, &b_iPhiEle2);
   fChain->SetBranchAddress("Hit_ES_Eta_Ele2", &Hit_ES_Eta_Ele2, &b_Hit_ES_Eta_Ele2);
   fChain->SetBranchAddress("Hit_ES_Phi_Ele2", &Hit_ES_Phi_Ele2, &b_Hit_ES_Phi_Ele2);
   fChain->SetBranchAddress("Hit_ES_X_Ele2", &Hit_ES_X_Ele2, &b_Hit_ES_X_Ele2);
   fChain->SetBranchAddress("Hit_ES_Y_Ele2", &Hit_ES_Y_Ele2, &b_Hit_ES_Y_Ele2);
   fChain->SetBranchAddress("Hit_ES_Z_Ele2", &Hit_ES_Z_Ele2, &b_Hit_ES_Z_Ele2);
   fChain->SetBranchAddress("ES_RecHitEnEle2", &ES_RecHitEnEle2, &b_ES_RecHitEnEle2);
   fChain->SetBranchAddress("Hit_Eta_Ele2", &Hit_Eta_Ele2, &b_Hit_Eta_Ele2);
   fChain->SetBranchAddress("Hit_Phi_Ele2", &Hit_Phi_Ele2, &b_Hit_Phi_Ele2);
   fChain->SetBranchAddress("Hit_X_Ele2", &Hit_X_Ele2, &b_Hit_X_Ele2);
   fChain->SetBranchAddress("Hit_Y_Ele2", &Hit_Y_Ele2, &b_Hit_Y_Ele2);
   fChain->SetBranchAddress("Hit_Z_Ele2", &Hit_Z_Ele2, &b_Hit_Z_Ele2);
   fChain->SetBranchAddress("RecHitEnEle2", &RecHitEnEle2, &b_RecHitEnEle2);
   fChain->SetBranchAddress("RecHitFracEle2", &RecHitFracEle2, &b_RecHitFracEle2);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nEle);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("energy_ecal", &energy_ecal, &b_energy_ecal);
   fChain->SetBranchAddress("energy_ecal_mustache", &energy_ecal_mustache, &b_energy_ecal_mustache);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
   fChain->SetBranchAddress("passMVAMediumId", &passMVAMediumId, &b_passMVAMediumId);
   fChain->SetBranchAddress("Ele_R9", &Ele_R9, &b_Ele_R9);
   fChain->SetBranchAddress("Ele_S4", &Ele_S4, &b_Ele_S4);
   fChain->SetBranchAddress("Ele_SigIEIE", &Ele_SigIEIE, &b_Ele_SigIEIE);
   fChain->SetBranchAddress("Ele_SigIPhiIPhi", &Ele_SigIPhiIPhi, &b_Ele_SigIPhiIPhi);
   fChain->SetBranchAddress("Ele_SCEtaW", &Ele_SCEtaW, &b_Ele_SCEtaW);
   fChain->SetBranchAddress("Ele_SCPhiW", &Ele_SCPhiW, &b_Ele_SCPhiW);
   fChain->SetBranchAddress("Ele_CovIEtaIEta", &Ele_CovIEtaIEta, &b_Ele_CovIEtaIEta);
   fChain->SetBranchAddress("Ele_CovIEtaIPhi", &Ele_CovIEtaIPhi, &b_Ele_CovIEtaIPhi);
   fChain->SetBranchAddress("Ele_ESSigRR", &Ele_ESSigRR, &b_Ele_ESSigRR);
   fChain->SetBranchAddress("Ele_SCRawE", &Ele_SCRawE, &b_Ele_SCRawE);
   fChain->SetBranchAddress("Ele_SC_ESEnByRawE", &Ele_SC_ESEnByRawE, &b_Ele_SC_ESEnByRawE);
   fChain->SetBranchAddress("Ele_HadOverEm", &Ele_HadOverEm, &b_Ele_HadOverEm);
   fChain->SetBranchAddress("HFEMClust_pt", &HFEMClust_pt, &b_HFEMClust_pt);
   fChain->SetBranchAddress("HFEMClust_eta", &HFEMClust_eta, &b_HFEMClust_eta);
   fChain->SetBranchAddress("HFEMClust_phi", &HFEMClust_phi, &b_HFEMClust_phi);
   fChain->SetBranchAddress("HFEMClust_energy", &HFEMClust_energy, &b_HFEMClust_energy);
   fChain->SetBranchAddress("HFEMClust_eLong1x1", &HFEMClust_eLong1x1, &b_HFEMClust_eLong1x1);
   fChain->SetBranchAddress("HFEMClust_eShort1x1", &HFEMClust_eShort1x1, &b_HFEMClust_eShort1x1);
   fChain->SetBranchAddress("HFEMClust_eLong3x3", &HFEMClust_eLong3x3, &b_HFEMClust_eLong3x3);
   fChain->SetBranchAddress("HFEMClust_eShort3x3", &HFEMClust_eShort3x3, &b_HFEMClust_eShort3x3);
   fChain->SetBranchAddress("HFEMClust_eLong5x5", &HFEMClust_eLong5x5, &b_HFEMClust_eLong5x5);
   fChain->SetBranchAddress("HFEMClust_eShort5x5", &HFEMClust_eShort5x5, &b_HFEMClust_eShort5x5);
   fChain->SetBranchAddress("HFEMClust_e1x1", &HFEMClust_e1x1, &b_HFEMClust_e1x1);
   fChain->SetBranchAddress("HFEMClust_e3x3", &HFEMClust_e3x3, &b_HFEMClust_e3x3);
   fChain->SetBranchAddress("HFEMClust_e5x5", &HFEMClust_e5x5, &b_HFEMClust_e5x5);
   fChain->SetBranchAddress("HFEMClust_eSeL", &HFEMClust_eSeL, &b_HFEMClust_eSeL);
   fChain->SetBranchAddress("HFEMClust_eCOREe9", &HFEMClust_eCOREe9, &b_HFEMClust_eCOREe9);
   fChain->SetBranchAddress("HFEMClust_e9e25", &HFEMClust_e9e25, &b_HFEMClust_e9e25);
   fChain->SetBranchAddress("HFEMClust_eCore", &HFEMClust_eCore, &b_HFEMClust_eCore);
   fChain->SetBranchAddress("Ele_Gen_Pt", &Ele_Gen_Pt, &b_Ele_Gen_Pt);
   fChain->SetBranchAddress("Ele_Gen_Eta", &Ele_Gen_Eta, &b_Ele_Gen_Eta);
   fChain->SetBranchAddress("Ele_Gen_Phi", &Ele_Gen_Phi, &b_Ele_Gen_Phi);
   fChain->SetBranchAddress("Ele_Gen_E", &Ele_Gen_E, &b_Ele_Gen_E);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NtupleVariables_cxx
