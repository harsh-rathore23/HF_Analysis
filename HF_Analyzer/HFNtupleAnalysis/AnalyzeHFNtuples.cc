#define AnalyzeHFNtuples_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHFNtuples.h"

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  AnalyzeHFNtuples zee(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  zee.EventLoop(data);

  return 0;
}

void AnalyzeHFNtuples::EventLoop(const char *data) {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    cout << "nentries " << nentries << endl;

    double evtWgt = 1.0;
    std::cout << "Dataset " << data << " evtWeight " << evtWgt << std::endl;

    Long64_t nbytes = 0, nb = 0;
    int decade = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        double progress = 10.0 * jentry / (1.0 * nentries);
        int k = int(progress);
        if (k > decade)
            cout << 10 * k << " %" << endl;
        decade = k;

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        hnElectrons.Fill(nElectrons);
        
        int nHFEMClust = HFEMClust_pt->size();
        h_nHFEMClust->Fill(nHFEMClust);

        for (size_t a = 0; a < pt->size(); a++) {
            hpt  -> Fill(pt->at(a));
            heta -> Fill(eta->at(a));
            hphi -> Fill(phi->at(a));
        }

     
        if (pt->size() > 0) {
            int pt1   =   pt->at(0);
            int pt2   =  (pt->size() > 1) ? pt->at(1) : 0;
            h_pt_e1  ->  Fill(pt1);
            h_pt_e2  ->  Fill(pt2);
        }
       
        if (nElectrons == 2) {
            float Zmz, Zpt, Zeta, Zphi, ZEz, R, pz, P ,Ratio1 ,Ratio2;
            TLorentzVector v1, v2, vz;
            v1.SetPtEtaPhiE(pt->at(0), eta->at(0), phi->at(0), energy->at(0));
            v2.SetPtEtaPhiE(pt->at(1), eta->at(1), phi->at(1), energy->at(1));
            vz            =   v1 + v2;
            Zmz           =   vz.M();
            Zpt           =   vz.Pt();
            Zeta          =   vz.Eta();
            Zphi          =   vz.Phi();
            ZEz           =   vz.E();
            pz            =   vz.Pz();
            P             =   vz.P();
            R             =   0.5*log((ZEz+pz)/(ZEz-pz));
            Ratio1        =   Zeta/R;
            Ratio2        =   P/ZEz;

            // boosting electrons in another frame
            float Theta;
            TVector3 b=vz.BoostVector();
            
            v1.Boost(-b);
            v2.Boost(-b);
            Theta=v1.Angle(v2.Vect());

            h_Z_mass     ->  Fill(Zmz);
            h_Z_pt       ->  Fill(Zpt);
            h_Z_eta      ->  Fill(Zeta);
            h_Z_phi      ->  Fill(Zphi);
            h_Z_Energy   ->  Fill(ZEz);
            h_Z_Rapi     ->  Fill(R);
            Rapi_vs_eta  ->  Fill(R,Zeta);
            h_test       ->  Fill(Ratio1,Ratio2);
            angle_e_Zf   ->  Fill(Theta);
        } 

        if (nElectrons == 1 && HFEMClust_pt->size() > 0) 
        {float mz;
         TLorentzVector v1, v2, vz;
         v1.SetPtEtaPhiE(pt->at(0), eta->at(0), phi->at(0), energy->at(0));
                for (int i = 0; i < HFEMClust_pt->size(); i++)
                {
                    v2.SetPtEtaPhiE(HFEMClust_pt->at(i), HFEMClust_eta->at(i), HFEMClust_phi->at(i), HFEMClust_energy->at(i));
                    vz = v1 + v2;
                    mz = vz.M();
                    h_Zmass_HFEMClust->Fill(mz);
                }
          }

        for (size_t a = 0; a < Ele_Gen_Pt->size(); a++) {
            h_Ele_Gen_Pt        ->  Fill(Ele_Gen_Pt->at(a));
            h_Ele_Gen_Phi       ->  Fill(Ele_Gen_Phi->at(a));
            h_Ele_Gen_E         ->  Fill(Ele_Gen_E->at(a));
            h_Ele_Gen_Eta       ->  Fill(Ele_Gen_Eta->at(a));
        }

          for(int i=0 ; i<HFEMClust_eLong3x3->size();i++)
          {
          h_HFEMClust_eLong3x3  ->  Fill(HFEMClust_eLong3x3->at(i));
          h_HFEMClust_eShort3x3 ->  Fill(HFEMClust_eShort3x3->at(i));  
          h_HFEMClust_eLong5x5  ->  Fill(HFEMClust_eLong5x5->at(i));  
          h_HFEMClust_eShort5x5 ->  Fill(HFEMClust_eShort5x5->at(i));  

          float w,q;
          w=(HFEMClust_eLong3x3->at(i))/(HFEMClust_eLong5x5->at(i));
          q=(HFEMClust_eShort3x3->at(i))/(HFEMClust_eShort5x5->at(i));
          
          h_test_long           ->  Fill(w);
          h_test_short          ->  Fill(q);
          R                     ->  Fill(((HFEMClust_eLong1x1->at(i))-(HFEMClust_eShort1x1->at(i)))/((HFEMClust_eLong1x1->at(i))+(HFEMClust_eShort1x1->at(i))));
          }
          
          for(int i=0 ; i<Z_Truth_e_Pt->size();i++)
          {     
                h_Z_Truth_e_Pt           -> Fill(Z_Truth_e_Pt->at(i));
                h_Z_Truth_e_Eta          -> Fill(Z_Truth_e_Eta->at(i));
                h_Z_Truth_e_Phi          -> Fill(Z_Truth_e_Phi->at(i));
                h_Z_Truth_e_E            -> Fill(Z_Truth_e_E->at(i));
          }
          for(int i=0 ; i<Z_Truth_all_Pt->size();i++)
         {      
                h_Z_Truth_all_Pt         -> Fill(Z_Truth_all_Pt->at(i));
                h_Z_Truth_all_Eta        -> Fill(Z_Truth_all_Eta->at(i));
                h_Z_Truth_all_Phi        -> Fill(Z_Truth_all_Phi->at(i));
                h_Z_Truth_all_E          -> Fill(Z_Truth_all_E->at(i));
         }
          for(int i=0 ; i<Z_Truth_no_pt_zero_Pt->size();i++)
          {
                h_Z_Truth_no_pt_zero_Pt  -> Fill(Z_Truth_no_pt_zero_Pt->at(i));
                h_Z_Truth_no_pt_zero_Eta -> Fill(Z_Truth_no_pt_zero_Eta->at(i));
                h_Z_Truth_no_pt_zero_Phi -> Fill(Z_Truth_no_pt_zero_Phi->at(i));
                h_Z_Truth_no_pt_zero_E   -> Fill(Z_Truth_no_pt_zero_E->at(i));
          }
    }
}