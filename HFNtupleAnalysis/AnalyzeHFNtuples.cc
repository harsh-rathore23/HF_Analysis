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

        // Check sizes before accessing elements
        int nHFEMClust = HFEMClust_pt->size();
        h_nHFEMClust->Fill(nHFEMClust);

        for (size_t a = 0; a < pt->size(); a++) {
            hpt->Fill(pt->at(a));
            heta->Fill(eta->at(a));
            hphi->Fill(phi->at(a));
        }

        if (pt->size() > 0) {
            int pt1 = pt->at(0);
            int pt2 = (pt->size() > 1) ? pt->at(1) : 0;
            h_pt_e1->Fill(pt1);
            h_pt_e2->Fill(pt2);
        }

        if (nElectrons == 2 && pt->size() > 1 && eta->size() > 1 && phi->size() > 1 && energy->size() > 1) {
            float emz, ept, eeta, ephi, eEz;
            TLorentzVector v1, v2, vz;
            v1.SetPtEtaPhiE(pt->at(0), eta->at(0), phi->at(0), energy->at(0));
            v2.SetPtEtaPhiE(pt->at(1), eta->at(1), phi->at(1), energy->at(1));
            vz = v1 + v2;
            emz = vz.M();
            ept = vz.Pt();
            eeta = vz.Eta();
            ephi = vz.Phi();
            eEz = vz.E();

            h_Z_mass->Fill(emz);
            h_Z_pt->Fill(ept);
            h_Z_eta->Fill(eeta);
            h_Z_phi->Fill(ephi);
            h_Z_Energy->Fill(eEz);
        } 

        if (nElectrons == 1 && HFEMClust_pt->size() > 0) {
            if (pt->size() > 0 && eta->size() > 0 && phi->size() > 0 && energy->size() > 0 &&
                HFEMClust_eta->size() == HFEMClust_pt->size() &&
                HFEMClust_phi->size() == HFEMClust_pt->size() &&
                HFEMClust_energy->size() == HFEMClust_pt->size()) {
                
                float mz;
                TLorentzVector v1, v2, vz;
                v1.SetPtEtaPhiE(pt->at(0), eta->at(0), phi->at(0), energy->at(0));

                for (int i = 0; i < HFEMClust_pt->size(); i++) {
                    v2.SetPtEtaPhiE(HFEMClust_pt->at(i), HFEMClust_eta->at(i), HFEMClust_phi->at(i), HFEMClust_energy->at(i));
                    vz = v1 + v2;
                    mz = vz.M();
                    h_Zmass_HFEMClust->Fill(mz);
                }
            } 
        }
    }
}
