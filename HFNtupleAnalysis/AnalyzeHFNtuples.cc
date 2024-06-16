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

  //  double weight = 1.0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;
    
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    
    int nHFEMClust = HFEMClust_pt->size();
    h_nHFEMClust->Fill(nHFEMClust);
      
   
  } // loop over ntuple entries
  
}

