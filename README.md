# HF_Analysis
HF electron analysis

Originally from Radjdeep
cmsrel CMSSW_13_0_7
cd CMSSW_13_0_7/src
cmsenv
git cms-init
git clone git@github.com:seemasharmafnal/HF_Analysis.git

================================
For making ntuples
cd HF_Analyzer
scram b -j32
cd test
cmsRun ZEE_RecHit_AOD_cfg.py 

Dataset to be used
/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4/AODSIM

=================================
For analyzing ntuples
cd HFNtupleAnalysis
make
./analyzeHFNtuples 
Please give 3 arguments runList  outputFileName dataset
./analyzeHFNtuples filelistRChatter-short.txt out.root HFe

(Jun 2024) Current list of files from Rajdeep are here
/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/rchatter/HF_Electron/MC/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_20240111_184939/240111_174946/0000/
