import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(4000)

options = VarParsing.VarParsing('standard')

options.register('inputFile',
        "~/",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File containing a list of the EXACT location of the output file  (default = ~/)"
        )
 #print(" it is printing ")

options.parseArguments()
options.inputFile = 'root://eoscms//' + options.inputFile
print(options.inputFile)
#Give input file name here###################
#input_file1 ="root://cms-xrd-global.cern.ch///store/mc/Run3Summer22EEDRPremix/DYto2L-2Jes_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postE_v1-v4/70000/cb8d5615-9998-49c3-834e-5e23e1442d2f.root"
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#          'root://cms-xrd-global.cern.ch///store/data/Run2022E/EGamma/AOD/27Jun2023-v1/25330005/d173e4c1-07a5-49d5-b2c4-605990493520.root'
#	input_file1
#          'root://cms-xrd-global.cern.ch///store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/70000/cb8d5615-9998-49c3-834e-5e23e1442d2f.root'
        'root://cms-xrd-global.cern.ch///store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/70003/18f26095-098b-4b45-92e9-449fd821a44a.root',
#         'root://cms-xrd-global.cern.ch///store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/70002/da0c25f3-6401-4bf0-a992-0f5382538ef7.root',
#          'root://cms-xrd-global.cern.ch///store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/70001/7d00b6a1-301c-447e-b4f2-ed572e74a70a.root'
#	 
#        'root://cms-xrd-global.cern.ch///store/mc/RunIIWinter19PFCalibDR/DoubleElectron_FlatPt-1To300/AODSIM/2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/40001/C6611D6B-CD06-C346-9828-BFA9EACE213E.root'
#                options.inputFile
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_postEE_v1', '')


#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10', '')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)
'''
from CondCore.DBCommon.CondDBSetup_cfi import *

process.GlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                               globaltag = cms.string('106X_upgrade2018_realistic_v11_L1v1'),
                               toGet = cms.VPSet(

                                        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EBCorrection_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EBCorrection_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),

        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EBUncertainty_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EBUncertainty_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),


        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EECorrection_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EECorrection_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),


        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EEUncertainty_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EEUncertainty_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),

                                )
                              )
'''
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']

#Checking if the given file coresponds to real data or simulation
#if (input_file1.find("SIM")!=-1):
#	checksim=False
#else:
#	checksim=True

process.nTuplelize = cms.EDAnalyzer('ZEE_RecHit_NTuplizer',
        electrons = cms.InputTag("gedGsfElectrons"),
#	isMC=cms.bool(checksim)
	genParticles = cms.InputTag("genParticles"),
	isMC=cms.bool(True)
       
	)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string("nTuple_MC.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.nTuplelize)

