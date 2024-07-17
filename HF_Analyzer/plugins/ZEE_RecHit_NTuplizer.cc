// -*- C++ -*-
//
// Package:    Electron_GNN_Regression/ZEE_RecHit_NTuplizer
// Class:      ZEE_RecHit_NTuplizer
//
/**\class ZEE_RecHit_NTuplizer ZEE_RecHit_NTuplizer.cc Electron_GNN_Regression/ZEE_RecHit_NTuplizer/plugins/ZEE_RecHit_NTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Fri, 21 Feb 2020 11:38:58 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"

//////////////////////Added by Harsh////////////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
////////////////////////////////////////////////////////////////////////

//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace std;
using namespace reco;

using reco::TrackCollection;

class ZEE_RecHit_NTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:

   
      explicit ZEE_RecHit_NTuplizer(const edm::ParameterSet&);
      ~ZEE_RecHit_NTuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     std::vector<float> Hit_ES_Eta[2];
     std::vector<float> Hit_ES_Phi[2];
     std::vector<float> Hit_ES_X[2];
     std::vector<float> Hit_ES_Y[2];
     std::vector<float> Hit_ES_Z[2];
     std::vector<float> ES_RecHitEn[2];


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::ofstream outFile;
//   cluster tools
      EcalClusterLazyTools *clustertools;
      noZS::EcalClusterLazyTools *clustertools_NoZS;

//   Identify if the SC lies in EB OR EE based on its seed
     bool isEB = 0;
     bool isEE = 0; // !isEB not sufficient since later will try to include the preshower as well

// Get the hits from the ES
//     std::vector<GlobalPoint> GetESPlaneRecHits(const reco::SuperCluster& sc, unsigned int planeIndex) const;
     void GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int elenum, unsigned int planeIndex);

//   clear the vectors 
     void ClearTreeVectors();

// ----------member data ---------------------------
     TTree* T;
     // Variables for Run info.
     int run;
     int event;
     int lumi;
     int n=0;
     // Electron variables
     int nElectrons_;
     Float_t rho;

     std::vector<float> iEta[2];
     std::vector<float> iPhi[2];
     std::vector<float> Hit_Eta[2];
     std::vector<float> Hit_Phi[2];
     std::vector<float> Hit_X[2];
     std::vector<float> Hit_Y[2];
     std::vector<float> Hit_Z[2];


     std::vector<float> RecHitFrac[2];
     std::vector<float> RecHitEn[2];

     std::vector<float> Ele_pt_;
     std::vector<float> Ele_eta_;
     std::vector<float> Ele_phi_;
     std::vector<float> Ele_energy_;
     std::vector<float> Ele_ecal_energy_;
     std::vector<float> Ele_ecal_mustache_energy_;

     std::vector<float> Ele_R9;
     std::vector<float> Ele_S4;
     std::vector<float> Ele_SigIEIE;
     std::vector<float> Ele_SigIPhiIPhi;
     std::vector<float> Ele_SCEtaW;
     std::vector<float> Ele_SCPhiW;
     std::vector<float> Ele_CovIEtaIEta;
     std::vector<float> Ele_CovIEtaIPhi;
     std::vector<float> Ele_ESSigRR;
     std::vector<float> Ele_SCRawE;
     std::vector<float> Ele_SC_ESEnByRawE;
     std::vector<float> Ele_HadOverEm;

     std::vector<float> HFEMClust_pt_;
     std::vector<float> HFEMClust_eta_;
     std::vector<float> HFEMClust_phi_;
     std::vector<float> HFEMClust_energy_;
    //energy in long or short fibers various cluster sizes
     std::vector<float> HFEMClust_eLong1x1_;
     std::vector<float> HFEMClust_eShort1x1_;
     std::vector<float> HFEMClust_eLong3x3_;
     std::vector<float> HFEMClust_eShort3x3_;
     std::vector<float> HFEMClust_eLong5x5_;
     std::vector<float> HFEMClust_eShort5x5_;
    //total energy in various clusters
     std::vector<float> HFEMClust_e1x1_;
     std::vector<float> HFEMClust_e3x3_;
     std::vector<float> HFEMClust_e5x5_;
    //Identification Variables
     std::vector<float> HFEMClust_eSeL_; //Longitudinal variable: E(3x3,short fibers)/E(3x3,long fibers)
     std::vector<float> HFEMClust_eCOREe9_; // Transverse Variable: E(Core of cluster)/E(3x3)
     std::vector<float> HFEMClust_e9e25_; // Shower Exclusion Variable: E(3x3)/E(5x5)
     std::vector<float> HFEMClust_eCore_; // energy in central highest energy cells (at least 50% energy of previous total energy startign with seed cell)


	
	

     std::vector<float> Ele_Gen_Pt;
     std::vector<float> Ele_Gen_Eta;
     std::vector<float> Ele_Gen_Phi;
     std::vector<float> Ele_Gen_E;


     std::vector<int> passLooseId_;
     std::vector<int> passMediumId_;
     std::vector<int> passTightId_;
     std::vector<int> passMVAMediumId_;

     std::vector<int> isTrue_;

// added by Harsh
////////Defining variables for property of Z Boson////////


//////////////////////////////////Z boson that produce e-////////////
     std::vector<float> Z_Truth_e_Pt;
     std::vector<float> Z_Truth_e_Eta;
     std::vector<float> Z_Truth_e_Phi;
     std::vector<float> Z_Truth_e_E;
//////////////////////////////////all the Z boson////////////
     std::vector<float> Z_Truth_all_Pt;      
     std::vector<float> Z_Truth_all_Eta;     
     std::vector<float> Z_Truth_all_Phi;       
     std::vector<float> Z_Truth_all_E;           
//////////////////////////////////all Z boson except pt=0//////
     std::vector<float> Z_Truth_no_pt_zero_Pt;            
     std::vector<float> Z_Truth_no_pt_zero_Eta;            
     std::vector<float> Z_Truth_no_pt_zero_Phi;       
     std::vector<float> Z_Truth_no_pt_zero_E;        



//////////////////////////////////////////////////////////

      // -----------------Handles--------------------------
      edm::Handle<EcalRecHitCollection> EBRechitsHandle;
      edm::Handle<EcalRecHitCollection> EERechitsHandle;
      edm::Handle<EcalRecHitCollection> ESRechitsHandle;
      edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
      edm::Handle<reco::SuperClusterCollection> HFSC;
      edm::Handle<reco::HFEMClusterShapeAssociationCollection> HFClusterCollection;
      edm::Handle<edm::View<reco::GsfElectron> > electrons;  


      ///////////////////////////////////////////////////////////////////////////////////////added by Harsh
      edm::Handle<GenEventInfoProduct> xhandle;
      edm::Handle<edm::View<reco::GenParticle> >gphandle;
      ////////////////////////////////////////////////////////////////////////////////////////////////



//      edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
//      edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
//      edm::Handle<edm::View<reco::GenParticle> > genParticles;
      //---------------- Input Tags-----------------------
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEBToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEEToken_;
      edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionESToken_;
      edm::EDGetTokenT<reco::RecoEcalCandidateCollection> HFElectronsToken_;
      edm::EDGetTokenT<reco::SuperClusterCollection> HFSCToken_;
      edm::EDGetTokenT<reco::HFEMClusterShapeAssociationCollection> HFClusterToken_;
      edm::EDGetToken electronsToken_;

      ////////////////added by Harsh///////////////////
      edm::EDGetTokenT<GenEventInfoProduct> xtoken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > gptoken_;
      edm::EDGetTokenT<reco::GenParticleCollection> gpctoken_;
      /////////////////////////////////////////////////


//      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
//      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
//      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//



ZEE_RecHit_NTuplizer::ZEE_RecHit_NTuplizer(const edm::ParameterSet& iConfig):
   recHitCollectionEBToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEB"))),
   recHitCollectionEEToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEE"))),
   recHitCollectionESToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsES"))),
   HFElectronsToken_(consumes<reco::RecoEcalCandidateCollection>(edm::InputTag("hfRecoEcalCandidate"))),
   HFSCToken_(consumes<reco::SuperClusterCollection>(edm::InputTag("hfEMClusters"))),
   HFClusterToken_(consumes<reco::HFEMClusterShapeAssociationCollection>(edm::InputTag("hfEMClusters"))),
   ///////added by Harsh/////
   xtoken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   gpctoken_(consumes<reco::GenParticleCollection>(edm::InputTag("genParticles")))
   //////////////////////////

   //now do what ever initialization is needed

   //added by Harsh
 {  gptoken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   ///////////////////////////////////////////


   electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"));
//   genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   usesResource("TFileService");

    outFile.open("zee_sem_events.txt");
    if (!outFile) {
        std::cerr << "Error opening file output.txt for writing!" << std::endl;
    }

 }


ZEE_RecHit_NTuplizer::~ZEE_RecHit_NTuplizer()
{
    if (outFile.is_open()) {
            outFile.close();}
   
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)




}


//
// member functions
//

// ------------ method called for each event  ------------


int number=0;
int p=0;
bool v=false;



void ZEE_RecHit_NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   iEvent.getByToken(recHitCollectionEBToken_, EBRechitsHandle);
   iEvent.getByToken(recHitCollectionEEToken_, EERechitsHandle);
   iEvent.getByToken(recHitCollectionESToken_, ESRechitsHandle);
   iEvent.getByToken(HFElectronsToken_, HFElectrons);
   iEvent.getByToken(HFSCToken_, HFSC);
   iEvent.getByToken(HFClusterToken_, HFClusterCollection);
   iEvent.getByToken(electronsToken_, electrons);



   ///////added by Harsh///////
   	//std::cout<<"This is simulated data. Finding the Generator level information"<<std::endl;
 
   	iEvent.getByToken(xtoken_, xhandle);
   	iEvent.getByToken(gptoken_, gphandle);
   ////////////////////////////



//Clear all vectors to be written to the tree
   ClearTreeVectors();
   run=0;  event=0;  lumi=0;

///////////////////////////Fill Electron/Photon related stuff/////////////////////////////////////////////////////
   nElectrons_ = 0;

   for (size_t i = 0; i < electrons->size(); ++i){
	if(electrons->size() > 2) break;
	const auto ele = electrons->ptrAt(i);
        if( ele->pt() < 1. ) continue;
	if(!(ele->ecalDrivenSeed())) continue;
	if(ele->parentSuperCluster().isNull()) continue;

        nElectrons_++;
        Ele_pt_.push_back( ele->pt() );
        Ele_eta_.push_back( ele->superCluster()->eta() );
        Ele_phi_.push_back( ele->superCluster()->phi() );
        Ele_energy_.push_back( ele->energy() );
	Ele_ecal_energy_.push_back( ele->correctedEcalEnergy() );
        Ele_R9.push_back(ele->full5x5_r9());
        Ele_SigIEIE.push_back(ele->full5x5_sigmaIetaIeta());
        Ele_SigIPhiIPhi.push_back(ele->full5x5_sigmaIphiIphi());
        Ele_SCEtaW.push_back(ele->superCluster()->etaWidth());
        Ele_SCPhiW.push_back(ele->superCluster()->phiWidth());
	Ele_HadOverEm.push_back(ele->hadronicOverEm());
        const CaloClusterPtr seed_clu = ele->superCluster()->seed();
	Ele_SCRawE.push_back(ele->superCluster()->rawEnergy());
        Ele_SC_ESEnByRawE.push_back( (ele->superCluster()->preshowerEnergy())/(ele->superCluster()->rawEnergy()) );

   }         

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////HFCluster variables/////////////////////////////////////////////

//   std::cout<<std::endl<<" HFClusters = "<<HFClusterCollection->size()<<" HFElectrons = "<<HFSC->size()<<std::endl;
   for (unsigned int iii = 0; iii < HFSC->size(); ++iii){
//	std::cout<<std::endl<<" Iteration number = "<<iii++<<std::endl;
	const SuperCluster& supClus = (*HFSC)[iii];	
	reco::SuperClusterRef theClusRef = edm::Ref<reco::SuperClusterCollection>(HFSC, iii);	
	const HFEMClusterShapeRef clusShapeRef = HFClusterCollection->find(theClusRef)->val;
	const HFEMClusterShape& clusShape=*clusShapeRef;

        HFEMClust_pt_.push_back(supClus.energy()/cosh(supClus.eta()) );
        HFEMClust_eta_.push_back(supClus.eta());
        HFEMClust_phi_.push_back(supClus.phi());
        HFEMClust_energy_.push_back(supClus.energy());
        HFEMClust_eLong1x1_.push_back(clusShape.eLong1x1());
        HFEMClust_eShort1x1_.push_back(clusShape.eShort1x1());
        HFEMClust_eLong3x3_.push_back(clusShape.eLong3x3());
        HFEMClust_eShort3x3_.push_back(clusShape.eShort3x3());
        HFEMClust_eLong5x5_.push_back(clusShape.eLong5x5());
        HFEMClust_eShort5x5_.push_back(clusShape.eShort5x5());
        HFEMClust_e1x1_.push_back(clusShape.e1x1());
        HFEMClust_e3x3_.push_back(clusShape.e3x3());
        HFEMClust_e5x5_.push_back(clusShape.e5x5());
        HFEMClust_eSeL_.push_back(clusShape.eSeL());
        HFEMClust_eCOREe9_.push_back(clusShape.eCOREe9());
        HFEMClust_e9e25_.push_back(clusShape.e9e25());
        HFEMClust_eCore_.push_back(clusShape.eCore());

   }


   //////////////////////////////////////Added by Harsh//////////////////////////////////////////////////////////////////////////////////
 //const reco::GenParticle* z_boson=nullptr;
 /*{outFile<<"S no."<<"\t"<<"Entry"
                <<"\t"<<"pdgId" 
                <<"\t"<<"status"              
                <<"\t"<<"pt"
                <<"\t\t"<<"eta"
                <<"\t\t"<<"phi"
                <<"\t\t"<<"energy"
                <<"\t\t"<<"Mother Id"
                <<"\t\t"<<"Daughter Id"
                <<endl;
//            <<"  event  :    "<<iEvent.id().event()<<endl;
                  }*/
                 
 
for(edm::View<GenParticle>::const_iterator part = gphandle->begin(); part!=gphandle->end(); ++part){
		if(abs(part ->pdgId())==11 && part ->mother()->pdgId()==23){
			Ele_Gen_Pt.push_back(part->pt());
			Ele_Gen_Eta.push_back(part->eta());
			Ele_Gen_Phi.push_back(part->phi());
			Ele_Gen_E.push_back(part->energy());
                 //       cout<<"Electron "<<"Pt: "<<part->pt()<<"    MotherID: "<<part->mother()->pdgId()<<"    Status:"<< part-> status()<<endl;
			}

               //////////////////////////////////////////////////////////////// filling of Z bosons that decay into e- and e+  /////////////////////////////////////////////////////
               if(part ->pdgId()==11 && part->mother()->pdgId()==23)
                {  Z_Truth_e_Pt.push_back(part->mother()->pt());
                   Z_Truth_e_Eta.push_back(part->mother()->eta());
                   Z_Truth_e_Phi.push_back(part->mother()->phi());
                   Z_Truth_e_E.push_back(part->mother()->energy());
                }
////////////////////////////////////////////////////////////////// Testing Z that decay to e+e- ////////////////////////////////////////////////////////////////////////////////////
                int nDZs = part->numberOfDaughters();
                if(part->pdgId()==23 && nDZs==2 && ((part->daughter(0)->pdgId()==11 && part->daughter(1)->pdgId()==-11)||(part->daughter(0)->pdgId()==-11 && part->daughter(1)->pdgId()==11))){ 
                         
        //        if((part->daughter(0)->pdgId()==11 && part->daughter(1)->pdgId()==-11)||(part->daughter(0)->pdgId()==-11 && part->daughter(1)->pdgId()==11))
                if(v==false && p<15)
                {outFile<<"S no."<<"\t"<<"Entry"
                       <<"\t"<<"pdgId" 
                       <<"\t"<<"status"              
                       <<"\t"<<"pt"
                       <<"\t\t"<<"eta"
                       <<"\t\t"<<"phi"
                       <<"\t\t"<<"energy"
                       <<"\t\t"<<"Mother Id"
                       <<"\t\t"<<"Daughter Id"
                       <<endl;
//                     <<"  event  :    "<<iEvent.id().event()<<endl;
                       v=true;
                }}
                 if(v==true){
                        outFile<<number
                        <<"\t"<<setprecision(2)<<p
                        <<"\t"<<setprecision(2)<<part->pdgId()
                        <<"\t"<<setprecision(2)<<part->status()
                        <<"\t"<<setprecision(2)<<part->pt()
                        <<"\t\t"<<setprecision(2)<<part->eta()
                        <<"\t\t"<<setprecision(2)<<part->phi()
                        <<"\t\t"<<setprecision(2)<<part->energy()<<"\t";
  //                    <<"  event  :    "<<iEvent.id().event()<<endl;
        
                        int nMothers = part->numberOfMothers(); 
                                for (int i = 0; i < nMothers; i++)
                                              {
                                                const auto* mother = part->mother(i);
                                                if (mother != nullptr) {
                                                outFile<<"\t"<<mother->pdgId();
                                                                        }
                                              }
                        outFile<<"\t\t";
        
                        int nDaughters = part->numberOfDaughters(); 
                                for (int i = 0; i < nDaughters; i++)
                                             {
                                              const auto* daughter = part->daughter(i); 
                                              if (daughter != nullptr) {
                                               outFile<<"\t"<<daughter->pdgId();
                                                                       }
                                              }
                        outFile<<endl;

                    number = number + 1;
         //           p      = p+1;
                    }
                
                
                ///////////////////////////////////////////////////////////////Testing the Z boson that decay into e+ and e-////////////////////////////////////////////////////////
   /*             if(part->pdgId()==23 && ((part->daughter(1)->pdgId()==11 && part->daughter(0)->pdgId()==-11)||(part->daughter(0)->pdgId()==11 && part->daughter(1)->pdgId()==-11)))
                { 
                int nMothers = part->numberOfMothers(); 
                for (int i = 0; i < nMothers; i++)
                              {
                                const auto* mother = part->mother(i);
                                if (mother != nullptr) {
                                cout << "mother " << i + 1 << " Id         :     "<<mother->pdgId()
                                     <<"     pt  :    "<<mother->pt() 
                                     <<"     eta  :    "<<mother->eta()
                                     <<"     Phi  :    "<<mother->phi()
                                     <<"     Energy  :    "<<mother->energy()
                                     <<"     Status  :    "<<mother->status() << endl;
                                                         }
                               }              
                    
                cout<<"Z boson Particle Id :     "<<part->pdgId()
                    <<"     pt  :    "<<part->pt()
                    <<"     eta  :    "<<part->eta()
                    <<"     Phi  :    "<<part->phi()
                    <<"     Energy  :    "<<part->energy()
                    <<"     Status  :    "<<part->status()<<endl;          
                                           
                int nDaughters = part->numberOfDaughters(); 
                for (int i = 0; i < nDaughters; i++)
                             {
                              const auto* daughter = part->daughter(i); 
                              if (daughter != nullptr) {
                              cout << "daughter " << i + 1 << " Id       :     " << daughter->pdgId()
                                   <<"     pt  :    "<<daughter->pt()
                                   <<"     eta  :    "<<daughter->eta()
                                   <<"     Phi  :    "<<daughter->phi()
                                   <<"     Energy  :    "<<daughter->energy()
                                   <<"     Status  :    "<<daughter->status() << endl;
                                                         }
                             }
                cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;      
                }                  
   */           
                ////////////////////////////////////////////////////////////////filling and testing of all the Z bosons///////////////////////////////////////////////////////////
                if(part->pdgId()==23){
                        Z_Truth_all_Pt.push_back(part->pt());
                        Z_Truth_all_Eta.push_back(part->eta());   
                        Z_Truth_all_Phi.push_back(part->phi());      
                        Z_Truth_all_E.push_back(part->energy());    
                        
         /*               if(part->pt()==0)
                        {int nMothers = part->numberOfMothers(); 
                        for (int i = 0; i < nMothers; i++)
                                      {
                                        const auto* mother = part->mother(i);
                                        if (mother != nullptr) {
                                        cout << "mother " << i + 1 << " Id         :     "<<mother->pdgId()
                                             <<"     pt  :    "<<mother->pt()
                                             <<"     eta  :    "<<mother->eta()
                                             <<"     Phi  :    "<<mother->phi()
                                             <<"     Energy  :    "<<mother->energy()
                                             <<"     Status  :    "<<mother->status() << endl;
                                                                 }
                                       }              

                        cout<<"Z boson Particle Id :     "<<part->pdgId()
                            <<"     pt  :    "<<part->pt()
                            <<"     eta  :    "<<part->eta()
                            <<"     Phi  :    "<<part->phi()
                            <<"     Energy  :    "<<part->energy()
                            <<"     Status  :    "<<part->status()<<endl;          
                        
                        int nDaughters = part->numberOfDaughters(); 
                        for (int i = 0; i < nDaughters; i++)
                                     {
                                      const auto* daughter = part->daughter(i); 
                                      if (daughter != nullptr) {
                                      cout << "daughter " << i + 1 << " Id       :     "<< daughter->pdgId()
                                           <<"     pt  :    "<<daughter->pt()
                                           <<"     eta  :    "<<daughter->eta()
                                           <<"     Phi  :    "<<daughter->phi()
                                           <<"     Energy  :    "<<daughter->energy()
                                           <<"     Status  :    "<<daughter->status() << endl;
                                                                 }
                                     }

        
                        cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
                 }  */
                }
                if(part->pdgId()==23 && part->pt()!=0){
                       Z_Truth_no_pt_zero_Pt.push_back(part->pt()); 
                       Z_Truth_no_pt_zero_Eta.push_back(part->eta());
                       Z_Truth_no_pt_zero_Phi.push_back(part->phi());
                       Z_Truth_no_pt_zero_E.push_back(part->energy());  

     /*                   int nMothers = part->numberOfMothers(); 
                        for (int i = 0; i < nMothers; i++)
                                      {
                                        const auto* mother = part->mother(i);
                                        if (mother != nullptr) {
                                        cout << "mother " << i + 1 << " Id         :     "<<mother->pdgId()
                                             <<"     pt  :    "<<mother->pt()
                                             <<"     eta  :    "<<mother->eta()
                                             <<"     Phi  :    "<<mother->phi()
                                             <<"     Energy  :    "<<mother->energy()
                                             <<"     Status  :    "<<mother->status() << endl;
                                                                 }
                                       }              

                        cout<<"Z boson Particle Id :     "<<part->pdgId()
                            <<"     pt  :    "<<part->pt()
                            <<"     eta  :    "<<part->eta()
                            <<"     Phi  :    "<<part->phi()
                            <<"     Energy  :    "<<part->energy()
                            <<"     Status  :    "<<part->status()<<endl;          
                        
                        int nDaughters = part->numberOfDaughters(); 
                        for (int i = 0; i < nDaughters; i++)
                                     {
                                      const auto* daughter = part->daughter(i); 
                                      if (daughter != nullptr) {
                                      cout << "daughter " << i + 1 << " Id       :     " << daughter->pdgId()
                                           <<"     pt  :    "<<daughter->pt()
                                           <<"     eta  :    "<<daughter->eta()
                                           <<"     Phi  :    "<<daughter->phi()
                                           <<"     Energy  :    "<<daughter->energy()
                                           <<"     Status  :    "<<daughter->status() << endl;
                                                                 }
                                     }

        
                        cout<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
        */     }
                
                
        /*        if(part->pdgId()==1||part->pdgId()==2||part->pdgId()==3||part->pdgId()==4||part->pdgId()==5||part->pdgId()==6||part->pdgId()==7||part->pdgId()==8)
                { int i = 0;
                  int nDaughters = part->numberOfDaughters();
                  for (int i = 0; i < nDaughters; i++)
                                     {
                                      const auto* daughter = part->daughter(i); 
                                      if (daughter != nullptr) {
                                      cout <<" daughter " << i + 1 << " Id       :     " << daughter->pdgId()
                                           <<"     pt  :    "<<daughter->pt()
                                           <<"     eta  :    "<<daughter->eta()
                                           <<"     Phi  :    "<<daughter->phi()
                                           <<"     Energy  :    "<<daughter->energy()
                                           <<"     Status  :    "<<daughter->status() << endl;
                                                                 }
                                     }
                }        
        */
        
        
}
if(v==true){
        outFile<<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
        v=false;
        p=p+1;
}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

/////////////////////////Run, event, lumi//////////////////////////////////
   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();
///////////////////////////////////////////////////////////////////////////



   T->Fill(); // Write out the events


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
ZEE_RecHit_NTuplizer::beginJob()
{
        edm::Service<TFileService> fs;
        T=fs->make<TTree>("T","MyTuple");
	T->Branch("iEtaEle1"  ,  &(iEta[0]));
        T->Branch("iPhiEle1"  ,  &(iPhi[0]));
	T->Branch("Hit_ES_Eta_Ele1"  ,  &(Hit_ES_Eta[0]));
        T->Branch("Hit_ES_Phi_Ele1"  ,  &(Hit_ES_Phi[0]));
        T->Branch("Hit_ES_X_Ele1"  ,  &(Hit_ES_X[0]));
        T->Branch("Hit_ES_Y_Ele1"  ,  &(Hit_ES_Y[0]));
        T->Branch("Hit_ES_Z_Ele1"  ,  &(Hit_ES_Z[0]));
        T->Branch("ES_RecHitEnEle1"  ,  &(ES_RecHitEn[0]));

	T->Branch("Hit_Eta_Ele1"  ,  &(Hit_Eta[0]));
        T->Branch("Hit_Phi_Ele1"  ,  &(Hit_Phi[0]));
	T->Branch("Hit_X_Ele1"  ,  &(Hit_X[0]));
        T->Branch("Hit_Y_Ele1"  ,  &(Hit_Y[0]));
	T->Branch("Hit_Z_Ele1"  ,  &(Hit_Z[0]));
        T->Branch("RecHitEnEle1"  ,  &(RecHitEn[0]));
	T->Branch("RecHitFracEle1"  ,  &(RecHitFrac[0]));
	T->Branch("iEtaEle2"  ,  &(iEta[1]));
        T->Branch("iPhiEle2"  ,  &(iPhi[1]));
	T->Branch("Hit_ES_Eta_Ele2"  ,  &(Hit_ES_Eta[1]));
        T->Branch("Hit_ES_Phi_Ele2"  ,  &(Hit_ES_Phi[1]));
        T->Branch("Hit_ES_X_Ele2"  ,  &(Hit_ES_X[1]));
        T->Branch("Hit_ES_Y_Ele2"  ,  &(Hit_ES_Y[1]));
        T->Branch("Hit_ES_Z_Ele2"  ,  &(Hit_ES_Z[1]));
	T->Branch("ES_RecHitEnEle2"  ,  &(ES_RecHitEn[1]));

	T->Branch("Hit_Eta_Ele2"  ,  &(Hit_Eta[1]));
        T->Branch("Hit_Phi_Ele2"  ,  &(Hit_Phi[1]));
	T->Branch("Hit_X_Ele2"  ,  &(Hit_X[1]));
        T->Branch("Hit_Y_Ele2"  ,  &(Hit_Y[1]));
        T->Branch("Hit_Z_Ele2"  ,  &(Hit_Z[1]));
        T->Branch("RecHitEnEle2"  ,  &(RecHitEn[1]));
        T->Branch("RecHitFracEle2"  ,  &(RecHitFrac[1]));

        T->Branch("nElectrons",  &nElectrons_ , "nEle/I");
        T->Branch("pt"  ,  &Ele_pt_);
        T->Branch("eta" ,  &Ele_eta_ );
        T->Branch("phi" ,  &Ele_phi_ );
        T->Branch("energy", &Ele_energy_);
	T->Branch("energy_ecal", &Ele_ecal_energy_);
	T->Branch("energy_ecal_mustache", &Ele_ecal_mustache_energy_);

        T->Branch("passMediumId" ,  &passMediumId_ );
        T->Branch("passTightId"  ,  &passTightId_ );
        T->Branch("passMVAMediumId", &passMVAMediumId_);

        T->Branch("Ele_R9"  ,  &Ele_R9);
        T->Branch("Ele_S4"  ,  &Ele_S4);
        T->Branch("Ele_SigIEIE"  ,  &Ele_SigIEIE);
        T->Branch("Ele_SigIPhiIPhi" , &Ele_SigIPhiIPhi);
        T->Branch("Ele_SCEtaW"  ,  &Ele_SCEtaW);
        T->Branch("Ele_SCPhiW"  ,  &Ele_SCPhiW);
        T->Branch("Ele_CovIEtaIEta"  ,  &Ele_CovIEtaIEta);
        T->Branch("Ele_CovIEtaIPhi"  ,  &Ele_CovIEtaIPhi);
        T->Branch("Ele_ESSigRR"  ,  &Ele_ESSigRR);
        T->Branch("Ele_SCRawE"  ,  &Ele_SCRawE);
        T->Branch("Ele_SC_ESEnByRawE"  ,  &Ele_SC_ESEnByRawE);
	T->Branch("Ele_HadOverEm"  ,  &Ele_HadOverEm);
	
	T->Branch("HFEMClust_pt", &HFEMClust_pt_);
	T->Branch("HFEMClust_eta", &HFEMClust_eta_);
        T->Branch("HFEMClust_phi", &HFEMClust_phi_);
        T->Branch("HFEMClust_energy", &HFEMClust_energy_);
        T->Branch("HFEMClust_eLong1x1", &HFEMClust_eLong1x1_);
        T->Branch("HFEMClust_eShort1x1", &HFEMClust_eShort1x1_);
        T->Branch("HFEMClust_eLong3x3", &HFEMClust_eLong3x3_);
        T->Branch("HFEMClust_eShort3x3", &HFEMClust_eShort3x3_);
        T->Branch("HFEMClust_eLong5x5", &HFEMClust_eLong5x5_);
        T->Branch("HFEMClust_eShort5x5", &HFEMClust_eShort5x5_);
        T->Branch("HFEMClust_e1x1", &HFEMClust_e1x1_);
        T->Branch("HFEMClust_e3x3", &HFEMClust_e3x3_);
        T->Branch("HFEMClust_e5x5", &HFEMClust_e5x5_);
        T->Branch("HFEMClust_eSeL", &HFEMClust_eSeL_);
        T->Branch("HFEMClust_eCOREe9", &HFEMClust_eCOREe9_);
        T->Branch("HFEMClust_e9e25", &HFEMClust_e9e25_);
        T->Branch("HFEMClust_eCore", &HFEMClust_eCore_);
/*	
        T->Branch("Ele_Gen_Pt" , &Ele_Gen_Pt);
        T->Branch("Ele_Gen_Eta" , &Ele_Gen_Eta);
        T->Branch("Ele_Gen_Phi" , &Ele_Gen_Phi);
        T->Branch("Ele_Gen_E" , &Ele_Gen_E);
  */
        T->Branch("rho", &rho, "rho/F");

        T->Branch("run",&run,"run/I");
        T->Branch("event",&event,"event/I");
        T->Branch("lumi",&lumi,"lumi/I");

        ////////////added by Harsh////////////
                T->Branch("Ele_Gen_Pt" , &Ele_Gen_Pt);
                T->Branch("Ele_Gen_Eta" , &Ele_Gen_Eta);
                T->Branch("Ele_Gen_Phi" , &Ele_Gen_Phi);
                T->Branch("Ele_Gen_E" , &Ele_Gen_E);

                T->Branch("ZTruth_e_Pt", &Z_Truth_e_Pt);
                T->Branch("ZTruth_e_Eta", &Z_Truth_e_Eta);
                T->Branch("ZTruth_e_Phi", &Z_Truth_e_Phi);
                T->Branch("ZTruth_e_E", &Z_Truth_e_E);

                T->Branch("Z_Truth_all_Pt", &Z_Truth_all_Pt); 
                T->Branch("Z_Truth_all_Eta", &Z_Truth_all_Eta);
                T->Branch("Z_Truth_all_Phi", &Z_Truth_all_Phi);
                T->Branch("Z_Truth_all_E", &Z_Truth_all_E);

                T->Branch("Z_Truth_no_pt_zero_Pt",&Z_Truth_no_pt_zero_Pt);
                T->Branch("Z_Truth_no_pt_zero_Eta",&Z_Truth_no_pt_zero_Eta);
                T->Branch("Z_Truth_no_pt_zero_Phi",&Z_Truth_no_pt_zero_Phi);
                T->Branch("Z_Truth_no_pt_zero_E",&Z_Truth_no_pt_zero_E);


}

// ------------ method called once each job just after ending the event loop  ------------
void
ZEE_RecHit_NTuplizer::endJob()
{
}


// Extract the rechits of the SC from the ES layers
void ZEE_RecHit_NTuplizer::GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int elenum, unsigned int planeIndex) {
        double RawenergyPlane = 0.;
        double pfRawenergyPlane = 0.;
//      if(!_ESRechitsHandle.isValid())
//              return RawenergyPlane;
              
//        if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
//                break;

	int NumHits = 0;

	const CaloSubdetectorGeometry* ecalESGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalPreshower));


        for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
                const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
                for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
                        //      std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
                        //      ++iCount;
                        for(ESRecHitCollection::const_iterator esItr = ESRechitsHandle->begin(); esItr != ESRechitsHandle->end(); ++esItr) {//loop over ES rechits to find the one in the cluster
                                ESDetId rhid = ESDetId(esItr->id());
                                if(rhid == (*rh).first) { // found ESrechit
//                                        std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
                                        if((int) rhid.plane() == (int) planeIndex) {
						std::shared_ptr<const CaloCellGeometry> geom = ecalESGeom->getGeometry(rhid);
						Hit_ES_Eta[elenum].push_back( geom->etaPos() );
                                                Hit_ES_Phi[elenum].push_back( geom->phiPos() );
						Hit_ES_X[elenum].push_back( geom->getPosition().x() );
						Hit_ES_Y[elenum].push_back( geom->getPosition().y() );
						Hit_ES_Z[elenum].push_back( geom->getPosition().z() ) ;
						ES_RecHitEn[elenum].push_back(esItr->energy());
//						std::cout << "Preshower" <<std::setprecision(4) << " Eta = " <<geom->etaPos() << " : " <<" Phi = "<< geom->phiPos() << " 3D position" << geom->getPosition().z() << std::endl;
                                                RawenergyPlane += esItr->energy();
                                                pfRawenergyPlane += rh->second;
						NumHits++;
                                        }
                                        break;
                                }
                        }
                }

        }

 //       return RawenergyPlane;
}





//Clear tree vectors each time analyze method is called
void ZEE_RecHit_NTuplizer::ClearTreeVectors()
{
	nElectrons_ = 0;
	iEta[0].clear();
	iPhi[0].clear();

	
	Hit_ES_Eta[0].clear();
        Hit_ES_Phi[0].clear();
        Hit_ES_X[0].clear();
        Hit_ES_Y[0].clear();
        Hit_ES_Z[0].clear();
	ES_RecHitEn[0].clear();


        Hit_Eta[0].clear();
	Hit_Phi[0].clear();
	Hit_X[0].clear();
	Hit_Y[0].clear();
	Hit_Z[0].clear();
	RecHitEn[0].clear();
	RecHitFrac[0].clear();
	iEta[1].clear();
        iPhi[1].clear();

	Hit_ES_Eta[1].clear();
        Hit_ES_Phi[1].clear();
        Hit_ES_X[1].clear();
        Hit_ES_Y[1].clear();
        Hit_ES_Z[1].clear();
	ES_RecHitEn[1].clear();

	Hit_Eta[1].clear();
        Hit_Phi[1].clear();
        Hit_X[1].clear();
        Hit_Y[1].clear();
        Hit_Z[1].clear();
        RecHitEn[1].clear();
        RecHitFrac[1].clear();
	Ele_pt_.clear();
	Ele_eta_.clear();
	Ele_phi_.clear();
	Ele_energy_.clear();
	Ele_ecal_energy_.clear();
	Ele_ecal_mustache_energy_.clear();		
	Ele_R9.clear();
	Ele_S4.clear();
	Ele_SigIEIE.clear();
	Ele_SigIPhiIPhi.clear();
	Ele_SCEtaW.clear();
	Ele_SCPhiW.clear();
	Ele_CovIEtaIEta.clear();
	Ele_CovIEtaIPhi.clear();
	Ele_ESSigRR.clear();
	Ele_SCRawE.clear();
	Ele_SC_ESEnByRawE.clear();
	Ele_HadOverEm.clear();

	HFEMClust_pt_.clear();
	HFEMClust_eta_.clear();
	HFEMClust_phi_.clear();
	HFEMClust_energy_.clear();
	HFEMClust_eLong1x1_.clear();
	HFEMClust_eShort1x1_.clear();
	HFEMClust_eLong3x3_.clear();
	HFEMClust_eShort3x3_.clear();
	HFEMClust_eLong5x5_.clear();
	HFEMClust_eShort5x5_.clear();
	HFEMClust_e1x1_.clear();
	HFEMClust_e3x3_.clear();
	HFEMClust_e5x5_.clear();
	HFEMClust_eSeL_.clear();
	HFEMClust_eCOREe9_.clear();
	HFEMClust_e9e25_.clear(); 
	HFEMClust_eCore_.clear();

/*
        Ele_Gen_Pt.clear();
	Ele_Gen_Eta.clear();
	Ele_Gen_Phi.clear();
	Ele_Gen_E.clear();
*/

/////////////////////added by Harsh////////////////////
	Ele_Gen_Pt.clear();
	Ele_Gen_Eta.clear();
	Ele_Gen_Phi.clear();
	Ele_Gen_E.clear();

        Z_Truth_e_Pt.clear();
        Z_Truth_e_Eta.clear();
        Z_Truth_e_Phi.clear();
        Z_Truth_e_E.clear();

        Z_Truth_all_Pt.clear(); 
        Z_Truth_all_Eta.clear();
        Z_Truth_all_Phi.clear();
        Z_Truth_all_E.clear();  

        Z_Truth_no_pt_zero_Pt.clear();  
        Z_Truth_no_pt_zero_Eta.clear(); 
        Z_Truth_no_pt_zero_Phi.clear(); 
        Z_Truth_no_pt_zero_E.clear();   
///////////////////////////////////////////////////////


	passMediumId_.clear();
	passTightId_ .clear();
	passMVAMediumId_.clear();

	isTrue_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZEE_RecHit_NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEE_RecHit_NTuplizer);