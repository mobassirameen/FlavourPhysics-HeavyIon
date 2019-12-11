//description on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BsJpsiPhi_AWG 
#include <cstddef>
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <cfloat>
#include <string>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/BsToJpsiPhiAnalysis.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
//#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerReport.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector.h"
#include "TLorentzRotation.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/KinematicFitInterface.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include <iostream>
#include <TMath.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;



BsToJpsiPhiAnalysis::BsToJpsiPhiAnalysis(const edm::ParameterSet& iConfig):theConfig_(iConfig),
nominalJpsiMass( 3.096916 ),
nominalPhiMass(1.019 ),
nominalElectronMass(0.00051099893),
nominalMuonMass(0.1056583),
nominalKaonMass(0.493677),
nominalPionMass(0.139570),
nominalKstarMass(0.892),
nominalBplusMass(5.2792)
{

  isMCstudy_ = iConfig.getParameter<bool>("isMCstudy");
//  TriggerTag = iConfig.getParameter<edm::InputTag>("TriggerTag");
//  TriggerTagTok = consumes<edm::TriggerResults>(TriggerTag);
  genParticlesLabel                 = iConfig.getParameter<InputTag>("genParticlesLabel");
  genParticlesTok                   = consumes<edm::View<reco::GenParticle>>(genParticlesLabel); 
  MuonTag                           = iConfig.getParameter<edm::InputTag>("MuonTag");
  MuonTagTok                        = consumes<edm::View<pat::Muon>>(MuonTag);
  JetCollection                     = iConfig.getParameter<edm::InputTag>("JetCollection");
  JetCollectionTok                  = consumes<edm::View<pat::Jet>>(JetCollection);
  PUInfo                            = iConfig.getParameter<InputTag>("PUInfo");
  PUInfoTok                         = consumes<edm::View<PileupSummaryInfo>>(PUInfo);
  vertexBeamSpot                    = iConfig.getParameter<edm::InputTag>("vertexBeamSpot");
  vertexBeamSpotTok                 = consumes<reco::BeamSpot>(vertexBeamSpot);
  primaryvertex                     = iConfig.getParameter<edm::InputTag>("primaryvertex");
  primaryvertexTok                  = consumes<edm::View<reco::Vertex>>(primaryvertex);
  triggerresults                    = iConfig.getParameter<edm::InputTag>("triggerresults");
  triggerresultsTok                 = consumes<edm::TriggerResults>(triggerresults);
  ElectronTag                       = iConfig.getParameter<edm::InputTag>("ElectronTag");
  ElectronTagTok                    = consumes<edm::View<pat::Electron>>(ElectronTag);
  track                             = iConfig.getParameter<edm::InputTag>("track");
  trackLabelK                       = consumes<edm::View<pat::PackedCandidate>>(track);
  isotrack                          = iConfig.getParameter<edm::InputTag>("isotrack");
  isotrackTok                       = consumes<edm::View<pat::IsolatedTrack>>(isotrack); 

//=================================================================================================================
//=================================================================================================================
  StoreDeDxInfo_ = iConfig.getParameter<bool>("StoreDeDxInfo");
  JpsiMassWindowBeforeFit_ = iConfig.getParameter<double>("JpsiMassWindowBeforeFit");

  BsLowerMassCutBeforeFit_  = iConfig.getParameter<double>("BsLowerMassCutBeforeFit");
  BsUpperMassCutBeforeFit_  = iConfig.getParameter<double>("BsUpperMassCutBeforeFit");
  BsLowerMassCutAfterFit_  = iConfig.getParameter<double>("BsLowerMassCutAfterFit");
  BsUpperMassCutAfterFit_  = iConfig.getParameter<double>("BsUpperMassCutAfterFit");

  JpsiMassWindowAfterFit_ = iConfig.getParameter<double>("JpsiMassWindowAfterFit");
  JpsiPtCut_ =  iConfig.getParameter<double>("JpsiPtCut");
  KaonTrackPtCut_ = iConfig.getParameter<double>("KaonTrackPtCut");
  BdKaonTrackPtCut_ = iConfig.getParameter<double>("BdKaonTrackPtCut");
  PhiMassWindowBeforeFit_ = iConfig.getParameter<double>("PhiMassWindowBeforeFit");
  PhiMassWindowAfterFit_ = iConfig.getParameter<double>("PhiMassWindowAfterFit");

  KstarMassWindowBeforeFit_ = iConfig.getParameter<double>("KstarMassWindowBeforeFit");
  KstarMassWindowAfterFit_ = iConfig.getParameter<double>("KstarMassWindowAfterFit");
  BdLowerMassCutBeforeFit_ = iConfig.getParameter<double>("BdLowerMassCutBeforeFit");
  BdUpperMassCutBeforeFit_ = iConfig.getParameter<double>("BdUpperMassCutBeforeFit");

  BdLowerMassCutAfterFit_ = iConfig.getParameter<double>("BdLowerMassCutAfterFit");
  BdUpperMassCutAfterFit_ = iConfig.getParameter<double>("BdUpperMassCutAfterFit");

  BdPDGMass_ = iConfig.getParameter<double>("BdPDGMass");
  BpPDGMass_ = iConfig.getParameter<double>("BpPDGMass");
  BsPDGMass_ = iConfig.getParameter<double>("BsPDGMass");

  verbose_                = iConfig.getParameter<bool>("verbose");
  TestVerbose_            = iConfig.getParameter<bool>("TestVerbose");
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile");
  event_counter_ = 0;
  elecounter_    = 0;
  muoncounter_   = 0;
  jetcounter_    = 0;
  tagmucounter_ =0;

edm::LogInfo("RecoVertex/BsToJpsiPhiAnalysis")<< "Initializing Bs to Jpsi Phi analyser  - Output file: " << outputFile_ <<"\n";
}

BsToJpsiPhiAnalysis::~BsToJpsiPhiAnalysis() {}
void BsToJpsiPhiAnalysis::beginJob()
{
  bsRootTree_ = new BsToJpsiPhiRootTree();
  bsRootTree_->createTree(outputFile_);
}
void BsToJpsiPhiAnalysis::endJob()
{
  bsRootTree_->writeFile();
  delete bsRootTree_;
  cout << "Total number of Events          : " << event_counter_ << endl;
  cout << "Total number of Tagged muons    : " << muoncounter_   << endl;
  cout << "Total number of Tagged electrons: " << elecounter_    << endl;
  cout << "Total number of Tagged jets     : " << jetcounter_    << endl;
  cout << "Max amount of Tag muons		    : " << tagmucounter_ << endl;
}

void BsToJpsiPhiAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{
  event_counter_++;
  //cout << "NEW EVENT " << endl;
 // double bestVtxProbBplus = -1;
 // double MinBdVtxProb = -999;
  double minVtxP = -99.;   //KK hypothesis
  bsRootTree_->resetEntries();
//===========================================================TLorentzVectors of B mesons after the kin fit
  TLorentzVector TheBp;
  TLorentzVector TheBs;
  TLorentzVector TheBd;
//===========================================================clear Bu branches
  bsRootTree_->BpPVTrkPt_->clear();
  bsRootTree_->BpPVTrkCharge_->clear();
  bsRootTree_->BpPVTrkEta_->clear();
  bsRootTree_->BpPVTrkPhi_->clear();
  bsRootTree_->BpBJetTrkCharge_->clear();
  bsRootTree_->BpBJetTrkPt_->clear();
//=========================================================Clear PV and BJet Branches
  bsRootTree_->PVTrkPt_->clear();
  bsRootTree_->PVTrkCharge_->clear();
  bsRootTree_->PVTrkEta_->clear();
  bsRootTree_->PVTrkPhi_->clear();
  bsRootTree_->BJetTrkCharge_->clear();
  bsRootTree_->BJetTrkPt_->clear();
//========================================================Create objects needed for building the B candidate
  pat::CompositeCandidate BCand_best;
  TrackRef trk1Ref_best;
  TrackRef trk2Ref_best;
  TrackRef trkMu1Ref_best;
  TrackRef trkMu2Ref_best;
  RefCountedKinematicParticle bs_best;
  pat::Muon mu1_best;
  pat::Muon mu2_best;
//=========================================================Define Primary vertices(PVs)
  Vertex PVvtxCosTheta;
  Vertex BpPVvtxCosTheta;
  Vertex BdPVvtxCosTheta;
 
  int BsPVVtxInd=0;
  vector<TrackRef> BpTrkRefs;
  vector<TrackRef> BsTrkRefs;
  vector<TrackRef> BdTrkRefs;
//=======================================================Save PU information
if(isMCstudy_){
//===================================================Creat handle for PU information
    
    edm:: Handle<edm::View<PileupSummaryInfo> > PUinfo;
    iEvent.getByToken( PUInfoTok, PUinfo);
    edm::View<PileupSummaryInfo>::const_iterator PVI;
    int numInteraction = 0;
    for(PVI = PUinfo->begin(); PVI != PUinfo->end(); ++PVI)
    {
      if (PVI->getBunchCrossing()==0)
        numInteraction += PVI->getPU_NumInteractions();
    }
    bsRootTree_->PUinteraction_ = numInteraction; /// SaveToTree -> total number of interactions (with PU)
    //cout<<numInteraction<<"\n";
   }

//================================================================Jets
 // edm::Handle<edm::View<pat::Jet>> myjets;
 // iEvent.getByToken(JetCollectionTok,myjets);
 // const edm::View<pat::Jet> & jets = *myjets;
//==========================================================Get Primary Vertices
 int    VtxIndex    = -99;
//====================================================Beam Spot
  double BSx         = -9999999.;
  double BSy         = -9999999.;
  double BSz         = -9999999.;
  double BSdx        = -9999999.;
  double BSdy        = -9999999.;
  double BSdz        = -9999999.;
  double BSdxdz      = -9999999.;
  double BSdydz      = -9999999.;
  double BSsigmaZ    = -9999999.;
  double BSdsigmaZ   = -9999999.;
//========================================================PV
  double PVx         = -9999999.;
  double PVy         = -9999999.;
  double PVz         = -9999999.;
  double PVerrx      = -9999999.;
  double PVerry      = -9999999.;
  double PVerrz      = -9999999.;

TLorentzVector kaontrack1, kaontrack2;
//==================================================Handle and assign beam spot parameters
   edm::Handle<reco::BeamSpot> vertexBeamSpot ;
   iEvent.getByToken(vertexBeamSpotTok,vertexBeamSpot);

         BSx = vertexBeamSpot->x0(); 
         BSy = vertexBeamSpot->y0();
         BSz = vertexBeamSpot->z0(); 
      //   BSdx = vertexBeamSpot->x0Error();
       //  BSdy = vertexBeamSpot->y0Error(); 
       //  BSdz = vertexBeamSpot->z0Error(); 
         BSdxdz = vertexBeamSpot->dxdz();
         BSdydz = vertexBeamSpot->dydz(); 
        // BSsigmaZ = vertexBeamSpot->sigmaZ(); 
        // BSdsigmaZ = vertexBeamSpot->sigmaZ0Error();



   edm::Handle<edm::View<reco::Vertex> > recVtxs;
   iEvent.getByToken(primaryvertexTok, recVtxs);
   
  bsRootTree_->NVertices_ = recVtxs->size();
  //double MinPtVertex = 0.;
//cout<<recVtxs<<"\n";
//==============================================================================================VERTEX 
  for(size_t iVtx = 0; iVtx < recVtxs->size(); ++ iVtx)
 { 
           VtxIndex = iVtx;
         //  cout<<iVtx<<"\n";
         // const Vertex &RecVtx = (*recVtxs)[iVtx];
         /*   
         //cout<<vtx<<"\n";
           double PtSumVertex = 0.;
           for (reco::Vertex::trackRef_iterator trackvertex = vtx.tracks_begin(); trackvertex != vtx.tracks_end(); trackvertex++)
           {//  cout<<trackvertex<<"\n";
          cout<<"pt: "<<(**trackvertex).pt()<<"\n";
          const reco::Track & VtxTrack = *(trackvertex->get());
          //const pat::PackedCandidate & VtxTrack = *(trackvertex->get());
          PtSumVertex = PtSumVertex + abs(VtxTrack.pt());
          // cout<<"ptsumvertex: "<<PtSumVertex<<"\n";
           }
           if(PtSumVertex > MinPtVertex)
          {int VtxIndex = -99;
           VtxIndex = iVtx;
           MinPtVertex = PtSumVertex;
           cout<<"minvertex"<<MinPtVertex<<"\n";
           }*/
 }//-----------------------------------------------------VERTEX LOOP ENDS

//===================================================================Save as PV the one with the highest Pt (called MinPtVertex)
      const Vertex &RecVtx = (*recVtxs)[VtxIndex];


  if(VtxIndex!=-99) /// Use PV
    {
      bsRootTree_->isPV_ = 1;
      PVx = RecVtx.x();
      PVy= RecVtx.y();
      PVz= RecVtx.z();
      PVerrx=RecVtx.xError();
      PVerry=RecVtx.yError();
      PVerrz=RecVtx.zError();
    }
  else {  /// Use BS
            bsRootTree_->isBS_ = 1;
            PVx=BSx;
            PVy=BSy;
            PVz=BSz;
            PVerrx=BSdx;
            PVerry=BSdy;
            PVerrz=BSdz;
            std:: cout<<PVerrz<<"\n";
       }

             bsRootTree_->getVtx(BSx,BSy,BSz,PVx,PVy,PVz,PVerrx,PVerry,PVerrz);
             
             bsRootTree_->BSdx_ = BSdx;
             bsRootTree_->BSdy_ = BSdy;
             bsRootTree_->BSdz_ = BSdz;
             bsRootTree_->BSsigmaZ_ = BSsigmaZ;
             bsRootTree_->BSdsigmaZ_ = BSdsigmaZ;
             bsRootTree_->BSdxdz_ = vertexBeamSpot->dxdz();
             bsRootTree_->BSdydz_ = vertexBeamSpot->dydz();

 if(verbose_ == true){
                         std::cout<<"BeamSpot   (x,y,z) = ("<< BSx << ", " << BSy << ", "<< BSz << ")\n";
                         std::cout<<"PrimaryVtx (x,y,z) = ("<< PVx <<" , " << PVy << ", "<< PVz << ")\n";
                     }

//========================================================================================================Run/Event/Lumiblock
         bsRootTree_->runNumber_ = iEvent.id().run();
        //  cout<<iEvent.id().run()<<"\n";
        bsRootTree_->eventNumber_ = (unsigned int)iEvent.id().event();
        bsRootTree_->lumiSection_ = iEvent.luminosityBlock();
//=============================================================================================GEN_PARTICLES
   edm::Handle<edm::View<reco::GenParticle> > genParticles;
   if(isMCstudy_)
              {
                     iEvent.getByToken(genParticlesTok, genParticles);
                     //std::cout<<"genparticles:"<<genParticles->size()<<"\n";
                     fillMCInfo(genParticles);
               }

           //vector<pat::IsolatedTrack>            "isolatedTracks"   
           //Create Handle to dE/dx
           //m:: Handle<DeDxDataValueMap> energyLossHandle;
           //edm::Handle<edm::View<pat::IsolatedTrack>> energyLossHandle;
           //if(StoreDeDxInfo_) iEvent.getByToken(isotrackTok, energyLossHandle);

//===================================================================================================================HANDLE OVER TRIGGER
   edm::Handle<edm::TriggerResults> hltresults;
   iEvent.getByToken(triggerresultsTok, hltresults);
                           //std::cout<<"hlt:  "<<hltresults->size()<<"\n";

  
   const  edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
   int ntrigs = hltresults->size();
   for (int itrig = 0; itrig != ntrigs; ++itrig)
        {
                                                  TString trigName = triggerNames_.triggerName(itrig);
                                                   //std::cout<<"triggernames:"<<itrig<<"::"<<trigName<<"\n";

                if (trigName=="HLT_DoubleMu4_JpsiTrk_Displaced_v12")       bsRootTree_->triggerbit_HLTmu4TkDis_                = hltresults->accept(itrig);
                if (trigName=="HLT_DoubleMu4_JpsiTrkTrk_Displaced_v4")     bsRootTree_->triggerbit_HLTmu4TkTkDis_              = hltresults->accept(itrig);
                if (trigName=="HLT_Dimuon0_Jpsi3p5_Muon2_v4")              bsRootTree_->triggerbit_HLTDimuon0JpsiMuon_         = hltresults->accept(itrig);
//HLT_DoubleMu4_JpsiTrkTrk_Displaced_v4
//HLT_Dimuon0_Jpsi3p5_Muon2_v4
//HLT_DoubleMu4_JpsiTrk_Displaced_v12
                string str = (string) trigName  ;
     
               if (str.compare(0,18,"HLT_DoubleMu4_3_Jpsi_Displaced_v12") == 0)
                   {
                        cout << "trigger found! " << str << endl;
                        cout <<"dimuon4_jpsidis  " << bsRootTree_->triggerbit_HLTmu4TkDis_  << endl;
                        cout <<"dimuon4_jpsitktkdis  " << bsRootTree_->triggerbit_HLTmu4TkTkDis_ << endl;
                    }

         }//Trigger loop ends
//===============================================================================================================MUONS_ELLECTRONS
            edm::Handle< View<pat::Muon> > allmuons;
            iEvent.getByToken(MuonTagTok, allmuons);
                          //if(allmuons->size()>0)std::cout<<"muonmultiplicity"<<allmuons->size()<<"\n";
            edm::Handle< View<pat::Electron> > allelectrons;
            iEvent.getByToken(ElectronTagTok, allelectrons);
                          //if(allelectrons->size()>0)std::cout<<"electronmultiplicity"<<allelectrons->size()<<"\n";

                     if(verbose_==true)
                      {
                          if(allmuons->size()>0)     std::cout<<"******found number of muons    = "<< allmuons->size() << std::endl;
                          if(allelectrons->size()>0) std::cout<<"******found number of electrons= "<< allelectrons->size() << std::endl;
                       }
             

                      

                        bsRootTree_->MuonMultiplicity_ = allmuons->size();
                        bsRootTree_->ElectronMultiplicity_ = allelectrons->size();
                        if(allelectrons->size() > 2 )
                        {tagmucounter_++;}
                        if (allmuons->size() > 25) cout << "WARNING! muon list oversize:" << allmuons->size() << endl;
     

//==================================================================================================================Muon LOOP
//-------------------------------------------------------First Muon
      for(size_t i=0; i < allmuons->size(); ++i)
          {
                    const pat::Muon & mu1 = (*allmuons)[i];
                    if (mu1.innerTrack().isNull()){continue;}
                    //cout<<"mu1pt: "<<mu1.pt()<<"\n";
                    if(verbose_ == true)
                        {
                                std::cout<<"Got one muon "<<mu1.pt()<<std::endl;
                        }

//------------------------------------------------------Loop over 2nd muon
                   for (size_t j=i+1; j < allmuons->size(); ++j)
                        {
                              const pat::Muon & mu2 = (*allmuons)[j];
                              if (mu2.innerTrack().isNull()){continue;}	   
                              //  cout<<"mu2pt: "<<mu2.pt()<<"\n";
                              if(verbose_ == true)
                                  {
                                         std::cout<<"Got the second muon "<<mu2.pt()<<std::endl;
                                   }

                     /*if (!isMCstudy_) {

                           if ( 
				(mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty() ||
                                mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty()  ) &&
				(mu1.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty() ||
                                  mu2.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty()    ) 
	                     	)  continue;
                 		else {cout << "dimu8_2 "<<mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty()<<"dimu8_1 "<<mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty() << endl;}
                            }
                         */
                  // if(not(mu1.passed(reco::Muon::SoftCutBasedId)))continue;
                  // if(not(mu2.passed(reco::Muon::SoftCutBasedId)))continue;
  
                   bsRootTree_->Mu1SoftID_ = mu1.passed(reco::Muon::SoftCutBasedId);
                   bsRootTree_->Mu2SoftID_ = mu2.passed(reco::Muon::SoftCutBasedId);
                   if(!mu1.isGlobalMuon() && !mu1.isTrackerMuon()) continue; // skip if mu1 is not GLB or TRK (PASS IF LooseMuId)
                   if(!mu2.isGlobalMuon() && !mu2.isTrackerMuon()) continue; // skip if mu2 is not GLB or TRK (PASS IF LooseMuId)

                             if(verbose_==true) {
                                          std::cout << "******mu1.isGlobalMuon() == "<<mu1.isGlobalMuon() << "\n";
                                          std::cout << "      mu1.isTrackerMuon()== "<<mu1.isTrackerMuon()<< "\n";
                                          std::cout << "      mu2.isGlobalMuon() == "<<mu2.isGlobalMuon() << "\n";
                                          std::cout << "      mu2.isTrackerMuon()== "<<mu2.isTrackerMuon()<< "\n";
                                       }

                  // if( mu1.isTrackerMuon() && !muon::isGoodMuon(mu1, muon::TrackerMuonArbitrated) ) continue; // Skip if mu1 is not TrackerMuonArbitrated
                  // if( mu2.isTrackerMuon() && !muon::isGoodMuon(mu2, muon::TrackerMuonArbitrated) ) continue; // Skip if mu2 is not TrackerMuonArbitrated



                   bsRootTree_->ihaveajpsi_=1;

                   if(mu1.charge()==mu2.charge()) continue; // Skip iif mu1 e mu2 have the same charge
                          if(verbose_==true) {
                                                  std::cout<<"******MUONS HAVE OPPOSITE CHARGE: mu1.charge() = " <<mu1.charge()<<" , mu2.charge() = "<<mu2.charge()<<std::endl;
                                             }

                                                if(bsRootTree_->iPassedCutIdent_   < 1 )   bsRootTree_->iPassedCutIdent_ = 1 ;
                                                if(bsRootTree_->iPassedCutIdentBd_   < 1 ) bsRootTree_->iPassedCutIdentBd_ = 1 ;

                    pat::CompositeCandidate Jpsi;
                    Jpsi.addDaughter(mu1);
                    Jpsi.addDaughter(mu2);
                    AddFourMomenta addP4;
                    addP4.set(Jpsi);

                          if(verbose_==true) 
                                          {
                                                 std::cout<<"******Di-Muon Mass= " <<Jpsi.mass()<<std::endl;
                                           }
                    if ( abs(Jpsi.mass() - nominalJpsiMass ) > JpsiMassWindowBeforeFit_ ) continue; // skip if mu1-mu2 combination mass is far from JPsi
                    if ( Jpsi.pt() < JpsiPtCut_) continue;                                          // skip if mu1-mu2 combination pt is less than JPsi Pt cut

                    //Jpsi window validity
 
                   if(bsRootTree_->iPassedCutIdent_   < 2 )   bsRootTree_->iPassedCutIdent_   = 2 ;
                   if(bsRootTree_->iPassedCutIdentBd_   < 2 ) bsRootTree_->iPassedCutIdentBd_ = 2 ;

                  edm::ESHandle<TransientTrackBuilder> theB;
                  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
                  TrackRef trkMu1Ref = mu1.get<TrackRef>();
                  TrackRef trkMu2Ref = mu2.get<TrackRef>();
                  TrackRef muonTrkP = mu1.track();
                  TrackRef muonTrkM = mu2.track();
                  vector<TransientTrack> trk_all;
                  TransientTrack mu1TT=(*theB).build(&trkMu1Ref);
                  TransientTrack mu2TT=(*theB).build(&trkMu2Ref);
                  trk_all.push_back(mu1TT);
                  trk_all.push_back(mu2TT);

            KalmanVertexFitter kvf(true);
            TransientVertex tv = kvf.vertex(trk_all);

    
            if(verbose_==true) 
                 {
                             std::cout<<"****** MUONS HAVE VALID VERTEX FIT"<< std::endl;
                 }



             // vertex validity
             if(bsRootTree_->iPassedCutIdent_   < 3 )   bsRootTree_->iPassedCutIdent_   = 3 ;     
             if(bsRootTree_->iPassedCutIdentBd_   < 3 ) bsRootTree_->iPassedCutIdentBd_ = 3 ;

             int isCowboy=0;
             if (mu1.charge()==1) {
             float mupPhi = atan(mu1.py()/mu1.px());
             if ( mu1.px() < 0 && mu1.py() < 0 ) mupPhi -= TMath::Pi();
             if ( mu1.px() < 0 && mu1.py() > 0 ) mupPhi += TMath::Pi();
             float mumPhi = atan(mu2.py()/mu2.px());
             if ( mu2.px() < 0 && mu2.py() < 0 ) mumPhi -= TMath::Pi();
             if ( mu2.px() < 0 && mu2.py() > 0 ) mumPhi += TMath::Pi();
             //std::cout<<"mumphi"<<mumPhi<<"\n";
             if ( (mupPhi - mumPhi)>0 ) isCowboy=1; } 
             else {
                   float mupPhi = atan(mu2.py()/mu2.px());
                   if ( mu2.px() < 0 && mu2.py() < 0 ) mupPhi -= TMath::Pi();
                   if ( mu2.px() < 0 && mu2.py() > 0 ) mupPhi += TMath::Pi();
                   float mumPhi = atan(mu1.py()/mu1.px());
                   if ( mu1.px() < 0 && mu1.py() < 0 ) mumPhi -= TMath::Pi();
                   if ( mu1.px() < 0 && mu1.py() > 0 ) mumPhi += TMath::Pi();
                   if ( (mupPhi - mumPhi)>0 ) isCowboy=1;
                  }

//int lalau =pow(isCowboy, 2);

//==========================================================================ceate vertex for mu1-mu2 combination
 Vertex vertex = tv;
//calculate variable in the cloasest way to trigger
      double  vtxProb_Jpsi = TMath::Prob(vertex.chi2(),(int)vertex.ndof());
      math::XYZVector      pperp(mu1.px() + mu2.px(), mu1.py() + mu2.py(), 0.);
      reco::Vertex::Point  vpoint=vertex.position();
//std::cout<<"vpoint"<<vpoint<<"\n";

//============================================================================================Translate into global point
      GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());
      GlobalPoint displacementFromBeamspot( -1*((BSx -  secondaryVertex.x()) +  (secondaryVertex.z() - BSz) * BSdxdz),-1*((BSy - secondaryVertex.y())+  (secondaryVertex.z() - BSz) * BSdydz), 0);
      reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
      double CosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
  //   std::cout<<"opening angle: "<<CosAlpha<<"\n";
       double MuonsDCA=999;
      TrajectoryStateClosestToPoint mu1TS = mu1TT.impactPointTSCP();
      TrajectoryStateClosestToPoint mu2TS = mu2TT.impactPointTSCP();
      if (mu1TS.isValid() && mu2TS.isValid()) {
       
       ClosestApproachInRPhi cApp;
       cApp.calculate(mu1TS.theState(), mu2TS.theState());
      MuonsDCA=cApp.distance();
      }
      double max_Dr1=fabs( (- (mu1.vx()-BSx) * mu1.py() + (mu1.vy()-BSy) * mu1.px() ) / mu1.pt() );
      double max_Dr2=fabs( (- (mu2.vx()-BSx) * mu2.py() + (mu2.vy()-BSy) * mu2.px() ) / mu2.pt() );
    // std::cout<<"maxDR1: "<<max_Dr1<<"\n";
//========================================================================================================================muon overlaping remover
      if ( muon::overlap(mu1,mu2,1,1,true) ) continue; /// Skip the mu-mu combination if the two muons overlap
      reco::Vertex::Error verr = vertex.error();
      //std::cout<<"Verror: "<<verr<<"\n";
//==================================================================================================================Translate to global error
      GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2) );
      float lxy = displacementFromBeamspot.perp();
      float lxyerr = sqrt(err.rerr(displacementFromBeamspot));
      //std::cout<<"LXYerror: "<<lxyerr<<"\t"<<"LXY: "<<lxy<<"\n";

//====================================================================================

//if( (!mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty()  == 1 && !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty()  == 1 && ( bsRootTree_->triggerbit_HLTmu3Tk_== 1 || bsRootTree_->triggerbit_HLTmu5_== 1 ) ) || (!mu2.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty()  == 1 && !mu1.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty()  == 1 && ( bsRootTree_->triggerbit_HLTmu3Tk_== 1 || bsRootTree_->triggerbit_HLTmu5_== 1 ) ) )


      bsRootTree_->JpsiNumberOfCandidates_++;
                   
        //bsRootTree_->JpsiMuon1MatchDimuon0_Novtx_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty() ;
        //bsRootTree_ ->JpsiMuon2MatchDimuon0_Novtx_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu43Jpsi").empty();
        //bsRootTree_->JpsiMuonMatch41_ = !mu1.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty() ;
        //bsRootTree_ ->JpsiMuonMatch42_ = !mu2.triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar").empty();
      
      bsRootTree_->JpsiM_alone_ = Jpsi.mass();
   //   std::cout<<"jpsimass:"<<Jpsi.mass()<<"\n";
      bsRootTree_->JpsiPhi_alone_ = Jpsi.phi();
      bsRootTree_->JpsiEta_alone_ = Jpsi.eta();
      bsRootTree_->JpsiPt_alone_ = Jpsi.pt();
     if( mu1.charge() == -1 )
       { 

                      bsRootTree_->Mu1Pt_beffit_   = mu1.pt();
                      bsRootTree_->Mu1Pz_beffit_   = mu1.pz();
                      bsRootTree_->Mu1Eta_beffit_  = mu1.eta();
                      bsRootTree_->Mu1Phi_beffit_  = mu1.phi();
                     // std::cout<<"mu1phi:"<<mu1.phi();
                      bsRootTree_->Mu2Phi_beffit_  = mu2.phi();
                      bsRootTree_->Mu2Pt_beffit_  = mu2.pt();
                      bsRootTree_->Mu2Pz_beffit_   = mu2.pz();
                      bsRootTree_->Mu2Eta_beffit_  = mu2.eta();
       }

      else 
      {
                            bsRootTree_->Mu2Pt_beffit_   = mu1.pt();
                            bsRootTree_->Mu2Pz_beffit_   = mu1.pz();
                            bsRootTree_->Mu2Eta_beffit_  = mu1.eta();
                            bsRootTree_->Mu2Phi_beffit_  = mu1.phi();
                            //std::cout<<"mu1phi"<<mu1.phi();
                            bsRootTree_->Mu1Phi_beffit_  = mu2.phi();
                            bsRootTree_->Mu1Pt_beffit_  = mu2.pt();
                            bsRootTree_->Mu1Pz_beffit_   = mu2.pz();
                            bsRootTree_->Mu1Eta_beffit_  = mu2.eta();
                     
      }

                bsRootTree_->JpsiMuMuDCA_beffit_  = MuonsDCA;
 		bsRootTree_->JpsiCosDeltaAlpha_  = CosAlpha;
                // std::cout<<"cosalpha"<<CosAlpha<<"\n";
 		bsRootTree_->JpsiLxySigma_   = lxyerr;
		bsRootTree_->JpsiLxy_  = lxy;
		bsRootTree_->JpsiLxyOverPt_ =	( displacementFromBeamspot.x() * Jpsi.px() + 	displacementFromBeamspot.y() * Jpsi.py() ) /	( Jpsi.pt()* Jpsi.pt()) ;
                //std::cout<<"jpsilxyoverpt: "<<( displacementFromBeamspot.x() * Jpsi.px() +    displacementFromBeamspot.y() * Jpsi.py() ) /    ( Jpsi.pt()* Jpsi.pt())<<"\n";
                 if( tv.isValid() )
                {
	          	 bsRootTree_->JpsiTrigVtxProb_  =  vtxProb_Jpsi;
                        // std::cout<<"prob:"<<vtxProb_Jpsi;
		}
  
                bool JpsiMu1Truth= MCmatchingJpsi( mu1, genParticles, bsRootTree_->JpsiMu1mcId_,bsRootTree_->JpsiMu1momId_, 443);
                bool JpsiMu2Truth= MCmatchingJpsi( mu2, genParticles, bsRootTree_->JpsiMu2mcId_,bsRootTree_->JpsiMu2momId_, 443);

                //std::cout << "Matching  " << bsRootTree_->isMatchedJpsi_<<  " JpsiCosDeltaAlpha " << bsRootTree_->JpsiCosDeltaAlpha_ << " JpsiLxySigma " << bsRootTree_->JpsiLxySigma_ << " JpsiTrigVtxProb_ " << bsRootTree_->JpsiTrigVtxProb_ << endl;
                if(JpsiMu1Truth==true && JpsiMu2Truth==true ) { bsRootTree_->isMatchedJpsi_ = 1; }
                if( mu1.isTrackerMuon() && muon::isGoodMuon(mu1, muon::TrackerMuonArbitrated) ) { bsRootTree_->JpsiMu1TrackerMuonArbitrated_ = 1 ;}		
		if( mu2.isTrackerMuon() && muon::isGoodMuon(mu2, muon::TrackerMuonArbitrated) ) { bsRootTree_->JpsiMu2TrackerMuonArbitrated_ = 1 ;}

//=============================MUON CHECK


      if (mu1.isTrackerMuon() && !mu1.isGlobalMuon())       bsRootTree_->JpsiMuon1Cat_alone_ = 1;
      else if (!mu1.isTrackerMuon() && mu1.isGlobalMuon())  bsRootTree_->JpsiMuon1Cat_alone_ = 2;
      else if (mu1.isTrackerMuon() && mu1.isGlobalMuon())   bsRootTree_->JpsiMuon1Cat_alone_ = 3;
      else if (!mu1.isTrackerMuon() && !mu1.isGlobalMuon()) bsRootTree_->JpsiMuon1Cat_alone_ = 4;
      if (mu1.isPFMuon())       bsRootTree_->JpsiMuonCat1_ = 1;

      if (mu2.isTrackerMuon() && !mu2.isGlobalMuon())       bsRootTree_->JpsiMuon2Cat_alone_ = 1;
      else if (!mu2.isTrackerMuon() && mu2.isGlobalMuon())  bsRootTree_->JpsiMuon2Cat_alone_ = 2;
      else if (mu2.isTrackerMuon() && mu2.isGlobalMuon())   bsRootTree_->JpsiMuon2Cat_alone_ = 3;
      else if (!mu2.isTrackerMuon() && !mu2.isGlobalMuon()) bsRootTree_->JpsiMuon2Cat_alone_ = 4;
      if (mu2.isPFMuon())       bsRootTree_->JpsiMuonCat2_ = 1;

     //numberOfHits changed to numberOfValidHits

      int pixhits1 = 0;
      const reco::HitPattern& pp1 = trkMu1Ref.get()->hitPattern();
      for (int iter=0; iter<pp1.numberOfAllHits(reco::HitPattern::TRACK_HITS); iter++) {
			uint32_t hit = pp1.getHitPattern(reco::HitPattern::TRACK_HITS,iter);
			if (pp1.validHitFilter(hit) && pp1.pixelBarrelHitFilter(hit)) pixhits1++;
			if (pp1.validHitFilter(hit) && pp1.pixelEndcapHitFilter(hit)) pixhits1++;
      }
      bsRootTree_->JpsiMu1nPixHits_alone_   = pixhits1;

      int pixhits2 = 0;
      const reco::HitPattern& pp2 = trkMu2Ref.get()->hitPattern();
      for (int iter=0; iter<pp2.numberOfAllHits(reco::HitPattern::TRACK_HITS); iter++) 
     {
			uint32_t hit = pp2.getHitPattern(reco::HitPattern::TRACK_HITS,iter);
                        //std::cout<<"hits"<<hit<<"\n";
			if (pp2.validHitFilter(hit) && pp2.pixelBarrelHitFilter(hit)) pixhits2++;
			if (pp2.validHitFilter(hit) && pp2.pixelEndcapHitFilter(hit)) pixhits2++;
      }
     
      bsRootTree_->JpsiMu2nPixHits_alone_   = pixhits2;
     


     //if (!tv.isValid()) continue;//skip if the mu-mu common vertex is not valid
     

      edm::Handle<View<pat::PackedCandidate>> allTracks;
      iEvent.getByToken(trackLabelK, allTracks); 
      //edm::Handle<edm::View<pat::IsolatedTrack>> allTracks;
      //iEvent.getByToken(isotrackTok, allTracks);  

     

//==================================================================================================================PHI RECONSTRUCTION
          //try{

            for (size_t k=0; k< allTracks->size(); ++k)
                  {
                           const pat::PackedCandidate & track1 = (*allTracks)[k];
                             if (!track1.hasTrackDetails())continue;
                             
                             if (track1.charge()<0)continue;                          
                             if (track1.pt() < KaonTrackPtCut_) continue;
                             if (track1.numberOfHits() < 5)continue;
                             if(!track1.trackHighPurity()) continue;
                             double DeltaRKaon1Jpsi = 999;
                             DeltaRKaon1Jpsi = deltaR(Jpsi.eta(), Jpsi.phi(), track1.eta(), track1.phi());
                             if (DeltaRKaon1Jpsi > 2.2) continue;
                             //std::cout<<" ********************DeltaRof jpsikaon1"<<DeltaRKaon1Jpsi<<"\n";
                             const reco::Track &  rtrk1 = (*allTracks)[k].pseudoTrack();
                             if (rtrk1.charge()<0) continue;
                             TransientTrack KPTT=(*theB).build(&rtrk1);
                             TrajectoryStateClosestToPoint KPTS = KPTT.impactPointTSCP();                            
                             if(!KPTS.isValid())continue;
                             if (!track1.clone()->hasTrackDetails())continue;
                             pat::PackedCandidate *track11 = track1.clone();


		     for (size_t l=k+1; l< allTracks->size(); ++l)
                      {
                             double kaonmass = 0.493677;
	                     double pionmass = 0.139570;           
	                     const pat::PackedCandidate & track2 = (*allTracks)[l];
	                     if ( !track2.hasTrackDetails() )continue;
                             if (track2.charge()>0) continue;                             
                             if (track2.pt() < KaonTrackPtCut_) continue;
                             if ( track2.numberOfHits()<5) continue;
                             if(!track2.trackHighPurity()) continue;
                             double DeltaRKaon2Jpsi = 999;
                             DeltaRKaon2Jpsi = deltaR(Jpsi.eta(), Jpsi.phi(), track2.eta(), track2.phi());
                           //  std::cout<<" ********************DeltaRof jpsikaon2"<<DeltaRKaon2Jpsi<<"\n";
                             if (DeltaRKaon2Jpsi >2.2) continue; 
                             const reco::Track &  rtrk2 = (*allTracks)[l].pseudoTrack();
                             if (rtrk2.charge()>0) continue;                            
                             TransientTrack KMTT=(*theB).build(&rtrk2);                           
                             TrajectoryStateClosestToPoint KMTS = KMTT.impactPointTSCP();                       
                             if(!KMTS.isValid())continue;
                             double KaonsDCA=999; 
                             if (KPTS.isValid() && KMTS.isValid()) {                            
                             ClosestApproachInRPhi cAppK;
                             cAppK.calculate(KPTS.theState(), KMTS.theState());
                             KaonsDCA=cAppK.distance();}
                            // std::cout<<"*****************************************///////////////kaonsDCA::"<<KaonsDCA<<"\n";
                             if(KaonsDCA > 0.5)continue;
                             if (!track2.clone()->hasTrackDetails())continue;
                             pat::PackedCandidate *track22 = track2.clone();
                             //bsRootTree_->KaonsDCA_   = KaonsDCA;            
                             if(bsRootTree_->iPassedCutIdent_ < 4 ) bsRootTree_->iPassedCutIdent_ = 4 ;
                             pat::CompositeCandidate PhiCand;
                             track11->setMass(kaonmass);
                             PhiCand.addDaughter(*track11);
                             track22->setMass(kaonmass);
                             PhiCand.addDaughter(*track22);
                             AddFourMomenta ad;
                             ad.set(PhiCand);
                             if (abs(PhiCand.mass()- nominalPhiMass) > PhiMassWindowBeforeFit_) continue;
                             std::cout<<"PhiCandMass"<<PhiCand.mass()<<"\n";
                              /* TLorentzVector phi4V;
                               kaontrack1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), kaonmass);
                               kaontrack2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), kaonmass);
                               phi4V =kaontrack1+kaontrack2;
                               if(phi4V.M()<0.950 || phi4V.M()>1.09) continue;
                               bsRootTree_->PhiMassFV_ = phi4V.M();
                                       
                               */
                             //std::cout<<"track1pt"<<track1.pt()<<"\n";
                             // std::cout<<"track2pt"<<track2.pt()<<"\n";
                             //std::cout<<"c***********"<<track1.hasTrackDetails()<<"\n";
	                     //std:std::cout<<"d***********"<<track2.hasTrackDetails()<<"\n";
	                     //std:std::cout<<"electron1**********"<<track1.isElectron()<<"\n";
                             //std:std::cout<<"electron2**********"<<track2.isElectron()<<"\n";
                             //std:std::cout<<"muon1**********"<<track1.isMuon()<<"\n";
                             //std:std::cout<<"muon2**********"<<track2.isMuon()<<"\n";
                             //std:std::cout<<"SAMuon1**********"<<track1.isStandAloneMuon()<<"\n";
	                     //	std:std::cout<<"SAMuon2**********"<<track2.isStandAloneMuon()<<"\n";
	                     //	std:std::cout<<"TMuon1**********"<<track1.isTrackerMuon()<<"\n";
	                     //	std:std::cout<<"TMuon2**********"<<track2.isTrackerMuon()<<"\n";
	                     //	std:std::cout<<"CMuon1**********"<<track1.isCaloMuon()<<"\n";
	                     //	std:std::cout<<"CMuon2**********"<<track2.isCaloMuon()<<"\n";
	                     //	std:std::cout<<"GMuon1**********"<<track1.isGlobalMuon()<<"\n";
	                     //	std:std::cout<<"GMuon2**********"<<track2.isGlobalMuon()<<"\n";
	                     //	std:std::cout<<"photon1**********"<<track1.isPhoton()<<"\n";
	                     //	std:std::cout<<"photon2**********"<<track2.isPhoton()<<"\n";
	                     //	std:std::cout<<"convertedphoton1**********"<<track1.isConvertedPhoton()<<"\n";
	                     //	std:std::cout<<"convertedphoton2**********"<<track2.isConvertedPhoton()<<"\n";
	                     //	std:std::cout<<"jet1**********"<<track1.isJet()<<"\n";
          		     // const reco::Candidate & track1 = (*allTracks)[k];
                             // const reco::Candidate & track2 = (*allTracks)[l];  	
                             //if (!track1 || !track2) continue;
                             //const reco::Track & rtrk1 = track1.pseudoTrack();
                             //const reco::Track & rtrk2 = track2.pseudoTrack(); 
                             //const reco::Candidate & rtrk1_ref = rtrk1;
                             //const reco::Candidate & rtrk2_ref = rtrk2; 
	  	             


                             //delete track11; delete track22;
                             //TrackRef trk2Ref = rtrk2
                             //TrackRef rtrk1Ref = rtrk1.get<TrackRef>();
                             //TrackRef rtrk2Ref = rtrk2.get<TrackRef>();                      
                             //select tracks under the assumpsion that they are kaons
                             //if (trk1Ref.isNull() || trk2Ref.isNull()){cout<<"Null trk  ref"<<endl; continue;}
                             if(track1.isElectron() >0){std::cout<<"electrontrack1:"<<track1.isElectron();}
                             if(track2.isElectron() >0){std::cout<<"electrontrack2:"<<track2.isElectron();}
                             if(track1.isMuon() >0){std::cout<<"Muontrack1:"<<track1.isMuon();}
                             if(track2.isMuon() >0){std::cout<<"Muontrack2:"<<track2.isMuon();}
                             if(track1.isPhoton() >0){std::cout<<"Photontrack1:"<<track1.isPhoton();}
                             if(track2.isPhoton() >0){std::cout<<"Photontrack2:"<<track2.isPhoton();}
                                                       
                              //met.setP4(reco::Candidate::LorentzVector(metPx, metPy, 0., TMath::Sqrt(metPx*metPx + metPy*metPy)));
                              // std::cout<<"ptetaphi1"<<track1.pt()<<track1.eta()<<track1.phi()<<"\n";
                              // std::cout<<"ptetaphi2"<<track2.pt()<<track2.eta()<<track2.phi()<<"\n";
       
       
       
                               //math::PtEtaPhiMLorentzVectorD* kaontrackI = new math::PtEtaPhiMLorentzVectorD(track1.pt(), track1.eta(), track1.phi(), kaonmass);
                               //math::PtEtaPhiMLorentzVectorD* kaontrackII = new math::PtEtaPhiMLorentzVectorD(track2.pt(), track2.eta(), track2.phi(), kaonmass);                       
       
                               // pat::PackedCandidate  &  kaoncand1 = kaoncand1.setP4(try1);
                                //pat::PackedCandidate  &  kaoncand2 = kaoncand2.setP4(kaontrackII);
                               //kaoncand1.setP4(kaontrackI);
                               //const Candidate  & kaoncand2.setP4(kaontrack2);
                               
                        
       	  	               if (abs(PhiCand.mass()- nominalPhiMass) > PhiMassWindowBeforeFit_) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 5 ) bsRootTree_->iPassedCutIdent_ = 5 ;
       	       	               bsRootTree_->PhiNumberOfCandidatesBeforeFit_++;
                                
       
                              // std::cout<<"alll:"<<4<<"\n";
                               if(bsRootTree_->iPassedCutIdent_   < 7 ) bsRootTree_->iPassedCutIdent_ = 7 ;
                               pat::CompositeCandidate BCand;
       	  		       BCand.addDaughter(mu1);
       	  		       BCand.addDaughter(mu2);
       	  	               BCand.addDaughter(*track11);
       	  		       BCand.addDaughter(*track22);
       	  		       AddFourMomenta add4mom;
       	  		       add4mom.set(BCand);
                               std::cout<<"mass B: "<<BCand.mass()<<"\n";
       	  		       if (BCand.mass() < BsLowerMassCutBeforeFit_ || BCand.mass() > BsUpperMassCutBeforeFit_) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 8 ) bsRootTree_->iPassedCutIdent_ = 8 ;
       
      /* 
       
                               reco::Track trk1Ref_best = rtrk1;//trk1Ref;
                               reco::Track trk2Ref_best = rtrk2;//trk2Ref;
                              trkMu1Ref_best=trkMu1Ref;
                              trkMu2Ref_best=trkMu2Ref;
                              mu1_best=mu1;
                              mu2_best=mu2;
       */
       	  		vector<TransientTrack> t_tracks;
       	  		t_tracks.push_back((*theB).build(&trkMu1Ref));
       	  		t_tracks.push_back((*theB).build(&trkMu2Ref));
       	  		t_tracks.push_back((*theB).build(&rtrk1));
       	  		t_tracks.push_back((*theB).build(&rtrk2));
       
       	  		if (!trkMu1Ref.isNonnull() || !trkMu2Ref.isNonnull() )continue;//|| !trk1Ref.isNonnull() || !trk2Ref.isNonnull()) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 9 ) bsRootTree_->iPassedCutIdent_ = 9 ;
       	  		       bsRootTree_->BsNumberOfCandidatesBeforeFit_++;
                               vector<TransientTrack> phi_tracks;
                               phi_tracks.push_back((*theB).build(&rtrk1));//trk1Ref
                               phi_tracks.push_back((*theB).build(&rtrk2));//trk2Ref
                               KalmanVertexFitter kvfphi;
                               TransientVertex tvphi = kvfphi.vertex(phi_tracks);
                               if (!tvphi.isValid()) continue;
       
                               if(verbose_==true) {
                                                     std::cout<<"****** KAONS HAVE VALID VERTEX FIT"<< std::endl;
                                                }
       			
                               Vertex vertexphi = tvphi;
                               double vtxProb_Phi = TMath::Prob(vertexphi.chi2(),(int)vertexphi.ndof());
                                //std::cout<<"phivtsprob"<<vtxProb_Phi<<"\n";
       
                        if(verbose_==true) 
                         {
                                     std::cout<<"Phi vertex probability:"<<vtxProb_Phi<< std::endl;
                        }

                        KalmanVertexFitter kvfbs;
                        TransientVertex kvfbsvertex = kvfbs.vertex(t_tracks);
                        Vertex vertexbskalman = kvfbsvertex;
                        if (!kvfbsvertex.isValid()) continue;
                        GlobalError gigi=kvfbsvertex.positionError();
                        // std::cout<<"error"<<gigi<<"\n";
                        bsRootTree_->K1Pt_beffit_   = track1.pt();
	 	 	bsRootTree_->K1Pz_beffit_   = track1.pz();
	  		bsRootTree_->K1Eta_beffit_  = track1.eta();
	  		bsRootTree_->K1Phi_beffit_  = track1.phi();
	 	 	bsRootTree_->K2Pt_beffit_   = track2.pt();
	  		bsRootTree_->K2Pz_beffit_   = track2.pz();
	  		bsRootTree_->K2Eta_beffit_  = track2.eta();
	  		bsRootTree_->K2Phi_beffit_  = track2.phi();

                        KinematicFitInterface Kfitter;
	  		bool fitSuccess = Kfitter.doFit(t_tracks, nominalMuonMass,  nominalKaonMass, nominalKaonMass);
                        if(fitSuccess != 1) continue;

                        if(bsRootTree_->iPassedCutIdent_   < 10 ) bsRootTree_->iPassedCutIdent_ = 10 ;
                        double vtxprob_Bs = TMath::Prob(vertexbskalman.chi2(),(int)vertexbskalman.ndof());
                       // std::cout<<"chi2for bs:"<<vertexbskalman.chi2()<<"\n";
                       // std::cout<<"nof for bs:"<<vertexbskalman.ndof()<<"\n";
                       // if(vtxprob_Bs >0.01){ std::cout<<"bsvtxprob:***********"<<"\t"<<vtxprob_Bs<<"\n";}
                 	RefCountedKinematicParticle bs = Kfitter.getParticle();
	  		RefCountedKinematicVertex bVertex = Kfitter.getVertex();
	  		AlgebraicVector7 b_par = bs->currentState().kinematicParameters().vector();
	  		AlgebraicSymMatrix77 bs_er = bs->currentState().kinematicParametersError().matrix();
                        
                        AlgebraicMatrix33 BVError(bVertex->error().matrix());
	  		double fittedBsMass = b_par[6];



                        if(!bVertex->vertexIsValid()) continue;


                        if(verbose_ == true){
                                                   std::cout<<"Good kinematic vertices"<<std::endl;
                                            }




                        TMatrix cova(2,2);
                        cova.IsSymmetric();
                        cova(0,0)=gigi.cxx();
                        cova(1,1)=gigi.cyy();
                        cova(0,1)=gigi.cyx();
                        cova(1,0)=gigi.cyx();


                      if(abs(Jpsi.mass() - nominalJpsiMass) < JpsiMassWindowAfterFit_  && Jpsi.pt() > JpsiPtCut_   && abs(PhiCand.mass() - nominalPhiMass) < PhiMassWindowAfterFit_    &&  fittedBsMass > BsLowerMassCutAfterFit_  &&  fittedBsMass < BsUpperMassCutAfterFit_ )
                         {
                                bsRootTree_->BsNumberOfCandidatesAfterFit_++;
                          }

                         if(vtxprob_Bs > minVtxP)
                         {

                              reco::Track trk1Ref_best = rtrk1;//trk1Ref;
                              reco::Track trk2Ref_best = rtrk2;//trk2Ref;
                              trkMu1Ref_best=trkMu1Ref;
                              trkMu2Ref_best=trkMu2Ref;
                              mu1_best=mu1;
                              mu2_best=mu2;



                               TrackRef mu1trkref = mu1.get<TrackRef>();
                               TrackRef mu2trkref = mu2.get<TrackRef>();
                               TrackRef K1trkRef = track1.get<TrackRef>();
                               TrackRef K2trkRef = track2.get<TrackRef>();

                               BsTrkRefs.clear();
                               BsTrkRefs.push_back(mu1trkref);
                               BsTrkRefs.push_back(mu2trkref);
                               BsTrkRefs.push_back(K1trkRef);
                               BsTrkRefs.push_back(K2trkRef);

                               if (abs(Jpsi.mass() - nominalJpsiMass) > JpsiMassWindowAfterFit_ || Jpsi.pt() < JpsiPtCut_) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 11 ) bsRootTree_->iPassedCutIdent_ = 11 ;
                               if (abs(PhiCand.mass() - nominalPhiMass) > PhiMassWindowAfterFit_) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 12 ) bsRootTree_->iPassedCutIdent_ = 12 ;
                               if (fittedBsMass < BsLowerMassCutAfterFit_ || fittedBsMass > BsUpperMassCutAfterFit_) continue;
                               if(bsRootTree_->iPassedCutIdent_   < 13 ) bsRootTree_->iPassedCutIdent_ = 13 ;
                               bsRootTree_->BsNumberOfCandidatesAfterBestFit_++;
                               minVtxP = vtxprob_Bs;
                               BCand_best = BCand;
   

                  
                               /*

                                 bsRootTree_->matchL11_=!mu1.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
                                 bsRootTree_->matchL12_=!mu2.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
                                 bsRootTree_->match2mu01_=!mu1.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
                                 bsRootTree_->match2mu02_=!mu2.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
                                 bsRootTree_->match1mu01_=!mu1.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
                                 bsRootTree_->match1mu02_=!mu2.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
                                 bsRootTree_->matchDoubleMu31J_ = !mu1.triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
                                 bsRootTree_->matchDoubleMu32J_ = !mu2.triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
                                 bsRootTree_->matchDoubleMu31Q_ = !mu1.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
                                 bsRootTree_->matchDoubleMu32Q_ = !mu2.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
                                 bsRootTree_->matchDoubleMu71_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();
                                 bsRootTree_->matchDoubleMu72_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();

                                 if(isMCstudy_){
	                             		
	                             		  
	                             			bsRootTree_->matchDoubleMu41_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
                                 	bsRootTree_->matchDoubleMu42_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
                                     bsRootTree_->matchDoubleMu41_v12_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4JpsiV4").empty();
                                 	bsRootTree_->matchDoubleMu42_v12_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4JpsiV4").empty();

	                             		}
                                  else{       					
                                     bsRootTree_->matchDoubleMu41_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
                                 	bsRootTree_->matchDoubleMu42_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
                              	}

                                 bsRootTree_->matchDoubleMu51_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
                                 bsRootTree_->matchDoubleMu52_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
                                 bsRootTree_->matchDoubleMu81_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon8Jpsi").empty();
                                 bsRootTree_->matchDoubleMu82_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon8Jpsi").empty();
	                             		cout << "Bs mu1 triggered with  dimuon8 " <<  bsRootTree_->matchDoubleMu81_  << endl;
	                             		cout << "Bs mu2  triggered with  dimuon8 " <<  bsRootTree_->matchDoubleMu82_  << endl;
	                             	// trigger matching for HLT_Dimuon0_Jpsi
	                             	
	                             	                         bsRootTree_->matchDoubleMu01_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
	                             	                         bsRootTree_->matchDoubleMu02_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();


                                 bsRootTree_->matchDoubleMu01_JpsiMuon_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
                                 bsRootTree_->matchDoubleMu02_JpsiMuon_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
                                 bsRootTree_->matchDoubleMu101_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
                                 bsRootTree_->matchDoubleMu102_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
                                 bsRootTree_->matchDoubleMu131_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();
                                 bsRootTree_->matchDoubleMu132_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();


//====================================================================================================================================

                                 bool matchedMu = false, matchedTrack = false;

                                 pat::TriggerObjectStandAloneCollection  mu0tkmuMatch = mu1.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
                                 for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
                                   if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;

                                 if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;

                                 }

                                 if (matchedMu) bsRootTree_->match2mu31_=1;
                                 else if (matchedTrack) bsRootTree_->match2mu31_=1;
                                 else bsRootTree_->match2mu31_=0;

                                 matchedMu = false; matchedTrack = false;
                                 mu0tkmuMatch = mu2.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
                                 for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
                                   if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
                                   if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
                                 }
                                 if (matchedMu) bsRootTree_->match2mu32_=1;
                                 else if (matchedTrack) bsRootTree_->match2mu32_=1;
                                 else bsRootTree_->match2mu32_=0;
//=====================================================================================
                                 matchedMu = false, matchedTrack = false;
                                 mu0tkmuMatch = mu1.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFilteredTight");
                                 for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
                                   if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
                                   if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
                                 }
                                 if (matchedMu) bsRootTree_->matchmu0tk01_=1;
                                 else if (matchedTrack) bsRootTree_->matchmu0tk01_=1;
                                 else bsRootTree_->matchmu0tk01_=0;

                                 matchedMu = false, matchedTrack = false;
                                 mu0tkmuMatch = mu2.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFilteredTight");
                                 for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
                                   if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
                                   if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
                                 }
                                 if (matchedMu) bsRootTree_->matchmu0tk02_=1;
                                 else if (matchedTrack) bsRootTree_->matchmu0tk02_=1;
                                 else bsRootTree_->matchmu0tk02_=0;
                               */
//================================================================================================================================
                             reco::Vertex reFitVertex;
                                 
                                 bsRootTree_->PVx_refit_ = reFitVertex.x();
                                 bsRootTree_->PVy_refit_ = reFitVertex.y();
                                 bsRootTree_->PVz_refit_ = reFitVertex.z();
                     
                                 bsRootTree_->PVerrx_refit_ = reFitVertex.xError();
                                 bsRootTree_->PVerry_refit_ = reFitVertex.yError();
                                 bsRootTree_->PVerrz_refit_ = reFitVertex.zError();
                                 bsRootTree_->BsFitChi2_  = bs->chiSquared();
                                 bsRootTree_->BsFitNdof_   =(int)bs->degreesOfFreedom();
                                 bsRootTree_->BsFitVtxProb_ = vtxprob_Bs;
                                 bsRootTree_->BsPhiVtxProb_ = vtxProb_Phi;
                                 //Save Jpsi Vertex probability
                                 bsRootTree_->JpsiVtxProb_ = vtxProb_Jpsi;
                                 bsRootTree_->CosDeltaAlpha_ = CosAlpha;
                                 bsRootTree_->MuMuDCA_ = MuonsDCA;
                                 bsRootTree_->MuMuDistance_ = lxy;
                                 bsRootTree_->MuMuDistanceSigma_ = lxyerr;
                                 bsRootTree_->MuDr1_ = max_Dr1;
                                 bsRootTree_->MuDr2_ = max_Dr2;
                                 //momentum vector for bs
                                 GlobalVector Bsvec(b_par[3], b_par[4], b_par[5]); // the fitted momentum vector of the Bs
                                 bsRootTree_->BsFitM_ = fittedBsMass;
                                 bsRootTree_->BsFitEta_ = Bsvec.eta();
                                 bsRootTree_->BsFitPt_  = Bsvec.perp();
                                 bsRootTree_->BsFitPz_  = Bsvec.z();
                                 bsRootTree_->BsFitPhi_ = Bsvec.phi();
                     
                                 TheBs.SetPtEtaPhiM(Bsvec.perp(),Bsvec.eta(),Bsvec.phi(),fittedBsMass);
                     
                                 RefCountedKinematicTree reftree = Kfitter.getTree();
                                 setFitParKK(reftree);
                                 RefCountedKinematicTree jpsitree = Kfitter.getJpsiTree();
                     
                     
                     
                                 vector< RefCountedKinematicParticle > bs_children = reftree->finalStateParticles();
                                 AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
                                 AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();
                     
                     
                                 vector< RefCountedKinematicParticle > jpsi_children = jpsitree->finalStateParticles();
                                 AlgebraicVector7 bs_par3 = jpsi_children[0]->currentState().kinematicParameters().vector();
                                 AlgebraicVector7 bs_par4 = jpsi_children[1]->currentState().kinematicParameters().vector();
                                 double pt1 = sqrt(bs_par1[3]*bs_par1[3]+bs_par1[4]*bs_par1[4]);
                                 double pt2 = sqrt(bs_par2[3]*bs_par2[3]+bs_par2[4]*bs_par2[4]);
                                 bsRootTree_->K1Pt_fit_ = pt1;
                                 bsRootTree_->K2Pt_fit_ = pt2;
                     
                     
                     
                     
                                 TLorentzVector pK1;
                                 double en1 = sqrt(bs_par1[3]*bs_par1[3]+bs_par1[4]*bs_par1[4]+bs_par1[5]*bs_par1[5]+bs_par1[6]*bs_par1[6]);
                                 pK1.SetPxPyPzE(bs_par1[3],bs_par1[4],bs_par1[5],en1);
                                 TLorentzVector pK2;
                                 double en2 = sqrt(bs_par2[3]*bs_par2[3]+bs_par2[4]*bs_par2[4]+bs_par2[5]*bs_par2[5]+bs_par2[6]*bs_par2[6]);
                                 pK2.SetPxPyPzE(bs_par2[3],bs_par2[4],bs_par2[5],en2);
                                 TLorentzVector PhiFit = pK1 + pK2;
                                 bsRootTree_->PhiM_fit_ = PhiFit.M();
                     
                     
                     
                                                      /*
                                                      if (mu1.isGlobalMuon()) bsRootTree_->BsMu1QualityG_=selGlobalMuon(mu1,RecVtx.position());
                                                      if (mu2.isGlobalMuon()) bsRootTree_->BsMu2QualityG_=selGlobalMuon(mu2,RecVtx.position());
                                                      if (mu1.isTrackerMuon()) bsRootTree_->BsMu1QualityT_=selTrackerMuon(mu1,RecVtx.position());
                                                      if (mu2.isTrackerMuon()) bsRootTree_->BsMu2QualityT_=selTrackerMuon(mu2,RecVtx.position());
                                                      */
                                                      bsRootTree_->BsFitVtx_x_ = bVertex->position().x();
                                                      bsRootTree_->BsFitVtx_y_ = bVertex->position().y();
                                                      bsRootTree_->BsFitVtx_z_ = bVertex->position().z();
                                                     //  bsRootTree_->KaonsDCA_   = KaonsDCA;            
                                                      bsRootTree_->BsM_nofit_ = BCand.mass();
                                                      bsRootTree_->BsPt_nofit_ = BCand.pt();
                                                      bsRootTree_->BsPz_nofit_ = BCand.pz();
                                                      bsRootTree_->BsPhi_nofit_ = BCand.phi();
                                                      bsRootTree_->BsEta_nofit_ = BCand.eta();
                                                      bsRootTree_->JpsiM_nofit_ = Jpsi.mass();
                                                      bsRootTree_->JpsiPhi_nofit_ = Jpsi.phi();
                                                      bsRootTree_->JpsiEta_nofit_ = Jpsi.eta();
                                                      bsRootTree_->JpsiPt_nofit_ = Jpsi.pt();
                                                      bsRootTree_->JpsiPz_nofit_ = Jpsi.pz();
                                                       bsRootTree_->deltaRkaon1Jpsi_ = DeltaRKaon1Jpsi;
                                                       bsRootTree_->deltaRkaon2Jpsi_ = DeltaRKaon2Jpsi;
 
                                                      bsRootTree_->PhiM_nofit_ = PhiCand.mass();
                                                      bsRootTree_->PhiPhi_nofit_ = PhiCand.phi();
                                                      bsRootTree_->PhiEta_nofit_ = PhiCand.eta();
                                                      bsRootTree_->PhiPt_nofit_ = PhiCand.pt();
                                                      bsRootTree_->PhiPz_nofit_ = PhiCand.pz();
                                          
                                                      bsRootTree_->K1Pt_nofit_   = track1.pt();
                                                      bsRootTree_->K1Pz_nofit_   = track1.pz();
                                                      bsRootTree_->K1Eta_nofit_  = track1.eta();
                                                      bsRootTree_->K1Phi_nofit_  = track1.phi();
                                                      //bsRootTree_->K1Key_nofit_  = trk1Ref.key();
                                                      bsRootTree_->K2Pt_nofit_   = track2.pt();
                                                      bsRootTree_->K2Pz_nofit_   = track2.pz();
                                                      bsRootTree_->K2Eta_nofit_  = track2.eta();
                                                      bsRootTree_->K2Phi_nofit_  = track2.phi();
                                                      //bsRootTree_->K2Key_nofit_  = trk2Ref.key();
                                          
                                          
                                          
                                          
                                          
                                          				//bsRootTree_->K1HighPurityTrack_ = trk1Ref->quality(reco::TrackBase::highPurity);
                                          				//bsRootTree_->K2HighPurityTrack_ = trk2Ref->quality(reco::TrackBase::highPurity);
                                          
                                          				//cout << "trk1Ref->quality(reco::TrackBase::highPurity) " << trk1Ref->quality(reco::TrackBase::highPurity) << endl; 
                                          				//cout << "trk2Ref->quality(reco::TrackBase::highPurity) " << trk2Ref->quality(reco::TrackBase::highPurity) << endl; 
                                          
                                                      //bsRootTree_->K1Chi2_ = trk1Ref.get()->normalizedChi2();
                                                      //bsRootTree_->K1nHits_= trk1Ref.get()->numberOfValidHits();
                                                      //bsRootTree_->K2Chi2_ = trk2Ref.get()->normalizedChi2();
                                                      //bsRootTree_->K2nHits_= trk2Ref.get()->numberOfValidHits();
                                                      bsRootTree_->Mu1Chi2_ = trkMu1Ref.get()->normalizedChi2();
                                                      bsRootTree_->Mu1nHits_= trkMu1Ref.get()->numberOfValidHits();
                                                      bsRootTree_->Mu2Chi2_ = trkMu2Ref.get()->normalizedChi2();
                                                      bsRootTree_->Mu2nHits_ =trkMu2Ref.get()->numberOfValidHits();
                                          
                                                      bsRootTree_->BsCowboy_=isCowboy;
                                                      bsRootTree_->Mu1d0_ = trkMu1Ref->d0();
                                                      bsRootTree_->Mu2d0_ = trkMu1Ref->d0();
                                                      bsRootTree_->Mu1dz_ = trkMu1Ref->dz();
                                                      bsRootTree_->Mu2dz_ = trkMu1Ref->dz();
                                          

				TLorentzVector Mu1, Mu2;
 				double MuMass =  0.1056583715;

 			 	Mu1.SetPtEtaPhiM( mu1.pt(), mu1.eta(), mu1.phi(), MuMass);
 				Mu2.SetPtEtaPhiM( mu2.pt(), mu2.eta(), mu2.phi(), MuMass);
  
 			 	bsRootTree_->MuonPairDR_ = Mu1.DeltaR(Mu2) ;

				if( mu1.innerTrack()->quality(reco::TrackBase::highPurity) )	{ bsRootTree_->Mu1InnerTrkHighQuality_ = 1 ;}
 	 			if( muon::isGoodMuon(mu1, muon::TMOneStationTight) )	{ bsRootTree_->Mu1isGood_ = 1 ;}
		
	 			bsRootTree_->Mu1TrkBSDxy_ = mu1.innerTrack()->dxy(  vertexBeamSpot->position() );
	 			bsRootTree_->Mu1TrkBSDz_ = mu1.innerTrack()->dz(  vertexBeamSpot->position() ) ;
	 			bsRootTree_->Mu1PixelHits_ = mu1.innerTrack()->hitPattern().pixelLayersWithMeasurement() ;
	 			bsRootTree_->Mu1TrackerHits_ = mu1.innerTrack()->hitPattern().trackerLayersWithMeasurement() ;
	   
				if( mu2.innerTrack()->quality(reco::TrackBase::highPurity) )	{ bsRootTree_->Mu2InnerTrkHighQuality_ = 1 ;}
 				if( muon::isGoodMuon(mu2, muon::TMOneStationTight) )	{ bsRootTree_->Mu2isGood_ = 1 ;}
	   	        	bsRootTree_->Mu2TrkBSDxy_ = mu2.innerTrack()->dxy(  vertexBeamSpot->position() );
	   	        	bsRootTree_->Mu2TrkBSDz_ = mu2.innerTrack()->dz(  vertexBeamSpot->position() ) ;
	   	        	bsRootTree_->Mu2PixelHits_ = mu2.innerTrack()->hitPattern().pixelLayersWithMeasurement() ;
	   	        	bsRootTree_->Mu2TrackerHits_ = mu2.innerTrack()->hitPattern().trackerLayersWithMeasurement();
                   	        if(mu1.charge() == -1)
                                {
                   
                                    bsRootTree_->JpsiMu1Pt_alone_ = mu1.pt();
                                    bsRootTree_->JpsiMu2Pt_alone_ = mu2.pt();
                                    bsRootTree_->JpsiMu1Phi_alone_ = mu1.phi();
                                    bsRootTree_->JpsiMu2Phi_alone_ = mu2.phi();
                                    bsRootTree_->JpsiMu1Eta_alone_ = mu1.eta();
                                    bsRootTree_->JpsiMu2Eta_alone_ = mu2.eta();
                                } 
                                  else
                                { 
                   					
                                    bsRootTree_->JpsiMu1Pt_alone_ = mu2.pt();
                                    bsRootTree_->JpsiMu2Pt_alone_ = mu1.pt();
                                    bsRootTree_->JpsiMu1Phi_alone_ = mu2.phi();
                                    bsRootTree_->JpsiMu2Phi_alone_ = mu1.phi();
                                    bsRootTree_->JpsiMu1Eta_alone_ = mu2.eta();
                                    bsRootTree_->JpsiMu2Eta_alone_ = mu1.eta();
                       
                               }
                   
                   
                               bsRootTree_->JpsiMu1d0_alone_ = trkMu1Ref->d0();
                               bsRootTree_->JpsiMu2d0_alone_ = trkMu1Ref->d0();
                               bsRootTree_->JpsiMu1dz_alone_ = trkMu1Ref->dz();
                               bsRootTree_->JpsiMu2dz_alone_ = trkMu1Ref->dz();
                               bsRootTree_->JpsiMu1chi2_alone_ = trkMu1Ref->chi2();
                               bsRootTree_->JpsiMu2chi2_alone_ = trkMu1Ref->chi2();
                               bsRootTree_->JpsiMu1ndof_alone_ = trkMu1Ref->ndof();
                               bsRootTree_->JpsiMu2ndof_alone_ = trkMu1Ref->ndof();
                               bsRootTree_->JpsiMu1nHits_alone_ =  trkMu1Ref->numberOfValidHits();
                               bsRootTree_->JpsiMu2nHits_alone_ =  trkMu2Ref->numberOfValidHits();
                   
                   
                   
                               
                                Vertex PVvtxHightestPt = (*recVtxs)[VtxIndex];//reVertex(iSetup, vertexBeamSpot, RecVtx, mu1, mu2, trk1Ref, trk2Ref);
                                bsRootTree_->PVx_refit_ = PVvtxHightestPt.x();std::cout<<"Primary vertex HightestPt"<<PVvtxHightestPt.x()<<"\n";
                                bsRootTree_->PVy_refit_ = PVvtxHightestPt.y();
                                bsRootTree_->PVz_refit_ = PVvtxHightestPt.z();

                               //GlobalVector Bsvec(b_par[3], b_par[4], b_par[5]);

                              bsRootTree_->BsCt3D_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxHightestPt.x())*Bsvec.x()+
                                                             (kvfbsvertex.position().y()-PVvtxHightestPt.y())*Bsvec.y()+
                                                             (kvfbsvertex.position().z()-PVvtxHightestPt.z())*Bsvec.z())/
                                                             (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());


                              bsRootTree_->BsCt2D_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxHightestPt.x())*Bsvec.x()+
                                                             (kvfbsvertex.position().y()-PVvtxHightestPt.y())*Bsvec.y())/
                                                             (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());
                              
                           


                           bsRootTree_->BsCt2DBS_ = BsPDGMass_*( (kvfbsvertex.position().x()-BSx)*Bsvec.x()+
                                                               (kvfbsvertex.position().y()-BSy)*Bsvec.y())/
                                                               (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());


                            Int_t PVCosThetaIndex = -1;
                            double MinDistance = 10000000;
	                        		double MinDistanceZ = 10000000;
	                        		Int_t PVClosestZIndex = -1;
                            if(verbose_ == true){
                              std::cout<<"Try to find the min Costheta PV"<<std::endl;
                            }


                           vector<TransientTrack> BsJpsitrks;
                           TrackRef Jpsimu1trkref = mu1.get<TrackRef>();
                           TrackRef Jpsimu2trkref = mu2.get<TrackRef>();

      		           TransientTrack BsJpsimu1TT= (*theB).build(&Jpsimu1trkref);
			   TransientTrack BsJpsimu2TT= (*theB).build(&Jpsimu2trkref);
    
      	           	BsJpsitrks.push_back(BsJpsimu1TT);
      	         	BsJpsitrks.push_back(BsJpsimu2TT);

      	         	KalmanVertexFitter BsJpsiVtxkvf(true);
      		        TransientVertex BsJpsiVtx = BsJpsiVtxkvf.vertex(BsJpsitrks);
 		
				
				

			 	bsRootTree_->BsSVx_ = kvfbsvertex.position().x();
			 	bsRootTree_->BsSVy_ = kvfbsvertex.position().y();
			 	bsRootTree_->BsSVz_ = kvfbsvertex.position().z();


                    for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++){
					 
                                        const Vertex &vtx = (*recVtxs)[VtxInd];
					if(  MinDistanceZ > abs( kvfbsvertex.position().z()- vtx.z() )  ){
					MinDistanceZ = abs( kvfbsvertex.position().z()- vtx.z() );
                                        //cout << "*****************MinDistance in z " <<   MinDistanceZ << " **************vtx.z "<< vtx.z() << endl;
					PVClosestZIndex = VtxInd;
				}		
			
Double_t PVSVvecDotBsPvec=(kvfbsvertex.position().x()-vtx.x())*Bsvec.x()+(kvfbsvertex.position().y()-vtx.y())*Bsvec.y()+(kvfbsvertex.position().z()-vtx.z())*Bsvec.z();
Double_t PVSVlength = TMath::Sqrt( pow((kvfbsvertex.position().x()- vtx.x()), 2.0) + pow((kvfbsvertex.position().y()-vtx.y()), 2.0) + pow((kvfbsvertex.position().z()- vtx.z()), 2.0) );

              Double_t BsPlength = TMath::Sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z());
              Double_t BsCosTheta = PVSVvecDotBsPvec / (BsPlength * PVSVlength);

              Double_t distance = 1-BsCosTheta;

 
              if(distance < MinDistance){
                 MinDistance = distance;
                 cout << "mindist***************************** "<< MinDistance << endl;
                 PVCosThetaIndex = VtxInd;
                 std::cout<<"costhetaindex"<<PVCosThetaIndex<<"\n";
                }
            }

            if(verbose_ == true){std::cout<<"Found the min Costheta PV "<<Bsvec<<std::endl;}
            if (PVCosThetaIndex==-1) continue;
	    if (PVClosestZIndex==-1) continue;

            BsPVVtxInd = PVCosThetaIndex;
            const Vertex &Bsvtx = (*recVtxs)[BsPVVtxInd];
            bsRootTree_->TrackMultiplicity_ = Bsvtx.nTracks();





             
            PVvtxCosTheta = (*recVtxs)[PVCosThetaIndex];//reVertex(iSetup, vertexBeamSpot,(*recVtxs)[PVCosThetaIndex], mu1, mu2, trk1Ref, trk2Ref);
            Vertex PVvtxClosestZ = (*recVtxs)[PVClosestZIndex];//reVertex(iSetup, vertexBeamSpot, (*recVtxs)[PVClosestZIndex], mu1, mu2, trk1Ref, trk2Ref);		
            

            bsRootTree_->PVx_refit_cosTheta_ =  PVvtxCosTheta.x();std::cout<<"costhetaX"<<PVvtxCosTheta.x()<<"\n";
            bsRootTree_->PVy_refit_cosTheta_ =  PVvtxCosTheta.y();
            bsRootTree_->PVz_refit_cosTheta_ =  PVvtxCosTheta.z();
	    bsRootTree_->PVx_refit_closestZ_ =  PVvtxClosestZ.x();std::cout<<"costhetaX"<<PVvtxClosestZ.x()<<"\n";
            bsRootTree_->PVy_refit_closestZ_ =  PVvtxClosestZ.y();
            bsRootTree_->PVz_refit_closestZ_ =  PVvtxClosestZ.z();
	    bsRootTree_->BsSoftMuon1_=mu1.isSoftMuon(PVvtxCosTheta);
	    bsRootTree_->BsSoftMuon2_=mu2.isSoftMuon(PVvtxCosTheta);


bsRootTree_->BsCt3DPVClosestZ_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxClosestZ.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxClosestZ.y())*Bsvec.y() + (kvfbsvertex.position().z()-PVvtxClosestZ.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );

	bsRootTree_->BsCt2DPVClosestZ_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxClosestZ.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxClosestZ.y())*Bsvec.y()  )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y()  );
				
				 //3D ctau with pointing angle PV selection

            bsRootTree_->BsCt3DPVCosTheta_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*Bsvec.y() + (kvfbsvertex.position().z()-PVvtxCosTheta.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );

           
				// 2D ctau with pointing angle PV selection
        bsRootTree_->BsCt2DPVCosTheta_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() );

     	bsRootTree_->BsCt2DPVCosThetaOld_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*BCand.px() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*BCand.py() )/( BCand.px()*BCand.px() + BCand.py()*BCand.py() );
	bsRootTree_->BsCt2DBSOld_ = BsPDGMass_*( (kvfbsvertex.position().x()-BSx)*BCand.px() + (kvfbsvertex.position().y()-BSy)*BCand.py() )/( BCand.px()*BCand.px() + BCand.py()*BCand.py() );
	bsRootTree_->BsCt2DPVClosestZOld_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxClosestZ.x())*BCand.px() + (kvfbsvertex.position().y()-PVvtxClosestZ.y())*BCand.py() )/( BCand.px()*BCand.px() + BCand.py()*BCand.py() );
	bsRootTree_->BsCt2DOld_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxHightestPt.x())*BCand.px() + (kvfbsvertex.position().y()-PVvtxHightestPt.y())*BCand.py() )/( BCand.px()*BCand.px() + BCand.py()*BCand.py() );

                   cout << "Ct 3D cosTheta " << bsRootTree_->BsCt3DPVCosTheta_ << endl;//cout to CO
                   bsRootTree_->BsCtErr2DBS_ = CalculateCtErrbeamspot(iSetup, vertexBeamSpot, kvfbsvertex, cova, Bsvec, BsPDGMass_) ;
           
				bsRootTree_->BsCtErr2D_ = CalculateCtErrvertex(iSetup,PVvtxHightestPt,kvfbsvertex,cova,Bsvec,BsPDGMass_); 
				bsRootTree_->BsCtErr2DClosestZ_ = CalculateCtErrvertex(iSetup, PVvtxClosestZ,kvfbsvertex,cova,Bsvec,BsPDGMass_); 
                                bsRootTree_->BsCtErr2DCostheta_ = CalculateCtErrvertex(iSetup, PVvtxCosTheta,kvfbsvertex,cova,Bsvec,BsPDGMass_); 
                                cout<<" ct ErrCostheta"<<  bsRootTree_->BsCtErr2DCostheta_ <<"\n";
				if( BsJpsiVtx.isValid() ){
			 	  bsRootTree_->JpsiSVx_ = BsJpsiVtx.position().x();
			 	  bsRootTree_->JpsiSVy_ = BsJpsiVtx.position().y();
			 	  bsRootTree_->JpsiSVz_ = BsJpsiVtx.position().z();
				  GlobalError BsJpsiVtxUnc = BsJpsiVtx.positionError();
         	  TMatrix JpsiCova(2,2);
         	  JpsiCova.IsSymmetric();
         	  JpsiCova(0,0)= BsJpsiVtxUnc.cxx();
         	  JpsiCova(1,1)= BsJpsiVtxUnc.cyy();
         	  JpsiCova(0,1)= BsJpsiVtxUnc.cyx();
         	  JpsiCova(1,0)= BsJpsiVtxUnc.cyx();
					bsRootTree_->BsCtErr2DBS_JpsiVtx_ = CalculateCtErrbeamspot(iSetup,vertexBeamSpot,BsJpsiVtx,JpsiCova,Bsvec,BsPDGMass_) ;
				   bsRootTree_->BsCtErr2D_JpsiVtx_ = CalculateCtErrvertex(iSetup,PVvtxHightestPt,BsJpsiVtx,JpsiCova,Bsvec,BsPDGMass_); 
					bsRootTree_->BsCtErr2DClosestZ_JpsiVtx_ = CalculateCtErrvertex(iSetup, PVvtxClosestZ,BsJpsiVtx,JpsiCova,Bsvec,BsPDGMass_); 
            	bsRootTree_->BsCtErr2DCostheta_JpsiVtx_ = CalculateCtErrvertex(iSetup, PVvtxCosTheta,BsJpsiVtx,JpsiCova,Bsvec,BsPDGMass_); 
				}
 cout << "Cterr Jpsi vtx " <<  bsRootTree_->BsCtErr2DCostheta_JpsiVtx_ << endl;

            GlobalVector BsVecNonFitted(BCand.px(), BCand.py(), BCand.pz()); 
            bsRootTree_->BsCtErr2DCosthetaOld_ = CalculateCtErrvertex(iSetup,PVvtxCosTheta,kvfbsvertex,cova,BsVecNonFitted,BsPDGMass_); 
            bsRootTree_->BsCtErr2DBSOld_ = CalculateCtErrbeamspot(iSetup,vertexBeamSpot,kvfbsvertex,cova,BsVecNonFitted,BsPDGMass_); 
            bsRootTree_->BsCtErr2DClosestZOld_ = CalculateCtErrvertex(iSetup,PVvtxClosestZ,kvfbsvertex,cova,BsVecNonFitted,BsPDGMass_); 
            bsRootTree_->BsCtErr2DOld_ = CalculateCtErrvertex(iSetup,PVvtxHightestPt,kvfbsvertex,cova,BsVecNonFitted,BsPDGMass_);



            AlgebraicMatrix31 pB;
            pB(0,0) = bs->currentState().globalMomentum().x();
            pB(1,0) = bs->currentState().globalMomentum().y();
            pB(2,0) = bs->currentState().globalMomentum().z();



            AlgebraicMatrix13 pBT;
            pBT(0,0) = bs->currentState().globalMomentum().x();
            pBT(0,1) = bs->currentState().globalMomentum().y();
            pBT(0,2) = bs->currentState().globalMomentum().z();

            AlgebraicMatrix31 PV;
            PV(0,0) = PVx;
            PV(0,1) = PVy;
            PV(0,2) = PVz;
            
            AlgebraicMatrix31 BV;
            BV(0,0) = bVertex->position().x();
            BV(0,1) = bVertex->position().y();
            BV(0,2) = bVertex->position().z();
            AlgebraicMatrix31 lxyz = BV-PV;
            AlgebraicMatrix33 PVError(RecVtx.error());
            AlgebraicMatrix33 BVError(bVertex->error().matrix());
            AlgebraicMatrix33 lxyzError = PVError + BVError;
            lxyzError.Invert();

            AlgebraicMatrix11 a = pBT * lxyzError * pB ;
            AlgebraicMatrix11 b = pBT * lxyzError * lxyz;
            double num(b(0,0));
            double deno(a(0,0));
            bsRootTree_->BsCtMPV_ = (num*bs->currentState().mass())/(deno);

//error on ctau 3D


            GlobalPoint SVpos( bVertex->position().x(), bVertex->position().y(), bVertex->position().z());
            GlobalPoint PVpos( PVx, PVy, PVz);
            GlobalError SVerr( bVertex->error() );
            GlobalError PVerr( RecVtx.error() );
            VertexDistance3D d1;
            Measurement1D measurement = d1.distance(VertexState(SVpos,SVerr),VertexState(PVpos,PVerr));
            double error3D = measurement.error();
            double scale1 = ((bVertex->position().x() - PVx)*Bsvec.x()+
            (bVertex->position().y() - PVy)*Bsvec.y()+
            (bVertex->position().z() - PVz)*Bsvec.z())/
            (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z())*
            sqrt((bVertex->position().x() - PVx)*(bVertex->position().x() - PVx)+
            (bVertex->position().y() - PVy)*(bVertex->position().y() - PVy)+
            (bVertex->position().z() - PVz)*(bVertex->position().z() - PVz)));
           bsRootTree_->BsCtErr3D_ = BsPDGMass_*(error3D*abs(scale1))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());

              // error on ctau 2D - 2 (first approximation)
 
            bsRootTree_->BsCtErr2D2_ = sqrt((1/(Bsvec.perp()*Bsvec.perp()))*
            (bs_er(1,1)*Bsvec.x()*Bsvec.x()+
            bs_er(2,2)*Bsvec.y()*Bsvec.y()+
            bs_er(1,2)*Bsvec.x()*Bsvec.y()));


          //error on ctau 3D MPV
            bsRootTree_->BsCtErrMPV_ = BsPDGMass_/sqrt(deno);


           //reco::Vertex reFitVertex;
             //proper decay time and proper decay length with the refitted vertex
            // ctau 3D

            bsRootTree_->BsCt3Drefit_ = BsPDGMass_*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
            +(bVertex->position().y()-reFitVertex.y())*Bsvec.y()+
            +(bVertex->position().z()-reFitVertex.z())*Bsvec.z())/
            (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());

            // ctau 2d
            bsRootTree_->BsCt2Drefit_ = BsPDGMass_*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
            +(bVertex->position().y()-reFitVertex.y())*Bsvec.y())/
            (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());



            //ctau 3D MPV

            AlgebraicMatrix31 pB2;

            pB2(0,0) = bs->currentState().globalMomentum().x();
            pB2(1,0) = bs->currentState().globalMomentum().y();
            pB2(2,0) = bs->currentState().globalMomentum().z();


            AlgebraicMatrix13 pBT2;

            pBT2(0,0) = bs->currentState().globalMomentum().x();
            pBT2(0,1) = bs->currentState().globalMomentum().y();
            pBT2(0,2) = bs->currentState().globalMomentum().z();

            AlgebraicMatrix31 PV2;
            PV2(0,0) = reFitVertex.x();
            PV2(0,1) = reFitVertex.y();
            PV2(0,2) = reFitVertex.z();
            AlgebraicMatrix31 BV2;
            BV2(0,0) = bVertex->position().x();
            BV2(0,1) = bVertex->position().y();
            BV2(0,2) = bVertex->position().z();
            AlgebraicMatrix31 lxyz2 = BV2-PV2;
            AlgebraicMatrix33 PVError2(reFitVertex.error());
            AlgebraicMatrix33 BVError2(bVertex->error().matrix());
            AlgebraicMatrix33 lxyzError2 = PVError2 + BVError2;
            lxyzError2.Invert();

            AlgebraicMatrix11 a2 = pBT2 * lxyzError2 * pB2 ;
            AlgebraicMatrix11 b2 = pBT2 * lxyzError2 * lxyz2;
            double num2(b2(0,0));
            double deno2(a2(0,0));
            bsRootTree_->BsCtMPVrefit_ = (num2*bs->currentState().mass())/(deno2);


            //Error on ctau 3D


            GlobalPoint SVpos2( bVertex->position().x(), bVertex->position().y(), bVertex->position().z());
            GlobalPoint PVpos2( reFitVertex.x(), reFitVertex.y(), reFitVertex.z());
            GlobalError SVerr2( bVertex->error() );
            GlobalError PVerr2( reFitVertex.error() );
            VertexDistance3D d12;
            Measurement1D measurement12 = d12.distance(VertexState(SVpos2,SVerr2),VertexState(PVpos2,PVerr2));
            double error3D2 = measurement12.error();
            double scale12 = ((bVertex->position().x() - reFitVertex.x())*Bsvec.x()+
                             (bVertex->position().y() - reFitVertex.y())*Bsvec.y()+
                             (bVertex->position().z() - reFitVertex.z())*Bsvec.z())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z())*
               sqrt((bVertex->position().x() - reFitVertex.x())*(bVertex->position().x() - reFitVertex.x())+
                    (bVertex->position().y() - reFitVertex.y())*(bVertex->position().y() - reFitVertex.y())+
                    (bVertex->position().z() - reFitVertex.z())*(bVertex->position().z() - reFitVertex.z())));
            bsRootTree_->BsCtErr3Drefit_ = BsPDGMass_*(error3D2*abs(scale12))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());

           //Error on catu 2D


            VertexDistanceXY d22;
            Measurement1D measurement22 = d22.distance(reFitVertex,bVertex->vertexState());
            double error2D2 = measurement22.error();
            double scale22 = ((bVertex->position().x() - reFitVertex.x())*Bsvec.x()+
			      (bVertex->position().y() - reFitVertex.y())*Bsvec.y())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
               sqrt((bVertex->position().x() - reFitVertex.x())*(bVertex->position().x() - reFitVertex.x())+
		    (bVertex->position().y() - reFitVertex.y())*(bVertex->position().y() - reFitVertex.y())));
            bsRootTree_->BsCtErr2Drefit_ = BsPDGMass_*(error2D2*abs(scale22))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());



// error on ctau 3D MPV

            bsRootTree_->BsCtErrMPVrefit_ = BsPDGMass_/sqrt(deno2);


                    VertexDistanceXY vdist;
                    if(Bsvec.perp()!=0) {
              	bsRootTree_->BsLxy_    = vdist.distance( reFitVertex, bVertex->vertexState() ).value();
            	bsRootTree_->BsLxyErr_ = vdist.distance( reFitVertex, bVertex->vertexState() ).error();

	      		if (  (bVertex->position().x()- reFitVertex.x())*Bsvec.x()+(bVertex->position().y()-reFitVertex.y())*Bsvec.y() < 0  )
						bsRootTree_->BsLxy_ = -1.0 * bsRootTree_->BsLxy_;   // in case negative sign is necessary
	      			bsRootTree_->BsCt_     = bsRootTree_->BsLxy_     *  fittedBsMass/Bsvec.perp();
	      			bsRootTree_->BsCtErr_  = bsRootTree_->BsLxyErr_  *  fittedBsMass/Bsvec.perp();
	   		}
            
		 	   bsRootTree_->BsErrX_  = bs_er(1,1);
	    		bsRootTree_->BsErrY_  = bs_er(2,2);
	    		bsRootTree_->BsErrXY_ = bs_er(1,2);

	    		VertexDistance3D vdist3d;
	    		bsRootTree_->BsDist3d_    = vdist3d.distance(bVertex->vertexState(),reFitVertex).value();
	    		bsRootTree_->BsDist3dErr_ = vdist3d.distance(bVertex->vertexState(),reFitVertex).error();
	    		bsRootTree_->BsTime3d_    = bsRootTree_->BsDist3d_    * fittedBsMass/Bsvec.perp() * 100. /3.;
	    		bsRootTree_->BsTime3dErr_ = bsRootTree_->BsDist3dErr_ * BCand.mass()/Bsvec.perp() * 100. /3.;

	    		bsRootTree_->BsDist2d_     = vdist.distance(bVertex->vertexState(),reFitVertex).value();
	    		bsRootTree_->BsDist2dErr_ = vdist.distance(bVertex->vertexState(),reFitVertex).error();
	    		bsRootTree_->BsTime2d_     = bsRootTree_->BsDist2d_ * fittedBsMass/Bsvec.perp() *100. /3.;
	    		bsRootTree_->BsTime2dErr_  = bsRootTree_->BsDist2dErr_ * fittedBsMass/Bsvec.perp() * 100. /3.;

            if(verbose_ == true){
              std::cout<<"Calculate the angles of Bs"<<std::endl;
            }

//====================================================================================Transversity basis angles
                        

                        TLorentzVector pbs;
	    		pbs.SetPxPyPzE(BCand.px(),BCand.py(),BCand.pz(),BCand.energy());

       		TLorentzVector pmuplus;
	    		TLorentzVector pmuminus;
            if (jpsi_children[0]->currentState().particleCharge() == 1) 
                 {
	      		pmuplus.SetXYZM(bs_par3[3],bs_par3[4],bs_par3[5],bs_par3[6]);
	      		pmuminus.SetXYZM(bs_par4[3],bs_par4[4],bs_par4[5],bs_par4[6]);
                 } 
            else {
	      		pmuminus.SetXYZM(bs_par3[3],bs_par3[4],bs_par3[5],bs_par3[6]);
	      		pmuplus.SetXYZM(bs_par4[3],bs_par4[4],bs_par4[5],bs_par4[6]);
                 }

	    		TLorentzVector pkplus;
	    		TLorentzVector pkminus;

            if (bs_children[0]->currentState().particleCharge() == 1) {
			      pkplus.SetXYZM(bs_par1[3],bs_par1[4],bs_par1[5],bs_par1[6]);
			      pkminus.SetXYZM(bs_par2[3],bs_par2[4],bs_par2[5],bs_par2[6]);
            } else {
			      pkminus.SetXYZM(bs_par1[3],bs_par1[4],bs_par1[5],bs_par1[6]);
			      pkplus.SetXYZM(bs_par2[3],bs_par2[4],bs_par2[5],bs_par2[6]);
            }

                        //boosting in JPsi restframe

	    		TLorentzVector pjpsi;
	    		pjpsi = pmuplus + pmuminus;
	    		TLorentzVector pphi;
	    		pphi = pkplus + pkminus;

	   		  //the betas for the boost

	    		TVector3 p3_JPsi;
	    		p3_JPsi = pjpsi.Vect();
	    		p3_JPsi *= -1./pjpsi.E();

	    		 //the boost matrix

	   		TLorentzRotation boost_jpsi(p3_JPsi);
	    		TLorentzVector p_JPsi_JPsi;
	    		p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);

	    		 //the different momenta in the new frame

	    		TLorentzVector p_JPsi_muplus;
	    		TLorentzVector p_JPsi_Kplus;
	    		TLorentzVector p_JPsi_phi;
	    		p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
	    		p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);
	    		p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);

	    		 //the 3-momenta

	    		TVector3 p3_JPsi_muplus;
	    		p3_JPsi_muplus = p_JPsi_muplus.Vect();
	    		TVector3 p3_JPsi_Kplus;
	    		p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	    		TVector3 p3_JPsi_phi;
	    		p3_JPsi_phi = p_JPsi_phi.Vect();


                       //coordinate system

			TVector3 x,y,z;
			x = p3_JPsi_phi.Unit();
	    		y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	    		y = y.Unit();
	    		z = x.Cross(y);

	    		//Transversity Basis
	    		angle_costheta = p3_JPsi_muplus.Unit() * z;

	    		double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costheta*angle_costheta);
	    		double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costheta*angle_costheta);
	    		angle_phi = TMath::ACos(cos_phi);
	    		if (sin_phi < 0){
	      		angle_phi =  -angle_phi;
	    		}

                       // boosting in phi restframe
                       // the betas for the boost
                        TVector3 p3_phi;
			p3_phi = pphi.Vect();
	    		p3_phi *= -1./pphi.E();
            // the boost matrix

            TLorentzRotation boost_phi(p3_phi);
	    TLorentzVector p_phi_phi;
	    p_phi_phi = boost_phi.VectorMultiplication(pphi);

            // the different momenta in the new frame

            TLorentzVector p_phi_Kplus;
	    TLorentzVector p_phi_JPsi;
	    TLorentzVector p_phi_Bs;
	    p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	    p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	    p_phi_Bs = boost_phi.VectorMultiplication(pbs);

           // the 3-momenta

            TVector3 p3_phi_Kplus;
	    p3_phi_Kplus = p_phi_Kplus.Vect();
	    TVector3 p3_phi_JPsi;
	    p3_phi_JPsi = p_phi_JPsi.Vect();
	    TVector3 p3_phi_Bs;
	    p3_phi_Bs = p_phi_Bs.Vect();
	    angle_cospsi = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();

           // set cos of angle between bs momentum and decay length
                                  AngleBsDecayLength = ((bVertex->position().x()-PVx) * BCand.px() + (bVertex->position().y()-PVy) * BCand.py() +
				  (bVertex->position().z()-PVz) * BCand.pz()) / sqrt(((bVertex->position().x()-PVx) * (bVertex->position().x()-PVx) +
										      (bVertex->position().y()-PVy) * (bVertex->position().y()-PVy) +
										      (bVertex->position().z()-PVz) * (bVertex->position().z()-PVz)) *
										     (BCand.px()*BCand.px() + BCand.py()*BCand.py() +
										      BCand.pz()*BCand.pz()));

	    bsRootTree_->getAngles(angle_costheta,angle_phi,angle_cospsi,AngleBsDecayLength);
                       
                /*
                // number of pixel/tracker/muons hits kaons
                int pixhits1 = 0;
                // hit pattern of the track
                const reco::HitPattern& p = trk1Ref.get()->hitPattern();
                // loop over the hits of the track
                for (int iter=0; iter<p.numberOfHits(HitPattern::TRACK_HITS); iter++) 
              {
	        	uint32_t hit = p.getHitPattern(HitPattern::TRACK_HITS,iter);

                        // if the hit is valid and in pixel barrel & endcap, print out the layer
                        if (p.validHitFilter(hit) && p.pixelBarrelHitFilter(hit)) pixhits1++;
	        	if (p.validHitFilter(hit) && p.pixelEndcapHitFilter(hit)) pixhits1++;
               }
	 bsRootTree_->K1pixH_ = pixhits1;
         // count the number of valid tracker *** hits ***
         bsRootTree_->K1trkH_= p.numberOfValidTrackerHits();
           // count the number of tracker *** layers *** with measurement
            bsRootTree_->K1trkLay_ =p.trackerLayersWithMeasurement();
	    bsRootTree_->K1muDTh_  =p.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->K1muCSCh_ =p.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->K1muRPCh_ =p.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC

	    int pixhits2=0;
	    const reco::HitPattern& p2 = trk2Ref.get()->hitPattern();
	    for (int iter=0; iter<p2.numberOfHits(HitPattern::TRACK_HITS); iter++) {
	      uint32_t hit = p2.getHitPattern(HitPattern::TRACK_HITS,iter);
	      if (p2.validHitFilter(hit) && p2.pixelBarrelHitFilter(hit)) pixhits2++;
	      if (p2.validHitFilter(hit) && p2.pixelEndcapHitFilter(hit)) pixhits2++;
	    }
	    bsRootTree_->K2pixH_   = pixhits2;
	    bsRootTree_->K2trkH_   = p2.numberOfValidTrackerHits();
	    bsRootTree_->K2trkLay_ = p2.trackerLayersWithMeasurement();
	    bsRootTree_->K2muDTh_  = p2.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->K2muCSCh_ = p2.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->K2muRPCh_ = p2.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
           // number of pixel/tracker/muons hits muons

            int pixhits3 = 0;
	    const reco::HitPattern& p3 = trkMu1Ref.get()->hitPattern();
	    for (int iter=0; iter<p3.numberOfHits(HitPattern::TRACK_HITS); iter++) {
	      uint32_t hit = p3.getHitPattern(HitPattern::TRACK_HITS,iter);
	      if (p3.validHitFilter(hit) && p3.pixelBarrelHitFilter(hit)) pixhits3++;
		if (p3.validHitFilter(hit) && p3.pixelEndcapHitFilter(hit)) pixhits3++;
	    }
	    bsRootTree_->Mu1pixH_   = pixhits3;
	    bsRootTree_->Mu1trkH_   = p3.numberOfValidTrackerHits();
	    bsRootTree_->Mu1trkLay_ = p3.trackerLayersWithMeasurement();
	    bsRootTree_->Mu1muDTh_  = p3.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->Mu1muCSCh_ = p3.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->Mu1muRPCh_ = p3.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC

	    int pixhits4=0;
	    const reco::HitPattern& p4 = trkMu2Ref.get()->hitPattern();
	    for (int iter=0; iter<p4.numberOfHits(HitPattern::TRACK_HITS); iter++) {
	      uint32_t hit = p4.getHitPattern(HitPattern::TRACK_HITS,iter);
	      if (p4.validHitFilter(hit) && p4.pixelBarrelHitFilter(hit)) pixhits4++;
	      if (p4.validHitFilter(hit) && p4.pixelEndcapHitFilter(hit)) pixhits4++;
	    }
	    bsRootTree_->Mu2pixH_   = pixhits4;
	    bsRootTree_->Mu2trkH_   = p4.numberOfValidTrackerHits();
	    bsRootTree_->Mu2trkLay_ = p4.trackerLayersWithMeasurement();
	    bsRootTree_->Mu2muDTh_  = p4.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->Mu2muCSCh_ = p4.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->Mu2muRPCh_ = p4.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
        */
           // deltaR matching
            bool K1Truth = MCmatching( track1, genParticles, bsRootTree_->K1mcId_, bsRootTree_->K1momId_, bsRootTree_->K1gmomId_, 333, 531);
	    bool K2Truth = MCmatching( track2, genParticles, bsRootTree_->K2mcId_, bsRootTree_->K2momId_, bsRootTree_->K2gmomId_, 333, 531);
	    bool Mu1Truth= MCmatching( mu1,    genParticles, bsRootTree_->Mu1mcId_,bsRootTree_->Mu1momId_,bsRootTree_->Mu1gmomId_, 443, 531);
	    bool Mu2Truth= MCmatching( mu2,    genParticles, bsRootTree_->Mu2mcId_,bsRootTree_->Mu2momId_,bsRootTree_->Mu2gmomId_, 443, 531);
	    if (K1Truth==1 && K2Truth==1 && Mu1Truth==1 && Mu2Truth==1)  bsRootTree_->isMatched_ = 1;
	    else bsRootTree_->isMatched_ = 0;

                            }
 
//delete track11;
//delete track22;

           }// trk2 loop
}//trk1 loop
}//muons1
}//muons2

                          edm::Handle<View<pat::PackedCandidate>> allTracksPi;
                          iEvent.getByToken(trackLabelK, allTracksPi);
 
                    //      double BcAngStore=-999.;
                          double minBcProb=-999.;
                      //    double maxBcP=-999.;

                          if(minVtxP>0 && bsRootTree_->BsFitM_>5.00 && bsRootTree_->BsFitM_<5.43 && bsRootTree_->triggerbit_HLTmu4TkTkDis_ == 1)
                            {
                                         for (size_t i2=0; i2< allTracksPi->size(); ++i2)
                                                {         
   
                                                          double pimass = 0.139570;
                                                             const pat::PackedCandidate & Pitrack = (*allTracksPi)[i2];
                                                          if (!Pitrack.hasTrackDetails())continue;
                                                          const reco::Track &  Pirtrk = (*allTracksPi)[i2].pseudoTrack();
                                                          const reco::Track &  rtrk1bc = (*allTracksPi)[i2].pseudoTrack();
                                                          std::cout<<" =======Pt'strack"<<Pitrack.pt()<<"========Pt's rtrk"<<rtrk1bc.pt()<<"\n";
                                                          if (!Pitrack.clone()->hasTrackDetails())continue;
                                                          pat::PackedCandidate *Pitrack1 = Pitrack.clone();
                                                          std::cout<<"=============BCand_best eta phi"<<BCand_best.eta()<<"\t"<<BCand_best.phi()<<"\n";
                                                          double DeltaRBCandPion = 999;
                                                          DeltaRBCandPion = deltaR(BCand_best.eta(), BCand_best.phi(), Pitrack.eta(), Pitrack.phi());
                                                          std::cout<<"================bccandidate pion delta R"<<DeltaRBCandPion<<"\n";
                                                          if (DeltaRBCandPion > 2.2) continue;
                                                          pat::CompositeCandidate BcCand;
                                                          //nominalPionMass, nominalPionMass(0.139570)
                                                          Pitrack1->setMass(pimass);
                                                          BcCand.addDaughter(BCand_best);
                                                          BcCand.addDaughter(*Pitrack1);
                                                          AddFourMomenta add4mom;
                                                          add4mom.set(BcCand);
                                                          std::cout<<" BcCandMass"<<BcCand.mass()<<"\n";                                                                       
                                                          if (BcCand.mass() < 5.50 || BcCand.mass() > 7) continue;
                                                          //std::cout<<" BcCandMass"<<BcCand.mass()<<"\n";   
                                                 ESHandle<MagneticField> bFieldHandle;
                                                 iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
                                                 TransientTrack pion1TT(Pirtrk, &(*bFieldHandle) );
                                                          if (Pitrack.numberOfHits()<5) continue;
                                                          std::cout<<"=====================Number of Hits Pitrack"<<Pitrack.numberOfHits()<<"\n";
                                                          //if (Pirtrk == rtrk1 ){continue;}
                                                          //if (Pirtrk == rtrk2 ){continue;} 
                                                          //if (trkMu1Ref_best== Pirtrk) { continue;}
                                                          //if (trkMu2Ref_best== Pirtrk) { continue;}
 std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Tralala"<<"\n";                                                
 edm::ESHandle<TransientTrackBuilder> theB;
                                                 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
//BadRefCore RefCore: Request to resolve a null or invalid reference to a product of type 'std::vector<reco::Track>' has been detected.
//Please modify the calling code to test validity before dereferencing.                                                    
                                                if(trk1Ref_best == trk2Ref_best)continue;
                                                vector<TransientTrack> t_tracks;
       	  	                            	t_tracks.push_back((*theB).build(&trkMu1Ref_best));
                          	  		t_tracks.push_back((*theB).build(&trkMu2Ref_best));
       	                         		t_tracks.push_back((*theB).build(&trk1Ref_best));
                         	  		t_tracks.push_back((*theB).build(&trk1Ref_best));
                                     std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Tralala"<<"\n";
                                     KinematicFitInterface Kfitter2;
                                     bool fitSuccess = Kfitter2.doFit(t_tracks, nominalMuonMass,  nominalKaonMass, nominalKaonMass);
                                     if(fitSuccess != 1) {continue;}
                                     RefCountedKinematicParticle bs_best = Kfitter2.getParticle();
                                     RefCountedKinematicVertex BsVertex = Kfitter2.getVertex();
                                     RefCountedKinematicTree Bs4BcTree = Kfitter2.getTree();

                                     // if the fit fails, do not consider this as candidate
                                     //if(Bs4BcTree->isEmpty()) continue;
                                     KinematicParticleFitter constBsFitter;
                                     float BsMsigma = 0.0002;
                                     KinematicConstraint * bs_const = new MassKinematicConstraint( BsPDGMass_, BsMsigma);
                                     KinematicParticleFitter constFitter;
                                     RefCountedKinematicTree BsTreeN1 = constFitter.fit(bs_const,Bs4BcTree);
                                     //if(BsTreeN1->isEmpty()) continue;
                                      //BsTreeN1->movePointerToTheTop();
                                      RefCountedKinematicParticle bsConstr = BsTreeN1->currentParticle();
                                      KinematicParticleFactoryFromTransientTrack p_Factory;
                                     ParticleMass pion_mass = 0.13957018;
                                    float pion_sigma = pion_mass*1.e-6;
                                    float chi = 0.;
                                    float ndf = 0.;
                                    std::cout<<"++++++++++++++++++++++LLLLL"<<"\n";
                                    vector<RefCountedKinematicParticle> bcParticles;
                                    bcParticles.push_back(p_Factory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
                                    bcParticles.push_back(bsConstr);
                                    KinematicConstrainedVertexFitter Bcfitter;   
                                    RefCountedKinematicTree BcVertexFitTree = Bcfitter.fit(bcParticles);                           
                                    //if (BcVertexFitTree->isEmpty()) {continue;}
                                    if (!BcVertexFitTree->isValid()) {continue; }
                                   // BcVertexFitTree->movePointerToTheTop();

                                        RefCountedKinematicParticle BcFromFit = BcVertexFitTree->currentParticle();
                                        RefCountedKinematicVertex BcVertex = BcVertexFitTree->currentDecayVertex();

                                int ClosestPVindex = -1;
	                        Double_t MinDistance = 10000000; 	
                                double BcPVx;	
                                double BcPVy;	
                  //              double BcPVz;
                               for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++)
                                   {
			                       const Vertex &vtx = (*recVtxs)[VtxInd];
	                          	       Double_t PVSVdistance = TMath::Abs(BcVertex->position().z()-vtx.z());
                                               std::cout<<"==========Primary vertex secondary vertex distance"<<PVSVdistance<<"\n";
		                               if(PVSVdistance < MinDistance)
                                                       {
				                              MinDistance = PVSVdistance;
				                              ClosestPVindex = VtxInd;	
			                               }			
                        		}
                        
                if (ClosestPVindex==-1){ continue;}
	        const reco::Vertex &BcRecVtx = (*recVtxs)[ClosestPVindex]; 
		BcPVx = BcRecVtx.x(); 
	        BcPVy = BcRecVtx.y();
	       /// BcPVz= BcRecVtx.z();
                                std::pair<bool,Measurement1D> ImpactPar3D = IPTools::absoluteImpactParameter3D(pion1TT, BcRecVtx);     
                                AlgebraicVector7 b_par = BcFromFit->currentState().kinematicParameters().vector();

                                                  GlobalVector Bcvec(b_par[3], b_par[4], b_par[5]);
                                                  math::XYZVector Bcpperp(b_par[3], b_par[4], b_par[5]);
                                                  reco::Vertex::Point Bcvperp(BcVertex->position().x(),BcVertex->position().y(),BcVertex->position().z());
                                            
                                                  double vtxprob_Bc = TMath::Prob(BcVertex->chiSquared(), (int)BcVertex->degreesOfFreedom());
                                                  std::cout<<"=========================================================================================vtxprobbc:"<<"\t"<<vtxprob_Bc<<"\n";
                                                  double AngBcBs = reco::deltaR(BCand_best.eta(),BCand_best.phi(),Bcvec.eta(),Bcvec.phi());
                                                  double AngBcPi = reco::deltaR(Bcvec.eta(),Bcvec.phi(),Pitrack.eta(),Pitrack.phi());
                                                  double BcCt;
                                            
                                                  VertexDistanceXY vdist;
                                                  VertexDistance3D dist3;
                                                  double BcLxy    = vdist.distance( BcRecVtx, BcVertex->vertexState() ).value();
                                                  double BsLxy    = vdist.distance( BcRecVtx, BsVertex->vertexState() ).value();
                                    //              double BcLxyz    = dist3.distance( BcRecVtx, BcVertex->vertexState() ).value();
                                      //            double BsLxyz    = dist3.distance( BcRecVtx, BsVertex->vertexState() ).value();
                                                  if(sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z())!=0) {
                                                    double BcL_ =(BcVertex->position().x()-BcPVx)*Bcvec.x()+(BcVertex->position().y()-BcPVy)*Bcvec.y();          
                                                   if ( BcL_ < 0 )BcL_ = -1.0 * BcL_;   // in case negative sign is necessary 
                                                                                 BcCt=BcL_*6.277/(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y());
                                                                       }
                                                                        
                                                     else {continue;}
                                                     if (BsLxy<BcLxy) continue;
                                                     if (AngBcPi>1.) continue;
                                                     if (TMath::Abs(Bcvperp.Dot(Bcpperp)/(Bcvperp.R()*Bcpperp.R()))<0.9) continue;
                                                     if (ImpactPar3D.second.value()>0.5) continue;
                                                     if (Pitrack.pt()<0.3) continue;
                                                  /*   if (energyLossHandle.isValid() && Pitrack.p()<1) 
                                                      {
                                                               const DeDxDataValueMap &  eloss  = *energyLossHandle;
                                                               double dedxTrk = eloss[rtrk1bc].dEdx();
                                                               double trackmass=TMath::Sqrt((dedxTrk-2.557)/2.579)*Pitrack.p();
                                                      }*/
                                                      if (vtxprob_Bc > minBcProb)
                                                             {
                                                                      minBcProb=vtxprob_Bc;
                                                                    std::cout<<" maxBcP:=============="<<sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z())<<"\n";
                                                                      //BcAngStore=AngBcBs;
                                                                      bsRootTree_->BcCt_=BcCt;
                                                                      bsRootTree_->BcCosAlpha_ = Bcvperp.Dot(Bcpperp)/(Bcvperp.R()*Bcpperp.R());
                                                                     bsRootTree_->BcIP3D_=ImpactPar3D.second.value();
                                                                            bsRootTree_->BcP_=sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y());
                                                                            bsRootTree_->BctrackPt_ = Pitrack.pt(); 
                                                                            bsRootTree_->BcM_ = BcFromFit->currentState().mass(); 
                                                                            bsRootTree_->BcAng_ = AngBcBs;
                                                                            bsRootTree_->BcAngPi_ = AngBcPi;
                                                                            bsRootTree_->BcProb_ = vtxprob_Bc;
                                                                            bsRootTree_->deltaRBCandPion_ = DeltaRBCandPion;
                                                                       } 
                                                                                                                                                                                            

                                                         
//}}}
  }}                                                                      



              bsRootTree_->fill();
              }//muon three braces
/*
//======================

//==========================================Brace<->main

if(verbose_ == true){ std::cout<<"Fill the TTree"<<std::endl;}
   bsRootTree_->fill();
}
*/
//////==============================================================================================================================================================
//////====================================================================================================================================================================================
//////=========================================================================================================================================================================================================
GlobalVector BsToJpsiPhiAnalysis::flightDirection(const reco::Vertex &pv, reco::Vertex &sv){
  GlobalVector res(sv.position().X() - pv.position().X(),
                    sv.position().Y() - pv.position().Y(),
                    sv.position().Z() - pv.position().Z());
  return res;
}


void BsToJpsiPhiAnalysis::fillMCInfo(edm::Handle<edm::View<reco::GenParticle>>& genParticles){
  int iNumberOfBdecays = 0;
  Int_t SVindex =0;
//B meson ids
std::set<int> listOfBmesonIds;
  listOfBmesonIds.insert(511 );   // Bd
  listOfBmesonIds.insert(521 );   // B+
  listOfBmesonIds.insert(10511 );    // B_0*0
  listOfBmesonIds.insert(10521 );    // B_0*+
  listOfBmesonIds.insert(513 );   // B*d
  listOfBmesonIds.insert(523 );   // B*d+
  listOfBmesonIds.insert(10513 );   // B1(L)0
  listOfBmesonIds.insert(10523 );   // B1(L)+
  listOfBmesonIds.insert(20513 );   // B1(H)0
  listOfBmesonIds.insert(20523 );   // B1(H)+
  listOfBmesonIds.insert(515 );    // B2*_0
  listOfBmesonIds.insert(525 );    // B2*_+
  listOfBmesonIds.insert(531 );   // Bs
  listOfBmesonIds.insert(10531 );    // B_s0*_0
  listOfBmesonIds.insert(533 );   // B*s
  listOfBmesonIds.insert(10533 );   // Bs1(L)0
  listOfBmesonIds.insert(20533 );   // Bs1(H)0
  listOfBmesonIds.insert(535 );    // Bs2*_0
  listOfBmesonIds.insert(541 );   // Bc+
  listOfBmesonIds.insert(10541 );   // B*c0+
  listOfBmesonIds.insert(543 );   // B*c+
  listOfBmesonIds.insert(10543 );   // Bc1(L)+
  listOfBmesonIds.insert(20543 );   // Bc1(H)+
  listOfBmesonIds.insert(545 );    // Bc2*_0
  listOfBmesonIds.insert(551 );   // etab(1S)
  listOfBmesonIds.insert(10551 );   // chib(1P)
  listOfBmesonIds.insert(100551 );   // etab(2S)
  listOfBmesonIds.insert(110551 );   // chib(2P)
  listOfBmesonIds.insert(200551 );   // etab(3S)
  listOfBmesonIds.insert(210551 );   // chib(3P)
  listOfBmesonIds.insert(553 );   // upsilon(1S)
  listOfBmesonIds.insert(10553 );   // hb(1P)
  listOfBmesonIds.insert(20553 );   // chib1(1P)
  listOfBmesonIds.insert(30553 );   // upsilon1(1D)
  listOfBmesonIds.insert(100553 );   // upsilon(2S)
  listOfBmesonIds.insert(110553 );   // hb(2P)
  listOfBmesonIds.insert(120553 );   // chib1(2P)
  listOfBmesonIds.insert(130553 );   // upsilon1(2D)
  listOfBmesonIds.insert(200553 );   // upsilon(3S)
  listOfBmesonIds.insert(210553 );   // hb(3P)
  listOfBmesonIds.insert(220553 );   // chib1(3P)
  listOfBmesonIds.insert(300553 );   // upsilon(4S)
  listOfBmesonIds.insert(9000553 );   // upsilon(10860)
  listOfBmesonIds.insert(9010553 );   // upsilon(11020)
  listOfBmesonIds.insert(555 );   // chib2(1P)
  listOfBmesonIds.insert(10555 );   // etab2(1D)
  listOfBmesonIds.insert(20555 );   // upsilon2(1D)
  listOfBmesonIds.insert(100555 );   // chib2(2P)
  listOfBmesonIds.insert(110555 );   // etab2(2D)
  listOfBmesonIds.insert(120555 );   // upsilon2(2D)
  listOfBmesonIds.insert(200555 );   // chib2(3P)
  listOfBmesonIds.insert(557 );   // upsilon3(1D)
  listOfBmesonIds.insert(100557 );   // upsilon3(2D)
  listOfBmesonIds.insert(5122 );   // lambda_b0
  listOfBmesonIds.insert(5112 );   // sigma_b-
  listOfBmesonIds.insert(5212 );   // sigma_b0
  listOfBmesonIds.insert(5222 );   // sigma_b+
  listOfBmesonIds.insert(5114 );   // sigma*_b-
  listOfBmesonIds.insert(5214 );   // sigma*_b0
  listOfBmesonIds.insert(5224 );   // sigma*_b+
  listOfBmesonIds.insert(5132 );   // Xi_b-
  listOfBmesonIds.insert(5232 );   // Xi_b0
  listOfBmesonIds.insert(5312 );   // Xi'_b-
  listOfBmesonIds.insert(5322 );   // Xi'_b0
  listOfBmesonIds.insert(5314 );   // Xi*_b-
  listOfBmesonIds.insert(5324 );   // Xi*_b0
  listOfBmesonIds.insert(5332 );   // Omega_b-
  listOfBmesonIds.insert(5334 );   // Omega*_b-
  listOfBmesonIds.insert(5142 );   // Xi_bc0
  listOfBmesonIds.insert(5242 );   // Xi_bc+
  listOfBmesonIds.insert(5412 );   // Xi'_bc0
  listOfBmesonIds.insert(5422 );   // Xi'_bc+
  listOfBmesonIds.insert(5414 );   // Xi*_bc0
  listOfBmesonIds.insert(5424 );   // Xi*_bc+
  listOfBmesonIds.insert(5342 );   // Omega_bc0
  listOfBmesonIds.insert(5432 );   // Omega'_bc0
  listOfBmesonIds.insert(5434 );   // Omega*_bc0
  listOfBmesonIds.insert(5442 );   // Omega_bcc+
  listOfBmesonIds.insert(5444 );   // Omega*_bcc+
  listOfBmesonIds.insert(5512 );   // Xi_bb-
  listOfBmesonIds.insert(5522 );   // Xi_bb0
  listOfBmesonIds.insert(5514 );   // Xi*_bb-
  listOfBmesonIds.insert(5524 );   // Xi*_bb0
  listOfBmesonIds.insert(5532 );   // Omega_bb-
  listOfBmesonIds.insert(5524 );   // Omega*_bb-
  listOfBmesonIds.insert(5542 );   // Omega_bbc0
  listOfBmesonIds.insert(5544 );   // Omega*_bbc0
  listOfBmesonIds.insert(554 );   // Omega_bbb-

  const Candidate * Jpsi = 0;
  const Candidate * Phi = 0;
  const Candidate * mup = 0;
  const Candidate * mum = 0;
  const Candidate * Kp = 0;
  const Candidate * Km = 0;

  if(verbose_ == true){
    cout<<"NewEvent-------"<<endl;
  }
//Begin loop over all particles

for( size_t i = 0; i < genParticles->size(); ++ i ) {
    const GenParticle & genBsCand = (*genParticles)[ i ];
    int MC_particleID=genBsCand.pdgId();
    int absMC_particleID = abs(MC_particleID);
	if(absMC_particleID ==443 ){
		unsigned int numJpsiDaus = genBsCand.numberOfDaughters();
		vector<unsigned int> JpsiMuIdx;
		vector<unsigned int> JpsiPhotoIdx;

		if (	numJpsiDaus > 1 && numJpsiDaus < 6){

 		 for(unsigned int a = 0; a < numJpsiDaus ; a++){
			if(abs(genBsCand.daughter(a)->pdgId() ) == 13){
				
				JpsiMuIdx.push_back(a);
			}
			if(abs(genBsCand.daughter(a)->pdgId()) == 22){
				JpsiPhotoIdx.push_back(a);
			
			}
			
 		 }	// ned of for loop
		} // end of if numJpsiDaus

//cout << "JpsiPhotoIdx.size() + JpsiMuIdx.size() = " <<  JpsiPhotoIdx.size() + JpsiMuIdx.size() << " = " << numJpsiDaus <<"\n";
if(JpsiMuIdx.size() == 2){
		  if( JpsiPhotoIdx.size() + JpsiMuIdx.size() == numJpsiDaus &&  genBsCand.daughter(JpsiMuIdx[0])->pdgId() == -genBsCand.daughter(JpsiMuIdx[1])->pdgId()  ) {

			bsRootTree_->JpsiGenNumberOfCandidates_++;
	 	   bsRootTree_->JpsiGenPVx_ = genBsCand.vx();
	      bsRootTree_->JpsiGenPVy_ = genBsCand.vy();
			bsRootTree_->JpsiGenPVz_ = genBsCand.vz();
			bsRootTree_->JpsiGenLxyOld_= sqrt( pow( (bsRootTree_->BSx_ - genBsCand.vx()) , 2) + pow( (bsRootTree_->BSy_ - genBsCand.vy() ),2 )  );
GlobalPoint GenDisplacementFromBeamspot( -1*( (bsRootTree_->BSx_ -  genBsCand.vx() ) + ( genBsCand.vz() - bsRootTree_->BSz_) * bsRootTree_->BSdxdz_ ), -1*( (bsRootTree_->BSy_ - genBsCand.vy() )+  (genBsCand.vz() - bsRootTree_->BSz_) * bsRootTree_->BSdydz_), 0 );
bsRootTree_->JpsiGenLxy_ = GenDisplacementFromBeamspot.perp();
bsRootTree_->JpsiGenLxyOverPt_ = ( GenDisplacementFromBeamspot.x()*genBsCand.px()  +  GenDisplacementFromBeamspot.y()*genBsCand.py()  ) /( genBsCand.px() * genBsCand.px() + genBsCand.py() * genBsCand.py() );
			bsRootTree_->JpsiGenPt_ = genBsCand.pt();

//cout << "JpsiGenLxy_ : " << bsRootTree_->JpsiGenLxy_ << " = " << GenDisplacementFromBeamspot.perp() << endl;
//cout << "JpsiGenLxy/Pt:  " << bsRootTree_->JpsiGenLxyOverPt_ << endl;

if( genBsCand.daughter(JpsiMuIdx[0])->pdgId() == 13 ){ 


bsRootTree_->JpsiMu1PtMC_ = genBsCand.daughter(JpsiMuIdx[0])->pt();
				bsRootTree_->JpsiMu1EtaMC_ =  genBsCand.daughter(JpsiMuIdx[0])->eta();
				bsRootTree_->JpsiMu1PhiMC_ =  genBsCand.daughter(JpsiMuIdx[0])->phi();

				bsRootTree_->JpsiMu2PtMC_ = genBsCand.daughter(JpsiMuIdx[1])->pt();
				bsRootTree_->JpsiMu2EtaMC_ = genBsCand.daughter(JpsiMuIdx[1])->eta();
				bsRootTree_->JpsiMu2PhiMC_ =  genBsCand.daughter(JpsiMuIdx[1])->phi();


//cout << "Gen Mu- " <<	genBsCand.daughter(JpsiMuIdx[0])->charge() << " Pt " << 	bsRootTree_->JpsiMu1PtMC_ << " eta " << bsRootTree_->JpsiMu1EtaMC_ << "phi "<< bsRootTree_->JpsiMu1PhiMC_ << endl;
//cout << "Gen Mu+ " <<  genBsCand.daughter(JpsiMuIdx[1])->charge() << " Pt " << 	bsRootTree_->JpsiMu2PtMC_  << " eta " << bsRootTree_->JpsiMu2EtaMC_ << "phi "<< bsRootTree_->JpsiMu2PhiMC_ << endl;
			}
else{

				bsRootTree_->JpsiMu1PtMC_ = genBsCand.daughter(JpsiMuIdx[1])->pt();
				bsRootTree_->JpsiMu1EtaMC_ =  genBsCand.daughter(JpsiMuIdx[1])->eta();
				bsRootTree_->JpsiMu1PhiMC_ =  genBsCand.daughter(JpsiMuIdx[1])->phi();

				bsRootTree_->JpsiMu2PtMC_ = genBsCand.daughter(JpsiMuIdx[0])->pt();
				bsRootTree_->JpsiMu2EtaMC_ = genBsCand.daughter(JpsiMuIdx[0])->eta();
				bsRootTree_->JpsiMu2PhiMC_ =  genBsCand.daughter(JpsiMuIdx[0])->phi();
			
//cout << "Gen else Mu- " <<  genBsCand.daughter(JpsiMuIdx[1])->charge() << " Pt " << 	bsRootTree_->JpsiMu1PtMC_<<  " eta " << bsRootTree_->JpsiMu1EtaMC_ << "phi "<< bsRootTree_->JpsiMu1PhiMC_ << endl;
}
math::XYZVector JpsiGenPtvec(genBsCand.px() , genBsCand.py(), 0);
		
 			 bsRootTree_->JpsiCosAlphaMC_ = ( GenDisplacementFromBeamspot.x()* genBsCand.px() +  genBsCand.py()* GenDisplacementFromBeamspot.y() )/( JpsiGenPtvec.R()* GenDisplacementFromBeamspot.perp() );
}//cout << "saving the jpsi  " << endl;
		  }
		
	} // end of if Jpsi
 bool isPosMu=0, isNegMu=0, isJpsi=0, isKplus=0, isNegPi=0, isPosPi=0, isKstar = 0, isNeutralPi=0;
if( listOfBmesonIds.find( absMC_particleID ) != listOfBmesonIds.end()){


bool hasBDaughter=0;
      int numBsDaughters = genBsCand.numberOfDaughters();
      for(int idau=0; idau < numBsDaughters; idau++)
        if( listOfBmesonIds.find( abs(genBsCand.daughter(idau)->pdgId())) != listOfBmesonIds.end() ) hasBDaughter=1;

if( abs(MC_particleID) == 521 ){
        int numBplusDaughters = genBsCand.numberOfDaughters();
        int jpsiIndex = -5;
        int kstarIndex = -5;
        int kpmIndex = -5;

//checking of decay channels B+ -> Jpsi(mu+,mu-) K+ and B+ -> Jpsi K*(K+, pi0)

if(numBplusDaughters == 2){ //check if Bplus has 2 daughters
          for(int k = 0; k < numBplusDaughters; k++ ){ //check if the daughter's are J/psi and (K+ or K*)
            if( genBsCand.daughter(k)->pdgId() == 443){
              isJpsi = 1;
              jpsiIndex = k;
            }
            else if( abs(genBsCand.daughter(k)->pdgId()) == 321 ){ isKplus = 1; kpmIndex = k; }
            else if( genBsCand.daughter(k)->pdgId() == 323 ){
              isKstar = 1;
              kstarIndex = k;
            }
          } // end of for-loop
          if( isJpsi == 1 && ( isKplus == 1 || isKstar == 1 ) ){
                                   
//check if J/psi has two muon daughters
const Candidate *JpsiCand = genBsCand.daughter(jpsiIndex);
            for(unsigned int j = 0; j < JpsiCand->numberOfDaughters(); j++ ){ //check if Jpsi daughters are mu+ and mu-
              const Candidate * Jpsidau = JpsiCand->daughter(j);
              if(Jpsidau->pdgId() == 13 ){ isNegMu = 1; }
              else if (Jpsidau->pdgId() == -13 ){ isPosMu = 1; }
            }  // end of for loop
if( isPosMu == 1 && isNegMu == 1 && isKplus == 1 ){
				double bpmomx=genBsCand.px();
		      double bpmomy=genBsCand.py();
	 		   double bppvx=genBsCand.vx();
	     		double bppvy=genBsCand.vy();
				double bppvz=genBsCand.vz();
		  		double bpsvx=genBsCand.daughter(0)->vx();
        		double bpsvy=genBsCand.daughter(0)->vy();
//Calculate gen ctau 2D
  double Lxy2D = sqrt(pow(bpsvx-bppvx,2)+pow(bpsvy-bppvy,2));
				bsRootTree_->BpCt2DMC_ = Lxy2D * BpPDGMass_/sqrt(bpmomx*bpmomx+bpmomy*bpmomy);
//cout << "BpCt 2D MC " <<  bsRootTree_->BpCt2DMC_ << endl;
bsRootTree_->genBpPV_z_ =  bppvz;
				 bsRootTree_->genBpPV_y_ =  bppvy;
				 bsRootTree_->genBpPV_x_ = bppvx;
 
              if(kpmIndex >=0) {
                const Candidate *KaonCand = genBsCand.daughter(kpmIndex);
                if( KaonCand->pdgId() == 321) {
                  bsRootTree_->BplusDecayChannel_ = 1;
                  if(verbose_ == true){
                    cout << "B+ DecayChannel = 1" <<endl;
                  }
                }
                if( KaonCand->pdgId() == -321) {
                  bsRootTree_->BplusDecayChannel_ = -1;
                  if(verbose_ == true){
                    cout << "B+ DecayChannel = -1" <<endl;
                  }
                }					

              }
            } // end of if( isPosMu == 1 && isNegMu == 1 && isKplus == 1 )

isKplus = 0;

            if(isKstar == 1){
              const Candidate *Kstar = genBsCand.daughter(kstarIndex);
              if (Kstar->numberOfDaughters() == 2){
                for(unsigned int j = 0; j < Kstar->numberOfDaughters(); j++ ){ //check if Kstar daughters are Kplus and neutral pion
                 const Candidate * Kstardau = Kstar->daughter(j);
                 if(Kstardau->pdgId() == 321 ){ isKplus = 1; }
                 else if (Kstardau->pdgId() == 111 ){ isNeutralPi = 1; }
                }
                if( isNegMu == 1 && isPosMu == 1 && isNeutralPi == 1 && isKplus == 1 ){
                  bsRootTree_->BplusDecayChannel_ = 3; // channel 3 = B+ -> Jpsi K*(K+ pi0)
                }
              }
            } // end of if(isKstar == 1)
          } // end of if(isJpsi ==1 && ( isKplus == 1 || isKstar == 1 ))

        } // end of  if(numBplusDaughters == 2)

//checking of decay channel B+ -> Jpsi K+ pi+ pi-
if(numBplusDaughters == 4){
          for(int j = 0; j < numBplusDaughters; j++ ){ //check if Bplus has 4 daughters
           if(genBsCand.daughter(j)->pdgId() == 443 ){ isJpsi = 1; } //check if the daughter's are J/psi, pi+, pi- and K+
           else if(genBsCand.daughter(j)->pdgId() == 321 ){ isKplus = 1; }
           else if(genBsCand.daughter(j)->pdgId() == 211 ){ isPosPi = 1; }
           else if(genBsCand.daughter(j)->pdgId() == -211 ){ isNegPi = 1; }
          }

          if(isJpsi == 1 && isKplus == 1 && isPosPi == 1 && isNegPi == 1 ){
            bsRootTree_->BplusDecayChannel_ = 2;
          }
        } // end of  B+ -> Jpsi K+ pi+ pi- channel search
      } // end of Bplus search

//End of Bp MC matching and the start of Bd MC matching

if(abs(MC_particleID) == 511 && abs(genBsCand.mother(0)->pdgId()) != 511 ) {
	  bool isBdJpsiMC=false;
	  bool isBdKstar=false;	
	  bool isBdKstarbar=false;	
	  bool isMixedBd = false;	
	  const Candidate *genBdCand;
	  double bdmomx=genBsCand.px();
     double bdmomy=genBsCand.py();
     double bdmomz=genBsCand.pz();
	  double bdpvx=genBsCand.vx();
     double bdpvy=genBsCand.vy();
     double bdpvz=genBsCand.vz();
	  double bdsvx = 0; 
	  double bdsvy = 0;
	  double bdsvz = 0;

     const Candidate *BdMuPlus=0;
     const Candidate *BdMuMinus=0;
     const Candidate *BdKPlus=0;
     const Candidate *BdKMinus=0;
if( genBsCand.numberOfDaughters() == 1 && abs( genBsCand.daughter(0)->pdgId() ) == 511 ){ 
			isMixedBd = true;
			genBdCand = genBsCand.daughter(0); 
                        //cout << "mixed "  << endl;
                        //const Candidate *genBdCand = genBsCand.daughter(0); 
			//const Candidate  &genBdCand = *(genBsCand.daughter(0)); 
				  }

 else {genBdCand = &genBsCand; }
if( genBdCand->numberOfDaughters() == 2 && abs( genBdCand->daughter(0)->pdgId() ) != 511 && abs( genBdCand->daughter(1)->pdgId() ) != 511 ){
		for(size_t l = 0; l < genBdCand->numberOfDaughters(); l++){
			const Candidate * genBdDau =genBdCand->daughter(l);
// Jpsi search
if( abs(genBdDau->pdgId() ) == 443 && genBdDau->numberOfDaughters() > 1 ){
				int genBdJpsiMuP_Id = 0;
				int genBdJpsiMuM_Id = 0;

				for(size_t m = 0; m < genBdDau->numberOfDaughters(); m++){
					if(genBdDau->daughter(m)->pdgId() == 13) {
						genBdJpsiMuP_Id = genBdDau->daughter(m)->pdgId(); 
						BdMuPlus = genBdDau->daughter(m);
					}
					if(genBdDau->daughter(m)->pdgId() == -13) {
						genBdJpsiMuM_Id = genBdDau->daughter(m)->pdgId();
						BdMuMinus = genBdDau->daughter(m); 
					}
				}			
				if( genBdJpsiMuP_Id == 13 && genBdJpsiMuM_Id == -13){
					isBdJpsiMC = true;
				}			
			} // end of Jpsi dau search

//K star search
if( abs(genBdDau->pdgId() ) == 313  ){
				int genBdPi = 0;
				int genBdK =0;			
			  for(size_t n=0; n < genBdDau->numberOfDaughters(); n++){
					int genBdKDauId = genBdDau->daughter(n)->pdgId();
					if(abs(genBdKDauId) == 211){
						genBdPi = genBdKDauId;		
						if(genBdKDauId < 0){BdKMinus = genBdDau->daughter(n); }
						else BdKPlus = genBdDau->daughter(n);		
						bsRootTree_->BdMCKstarPion_ = genBdKDauId;							
					}
					if(abs(genBdKDauId) == 321){
						genBdK = genBdKDauId;		
						if(genBdKDauId < 0){BdKMinus = genBdDau->daughter(n); }
						else BdKPlus = genBdDau->daughter(n);
						bsRootTree_->BdMCKstarKaon_ = genBdKDauId;										
					}
				}
				if( genBdK == 321 && genBdPi == -211) { //pi- and K+  
					isBdKstar = true;//cout << " K* found" << endl; 
				}
				if( genBdK == -321 && genBdPi == 211) { //pi+ and K-  
					isBdKstarbar = true;//cout << " K*bar found" << endl; 
				}
			} // end of K* dau search
		} // end of dau loop
	 } // end of if genBsCand.numberOfDaughters() == 2 

	 if(isBdJpsiMC == true && (isBdKstarbar == true || isBdKstar == true) ){
		bsRootTree_->BdIniFlavour_=(int)(MC_particleID/511);
		bsRootTree_->BdEndFlavour_=(int)(genBdCand->pdgId()/511);
		if(isMixedBd == true){ 
			//cout << "IniFlavour "  << 511*bsRootTree_->BdIniFlavour_ << endl;
			//cout << "FinalFlavour " << 511*bsRootTree_->BdEndFlavour_ << endl; 
		}

		int flavcheck =0;	
		if(isMixedBd == true){ flavcheck  = genBdCand->mother(0)->pdgId() ; }		
		else flavcheck = genBdCand->pdgId();
		if(511*bsRootTree_->BdIniFlavour_ != flavcheck ) { 
			//cout << "Odd flavour " << endl;
		}
		if(isBdKstarbar == true ){ bsRootTree_->BdChannelID_=2;}
		if(isBdKstar == true ){ bsRootTree_->BdChannelID_=1; }
	   bdsvx=genBdCand->daughter(0)->vx();
      bdsvy=genBdCand->daughter(0)->vy();
      bdsvz=genBdCand->daughter(0)->vz();

//calculate the gen ctau 2D
 double Lxy2D = sqrt(pow(bdsvx-bdpvx,2)+pow(bdsvy-bdpvy,2));
      bsRootTree_->BdCt2DMC_ = Lxy2D * BdPDGMass_/sqrt(bdmomx*bdmomx+bdmomy*bdmomy);

double Lxy3D = sqrt(pow(bdsvx-bdpvx,2)+pow(bdsvy-bdpvy,2)+pow(bdsvz-bdpvz,2));
      bsRootTree_->BdCt3DMC_ = Lxy3D * BdPDGMass_/sqrt(bdmomx*bdmomx+bdmomy*bdmomy+bdmomz*bdmomz);
	  // cout << "BdCt2DMC_ " << bsRootTree_->BdCt2DMC_ << endl;
TLorentzVector PMuPlus;
      TLorentzVector PMuMinus;	      
	   TLorentzVector PKPlus;
	   TLorentzVector PKMinus;
      PMuPlus.SetXYZM(BdMuPlus->px(),BdMuPlus->py(),BdMuPlus->pz(),BdMuPlus->mass());
      PMuMinus.SetXYZM(BdMuMinus->px(),BdMuMinus->py(),BdMuMinus->pz(),BdMuMinus->mass());
      PKPlus.SetXYZM(BdKPlus->px(),BdKPlus->py(),BdKPlus->pz(),BdKPlus->mass());
	   PKMinus.SetXYZM(BdKMinus->px(),BdKMinus->py(),BdKMinus->pz(),BdKMinus->mass());



//boosting in JPsi restframe
          TLorentzVector pjpsi;
	  pjpsi = PMuPlus + PMuMinus;
	  TLorentzVector pphi;
	  pphi = PKPlus + PKMinus;

//betas for the boost
TVector3 p3_JPsi;
	  p3_JPsi = pjpsi.Vect();
	  p3_JPsi *= -1./pjpsi.E();
//the boost matrix
TLorentzRotation boost_jpsi(p3_JPsi);
	  TLorentzVector p_JPsi_JPsi;
	  p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);

//different momenta in the new frame
TLorentzVector p_JPsi_muplus;
	  TLorentzVector p_JPsi_Kplus;
	  TLorentzVector p_JPsi_phi;
	  p_JPsi_muplus = boost_jpsi.VectorMultiplication(PMuPlus);
	  p_JPsi_Kplus = boost_jpsi.VectorMultiplication(PKPlus);
     p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);
//the three momenta
 TVector3 p3_JPsi_muplus;
	  p3_JPsi_muplus = p_JPsi_muplus.Vect();
	  TVector3 p3_JPsi_Kplus;
	  p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	  TVector3 p3_JPsi_phi;
	  p3_JPsi_phi = p_JPsi_phi.Vect();
//coordinate system
TVector3 x,y,z;
	  x = p3_JPsi_phi.Unit();
	  y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	  y = y.Unit();
	  z = x.Cross(y);
     double Bdangle_costheta=p3_JPsi_muplus.Unit() * z;   
	  bsRootTree_->BdcosthetaMC_ = Bdangle_costheta;
	  double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	  double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	  double Bdangle_phi = TMath::ACos(cos_phi);
	  if (sin_phi < 0){
	    Bdangle_phi =  -Bdangle_phi;
          }
	  bsRootTree_->BdphiMC_ = Bdangle_phi;

// boosting in phi restframe

	  TVector3 p3_phi;
	  p3_phi = pphi.Vect();
	  p3_phi *= -1./pphi.E();
	      
//	   the boost matrix
	  TLorentzRotation boost_phi(p3_phi);
	  TLorentzVector p_phi_phi;
	  p_phi_phi = boost_phi.VectorMultiplication(pphi);
	      
//	   the different momenta in the new frame
	  TLorentzVector p_phi_Kplus;
	  TLorentzVector p_phi_JPsi;
	  TLorentzVector p_phi_Bs;   
	  p_phi_Kplus = boost_phi.VectorMultiplication(PKPlus);
	  p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
//	  the 3-momenta
	  TVector3 p3_phi_Kplus;
	  p3_phi_Kplus = p_phi_Kplus.Vect();
	  TVector3 p3_phi_JPsi;
	  p3_phi_JPsi = p_phi_JPsi.Vect();
	  bsRootTree_->BdcospsiMC_ = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
	 }	// end of if(isBdJpsiMC == true && (isBdKstarbar == true || isBdKstar == true) )
  } // end of if(abs(MC_particleID) == 511 && abs(genBsCand.mother(0)->pdgId()) != 511) 

//////End of Bd MC id written by Terhi
//
//
//      //BEGIN Bs MC matching

if( abs(genBsCand.mother(0)->pdgId()) != 531 && abs(MC_particleID)==531) {
        vector<int> figlie;
        vector<int> request;
        vector<int> requestDauPhi;
        vector<int> requestDauPsi;
        vector<int> requestDauPsi2;
        vector<int> requestDauPsi3;
        vector<int> requestDauPsi4;
        request.push_back(333);
        request.push_back(443);
        requestDauPhi.push_back(321);
        requestDauPhi.push_back(321);
        requestDauPsi.push_back(13);
        requestDauPsi.push_back(13);
        requestDauPsi2.push_back(13);
        requestDauPsi2.push_back(13);
        requestDauPsi2.push_back(22);
        requestDauPsi3.push_back(13);
        requestDauPsi3.push_back(13);
        requestDauPsi3.push_back(22);
        requestDauPsi3.push_back(22);
        requestDauPsi4.push_back(13);
        requestDauPsi4.push_back(13);
        requestDauPsi4.push_back(22);
        requestDauPsi4.push_back(22);
        requestDauPsi4.push_back(22);
        bool MCpsi=false;
        bool MCphi=false;
        bool MCpsi2=false;
        bool MCpsi3=false;
        bool MCpsi4=false;
        bool MCbs=false;

        double bssvx=0;
        double bssvy=0;
        double bssvz=0;
        double bspvx=genBsCand.vx();
        double bspvy=genBsCand.vy();
        double bspvz=genBsCand.vz();
        double bsmomx=genBsCand.px();
        double bsmomy=genBsCand.py();
        double bsmomz=genBsCand.pz();
        TLorentzVector pmuplus;
        TLorentzVector pmuminus;
        TLorentzVector pkplus;
        TLorentzVector pkminus;
        int BsHasMixed=-1;
        if (MC_particleID==(-1)*genBsCand.daughter(0)->pdgId()) {
          BsHasMixed=1;
          bsRootTree_->BsEndFlavour_=(int)(genBsCand.daughter(0)->pdgId()/531);
          
         // cout<<"Mixing: "<<bsRootTree_->BsEndFlavour_<<endl;
const Candidate * genBsCand2 =genBsCand.daughter(0);
          if(verbose_ == true){
            cout<<"Mixed - First vertex:"<<genBsCand.vx()<<" "<<genBsCand.vy()<<" "<<genBsCand.vz()<<" Second Vertex:"<<genBsCand2->vx()<<" "<<genBsCand2->vy()<<" "<<genBsCand2->vz()<<endl;
          }
          int numBsDau = genBsCand2->numberOfDaughters();
          for (int ghepensimi=0; ghepensimi<numBsDau; ghepensimi++) {
            int numBsDauDau = genBsCand2->daughter(ghepensimi)->numberOfDaughters();
            figlie.push_back(genBsCand2->daughter(ghepensimi)->pdgId());

vector<int> figlieFiglie;
            for (int ghepensimi2=0; ghepensimi2<numBsDauDau; ghepensimi2++) {
              if(verbose_ == true){
                cout<<"             >"<<genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()<<endl;
              }
              figlieFiglie.push_back(abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()));
              if (abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->charge()==-1)  pmuminus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->charge()==1)  pmuplus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==321)  pkplus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==-321)  pkminus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());

            }
            sort(figlieFiglie.begin(),figlieFiglie.end());
            if (figlieFiglie.size()==requestDauPhi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPhi.begin())) MCphi=true;
            if (figlieFiglie.size()==requestDauPsi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi.begin())) MCpsi=true;
            if (figlieFiglie.size()==requestDauPsi2.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi2.begin())) MCpsi2=true;
            if (figlieFiglie.size()==requestDauPsi3.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi3.begin())) MCpsi3=true;
            if (figlieFiglie.size()==requestDauPsi4.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi4.begin())) MCpsi4=true;
            if (abs(genBsCand2->daughter(ghepensimi)->pdgId())==443) {
              if(verbose_ == true){
                cout<<"Last Vertex:"<<genBsCand2->daughter(ghepensimi)->vx()<<" "<<genBsCand2->daughter(ghepensimi)->vy()<<" "<<genBsCand2->daughter(ghepensimi)->vz()<<endl;
              }
              bssvx=genBsCand2->daughter(ghepensimi)->vx();
              bssvy=genBsCand2->daughter(ghepensimi)->vy();
              bssvz=genBsCand2->daughter(ghepensimi)->vz();

            }
          }

 sort(figlie.begin(),figlie.end());
          if (figlie.size()==request.size() && equal(figlie.begin(),figlie.end(),request.begin())) MCbs=true;
          if ((MCpsi || MCpsi2 || MCpsi3) && MCbs) {
            bsRootTree_->ChannelID_=0;
            if(verbose_ == true){ cout<<"Channel!"<<endl;}
          }
          if (MCphi && MCpsi && MCbs) {
            bsRootTree_->ChannelID_=1;
            if(verbose_ == true){ cout<<"Signal!"<<endl;}
          }
          if (MCphi && MCpsi2 && MCbs) {
            bsRootTree_->ChannelID_=2;
            if(verbose_ == true){ cout<<"Strange!"<<endl;}
          }
          if (MCphi && MCpsi3 && MCbs) {
            bsRootTree_->ChannelID_=3;
            if(verbose_ == true){ cout<<"Strange2!"<<endl; }
          }
          if (MCphi && MCpsi4 && MCbs) {
            bsRootTree_->ChannelID_=4;
           if(verbose_ == true){ cout<<"Strange3!"<<endl; }
          }
        }
else {
          bsRootTree_->BsEndFlavour_=(int)(MC_particleID/531);
         // cout<<"No mixing: "<<bsRootTree_->BsEndFlavour_<<endl;
          int numBsDau = genBsCand.numberOfDaughters();
          for (int ghepensimi=0; ghepensimi<numBsDau; ghepensimi++) {
            int numBsDauDau = genBsCand.daughter(ghepensimi)->numberOfDaughters();
            figlie.push_back(genBsCand.daughter(ghepensimi)->pdgId());
            if(verbose_ == true){
              cout<<"Figlie di Bs:"<<genBsCand.daughter(ghepensimi)->pdgId()<<endl;
            }
            vector<int> figlieFiglie;
            for (int ghepensimi2=0; ghepensimi2<numBsDauDau; ghepensimi2++) {
              if(verbose_ == true){
                cout<<"             >"<<genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()<<endl;
              }
              figlieFiglie.push_back(abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()));
              if (abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->charge()==-1)  pmuminus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->charge()==1)  pmuplus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==321)  pkplus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==-321)  pkminus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
            }
            sort(figlieFiglie.begin(),figlieFiglie.end());
            if (figlieFiglie.size()==requestDauPhi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPhi.begin())) MCphi=true;
            if (figlieFiglie.size()==requestDauPsi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi.begin())) MCpsi=true;
            if (figlieFiglie.size()==requestDauPsi2.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi2.begin())) MCpsi2=true;
            if (figlieFiglie.size()==requestDauPsi3.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi3.begin())) MCpsi3=true;
            if (figlieFiglie.size()==requestDauPsi4.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi4.begin())) MCpsi4=true;

            if (abs(genBsCand.daughter(ghepensimi)->pdgId())==443) {
											
               if(verbose_ == true){
                cout<<"Last Vertex:"<<genBsCand.daughter(ghepensimi)->vx()<<" "<<genBsCand.daughter(ghepensimi)->vy()<<" "<<genBsCand.daughter(ghepensimi)->vz()<<endl;
                cout<<"Bs Vertex:"<<genBsCand.vx()<<" "<<genBsCand.vy()<<" "<<genBsCand.vz()<<endl;
                cout<<"Vertex:"<<genBsCand.mother(0)->vx()<<" "<<genBsCand.mother(0)->vy()<<" "<<genBsCand.mother(0)->vz()<<endl;
              }

 bssvx=genBsCand.daughter(ghepensimi)->vx();
              bssvy=genBsCand.daughter(ghepensimi)->vy();
              bssvz=genBsCand.daughter(ghepensimi)->vz();
            }
          }
          sort(figlie.begin(),figlie.end());
          if (figlie.size()==request.size() && equal(figlie.begin(),figlie.end(),request.begin())) MCbs=true;
          if (MCpsi && MCbs) {
            bsRootTree_->ChannelID_=0;
            if(verbose_ == true){ cout<<"Channel!"<<endl; }
          }
          if (MCphi && MCpsi && MCbs) {
            bsRootTree_->ChannelID_=1;
            if(verbose_ == true){ 
				cout<<"Signal!"<<endl;							
				}
          }
          if (MCphi && MCpsi2 && MCbs) {
            bsRootTree_->ChannelID_=2;
            if(verbose_ == true){ cout<<"Strange!"<<endl; }
          }
          if (MCphi && MCpsi3 && MCbs) {
            bsRootTree_->ChannelID_=3;
            if(verbose_ == true){ cout<<"Strange2!"<<endl; }
          }
          if (MCphi && MCpsi4 && MCbs) {
            bsRootTree_->ChannelID_=4;
            if(verbose_ == true){ cout<<"Strange3!"<<endl;}
          }
        }
        if (MC_particleID == 531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) { bsRootTree_->BsIniFlavour_=1; bsRootTree_->BsEndFlavour_=(-1)*BsHasMixed;}
        if (MC_particleID == -531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) { bsRootTree_->BsIniFlavour_=-1; bsRootTree_->BsEndFlavour_=BsHasMixed;}
if (MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) {
          bsRootTree_->SVZpos_[SVindex] = Double_t(bssvz);
          bsRootTree_->FirstBsMCZpos_ = genBsCand.vz();
//vertex z pos where the first Bs generated: genBsCand.vz() i.e primary vertex z
          SVindex++;
          //calculate gen ctau 2D
         // double Lxy2D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy);
          double Lxy2D = sqrt(pow(bssvx-bspvx,2)+pow(bssvy-bspvy,2));
          bsRootTree_->BsLxy2DMC_ = Lxy2D;
          bsRootTree_->BsMassMC_ = genBsCand.mass();
          bsRootTree_->BsCt2DMC_ = Lxy2D * BsPDGMass_/sqrt(bsmomx*bsmomx+bsmomy*bsmomy);
          bsRootTree_->BsCt2DMC_TrueMass_ = Lxy2D * genBsCand.mass()/sqrt(bsmomx*bsmomx+bsmomy*bsmomy);
          bsRootTree_->BsPtMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy);
	//		 cout << "gen Bs Ct "<< bsRootTree_->BsCt2DMC_  << endl;
			 bsRootTree_->genBsPV_x_= genBsCand.vx();
          bsRootTree_->genBsPV_y_= genBsCand.vy();
          bsRootTree_->genBsPV_z_= genBsCand.vz();

          //calculate gen ctau 3D

          //double Lxy3D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy+(bssvz-bspvz)*bsmomz);
          double Lxy3D = sqrt(pow(bssvx-bspvx,2)+pow(bssvy-bspvy,2)+pow(bssvz-bspvz,2));
          bsRootTree_->BsCt3DMC_ = Lxy3D * BsPDGMass_/sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);
          bsRootTree_->BsLxy3DMC_ = Lxy3D;
          bsRootTree_->BsPMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);


//           boosting in JPsi restframe

          TLorentzVector pjpsi;
          pjpsi = pmuplus + pmuminus;
          TLorentzVector pphi;
          pphi = pkplus + pkminus;


			 bsRootTree_->BsJpsiMu1PtMC_ =  pmuminus.Pt();
			 bsRootTree_->BsJpsiMu2PtMC_ =  pmuplus.Pt();
			 bsRootTree_->BsJpsiMu1EtaMC_ =  pmuminus.Eta();
			 bsRootTree_->BsJpsiMu2EtaMC_ =  pmuplus.Eta();
			 bsRootTree_->BsJpsiMu1PhiMC_ =  pmuminus.Phi();
			 bsRootTree_->BsJpsiMu2PhiMC_ =  pmuplus.Phi();
			 bsRootTree_->BsJpsiPtMC_ =  pjpsi.Pt();
			 bsRootTree_->BsJpsiMassMC_ =  pjpsi.M();
			 bsRootTree_->BsPhiPtMC_ =  pphi.Pt();
			 bsRootTree_->BsPhiMassMC_ =  pphi.M();
			 bsRootTree_->BsPhiK1PtMC_ =  pkminus.Pt();
			 bsRootTree_->BsPhiK2PtMC_ =  pkplus.Pt();
			 reco::Vertex::Point BSvec( bssvx-genBsCand.vx(), bssvy-genBsCand.vy(),0 );
 			 math::XYZVector Jpsivec(pjpsi.Px() , pjpsi.Py(), 0);
 			 bsRootTree_->BsCosAlphaMC_ = BSvec.Dot(Jpsivec)/( BSvec.R()*Jpsivec.R());


        //  the betas for the boost

          TVector3 p3_JPsi;
         p3_JPsi = pjpsi.Vect();
          p3_JPsi *= -1./pjpsi.E();

          // the boost matrix

          TLorentzRotation boost_jpsi(p3_JPsi);
          TLorentzVector p_JPsi_JPsi;
          p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);

          // the different momenta in the new frame
          TLorentzVector p_JPsi_muplus;
          TLorentzVector p_JPsi_Kplus;
          TLorentzVector p_JPsi_phi;
          p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
          p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);
         p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);

           //the 3-momenta

          TVector3 p3_JPsi_muplus;
          p3_JPsi_muplus = p_JPsi_muplus.Vect();
          TVector3 p3_JPsi_Kplus;
          p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
          TVector3 p3_JPsi_phi;
        p3_JPsi_phi = p_JPsi_phi.Vect();



          // coordinate system

          TVector3 x,y,z;
          x = p3_JPsi_phi.Unit();
          y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
          y = y.Unit();
          z = x.Cross(y);

          double angle_costhetaMC=p3_JPsi_muplus.Unit() * z;
          bsRootTree_->BscosthetaMC_= angle_costhetaMC;
         //      cout<<"angle_costhetaMC"<<angle_costhetaMC<<endl;
          double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
          double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
          double angle_phiMC = TMath::ACos(cos_phi);
          if (sin_phi < 0){
            angle_phiMC =  -angle_phiMC;
          }
          bsRootTree_->BsphiMC_ = angle_phiMC;

          //boosting in phi restframe

          TVector3 p3_phi;
          p3_phi = pphi.Vect();
          p3_phi *= -1./pphi.E();

           //the boost matrix

          TLorentzRotation boost_phi(p3_phi);
          TLorentzVector p_phi_phi;
          p_phi_phi = boost_phi.VectorMultiplication(pphi);


          // the different momenta in the new frame

          TLorentzVector p_phi_Kplus;
          TLorentzVector p_phi_JPsi;
          TLorentzVector p_phi_Bs;
          p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
          p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);

          // the 3-momenta

          TVector3 p3_phi_Kplus;
          p3_phi_Kplus = p_phi_Kplus.Vect();
          TVector3 p3_phi_JPsi;
          p3_phi_JPsi = p_phi_JPsi.Vect();
          bsRootTree_->BscospsiMC_ = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
        }
      }      //END Bs MC matching

// if this is a real B decay (no B mesons as daughters
 if(hasBDaughter == 0){
//count the number of B decays, should be equal two, for bbbar events
 iNumberOfBdecays++;
        int arrayIndex = iNumberOfBdecays - 1; // array index starts at zero
// protect array bounds
if(arrayIndex>=9) break;
        if(MC_particleID==(-1)*genBsCand.mother(0)->pdgId()){
//	  cout << "Mixing! " << endl;
 bsRootTree_->BmesonsId_[arrayIndex] = genBsCand.mother(0)->pdgId();
//	  cout << "Original flavour " << genBsCand.mother(0)->pdgId() << endl;
 }
        else bsRootTree_->BmesonsId_[arrayIndex] = MC_particleID;
        bsRootTree_->GenNumberOfDaughters_[arrayIndex] = numBsDaughters;
// generator variables
        bsRootTree_->BMMC_[arrayIndex] = genBsCand.mass();
        bsRootTree_->BPtMC_[arrayIndex] = genBsCand.pt();
        bsRootTree_->BPxMC_[arrayIndex] = genBsCand.px();
        bsRootTree_->BPyMC_[arrayIndex] = genBsCand.py();
        bsRootTree_->BPzMC_[arrayIndex] = genBsCand.pz();
        bsRootTree_->BEtaMC_[arrayIndex] = genBsCand.eta();
        bsRootTree_->BPhiMC_[arrayIndex] = genBsCand.phi();
        bsRootTree_->BVtxMC_x_[arrayIndex] = genBsCand.mother(0)->vx();
        bsRootTree_->BVtxMC_y_[arrayIndex] = genBsCand.mother(0)->vy();
        bsRootTree_->BVtxMC_z_[arrayIndex] = genBsCand.mother(0)->vz();

//generated primary vertex
 if (abs(MC_particleID)== 531){
          bsRootTree_->genBsVtx_x_= genBsCand.mother(0)->vx();
          bsRootTree_->genBsVtx_y_= genBsCand.mother(0)->vy();
          bsRootTree_->genBsVtx_z_= genBsCand.mother(0)->vz();
 }

        for(int j = 0; j < numBsDaughters; ++ j) {
          if(j>=14) break; // protect array bounds
	  const Candidate * Bsdau = genBsCand.daughter( j );
          if (abs(Bsdau->pdgId()) == 443) Jpsi = Bsdau;
          if (abs(Bsdau->pdgId()) == 333) Phi = Bsdau;

// generator variables
          bsRootTree_->BDauIdMC_[arrayIndex][j] = Bsdau->pdgId();
          bsRootTree_->BDauMMC_[arrayIndex][j] = Bsdau->mass();
          bsRootTree_->BDauPtMC_[arrayIndex][j] = Bsdau->pt();
          bsRootTree_->BDauPzMC_[arrayIndex][j] = Bsdau->pz();
          bsRootTree_->BDauEtaMC_[arrayIndex][j] = Bsdau->eta();
          bsRootTree_->BDauPhiMC_[arrayIndex][j] = Bsdau->phi();
          bsRootTree_->BSVtxMC_x_[arrayIndex]   = Bsdau->vx();
          bsRootTree_->BSVtxMC_y_[arrayIndex]   = Bsdau->vy();
          bsRootTree_->BSVtxMC_z_[arrayIndex]   = Bsdau->vz();
 //Generated secondary vertex
 if ( abs(Bsdau->pdgId())== 443){
            bsRootTree_->genBsSVtx_x_= Bsdau->vx();
            bsRootTree_->genBsSVtx_y_= Bsdau->vy();
            bsRootTree_->genBsSVtx_z_= Bsdau->vz();

 }
// daughter of daughter (muons, kaons in case of jpsi phi)
int numBsDaughtersDaughters = Bsdau->numberOfDaughters();
          bsRootTree_->GenNumberOfDaughtersDaughters_[arrayIndex][j] = numBsDaughtersDaughters;
          for(int k=0; k< numBsDaughtersDaughters; k++){
            if(k>=9) break; //protect array bounds
	    const Candidate * Bsdaudau = Bsdau->daughter(k);
            if ( Bsdaudau->pdgId() == -13) mup = Bsdaudau;
            if ( Bsdaudau->pdgId() == 13) mum = Bsdaudau;
            if ( Bsdaudau->pdgId() == 321) Kp = Bsdaudau;
            if ( Bsdaudau->pdgId() == -321) Km = Bsdaudau;
// generator variables
bsRootTree_->BDauDauIdMC_[arrayIndex][j][k] = Bsdaudau->pdgId();
            bsRootTree_->BDauDauMMC_[arrayIndex][j][k] = Bsdaudau->mass();
            bsRootTree_->BDauDauPtMC_[arrayIndex][j][k] = Bsdaudau->pt();
            bsRootTree_->BDauDauPzMC_[arrayIndex][j][k] = Bsdaudau->pz();
            bsRootTree_->BDauDauEtaMC_[arrayIndex][j][k] = Bsdaudau->eta();
            bsRootTree_->BDauDauPhiMC_[arrayIndex][j][k] = Bsdaudau->phi();

            if (abs(MC_particleID)== 531 && numBsDaughters == 2 && Jpsi && Phi && mup && mum && Kp && Km){
              TLorentzVector pmuplus;
              TLorentzVector pmuminus;
              TLorentzVector pkplus;
              TLorentzVector pkminus;

// extra check
              if (mup)
                pmuplus.SetXYZM(mup->px(),mup->py(),mup->pz(),mup->mass());
              if (mum)
                pmuminus.SetXYZM(mum->px(),mum->py(),mum->pz(),mum->mass());
              if (Kp)
                pkplus.SetXYZM(Kp->px(),Kp->py(),Kp->pz(),Kp->mass());
              if (Km)
                pkminus.SetXYZM(Km->px(),Km->py(),Km->pz(),Km->mass());


               //boosting in JPsi restframe
              TLorentzVector pjpsi;
              pjpsi = pmuplus + pmuminus;
              TLorentzVector pphi;
              pphi = pkplus + pkminus;



              //the betas for the boost
              TVector3 p3_JPsi;
              p3_JPsi = pjpsi.Vect();
              p3_JPsi *= -1./pjpsi.E();

              // the boost matrix
              TLorentzRotation boost_jpsi(p3_JPsi);
              TLorentzVector p_JPsi_JPsi;
              p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);

              // the different momenta in the new frame

              TLorentzVector p_JPsi_muplus;
              TLorentzVector p_JPsi_Kplus;
              TLorentzVector p_JPsi_phi;

              p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
              p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);
              p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);

               //the 3-momenta

              TVector3 p3_JPsi_muplus;
              p3_JPsi_muplus = p_JPsi_muplus.Vect();
              TVector3 p3_JPsi_Kplus;
              p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
              TVector3 p3_JPsi_phi;
              p3_JPsi_phi = p_JPsi_phi.Vect();

              //coordinate system

              TVector3 x,y,z;
              x = p3_JPsi_phi.Unit();
              y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
              y = y.Unit();
              z = x.Cross(y);

              double angle_costhetaMC=p3_JPsi_muplus.Unit() * z;
              bsRootTree_->costhetaMC_[arrayIndex] = angle_costhetaMC;
              double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
              double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
              double angle_phiMC = TMath::ACos(cos_phi);
              if (sin_phi < 0){
                angle_phiMC =  -angle_phiMC;
              }
              bsRootTree_->phiMC_[arrayIndex] = angle_phiMC;

            //boosting in phi rest frame
              TVector3 p3_phi;
              p3_phi = pphi.Vect();
              p3_phi *= -1./pphi.E();

              // the boost matrix
              TLorentzRotation boost_phi(p3_phi);
              TLorentzVector p_phi_phi;
              p_phi_phi = boost_phi.VectorMultiplication(pphi);

             // the different momenta in the new frame

              TLorentzVector p_phi_Kplus;
              TLorentzVector p_phi_JPsi;
              TLorentzVector p_phi_Bs;
              p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
              p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);

               //the 3-momenta

              TVector3 p3_phi_Kplus;
              p3_phi_Kplus = p_phi_Kplus.Vect();
              TVector3 p3_phi_JPsi;
              p3_phi_JPsi = p_phi_JPsi.Vect();
              bsRootTree_->cospsiMC_[arrayIndex] = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
//cout << "(costheta,phi,cospsi): "<< p3_JPsi_muplus.Unit() * z << "," << angle_phi << "," << -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit() << ")" << endl;
 }
          }// loop Bs daughters daughters
          } // loop Bs daughters

//calculate ctau 2D
double Lxy2D = ((bsRootTree_->BSVtxMC_x_[arrayIndex] - bsRootTree_->BVtxMC_x_[arrayIndex])*bsRootTree_->BPxMC_[arrayIndex]+
          (bsRootTree_->BSVtxMC_y_[arrayIndex] - bsRootTree_->BVtxMC_y_[arrayIndex])*bsRootTree_->BPyMC_[arrayIndex]);
          bsRootTree_->BCt_MC2D_[arrayIndex] = Lxy2D * bsRootTree_->BMMC_[arrayIndex]/(bsRootTree_->BPtMC_[arrayIndex]*bsRootTree_->BPtMC_[arrayIndex]);

//calculate gen ctau 3D
 double Lxy3D = ((bsRootTree_->BSVtxMC_x_[arrayIndex] - bsRootTree_->BVtxMC_x_[arrayIndex])*bsRootTree_->BPxMC_[arrayIndex]+
          (bsRootTree_->BSVtxMC_y_[arrayIndex] - bsRootTree_->BVtxMC_y_[arrayIndex])*bsRootTree_->BPyMC_[arrayIndex]+
          (bsRootTree_->BSVtxMC_z_[arrayIndex] - bsRootTree_->BVtxMC_z_[arrayIndex])*bsRootTree_->BPzMC_[arrayIndex]);
          bsRootTree_->BCt_MC3D_[arrayIndex] = Lxy3D * bsRootTree_->BMMC_[arrayIndex]/
          (bsRootTree_->BPtMC_[arrayIndex]*bsRootTree_->BPtMC_[arrayIndex]+
          bsRootTree_->BPzMC_[arrayIndex]*bsRootTree_->BPzMC_[arrayIndex]);

//calculate gen ctau
double deltaX =  bsRootTree_->BSVtxMC_x_[arrayIndex] - 	bsRootTree_->BVtxMC_x_[arrayIndex];
          double deltaY =  bsRootTree_->BSVtxMC_y_[arrayIndex] - 	bsRootTree_->BVtxMC_y_[arrayIndex];
          bsRootTree_->BLxy_MC_[arrayIndex] = sqrt( deltaX*deltaX + deltaY*deltaY);
          if(deltaX * genBsCand.px() + deltaY * genBsCand.py() < 0 )  bsRootTree_->BLxy_MC_[arrayIndex] = -1. *  bsRootTree_->BLxy_MC_[arrayIndex];
          bsRootTree_->BCt_MC_[arrayIndex] = bsRootTree_->BLxy_MC_[arrayIndex] * bsRootTree_->BMMC_[arrayIndex] / bsRootTree_->BPtMC_[arrayIndex];
          }
        }
// check if there is a Jpsi (prompt or non-prompt) in the event
 if(absMC_particleID == 443 ) bsRootTree_->isGenJpsiEvent_ = 1;
}
//END loop over all particles
bsRootTree_->GenNumberOfBdecays_ = iNumberOfBdecays;
}
/*

bool BsToJpsiPhiAnalysis::selGlobalMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx) {
  TrackRef iTrack = aMuon.innerTrack();
  const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(aMuon.originalObject());
  const reco::HitPattern& p = iTrack->hitPattern();
TrackRef gTrack = aMuon.globalTrack();
return (// isMuonInAccept(aMuon) &&
          aMuon.isPFMuon() &&
          gTrack->chi2()/gTrack->ndof() < 10.0 &&
          gTrack->hitPattern().numberOfValidMuonHits()>0 &&
          rmu1->numberOfMatchedStations()>1 &&
          p.numberOfValidMuonHits() > 0 &&
          p.trackerLayersWithMeasurement()>5);

}



bool BsToJpsiPhiAnalysis::selTrackerMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx) {
  TrackRef iTrack = aMuon.innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
bool trackOK = false;
 trackOK = (iTrack->found() > 10);
return (// isMuonInAccept(aMuon) &&
  	  trackOK &&
   	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
  	  aMuon.muonID("TrackerMuonArbitrated") &&
 	  aMuon.muonID("TMOneStationTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

*/

void BsToJpsiPhiAnalysis::setFitParKK(RefCountedKinematicTree& myTree)
{
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();

  for(int i=0; i<7; i++)  bsRootTree_->K1Fit_par_[i] = bs_par1[i];
  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->K1Fit_sigX_ = sqrt(bs_err1(0,0));
  bsRootTree_->K1Fit_sigY_ = sqrt(bs_err1(1,1));
  bsRootTree_->K1Fit_sigZ_ = sqrt(bs_err1(2,2));
  bsRootTree_->K1Fit_sigPX_ = sqrt(bs_err1(3,3));
  bsRootTree_->K1Fit_sigPY_ = sqrt(bs_err1(4,4));
  bsRootTree_->K1Fit_sigPZ_ = sqrt(bs_err1(5,5));
// first particle: kaon 2
 AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();

  for(int i=0; i<7; i++) bsRootTree_->K2Fit_par_[i] = bs_par2[i];

  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->K2Fit_sigX_ = sqrt(bs_err2(0,0));
  bsRootTree_->K2Fit_sigY_ = sqrt(bs_err2(1,1));
  bsRootTree_->K2Fit_sigZ_ = sqrt(bs_err2(2,2));
  bsRootTree_->K2Fit_sigPX_ = sqrt(bs_err2(3,3));
  bsRootTree_->K2Fit_sigPY_ = sqrt(bs_err2(4,4));
  bsRootTree_->K2Fit_sigPZ_ = sqrt(bs_err2(5,5));
}
/*
void BsToJpsiPhiAnalysis::setFitParHyp1(RefCountedKinematicTree& myTree)
{
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();

  for(int i=0; i<7; i++)  bsRootTree_->BdK1_kpi_par_Hyp1_[i] = bs_par1[i];

  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK1_kpi_sigX_Hyp1_ = sqrt(bs_err1(0,0));
  bsRootTree_->BdK1_kpi_sigY_Hyp1_ = sqrt(bs_err1(1,1));
  bsRootTree_->BdK1_kpi_sigZ_Hyp1_ = sqrt(bs_err1(2,2));
  bsRootTree_->BdK1_kpi_sigPX_Hyp1_ = sqrt(bs_err1(3,3));
  bsRootTree_->BdK1_kpi_sigPY_Hyp1_ = sqrt(bs_err1(4,4));
  bsRootTree_->BdK1_kpi_sigPZ_Hyp1_ = sqrt(bs_err1(5,5));
 AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();

  for(int i=0; i<7; i++) bsRootTree_->BdK2_kpi_par_Hyp1_[i] = bs_par2[i];

  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK2_kpi_sigX_Hyp1_ = sqrt(bs_err2(0,0));
  bsRootTree_->BdK2_kpi_sigY_Hyp1_ = sqrt(bs_err2(1,1));
  bsRootTree_->BdK2_kpi_sigZ_Hyp1_ = sqrt(bs_err2(2,2));
  bsRootTree_->BdK2_kpi_sigPX_Hyp1_ = sqrt(bs_err2(3,3));
  bsRootTree_->BdK2_kpi_sigPY_Hyp1_ = sqrt(bs_err2(4,4));
  bsRootTree_->BdK2_kpi_sigPZ_Hyp1_ = sqrt(bs_err2(5,5));

}
void BsToJpsiPhiAnalysis::setFitParHyp2(RefCountedKinematicTree& myTree)
{
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
 // first particle: kaon 1
AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++)  bsRootTree_->BdK1_kpi_par_Hyp2_[i] = bs_par1[i];

  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK1_kpi_sigX_Hyp2_ = sqrt(bs_err1(0,0));
  bsRootTree_->BdK1_kpi_sigY_Hyp2_ = sqrt(bs_err1(1,1));
  bsRootTree_->BdK1_kpi_sigZ_Hyp2_ = sqrt(bs_err1(2,2));
  bsRootTree_->BdK1_kpi_sigPX_Hyp2_ = sqrt(bs_err1(3,3));
  bsRootTree_->BdK1_kpi_sigPY_Hyp2_ = sqrt(bs_err1(4,4));
  bsRootTree_->BdK1_kpi_sigPZ_Hyp2_ = sqrt(bs_err1(5,5));
// first particle: kaon 2
AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();

  for(int i=0; i<7; i++) bsRootTree_->BdK2_kpi_par_Hyp2_[i] = bs_par2[i];

  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK2_kpi_sigX_Hyp2_ = sqrt(bs_err2(0,0));
  bsRootTree_->BdK2_kpi_sigY_Hyp2_ = sqrt(bs_err2(1,1));
  bsRootTree_->BdK2_kpi_sigZ_Hyp2_ = sqrt(bs_err2(2,2));
  bsRootTree_->BdK2_kpi_sigPX_Hyp2_ = sqrt(bs_err2(3,3));
  bsRootTree_->BdK2_kpi_sigPY_Hyp2_ = sqrt(bs_err2(4,4));
  bsRootTree_->BdK2_kpi_sigPZ_Hyp2_ = sqrt(bs_err2(5,5));
}
*/
bool  BsToJpsiPhiAnalysis::MCmatching(const Candidate & track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles, int &K1mcId, int &K1momId, int &K1gmomId, int condMom, int condGMom){

  if(!isMCstudy_ ) return 0;
  bool K1Truth = 0;

  double MinDRK=999.;
  K1mcId = -9999999;
  K1momId = -9999999;
  K1gmomId = -9999999;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    double DeltaRK1 = deltaR(p.eta(), p.phi(), track1.eta(), track1.phi() );
    if(DeltaRK1<MinDRK && DeltaRK1<0.05){
      MinDRK=DeltaRK1;
      K1mcId=p.pdgId();
      if(p.mother()!=0) K1momId=p.mother()->pdgId();
      if(p.mother()!=0 && p.mother()->mother()!=0) {K1gmomId=p.mother()->mother()->pdgId(); }
      if (abs(K1momId)==condMom && abs(K1gmomId)==condGMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }
cout << "Muon id " << K1mcId << endl;
return K1Truth;
}

bool  BsToJpsiPhiAnalysis::MCmatchingJpsi (const Candidate & track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles, int &K1mcId, int &K1momId, int condMom){
  if(!isMCstudy_ ) return 0;
  bool K1Truth = 0;
  double MinDRK=999.;

  K1mcId = -9999999;
  K1momId = -9999999;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    double DeltaRK1 = deltaR(p.eta(), p.phi(), track1.eta(), track1.phi() );

    if(DeltaRK1<MinDRK && DeltaRK1<0.05){
      MinDRK=DeltaRK1;
      K1mcId=p.pdgId();

      if(p.mother()!=0) K1momId=p.mother()->pdgId();
      if (abs(K1momId)==condMom ) K1Truth = 1;
      else K1Truth = 0;
    }
  }
 return K1Truth;
}
/*
bool BsToJpsiPhiAnalysis::MCmatchingBplusK(const Candidate & track1, edm::Handle<edm::View<reco::GenParticle>>& genParticles, int &K1mcId, int &K1momId, int condMom){
  if(!isMCstudy_ ) return 0;
  bool K1Truth = 0;
  double MinDRK=999.;

  K1mcId = -9999999;
  K1momId = -9999999;

  for(size_t i = 0; i < genParticles->size(); ++ i){
    const GenParticle & p = (*genParticles)[i];
    double DeltaRK1 = deltaR(p.eta(), p.phi(), track1.eta(), track1.phi() );

    if(DeltaRK1<MinDRK && DeltaRK1<0.05){
      MinDRK=DeltaRK1;
      K1mcId=p.pdgId();

      if(p.mother()!=0) {K1momId=p.mother()->pdgId(); }
      if (abs(K1momId)==condMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }

return K1Truth;
}*/
double BsToJpsiPhiAnalysis::CalculateCtErrvertex(const edm::EventSetup& iSetup, reco::Vertex PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass){


VertexDistanceXY d;

   const GlobalPoint myPV(PV.x(),PV.y(),PV.z());
   const GlobalPoint mySV(SV.position().x(),SV.position().y(),SV.position().z());
	const GlobalError myPVErr = PV.covariance();
	const GlobalError mySVErr = SV.positionError();
   Measurement1D VtxDist = d.distance(VertexState(mySV,mySVErr),VertexState(myPV,myPVErr));
   double VtxDistErr = VtxDist.error();
		double pt = sqrt( pow(Ptvec.y(),2.)  + pow(Ptvec.x(),2.)  );
   if( pt != 0 && VtxDist.value() != 0 ){

   double cos = ( ( SV.position().x() - PV.x() )*Ptvec.x() + (SV.position().y() - PV.y() )*Ptvec.y() ) / ( pt * VtxDist.value() );
                 
   TVector LengthVector(2);
       
   LengthVector(0)= SV.position().x() - PV.x();      
   LengthVector(1)= SV.position().y() - PV.y();

	double firstTerm2 = pow( (Bmass * VtxDistErr * abs(cos) / pt) , 2. ); // error proportional to Lxy
	double secondTerm2 = pow( (Bmass * abs(cos) / (pt*pt) ), 2. ) * SVcovariance2D.Similarity(LengthVector)  ; // error proportional to p_T
	double ctErr = sqrt(	firstTerm2 + secondTerm2 );
	return ctErr;	
 }   
  else {return -9999999; }
}

double BsToJpsiPhiAnalysis::CalculateCtErrbeamspot(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass){


VertexDistanceXY d;

   const GlobalPoint myPV(PV->x0(),PV->y0(),PV->z0());
   const GlobalPoint mySV(SV.position().x(),SV.position().y(),SV.position().z());
	const GlobalError myPVErr = PV->covariance3D();
	const GlobalError mySVErr = SV.positionError();
   Measurement1D VtxDist = d.distance(VertexState(mySV,mySVErr),VertexState(myPV,myPVErr));
   double VtxDistErr = VtxDist.error();
		double pt = sqrt( pow(Ptvec.y(),2.)  + pow(Ptvec.x(),2.)  );
   if( pt != 0 && VtxDist.value() != 0 ){

   double cos = ( ( SV.position().x() - PV->x0() )*Ptvec.x() + (SV.position().y() - PV->y0() )*Ptvec.y() ) / ( pt * VtxDist.value() );
                 
   TVector LengthVector(2);
        
   LengthVector(0)= SV.position().x() - PV->x0();      
   LengthVector(1)= SV.position().y() - PV->y0();

	double firstTerm2 = pow( (Bmass * VtxDistErr * abs(cos) / pt) , 2. ); // error proportional to Lxy
	double secondTerm2 = pow( (Bmass * abs(cos) / (pt*pt) ), 2. ) * SVcovariance2D.Similarity(LengthVector)  ; // error proportional to p_T
	double ctErr = sqrt(	firstTerm2 + secondTerm2 );
	return ctErr;	
 }   
  else {return -9999999; }
}
/*
reco::Vertex BsToJpsiPhiAnalysis::reVertex(const edm::EventSetup& iSetup, reco::BeamSpot vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, TrackRef trk3, TrackRef trk4){

        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
        std::vector<reco::TransientTrack> vrtxRefit;
        for (reco::Vertex::trackRef_iterator trackvertex=RecVtx.tracks_begin(); trackvertex!=RecVtx.tracks_end(); trackvertex++){
          
          TrackRef tref = trackvertex->castTo<TrackRef>();
          if (!tref.isNonnull()) continue;
          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mu1.originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mu2.originalObject());

			 	if( rmu1->track().key()!=trackvertex->key() && rmu2->track().key()!=trackvertex->key() && trk3.key()!=trackvertex->key() && trk4.key()!=trackvertex->key()  ){
            TransientTrack tTrk = (*theB).build(&tref); 
            if (tTrk.isValid())  vrtxRefit.push_back(tTrk);
          }
        }
        AdaptiveVertexFitter avf;
        TransientVertex newVtx = avf.vertex(vrtxRefit, vertexBeamSpot);
        if (newVtx.isValid()) {
			
          return reco::Vertex(newVtx);
        } else {
          return RecVtx;
        }
}


///

reco::Vertex BsToJpsiPhiAnalysis::reVertex(const edm::Event& theEvent, const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, reco::Track & rtrk1, reco::Track & rtrk2){


edm::Handle< View<reco::Vertex> > refvtx;
theEvent.getByToken(primaryvertexTok, refvtx);

        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
        std::vector<reco::TransientTrack> vrtxRefit;
//In veretx we dont trackref any more , instead of that we can have 
//Track refittedTrack(const TrackRef & track) const;
//But for this reffitted tracks we have to chek overlaps with muon tracks and have to returns vertex index
for(View<reco::Vertex>::const_iterator pvtx = refvtx->begin(); pvtx != refvtx->end(); pvtx++)
{     Track trkReffit = pvtx->refittedTrack();
      TransientTrack tTrk = (*theB).build(&trkReffit);
      if (tTrk.isValid())  vrtxRefit.push_back(tTrk);
}
         AdaptiveVertexFitter avf;
         TransientVertex newVtx = avf.vertex(vrtxRefit);       

if (newVtx.isValid()) {
                        cout << "vtxCL " <<  TMath::Prob( newVtx.totalChiSquared(),(int)( newVtx.degreesOfFreedom()))  << endl;
                        return reco::Vertex(newVtx);
                      }  
                      else {
                             return RecVtx;
                           }
}
       for (reco::Vertex::trackRef_iterator trackvertex=RecVtx.tracks_begin(); trackvertex!=RecVtx.tracks_end(); trackvertex++){
TrackRef tref = trackvertex->castTo<TrackRef>();
          if (!tref.isNonnull()) continue;
          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mu1.originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mu2.originalObject());
if( rmu1->track().key()!=trackvertex->key() && rmu2->track().key()!=trackvertex->key() && trk3.key()!=trackvertex->key() && trk4.key()!=trackvertex->key()  ){
            TransientTrack tTrk = (*theB).build(&tref); //  fpTTB->build(*(*itTBR));i
            if (tTrk.isValid())  vrtxRefit.push_back(tTrk);
          }
        }
//virtual CachingVertex<5> vertex( const std::vector<reco::TransientTrack> & ) const;
        AdaptiveVertexFitter avf;
         TransientVertex newVtx = avf.vertex(vrtxRefit);
       // TransientVertex newVtx = avf.vertex(vrtxRefit, vertexBeamSpot);
        if (newVtx.isValid()) {
 cout << "vtxCL " <<  TMath::Prob( newVtx.totalChiSquared(),(int)( newVtx.degreesOfFreedom()))	<< endl;
          return reco::Vertex(newVtx);
        } else {
          return RecVtx;
        }
}

reco::Vertex BsToJpsiPhiAnalysis::reVertex2(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, TrackRef trk3, double &vtxCL){
        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
        std::vector<reco::TransientTrack> vrtxRefit;
        for (reco::Vertex::trackRef_iterator trackvertex=RecVtx.tracks_begin(); trackvertex!=RecVtx.tracks_end(); trackvertex++){

TrackRef tref = trackvertex->castTo<TrackRef>();
          if (!tref.isNonnull()) continue;
          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mu1.originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mu2.originalObject());

			 	if( rmu1->track().key()!=trackvertex->key() && rmu2->track().key()!=trackvertex->key() && trk3.key()!=trackvertex->key() ){
            TransientTrack tTrk = (*theB).build(&tref); //  fpTTB->build(*(*itTBR));i
            if (tTrk.isValid())  vrtxRefit.push_back(tTrk);
          }
        }
        AdaptiveVertexFitter avf;
        TransientVertex newVtx = avf.vertex(vrtxRefit);
       // TransientVertex newVtx = avf.vertex(vrtxRefit, vertexBeamSpot);

if (newVtx.isValid()) {
			 vtxCL = TMath::Prob( newVtx.totalChiSquared(),(int)( newVtx.degreesOfFreedom()));
cout << "vtxCL " <<  TMath::Prob( newVtx.totalChiSquared(),(int)( newVtx.degreesOfFreedom()))	<< endl;
return reco::Vertex(newVtx);
        } else {
          return RecVtx;
        }
}


double BsToJpsiPhiAnalysis::MuonChargeCone(const edm::Event& theEvent, TrackRef muTrackRef, const double Dr, const double KExp, bool IncludingMuon)
{
  if( Dr < 0.)
  {
    std::cout << "E R R O R! DeltaR must be positive.\n";
    exit(1);
  }
  if( KExp < 0.)
  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }

  double num = 0.;
  double den = 0.;

//edm::Handle<View<pat::PackedCandidate>> allTracks;
 //  iEvent.getByToken(trackLabelK, allTracks);


   edm::Handle<View<pat::PackedCandidate>> recoTracks;
   theEvent.getByToken(trackLabelK, recoTracks);


//  edm::Handle<reco::TrackCollection> recoTracks;
//  theEvent.getByLabel("generalTracks", recoTracks);
  for(size_t itTrack = 0; itTrack < recoTracks->size(); ++itTrack)
  {
    reco::TrackRef trkRef(recoTracks, itTrack);
    if( trkRef->pt()        < 0.5 ) continue;
    if( fabs(trkRef->eta()) > 2.5 ) continue;
    if( trkRef->numberOfValidHits() <= 5 ) continue;
    if ( !IncludingMuon && ( trkRef == muTrackRef) ) continue;
    double DeltaR = deltaR(muTrackRef->eta(), muTrackRef->phi(), trkRef->eta(), trkRef->phi());
    if(DeltaR < Dr)
    {
      num = num + trkRef->charge()*pow(trkRef->pt(),KExp);
      den = den + pow(trkRef->pt(),KExp);
    }
  }
  return (den > 0.) ? num/den : -9999999;
}

double BsToJpsiPhiAnalysis::MuonChargeConeWrtPV(const edm::Event& theEvent, TrackRef muTrackRef, Vertex PVtx, const double Dr, const double KExp, bool IncludingMuon)

{
  if( Dr < 0.)
  {
    std::cout << "E R R O R! DeltaR must be positive.\n";
    exit(1);
  }
  if( KExp < 0.)
  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);

  }
  double num = 0.;
  double den = 0.;
  for (reco::Vertex::trackRef_iterator trackvertex = PVtx.tracks_begin(); trackvertex != PVtx.tracks_end(); trackvertex++)

  {

    const reco::Track & VtxTrack = *(trackvertex->get());
    if( VtxTrack.pt()        < 0.5 ) continue;
    if( fabs(VtxTrack.eta()) > 2.5 ) continue;
    if( VtxTrack.numberOfValidHits() <= 5 ) continue;

    if( !IncludingMuon                                             &&
        ( VtxTrack.charge()     - muTrackRef->charge()     == 0  ) &&
        ( fabs(VtxTrack.pt()    - muTrackRef->pt())       < 1e-6 ) &&
        ( fabs(VtxTrack.eta()   - muTrackRef->eta())      < 1e-6 ) &&
        ( fabs(VtxTrack.phi()   - muTrackRef->phi())      < 1e-6 ) ) continue;
    double DeltaR = deltaR(muTrackRef->eta(), muTrackRef->phi(), VtxTrack.eta(), VtxTrack.phi());
    if(DeltaR < Dr)
    {
      num = num + VtxTrack.charge()*pow(VtxTrack.pt(),KExp);
      den = den + pow(VtxTrack.pt(),KExp);
    }
  }
  return (den > 0.) ? num/den : -9999999;
}


double BsToJpsiPhiAnalysis::ElectronChargeConeWrtPV(const edm::Event& theEvent, reco::GsfTrackRef eleGsfTrackRef, Vertex PVtx, const double Dr, const double KExp, bool IncludingElectron)

{ if( Dr < 0.)
  {

    std::cout << "E R R O R! DeltaR must be positive.\n";
    exit(1);

  }
  if( KExp < 0.)

  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);

  }


  double num = 0.;
  double den = 0.;

  size_t ClosestTrack = 0;
  double minDeltaR    = 999.;
   edm::Handle<View<pat::PackedCandidate>> recoTracks;
   theEvent.getByToken(trackLabelK, recoTracks);

//  edm::Handle<reco::TrackCollection> recoTracks;
//  theEvent.getByLabel("generalTracks", recoTracks);
  for(size_t itTrack = 0; itTrack < recoTracks->size(); ++itTrack)
  {
    reco::TrackRef trkRef(recoTracks, itTrack);
    double DeltaR = deltaR(eleGsfTrackRef->eta(), eleGsfTrackRef->phi(), trkRef->eta(), trkRef->phi());
    if( DeltaR < minDeltaR && eleGsfTrackRef->charge() == trkRef->charge()){
      minDeltaR    = DeltaR;
      ClosestTrack = itTrack;
    }
  }
  reco::TrackRef ClosestTrkRef(recoTracks, ClosestTrack);
  for (reco::Vertex::trackRef_iterator trackvertex = PVtx.tracks_begin(); trackvertex != PVtx.tracks_end(); trackvertex++)
  {
    const reco::Track & VtxTrack = *(trackvertex->get());
    if( !IncludingElectron                                            &&
        ( VtxTrack.charge()     - ClosestTrkRef->charge()     == 0  ) &&
        ( fabs(VtxTrack.pt()    - ClosestTrkRef->pt())       < 1e-6 ) &&
        ( fabs(VtxTrack.eta()   - ClosestTrkRef->eta())      < 1e-6 ) &&
        ( fabs(VtxTrack.phi()   - ClosestTrkRef->phi())      < 1e-6 ) ) continue;

    if( VtxTrack.pt()        < 0.5 ) continue;
    if( fabs(VtxTrack.eta()) > 2.5 ) continue;
    if( VtxTrack.numberOfValidHits() <= 5 ) continue;
    double DeltaR = deltaR(eleGsfTrackRef->eta(), eleGsfTrackRef->phi(), VtxTrack.eta(), VtxTrack.phi());
    if(DeltaR < )
    {
      num = num + VtxTrack.charge()*pow(VtxTrack.pt(),KExp);
      den = den + pow(VtxTrack.pt(),KExp);
    }
  }
  return (den > 0.) ? num/den : -9999999;
}

double BsToJpsiPhiAnalysis::ElectronChargeCone(const edm::Event& theEvent, reco::GsfTrackRef eleGsfTrackRef, const double Dr, const double KExp, bool IncludingElectron)
{
  if( Dr < 0.)
  {
    std::cout << "E R R O R! DeltaR must be positive.\n";
    exit(1);
  }
  if( KExp < 0.)

  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }
  double num = 0.;
  double den = 0.;

  edm::Handle<reco::TrackCollection> recoTracks;
  theEvent.getByLabel("generalTracks", recoTracks);

  size_t ClosestTrack = 0;
  double minDeltaR    = 999.;
  std::vector<reco::TrackRef> TrackRefInCone;
  for(size_t itTrack = 0; itTrack < recoTracks->size(); ++itTrack)

  {
    reco::TrackRef trkRef(recoTracks, itTrack);
    double DeltaR = deltaR(eleGsfTrackRef->eta(), eleGsfTrackRef->phi(), trkRef->eta(), trkRef->phi());
    if( DeltaR < minDeltaR  && eleGsfTrackRef->charge() == trkRef->charge() ){
      minDeltaR    = DeltaR;
      ClosestTrack = itTrack;

    }
    if( trkRef->pt()        < 0.5 )             continue;
    if( fabs(trkRef->eta()) > 2.5 )             continue;
    if( trkRef->numberOfValidHits() <= 5 )      continue;
    TrackRefInCone.push_back(trkRef);
  }
  reco::TrackRef ClosestTrkRef(recoTracks, ClosestTrack);
  for(unsigned int iTrk = 0 ; iTrk < TrackRefInCone.size() ; iTrk++ ){
    if( !IncludingElectron && TrackRefInCone[iTrk] == ClosestTrkRef) continue;
    num = num + TrackRefInCone[iTrk]->charge()*pow(TrackRefInCone[iTrk]->pt(),KExp);
    den = den + pow(TrackRefInCone[iTrk]->pt(),KExp);
  }
  return (den > 0.) ? num/den : -9999999;

}

double BsToJpsiPhiAnalysis::LeptonChargeCone(const reco::PFCandidateCollection &PFCand, const reco::PFCandidate theLept, const double Dr, const double KExp, bool IncludingLepton)

{
  if( Dr < 0.)
  {
    std::cout << "E R R O R! DeltaR must be positive.\n";
    exit(1);
  }
  if( KExp < 0.)
  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }
  double num = 0.;
  double den = 0.;

  for(size_t pfc =0; pfc< PFCand.size(); pfc++){
    double DeltaR = deltaR(theLept.eta(), theLept.phi(), PFCand[pfc].eta(), PFCand[pfc].phi());
    if( PFCand[pfc].particleId() < 1 || PFCand[pfc].particleId() > 3 )  continue;
    if( PFCand[pfc].pt()        < 0.5 )                                 continue;
    if( fabs(PFCand[pfc].eta()) > 2.5 )                                 continue;
    if( PFCand[pfc].particleId() == 2){
      if( PFCand[pfc].gsfTrackRef().isNull() )                            continue;
      if( PFCand[pfc].gsfTrackRef()->numberOfValidHits() <= 5 )           continue;
    }
    else{
      if( PFCand[pfc].trackRef().isNull() )                               continue;
      if( PFCand[pfc].trackRef()->numberOfValidHits() <= 5 )              continue;
    }
    if( DeltaR > Dr )                                                   continue;

if( !IncludingLepton                                          &&
        ( theLept.particleId() - PFCand[pfc].particleId() == 0  ) &&
        ( theLept.charge()     - PFCand[pfc].charge()     == 0  ) &&
        ( fabs(theLept.pt()    - PFCand[pfc].pt())       < 1e-6 ) &&
        ( fabs(theLept.eta()   - PFCand[pfc].eta())      < 1e-6 ) &&
        ( fabs(theLept.phi()   - PFCand[pfc].phi())      < 1e-6 ) ) continue;

    num = num + PFCand[pfc].charge()*pow(PFCand[pfc].pt(),KExp);
    den = den + pow(PFCand[pfc].pt(),KExp);
  }
  return (den > 0.) ? num/den : -9999999;
}
*/
double BsToJpsiPhiAnalysis::JetChargeCone(const pat::Jet theJet, const double Dr, const double KExp)

{
  if( verbose_ && Dr < 0.)
  {
    std::cout << "W A R N I N G! DeltaR set negative. Using all the objects associated to the jet.\n";
  }

  if( KExp < 0.)

  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }

  double num = 0.;
  double den = 0.;
  for(size_t pfc =0; pfc< theJet.getPFConstituents().size(); pfc++){
    double DeltaR = deltaR(theJet.eta(), theJet.phi(), theJet.getPFConstituent(pfc)->eta(), theJet.getPFConstituent(pfc)->phi());
    if( theJet.getPFConstituent(pfc)->particleId() < 1 || theJet.getPFConstituent(pfc)->particleId() > 3 )      continue;
    if( theJet.getPFConstituent(pfc)->pt()        < 0.5 )                                                       continue;
    if( fabs(theJet.getPFConstituent(pfc)->eta()) > 2.5 )                                                       continue;
    if( theJet.getPFConstituent(pfc)->particleId() == 2){
      if( theJet.getPFConstituent(pfc)->gsfTrackRef().isNull() )                                                  continue;
      if( theJet.getPFConstituent(pfc)->gsfTrackRef()->numberOfValidHits() <= 5 )                                 continue;
    }
    else{
      if( theJet.getPFConstituent(pfc)->trackRef().isNull() )                                                     continue;
      if( theJet.getPFConstituent(pfc)->trackRef()->numberOfValidHits() <= 5 )                                    continue;
    }
    if( Dr > 0 && DeltaR > Dr )                                                                                 continue;
    num = num + theJet.getPFConstituent(pfc)->charge()*pow(theJet.getPFConstituent(pfc)->pt(),KExp);
    den = den + pow(theJet.getPFConstituent(pfc)->pt(),KExp);
  }
  return (den > 0.) ? num/den : -9999999;

}

double BsToJpsiPhiAnalysis::JetTrackChargeCone(const pat::Jet theJet, const double Dr, const double KExp)

{
  if( verbose_ && Dr < 0.)
  {
    std::cout << "W A R N I N G! DeltaR set negative. Using all the objects associated to the jet.\n";
  }
  if( KExp < 0.)
  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }
  double num = 0.;
  double den = 0.;
  for(size_t pfc =0; pfc< theJet.associatedTracks().size(); pfc++){
    double DeltaR = deltaR(theJet.eta(), theJet.phi(), theJet.associatedTracks().at(pfc)->eta(), theJet.associatedTracks().at(pfc)->phi());
    if( theJet.associatedTracks().at(pfc)->pt()        < 0.5 )          continue;
    if( fabs(theJet.associatedTracks().at(pfc)->eta()) > 2.5 )          continue;
    if( theJet.associatedTracks().at(pfc)->numberOfValidHits() <= 5 )   continue;
    if( Dr > 0 && DeltaR > Dr )                                         continue;
    num = num + theJet.associatedTracks().at(pfc)->charge()*pow(theJet.associatedTracks().at(pfc)->pt(),KExp);
    den = den + pow(theJet.associatedTracks().at(pfc)->pt(),KExp);
  }
  return (den > 0.) ? num/den : -9999999;
}


double BsToJpsiPhiAnalysis::SVTrackChargeCone(const reco::SecondaryVertexTagInfo * SVtag, const double KExp)

{
  if( KExp < 0.)
  {
    std::cout << "E R R O R! Exponent K must be positive.\n";
    exit(1);
  }
  double num = 0.;
  double den = 0.;

  for(size_t iTr = 0 ; iTr < SVtag->vertexTracks(0).size() ; iTr ++){
    num = num + SVtag->vertexTracks(0).at(iTr)->charge()*pow(SVtag->vertexTracks(0).at(iTr)->pt(),KExp);
    den = den + pow(SVtag->vertexTracks(0).at(iTr)->pt(),KExp);
  }
  return (den > 0.) ? num/den : -9999999;
}



double BsToJpsiPhiAnalysis::EvalPtRel(const reco::PFCandidate theLept, const pat::Jet theJet, bool PtIn)
{
  TLorentzVector pLept, pJet;
  pLept.SetPtEtaPhiE(theLept.pt(),theLept.eta(),theLept.phi(),theLept.energy());
  pJet.SetPtEtaPhiE(theJet.pt(),theJet.eta(),theJet.phi(),theJet.energy());
  bool LeptIntoJet = false;
  for(size_t pfc = 0; pfc < theJet.getPFConstituents().size(); pfc++)
  {

if( ( theLept.particleId() - theJet.getPFConstituent(pfc)->particleId() == 0  ) &&
        ( theLept.charge()     - theJet.getPFConstituent(pfc)->charge()     == 0  ) &&
        ( fabs(theLept.pt()    - theJet.getPFConstituent(pfc)->pt())       < 1e-6 ) &&
        ( fabs(theLept.eta()   - theJet.getPFConstituent(pfc)->eta())      < 1e-6 ) &&
        ( fabs(theLept.phi()   - theJet.getPFConstituent(pfc)->phi())      < 1e-6 ) ) continue;
    LeptIntoJet = true;
  }
if ( !PtIn ){
if( !LeptIntoJet ) {
      if(verbose_) std::cout << "W A R N I N G! Evaluating PtOut (Lepton removed from Jet), but Lepton is not in the list of Jet's constituents...\n";
    }
    pJet = pJet - pLept;
  }
  double PtRel = pLept.Pt(pJet.Vect());
  return PtRel;
}
bool BsToJpsiPhiAnalysis::LooseJetId(const pat::Jet theJet)

{

  if( !theJet.isPFJet() )                             return false;
//loose jet id
if( theJet.neutralHadronEnergyFraction()  >  0.99 )  return false;
  if( theJet.neutralEmEnergyFraction()      >  0.99 )  return false;
  if( theJet.getPFConstituents().size()     <=    1 )  return false;
  if( fabs(theJet.eta()) > 2.4)
  {
    if( theJet.chargedHadronEnergyFraction() <= 0.00 ) return false;
    if( theJet.chargedMultiplicity()         <=    0 ) return false;
    if( theJet.chargedEmEnergyFraction()     >  0.99 ) return false;
  }
  return true;
}

short int BsToJpsiPhiAnalysis::FindMuonMCSimpleCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles)

{
//edm::Handle<edm::View<reco::GenParticle>>& genParticles
short int mccode = FindMuonMCCode(theMu, genParticles);
  if (mccode < 1 || mccode > 2) mccode = 3;
  return mccode;
}


short int BsToJpsiPhiAnalysis::FindElectronMCSimpleCode(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles)

{
short int mccode = FindElectronMCCode(theEle, genParticles);
  if (mccode < 1 || mccode > 2) mccode = 3;
  return mccode;
}

short int BsToJpsiPhiAnalysis::FindMuonMCCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles)

{
int mccode   = -9999999;
double DR  = 999999.;
  double DPt = 999999.;
  size_t genPIndex=0;

for( size_t igenp = 0; igenp < genParticles->size(); ++ igenp )

  {
    const GenParticle & genP = (*genParticles)[igenp];
    double DeltaR  = deltaR(theMu.innerTrack()->eta(),theMu.innerTrack()->phi(),genP.eta(),genP.phi());
    double DeltaPt = ( theMu.innerTrack()->pt() - genP.pt() ) / genP.pt();
    if( DeltaR > 0.12 || fabs(DeltaPt) > 0.30 ) continue;
    if( theMu.charge() != genP.charge() )       continue;
    if (DeltaR < DR)
    {
      DR  = DeltaR;
      DPt = DeltaPt;
      genPIndex = igenp;
    }
  }
  if( DR < 999999. && DPt < 999999.)
  {
    GenParticle MuonGenP    = (*genParticles)[genPIndex];
    if( abs(MuonGenP.pdgId()) == 13 )
    {
       std::cout << "Matched GenP is aMu -> Look for MCtruth.\n";
mccode = LookForMotherString(MuonGenP);
    }
    else
    {
mccode = 4;
    }
  }
  else
  {
mccode = 5;
  }
return mccode;
}
short int BsToJpsiPhiAnalysis::FindElectronMCCode(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles)
{
int mccode   = -9999999;
double DR  = 999999.;
  double DPt = 999999.;
  size_t genPIndex=0;
 std::cout << "\nTagElectron : " << theEle.pdgId() << "\n";
  for( size_t igenp = 0; igenp < genParticles->size(); ++ igenp )
  {
    const GenParticle & genP = (*genParticles)[igenp];
    double DeltaR  = deltaR(theEle.gsfTrack()->eta(),theEle.gsfTrack()->phi(),genP.eta(),genP.phi());
    double DeltaPt = ( theEle.gsfTrack()->pt() - genP.pt() ) / genP.pt();

    if( DeltaR > 0.12 || fabs(DeltaPt) > 0.30 )         continue;
    if( theEle.gsfTrack()->charge() != genP.charge() )    continue;

    if (DeltaR < DR)
    {
      DR  = DeltaR;
      DPt = DeltaPt;
      genPIndex = igenp;
    }
  }
 if( DR < 999999. && DPt < 999999.)
  {
    GenParticle ElectronGenP    = (*genParticles)[genPIndex];
    if( abs(ElectronGenP.pdgId()) == 11 )
    {
mccode = LookForMotherString(ElectronGenP);
    }
    else
    {
       std::cout << "Matched GenP is not a Ele -> possibile DIF/PT.\n";
       std::cout << "Returning 4\n";
      mccode = 4;
    }
  }
  else
  {
 mccode = 5;
  }
return mccode;
}


int BsToJpsiPhiAnalysis::FindMuonAncestor(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles)

{
//edm::Handle<edm::View<reco::GenParticle>>& genParticles
  int PdgId   = -9999999;
double DR  = 999999.;
  double DPt = 999999.;
  size_t genPIndex = 0;

  for( size_t igenp = 0; igenp < genParticles->size(); ++ igenp )
  {
    const GenParticle & genP = (*genParticles)[igenp];
 double DeltaR  = deltaR(theMu.innerTrack()->eta(),theMu.innerTrack()->phi(),genP.eta(),genP.phi());
    double DeltaPt = ( theMu.innerTrack()->pt() - genP.pt() ) / genP.pt();

    if( DeltaR > 0.12 || fabs(DeltaPt) > 0.30 ) continue;
    if( theMu.charge() != genP.charge() )       continue;
    if (DeltaR < DR)

    {
      DR  = DeltaR;
      DPt = DeltaPt;
      genPIndex = igenp;
    }
  }

  if( DR < 999999. && DPt < 999999.)
  {
    GenParticle MuonGenP    = (*genParticles)[genPIndex];
    if( abs(MuonGenP.pdgId()) == 13 )
    {
      if(verbose_) std::cout << "Matched GenP is a Mu -> Look for MCtruth.\n";
      PdgId = LookForMotherStringId(MuonGenP);
    }
    else
    {
      if(verbose_){
        std::cout << "Matched GenP is not a Mu -> possibile DIF/PT.\n";
        std::cout << "Returning 4\n";
      }
      PdgId = 0;
    }
  }
  else
  {
    if(verbose_) {
      std::cout << "There is no GenP matched to the Mu -> possibile FAKE/DIF/PT.\n";
      std::cout << "Returning 5\n";
    }
    PdgId = -9999999;
  }

return PdgId;

}

int BsToJpsiPhiAnalysis::FindElectronAncestor(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles)
{
//edm::Handle<edm::View<reco::GenParticle>>& genParticles
  int PdgId   = -9999999;
double DR  = 999999.;
  double DPt = 999999.;
  size_t genPIndex = 0;

  for( size_t igenp = 0; igenp < genParticles->size(); ++ igenp )

  {

    const GenParticle & genP = (*genParticles)[igenp];
    double DeltaR  = deltaR(theEle.gsfTrack()->eta(),theEle.gsfTrack()->phi(),genP.eta(),genP.phi());
    double DeltaPt = ( theEle.gsfTrack()->pt() - genP.pt() ) / genP.pt();

    if( DeltaR > 0.12 || fabs(DeltaPt) > 0.30 )         continue;
    if( theEle.gsfTrack()->charge() != genP.charge() )    continue;
   if (DeltaR < DR)
    {
      DR  = DeltaR;
      DPt = DeltaPt;
      genPIndex = igenp;
    }
  }
  if( DR < 999999. && DPt < 999999.)
  {
    GenParticle ElectronGenP    = (*genParticles)[genPIndex];
    if( abs(ElectronGenP.pdgId()) == 11 )
    {
      if(verbose_) std::cout << "Matched GenP is a Ele -> Look for MCtruth.\n";
      PdgId = LookForMotherStringId(ElectronGenP);
    }
    else
    {
      if(verbose_) {
        std::cout << "Matched GenP is not a Ele -> possibile DIF/PT.\n";
        std::cout << "Returning 4\n";
      }
      PdgId = 0;
    }
  }
  else
  {
    if(verbose_) {
      std::cout << "There is no GenP matched to the Ele -> possibile FAKE/DIF/PT.\n";
      std::cout << "Returning 5\n";
    }
    PdgId = -9999999;
  }
return PdgId;
}
short int BsToJpsiPhiAnalysis::LookForMotherString(GenParticle theGenP)

{
  if (theGenP.numberOfMothers() != 1)
  {
 return 0;
  }
  const Candidate * iGenP    = theGenP.mother(0);
  const Candidate * iGenPDau = theGenP.mother(0);
  int   id = iGenP->pdgId();
  if(id > 80 && id < 101)
  {
return(-1);

  }
  while(id < 81 || id > 100)
  {
    unsigned short nMoms = iGenP->numberOfMothers();
    if(nMoms > 1)
    { /// more than 1 mother
      break;
    }
    else if(nMoms == 0)
    { /// no mothers
      if(id == 2212)
      { /// the genp is the colliding proton itself
        break;
      }
      else
      { /// genp has no mother but is not the genp
        std::cout << "E R R O R! Genp has no mothers and particle does not come from the proton!\n";
        std::cout << "           Exiting...\n";
        exit(1);
      }
    }
    else
    {
      iGenPDau = iGenP;
      iGenP    = iGenP->mother();
      id       = iGenP->pdgId();
}
  }

  int idGenMuon = theGenP.pdgId();
  short int muonSign = idGenMuon>0?1:-1; // Sign of muon : - for mu+(-13) / + for mu-(+13)
  int firstId = iGenPDau->pdgId();
  short int firstSign = firstId>0?1:-1; // Sign/CP of first hadron after the string : + for B+(+521) / - for B-(-521) / + for B0(+511) / - for antiB0 (-511) / ...
  int firstType = firstId%10000;
  int firstSubType;
  if(abs(firstType) < 1000)
    firstSubType = firstType/100;
  else
  {
    firstSubType = firstType/1000;
    firstSign    = -firstSign;
  }
  if(abs(firstId) == 130)
    firstSubType = firstId>0?3:-3;
  if(abs(firstType) < 1000 && (abs(firstSubType) == 3 || abs(firstSubType) == 5))
    firstSubType = -firstSubType;
  if(abs(firstSubType) == 5) // first particle after string is a B
  {
    if( firstSign == (- muonSign) )
    {
 return(1);
    }
    else if( firstSign == muonSign )
    {
 return(2);
    }
    else
    {
      std::cout << "E R R O R ! Muon and b sign can only be Same or Opposite\n";
      std::cout << "            Muon sign is " << muonSign << ", b sign is " << firstSign << "\n";
      exit(1);
    }
  }
  else
  {
return(3);
  }
  std::cout << "E R R O R! We should never be here!\n";
  exit(1);
}

int BsToJpsiPhiAnalysis::LookForMotherStringId(GenParticle theGenP)
{
  if (theGenP.numberOfMothers() != 1)
  {
 return 0;
  }
  const Candidate * iGenP    = theGenP.mother(0);
  const Candidate * iGenPDau = theGenP.mother(0);
  int   id = iGenP->pdgId();
  if(id > 80 && id < 101)
  {
return(-1);

  }

  while(id < 81 || id > 100)

  {
    unsigned short nMoms = iGenP->numberOfMothers();
    if(nMoms > 1)
    { /// more than 1 mother
      break;
    }
    else if(nMoms == 0)
    { /// no mothers
      if(id == 2212)
      { /// the genp is the colliding proton itself
        break;
      }
      else
      { /// genp has no mother but is not the genp
        std::cout << "E R R O R! Genp has no mothers and particle does not come from the proton!\n";
        std::cout << "           Exiting...\n";
        exit(1);
      }
    }
    else
    {
      iGenPDau = iGenP;
      iGenP    = iGenP->mother();
      id       = iGenP->pdgId();
std::cout << " " << iGenP->pdgId() << " -> " << iGenPDau->pdgId() << "\n";
    }
  }
if(abs(iGenPDau->pdgId())<1000)
  {
    if( abs(iGenPDau->pdgId())%10<2 )
    {
 return iGenPDau->pdgId();
    }
    else
    { for(size_t idau = 0; idau < iGenPDau->numberOfDaughters(); idau++ )
      {
        const Candidate * tempGenP    = iGenPDau->daughter(idau);
if( abs(tempGenP->pdgId())/100 == abs(iGenPDau->pdgId())/100 )
        {
 return tempGenP->pdgId();
        }
}
    }
  }
  else if(abs(iGenPDau->pdgId())>=1000)
  {
    if( abs(iGenPDau->pdgId())%10==2 || abs(iGenPDau->pdgId())==2224 || abs(iGenPDau->pdgId())==2214 || abs(iGenPDau->pdgId())==2114 || abs(iGenPDau->pdgId())==1114)
    {
return iGenPDau->pdgId();
    }
    else
    {
 for(size_t idau = 0; idau < iGenPDau->numberOfDaughters(); idau++ )
      {
        const Candidate * tempGenP    = iGenPDau->daughter(idau);
if( abs(tempGenP->pdgId())/1000 == abs(iGenPDau->pdgId())/1000 )
        {
return tempGenP->pdgId();
        }
 }
    }
  }
std::cout << "W A R N I N G! It should never get here.\n";
  std::cout << "               Returning +9999999\n";
  return +9999999;
}









