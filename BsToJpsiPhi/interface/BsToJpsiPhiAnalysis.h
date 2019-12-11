#ifndef HeavyFlavorAnalysis_BsToJpsiPhi_BsToJpsiPhiAnalysis_h
#define HeavyFlavorAnalysis_BsToJpsiPhi_BsToJpsiPhiAnalysis_h

// -*- C++ -*-
//
// Package:    BsToJpsiPhi
// Class:      BsToJpsiPhiAnalysis
//


// system include files
#include <memory>
#include <string>

// user include files
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
// #include "RecoVertex/KalmanVertexFit/test/SimpleVertexTree.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/BsToJpsiPhiRootTree.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// #include "RecoVertex/KinematicFitPrimitives/interface/"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include <TFile.h>
#include <TH1F.h>

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//class BPHRecoCandidate;

class BsToJpsiPhiAnalysis : public edm::EDAnalyzer {
public:
  explicit BsToJpsiPhiAnalysis(const edm::ParameterSet&);
  ~BsToJpsiPhiAnalysis();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();


  void fillMCInfo( edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  void setFitParKK(RefCountedKinematicTree& myTree);
  void setFitParHyp1(RefCountedKinematicTree& myTree);
  void setFitParHyp2(RefCountedKinematicTree& myTree);

private:
  bool MCmatching(const reco::Candidate & track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles,
                  int &K1mcId, int &K1momId, int &K1gmomId,
                  int condMom, int condGMom);

bool MCmatchingBplusK(const reco::Candidate & track1, edm::Handle<edm::View<reco::GenParticle>>& genParticles,
                        int &K1mcId, int &K1momId,
                        int condMom);

 bool MCmatchingJpsi(const reco::Candidate& track1,  edm::Handle<edm::View<reco::GenParticle>>& genParticles, int &K1mcId, int &K1momId,	int condMom);

reco::Vertex reVertex(const edm::Event& theEvent,  const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, reco::Track & rtrk1, reco::Track & rtrk2);
  //reco::Vertex reVertex(const edm::EventSetup& ,edm::Handle<reco::BeamSpot> , reco::Vertex , pat::Muon , pat::Muon , reco::TrackRef , reco::TrackRef );
//  reco::Vertex reVertex2(const edm::EventSetup& ,edm::Handle<reco::BeamSpot>, reco::Vertex , pat::Muon , pat::Muon , reco::TrackRef, double &vtxCL);
reco::Vertex reVertex2(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, reco::TrackRef trk3, double &vtxCL);


double CalculateCtErrvertex(const edm::EventSetup& iSetup, reco::Vertex PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass);
  double CalculateCtErrbeamspot(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass); // for Beam spot 


 // double CalculateCtErr(const edm::EventSetup& iSetup, edm::Handle<reco::BeamSpot>& PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass);
 // double CalculateCtErr(const edm::EventSetup& iSetup, reco::Vertex PV, TransientVertex SV, TMatrix SVcovariance2D ,GlobalVector Ptvec, double Bmass); // for Beam spot 
  void RecursivelyPrintMother( const reco::Candidate & genp );

  BsToJpsiPhiRootTree * bsRootTree_;

  GlobalVector flightDirection(const reco::Vertex &pv, reco::Vertex &sv);


  edm::ParameterSet theConfig_;

  bool selGlobalMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx);
  bool selTrackerMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx);
  bool selTightMuon(const pat::Muon aMuon, const reco::Vertex RecVtx);

  /// MuonChargeCone - TrackBased
  double MuonChargeCone(const edm::Event& theEvent, reco::TrackRef muTrackRef, const double Dr, const double KExp, bool IncludingMuon);
  double MuonChargeConeWrtPV(const edm::Event& theEvent, reco::TrackRef muTrackRef, reco::Vertex PVtx, const double Dr, const double KExp, bool IncludingMuon);

  /// ElectronChargeCone - TrackBased
  double ElectronChargeCone(const edm::Event& theEvent, reco::GsfTrackRef eleGsfTrackRef, const double Dr, const double KExp, bool IncludingElectron);
  double ElectronChargeConeWrtPV(const edm::Event& theEvent, reco::GsfTrackRef eleGsfTrackRef, reco::Vertex PVtx, const double Dr, const double KExp, bool IncludingElectron);

  /// LeptonChargeCone - PFcandidateBased
  double LeptonChargeCone(const reco::PFCandidateCollection & PFCand, const reco::PFCandidate theLept, const double Dr, const double KExp, bool IncludingLepton);

  /// Jet and SV ChargeCone
  double JetChargeCone(const pat::Jet theJet, const double Dr, const double KExp);
  double JetTrackChargeCone(const pat::Jet theJet, const double Dr, const double KExp);
  double SVTrackChargeCone(const reco::SecondaryVertexTagInfo * SVtag, const double KExp);

  /// PtRel
  double EvalPtRel(const reco::PFCandidate theLept, const pat::Jet theJet, bool PtIn);

  /// JetId
  bool LooseJetId(const pat::Jet theJet);

  /// MuonMC
  short int FindMuonMCCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  short int FindMuonMCSimpleCode(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  int FindMuonAncestor(const pat::Muon theMu, edm::Handle<edm::View<reco::GenParticle>>& genParticles);

  /// ElectronMC 
  short int FindElectronMCCode(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  short int FindElectronMCSimpleCode(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles);
  int FindElectronAncestor(const pat::Electron theEle, edm::Handle<edm::View<reco::GenParticle>>& genParticles);

 /* void dumpRecoCand( const std::string& name,
                                     const BPHRecoCandidate* cand );
*/

  /// GenericMC
  int LookForMotherStringId(reco::GenParticle theGenP);
  short int LookForMotherString(reco::GenParticle theGenP);

  const TrackerGeometry* m_tracker;

  bool isMCstudy_;
  //edm::InputTag TriggerTag;
  //edm::EDGetTokenT<edm::TriggerResults> TriggerTagTok;// only when vector edm::View<edm::variable>
  edm::InputTag genParticlesLabel;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticlesTok;  
  edm::InputTag MuonTag;
  edm::EDGetTokenT<edm::View<pat::Muon>> MuonTagTok;
  edm::InputTag JetCollection;
  edm::EDGetTokenT<edm::View<pat::Jet>> JetCollectionTok;
  edm::InputTag PUInfo;
  edm::EDGetTokenT<edm::View<PileupSummaryInfo>> PUInfoTok;
  edm::InputTag vertexBeamSpot;
  edm::EDGetTokenT<reco::BeamSpot> vertexBeamSpotTok;
  edm::InputTag primaryvertex;
  edm::EDGetTokenT<edm::View<reco::Vertex>> primaryvertexTok;
  edm::InputTag triggerresults;
  edm::EDGetTokenT<edm::TriggerResults> triggerresultsTok;
  edm::InputTag ElectronTag;
  edm::EDGetTokenT<edm::View<pat::Electron>> ElectronTagTok;
  edm::InputTag track;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackLabelK;
  edm::InputTag isotrack;
  edm::EDGetTokenT<edm::View<pat::IsolatedTrack>> isotrackTok;
  bool StoreDeDxInfo_;
  bool verbose_;
  bool TestVerbose_;

  const double nominalJpsiMass;
  const double nominalPhiMass;
  const double nominalElectronMass;
  const double nominalMuonMass;
  const double nominalKaonMass;
  const double nominalPionMass;
  const double nominalKstarMass;
  const double nominalBplusMass;

  double JpsiMassWindowBeforeFit_;
  double JpsiMassWindowAfterFit_;
  double JpsiPtCut_;
  double KaonTrackPtCut_;
  double BdKaonTrackPtCut_;
  double PhiMassWindowAfterFit_;
  double PhiMassWindowBeforeFit_;
  double BsLowerMassCutBeforeFit_;
  double BsUpperMassCutBeforeFit_;
  double BsLowerMassCutAfterFit_ ;
  double BsUpperMassCutAfterFit_ ;
  double KstarMassWindowBeforeFit_;
  double KstarMassWindowAfterFit_;
  double BdLowerMassCutBeforeFit_;
  double BdUpperMassCutBeforeFit_;
  double BdLowerMassCutAfterFit_;
  double BdUpperMassCutAfterFit_;

  double BsPDGMass_;
  double BdPDGMass_;
  double BpPDGMass_;

  std::string outputFile_; // output file

  int Mu1Truth;


  int match[15][10];
  int match2[15][10];
  int matching[15][10];
  int L1_mu_size0;
  int L1_mu_size1;
  int L1_mu_size;
  int L1_mu_size2;

  unsigned int tagmucounter_;
  unsigned int event_counter_;
  unsigned int elecounter_;
  unsigned int muoncounter_;
  unsigned int jetcounter_;

  double angle_costheta;
  double angle_phi;
  double angle_cospsi;
  double AngleBsDecayLength;

};
#endif
