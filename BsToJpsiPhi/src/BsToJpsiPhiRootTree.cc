#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/BsToJpsiPhiRootTree.h"

using namespace std;


BsToJpsiPhiRootTree::BsToJpsiPhiRootTree()
{
  resetEntries();
  bsTree_ = 0;
  bsFile_ = 0;
}


void BsToJpsiPhiRootTree::createTree(const std::string filename)
{
  // open root file
  bsFile_ = new TFile (filename.c_str(), "RECREATE" );
  int bufsize = 256000;
  // create tree structure
  bsTree_ = new TTree("BsTree","BsTree",bufsize);

  PVTrkCharge_  = new std::vector<int>();
  PVTrkPt_      = new std::vector<float>();
  PVTrkEta_     = new std::vector<float>();
  PVTrkPhi_     = new std::vector<float>();

  BpPVTrkCharge_= new std::vector<int>();
  BpPVTrkPt_    = new std::vector<float>();
  BpPVTrkEta_   = new std::vector<float>();
  BpPVTrkPhi_   = new std::vector<float>();

  BJetTrkCharge_= new std::vector<int>();
  BJetTrkPt_    = new std::vector<float>();

  BpBJetTrkCharge_ = new std::vector<int>();
  BpBJetTrkPt_ = new std::vector<float>();

  /// Multiplicities - BEGIN
  bsTree_->Branch("MuonMultiplicity", &MuonMultiplicity_, "MuonMultiplicity/I");
  bsTree_->Branch("ElectronMultiplicity", &ElectronMultiplicity_, "ElectronMultiplicity/I");
  bsTree_->Branch("TrackMultiplicity", &TrackMultiplicity_, "TrackMultiplicity/I");
  bsTree_->Branch("TrackMultiplicityBp", &TrackMultiplicityBp_, "TrackMultiplicityBp/I");
  bsTree_->Branch("TrackMultiplicityBd", &TrackMultiplicityBd_, "TrackMultiplicityBd/I");
  /// Multiplicities - END

  /// Jet Variables - BEGIN
  bsTree_->Branch("BJetTrkCharge","vector<int>",BJetTrkCharge_);
  bsTree_->Branch("BJetTrkPt","vector<float>",BJetTrkPt_);

  bsTree_->Branch("BpBJetTrkCharge","vector<int>",BpBJetTrkCharge_);
  bsTree_->Branch("BpBJetTrkPt","vector<float>",BpBJetTrkPt_);

  bsTree_->Branch("PVTrkCharge","vector<int>",PVTrkCharge_);
  bsTree_->Branch("PVTrkPt","vector<float>",PVTrkPt_);
  bsTree_->Branch("PVTrkEta","vector<float>",PVTrkEta_);
  bsTree_->Branch("PVTrkPhi","vector<float>",PVTrkPhi_);

  bsTree_->Branch("BpPVTrkCharge","vector<int>",BpPVTrkCharge_);
  bsTree_->Branch("BpPVTrkPt","vector<float>",BpPVTrkPt_);
  bsTree_->Branch("BpPVTrkEta","vector<float>",BpPVTrkEta_);
  bsTree_->Branch("BpPVTrkPhi","vector<float>",BpPVTrkPhi_);

  bsTree_->Branch("BJetParton", &BJetParton_, "BJetParton/I");
  bsTree_->Branch("BJetEta", &BJetEta_, "BJetEta/D");
  bsTree_->Branch("BJetPhi", &BJetPhi_, "BJetPhi/D");
  bsTree_->Branch("BJetPt", &BJetPt_, "BJetPt/D");
  bsTree_->Branch("JetBTagProb", &JetBTagProb_, "JetBTagProb/D");

  bsTree_->Branch("BpBJetParton", &BpBJetParton_, "BpBJetParton/I");
  bsTree_->Branch("BpBJetEta", &BpBJetEta_, "BpBJetEta/D");
  bsTree_->Branch("BpBJetPhi", &BpBJetPhi_, "BpBJetPhi/D");
  bsTree_->Branch("BpBJetPt", &BpBJetPt_, "BpBJetPt/D");
  bsTree_->Branch("BpJetBTagProb", &BpJetBTagProb_, "BpJetBTagProb/D");
  /// Jet Variables - END

  bsTree_->Branch("BplusCharge",&BplusCharge_, "BplusCharge/I");
  bsTree_->Branch("BplusMu1Eta",&BplusMu1Eta_, "BplusMu1Eta/D");
  bsTree_->Branch("BplusMu2Eta",&BplusMu2Eta_, "BplusMu2Eta/D");
  bsTree_->Branch("BpSoftMuon1",&BpSoftMuon1_,"BpSoftMuon1/I");
  bsTree_->Branch("BpSoftMuon2",&BpSoftMuon2_,"BpSoftMuon2/I"); 
 
 
  bsTree_->Branch("BdMuonCat1",&BdMuonCat1_, "BdMuonCat1/I");
  bsTree_->Branch("BdMuonCat2",&BdMuonCat2_, "BdMuonCat2/I");
  bsTree_->Branch("BplusM_fit",&BplusM_fit_, "BplusM_fit/D");
  bsTree_->Branch("BplusVtxProb",&BplusVtxProb_, "BplusVtxProb/D");
  bsTree_->Branch("BplusChi2",&BplusChi2_, "BplusChi2/D");
  bsTree_->Branch("BplusPt",&BplusPt_, "BplusPt/D");
  bsTree_->Branch("BplusPtot",&BplusPtot_, "BplusPtot/D");
  bsTree_->Branch("KplusPt",&KplusPt_, "KplusPt/D");
  bsTree_->Branch("KplusHighPurityTrack",&KplusHighPurityTrack_, "KplusHighPurityTrack/D");

  bsTree_->Branch("KplusPtot",&KplusPtot_, "KplusPtot/D");
  bsTree_->Branch("BplusMu1Pt",&BplusMu1Pt_, "BplusMu1Pt/D");
  bsTree_->Branch("BplusMu2Pt",&BplusMu2Pt_, "BplusMu2Pt/D");
  bsTree_->Branch("BpJpsiVtxProb",&BpJpsiVtxProb_, "BpJpsiVtxProb/D");
  bsTree_->Branch("BpCosDeltaAlpha",&BpCosDeltaAlpha_, "BpCosDeltaAlpha/D");
  bsTree_->Branch("BpMuMuDCA",&BpMuMuDCA_, "BpMuMuDCA/D");
  bsTree_->Branch("BpMuMuDistance",&BpMuMuDistance_, "BpMuMuDistance/D");
  bsTree_->Branch("BpMuMuDistanceSigma",&BpMuMuDistanceSigma_, "BpMuMuDistanceSigma/D");
  bsTree_->Branch("BpMuDr1",&BpMuDr1_, "BpMuDr1/D");
  bsTree_->Branch("BpMuDr2",&BpMuDr2_, "BpMuDr2/D");
  bsTree_->Branch("BpMuonCat1",&BpMuonCat1_, "BpMuonCat1/I");
  bsTree_->Branch("BpMuonCat2",&BpMuonCat2_, "BpMuonCat2/I");
  bsTree_->Branch("BplusMu1Ptot",&BplusMu1Ptot_, "BplusMu1Ptot/D");
  bsTree_->Branch("BplusMu2Ptot",&BplusMu2Ptot_, "BplusMu2Ptot/D");
  bsTree_->Branch("BplusEta",&BplusEta_, "BplusEta/D");
  bsTree_->Branch("BplusPhi",&BplusPhi_, "BplusPhi/D");
  bsTree_->Branch("JpsiMass_bplus",&JpsiMass_bplus_, "JpsiMass_bplus/D");
  bsTree_->Branch("JpsiPt_bplus",&JpsiPt_bplus_, "JpsiPt_bplus/D");
  bsTree_->Branch("BplusPVindex",&BplusPVindex_, "BplusPVindex/I");

  bsTree_->Branch("BpCt2DPVCosTheta",&BpCt2DPVCosTheta_, "BpCt2DPVCosTheta/D");
  bsTree_->Branch("BpCt3DPVCosTheta",&BpCt3DPVCosTheta_, "BpCt3DPVCosTheta/D");
  bsTree_->Branch("BpCt2DPVClosestZ",&BpCt2DPVClosestZ_, "BpCt2DPVClosestZ/D");
  bsTree_->Branch("BpCt3DPVClosestZ",&BpCt3DPVClosestZ_, "BpCt3DPVClosestZ/D");
  bsTree_->Branch("BpCt2D",&BpCt2D_, "BpCt2D/D");
  bsTree_->Branch("BpCt3D",&BpCt3D_, "BpCt3D/D");
  bsTree_->Branch("BpCt2DBS",&BpCt2DBS_, "BpCt2DBS/D");

  bsTree_->Branch("BpPx",&BpPx_,"BpPx/D");
  bsTree_->Branch("BpPy",&BpPy_,"BpPy/D");
  bsTree_->Branch("BpPz",&BpPz_,"BpPz/D");

  bsTree_->Branch("BpSVx",&BpSVx_, "BpSVx/D");
  bsTree_->Branch("BpSVy",&BpSVy_, "BpSVy/D");
  bsTree_->Branch("BpSVz",&BpSVz_, "BpSVz/D");
  
  bsTree_->Branch("BpPVx_refit",&BpPVx_refit_, "BpPVx_refit/D");
  bsTree_->Branch("BpPVy_refit",&BpPVy_refit_, "BpPVy_refit/D");
  bsTree_->Branch("BpPVz_refit",&BpPVz_refit_, "BpPVz_refit/D");

  bsTree_->Branch("BpPVx_refit_closestZ",&BpPVx_refit_closestZ_, "BpPVx_refit_closestZ/D");
  bsTree_->Branch("BpPVy_refit_closestZ",&BpPVy_refit_closestZ_, "BpPVy_refit_closestZ/D");
  bsTree_->Branch("BpPVz_refit_closestZ",&BpPVz_refit_closestZ_, "BpPVz_refit_closestZ/D");

  bsTree_->Branch("BpPVx_refit_cosTheta",&BpPVx_refit_cosTheta_, "BpPVx_refit_cosTheta/D");
  bsTree_->Branch("BpPVy_refit_cosTheta",&BpPVy_refit_cosTheta_, "BpPVy_refit_cosTheta/D");
  bsTree_->Branch("BpPVz_refit_cosTheta",&BpPVz_refit_cosTheta_, "BpPVz_refit_cosTheta/D");

  bsTree_->Branch("BpPVCL_ClosestZ",&BpPVCL_ClosestZ_, "BpPVCL_ClosestZ/D");
  bsTree_->Branch("BpPVCL_CosTheta",&BpPVCL_CosTheta_, "BpPVCL_CosTheta/D");
  bsTree_->Branch("BpPVCL",&BpPVCL_, "BpPVCL/D");

  bsTree_->Branch("BpJpsiVtxCL",&BpJpsiVtxCL_, "BpJpsiVtxCL/D");
  bsTree_->Branch("BpVtxCL",&BpVtxCL_,"BpVtxCL/D");

  bsTree_->Branch("BpdRKaonJpsi",&BpdRKaonJpsi_,"BpdRKaonJpsi/D");
  bsTree_->Branch("BpCtErr2DBS",&BpCtErr2DBS_, "BpCtErr2DBS/D");
  bsTree_->Branch("BpCtErr2DCostheta",&BpCtErr2DCostheta_, "BpCtErr2DCostheta/D");
  bsTree_->Branch("BpCtErr2D",&BpCtErr2D_, "BpCtErr2D/D");
  bsTree_->Branch("BpCtErr2DClosestZ",&BpCtErr2DClosestZ_, "BpCtErr2DClosestZ/D");
  bsTree_->Branch("Mu1SoftID", &Mu1SoftID_, "Mu1SoftID/I");
  bsTree_->Branch("Mu2SoftID", &Mu2SoftID_, "Mu2SoftID/I");
  bsTree_->Branch(  "triggerbit_HLTmu4TkDis"             , &triggerbit_HLTmu4TkDis_,                "triggerbit_HLTmu4TkDis_/I");
  bsTree_->Branch(  "triggerbit_HLTm4TkTk"		  , &triggerbit_HLTmu4TkTkDis_,                "triggerbit_HLTmu4TkTkDis_/I");
  bsTree_->Branch("triggerbit_JpsiMuon", &triggerbit_HLTDimuon0JpsiMuon_, "triggerbit_HLTDimuon0JpsiMuon_/I");
  bsTree_->Branch("IP3DKandJpsiVtx",&IP3DKandJpsiVtx_, "IP3DKandJpsiVtx/D");
  bsTree_->Branch("IP3DKandJpsiVtxErr",&IP3DKandJpsiVtxErr_, "IP3DKandJpsiVtxErr/D");
  bsTree_->Branch("BpmatchDoubleMu01",&BpmatchDoubleMu01_, "BpmatchDoubleMu01/I");
  bsTree_->Branch("BpmatchDoubleMu02",&BpmatchDoubleMu02_, "BpmatchDoubleMu02/I");
  bsTree_->Branch("BpmatchDoubleMu41",&BpmatchDoubleMu41_, "BpmatchDoubleMu41/I");
  bsTree_->Branch("BpmatchDoubleMu42",&BpmatchDoubleMu42_, "BpmatchDoubleMu42/I");
  bsTree_->Branch("BpmatchDoubleMu41_v12",&BpmatchDoubleMu41_v12_, "BpmatchDoubleMu41_v12/I");
  bsTree_->Branch("BpmatchDoubleMu42_v12",&BpmatchDoubleMu42_v12_, "BpmatchDoubleMu42_v12/I");

  bsTree_->Branch("BpmatchDoubleMu01DiMuon0",&BpmatchDoubleMu01DiMuon0_, "BpmatchDoubleMu01DiMuon0/I");
  bsTree_->Branch("BpmatchDoubleMu02DiMuon0",&BpmatchDoubleMu02DiMuon0_, "BpmatchDoubleMu02DiMuon0/I");
  bsTree_->Branch("BplusKmcId",&BplusKmcId_, "BplusKmcId/I");
  bsTree_->Branch("BplusKmomId",&BplusKmomId_, "BplusKmomId/I");
  bsTree_->Branch("BplusMu1mcId",&BplusMu1mcId_, "BplusMu1mcId/I");
  bsTree_->Branch("BplusMu1momId",&BplusMu1momId_, "BplusMu1momId/I");
  bsTree_->Branch("BplusMu1gmomId",&BplusMu1gmomId_, "BplusMu1gmomId/I");
  bsTree_->Branch("BplusMu2mcId",&BplusMu2mcId_, "BplusMu2mcId/I");
  bsTree_->Branch("BplusMu2momId",&BplusMu2momId_, "BplusMu2momId/I");
  bsTree_->Branch("BplusMu2gmomId",&BplusMu2gmomId_, "BplusMu2gmomId/I");
  bsTree_->Branch("isMatchedBplus",&isMatchedBplus_, "isMatchedBplus/I");

  bsTree_->Branch("BplusDecayChannel",&BplusDecayChannel_, "BplusDecayChannel/I");

  bsTree_->Branch( "BpMuonPairDR", &BpMuonPairDR_, "BpMuonPairDR/D" );

  bsTree_->Branch(  "BcP"             , &BcP_,                "BcP/D");
  bsTree_->Branch(  "BcCosAlpha"             , &BcCosAlpha_,                "BcCosAlpha/D");
  bsTree_->Branch(  "BcIP3D"             , &BcIP3D_,                "BcIP3D/D");
  bsTree_->Branch(  "BcCt"             , &BcCt_,                "BcCt/D");
  bsTree_->Branch(  "BcM"             , &BcM_,                "BcM/D");
  bsTree_->Branch(  "BcMC"             , &BcMC_,                "BcMC/I");
  bsTree_->Branch(  "BcAng"             , &BcAng_,                "BcAng/D");
  bsTree_->Branch(  "BcAngPi"             , &BcAngPi_,                "BcAngPi/D");
  bsTree_->Branch(  "BcProb"             , &BcProb_,                "BcProb/D");
  bsTree_->Branch(  "BctrackPt"             , &BctrackPt_,                "BctrackPt/D");

  bsTree_->Branch("BsLxy3DMC", &BsLxy3DMC_, "BsLxy3DMC/D");
  bsTree_->Branch("BsLxy2DMC", &BsLxy2DMC_, "BsLxy2DMC/D");
  bsTree_->Branch("BsPMC", &BsPMC_, "BsPMC/D");
  bsTree_->Branch(  "PVZpos"             , PVZpos_ ,"PVZpos_[30]/D");
  bsTree_->Branch(  "NTracksInPV"             , NTracksInPV_ ,"NTracksInPV_[30]/I" );
  bsTree_->Branch(  "PVAbsPt"             , PVAbsPt_ ,"PVAbsPt_[30]/D");
  bsTree_->Branch("SVZpos", SVZpos_, "SVZpos_[30]/D");
  bsTree_->Branch("BsPtMC", &BsPtMC_, "BsPtMC/D");


  bsTree_->Branch("BsJpsiMu2PtMC", &BsJpsiMu2PtMC_, "BsJpsiMu2PtMC/D");
  bsTree_->Branch("BsJpsiMu1PtMC", &BsJpsiMu1PtMC_, "BsJpsiMu1PtMC/D");
  bsTree_->Branch("BsJpsiMu1EtaMC", &BsJpsiMu1EtaMC_, "BsJpsiMu1EtaMC/D");
  bsTree_->Branch("BsJpsiMu2EtaMC", &BsJpsiMu2EtaMC_, "BsJpsiMu2EtaMC/D");
  bsTree_->Branch("BsJpsiMu1PhiMC", &BsJpsiMu1PhiMC_, "BsJpsiMu1PhiMC/D");
  bsTree_->Branch("BsJpsiMu2PhiMC", &BsJpsiMu2PhiMC_, "BsJpsiMu2PhiMC/D");

  bsTree_->Branch("BsJpsiPtMC", &BsJpsiPtMC_, "BsJpsiPtMC/D");
  bsTree_->Branch("BsJpsiMassMC", &BsJpsiMassMC_, "BsJpsiMassMC/D");
  bsTree_->Branch("BsPhiPtMC", &BsPhiPtMC_, "BsPhiPtMC/D");
  bsTree_->Branch("BsPhiMassMC", &BsPhiMassMC_, "BsPhiMassMC/D");

  bsTree_->Branch("BsPhiK1PtMC", &BsPhiK1PtMC_, "BsPhiK1PtMC/D");
  bsTree_->Branch("BsPhiK2PtMC", &BsPhiK2PtMC_, "BsPhiK2PtMC/D");

  bsTree_->Branch("BsCosAlphaMC", &BsCosAlphaMC_, "BsCosAlphaMC/D");
  bsTree_->Branch("BsMassMC", &BsMassMC_, "BsMassMC/D");
  bsTree_->Branch("BsCt2DMC_TrueMass", &BsCt2DMC_TrueMass_, "BsCt2DMC_TrueMass/D");


  bsTree_->Branch(  "ihaveajpsi"             , &ihaveajpsi_,                "ihaveajpsi/I");
  bsTree_->Branch(  "BsCowboy"             , &BsCowboy_,                "BsCowboy/I");
  bsTree_->Branch(  "BdCowboy"             , &BdCowboy_,                "BdCowboy/I");
  bsTree_->Branch(  "BsPhiVtxProb"             , &BsPhiVtxProb_,                "BsPhiVtxProb/D");
  bsTree_->Branch(  "BsMu1QualityG"             , &BsMu1QualityG_,                "BsMu1QualityG/I");
  bsTree_->Branch(  "BsMu2QualityG"             , &BsMu2QualityG_,                "BsMu2QualityG/I");
  bsTree_->Branch(  "BsMu1QualityT"             , &BsMu1QualityT_,                "BsMu1QualityT/I");
  bsTree_->Branch(  "BsMu2QualityT"             , &BsMu2QualityT_,                "BsMu2QualityT/I");
  bsTree_->Branch(  "BdMu1QualityG"             , &BdMu1QualityG_,                "BdMu1QualityG/I");
  bsTree_->Branch(  "BdMu2QualityG"             , &BdMu2QualityG_,                "BdMu2QualityG/I");
  bsTree_->Branch(  "BdMu1QualityT"             , &BdMu1QualityT_,                "BdMu1QualityT/I");
  bsTree_->Branch(  "BdMu2QualityT"             , &BdMu2QualityT_,                "BdMu2QualityT/I");
  bsTree_->Branch(  "NVertices"             , &NVertices_,                "NVertices/I");


 bsTree_->Branch(  "JpsiMuonMatch41"         , &JpsiMuonMatch41_, "JpsiMuonMatch41/I");
 bsTree_->Branch(  "JpsiMuonMatch42"         , &JpsiMuonMatch42_, "JpsiMuonMatch42/I");

 bsTree_->Branch(  "JpsiCosDeltaAlpha"         , &JpsiCosDeltaAlpha_, "JpsiCosDeltaAlpha/D");
 bsTree_->Branch(  "JpsiLxySigma"         , &JpsiLxySigma_, "JpsiLxySigma/D");
 bsTree_->Branch(  "JpsiLxy"         , &JpsiLxy_, "JpsiLxy/D");
 bsTree_->Branch(  "JpsiGenLxy"         , &JpsiGenLxy_, "JpsiGenLxy/D");
 bsTree_->Branch(  "JpsiGenPt"         , &JpsiGenPt_, "JpsiGenPt/D");

 bsTree_->Branch(  "JpsiGenLxyOld"         , &JpsiGenLxyOld_, "JpsiGenLxyOld/D");
 bsTree_->Branch(  "JpsiGenLxyOverPt"         , &JpsiGenLxyOverPt_, "JpsiGenLxyOverPt/D");
 bsTree_->Branch(  "JpsiLxyOverPt"         , &JpsiLxyOverPt_, "JpsiLxyOverPt/D");


 bsTree_->Branch(  "JpsiMu1EtaMC"         , &JpsiMu1EtaMC_, "JpsiMu1EtaMC/D");
 bsTree_->Branch(  "JpsiMu2EtaMC"         , &JpsiMu2EtaMC_, "JpsiMu2EtaMC/D");

 bsTree_->Branch(  "JpsiMu1PtMC"         , &JpsiMu1PtMC_, "JpsiMu1PtMC/D");
 bsTree_->Branch(  "JpsiMu2PtMC"         , &JpsiMu2PtMC_, "JpsiMu2PtMC/D");

 bsTree_->Branch(  "JpsiMu1PhiMC"         , &JpsiMu1PhiMC_, "JpsiMu1PhiMC/D");
 bsTree_->Branch(  "JpsiMu2PhiMC"         , &JpsiMu2PhiMC_, "JpsiMu2PhiMC/D");

 bsTree_->Branch(  "JpsiCosAlphaMC"         , &JpsiCosAlphaMC_, "JpsiCosAlphaMC/D");				


 bsTree_->Branch(  "JpsiGenPVy"         , &JpsiGenPVy_, "JpsiGenPVy/D");
 bsTree_->Branch(  "JpsiGenPVx"         , &JpsiGenPVx_ ,"JpsiGenPVx/D");
 bsTree_->Branch(  "JpsiGenPVz"         , &JpsiGenPVz_, "JpsiGenPVz/D");

 bsTree_->Branch(  "JpsiTrigVtxProb"         , &JpsiTrigVtxProb_,"JpsiTrigVtxProb/D"); 
 bsTree_->Branch(  "JpsiMu1TrackerMuonArbitrated"         , &JpsiMu1TrackerMuonArbitrated_,"JpsiMu1TrackerMuonArbitrated/I" ); 
 bsTree_->Branch(  "JpsiMu2TrackerMuonArbitrated"         , &JpsiMu2TrackerMuonArbitrated_ ,"JpsiMu2TrackerMuonArbitrated/I"); 


  bsTree_->Branch(  "BSx"				  , &BSx_,                              "BSx/D");
  bsTree_->Branch(  "BSy"				  , &BSy_,                              "BSy/D");
  bsTree_->Branch(  "BSz"				  , &BSz_,                              "BSz/D");
  bsTree_->Branch(  "BSdx"                           , &BSdx_,                              "BSdx/D");
  bsTree_->Branch(  "BSdy"                           , &BSdy_,                              "BSdy/D");
  bsTree_->Branch(  "BSdz"                           , &BSdz_,                              "BSdz/D");

  bsTree_->Branch(  "BSdxdz"                           , &BSdxdz_,                              "BSdxdz/D");
  bsTree_->Branch(  "BSdydz"                           , &BSdydz_,                              "BSdydz/D");




  bsTree_->Branch(  "BSsigmaZ"                           , &BSsigmaZ_,                              "BSsigmaZ/D");
  bsTree_->Branch(  "BSdsigmaZ"                           , &BSdsigmaZ_,                              "BSdsigmaZ/D");


  bsTree_->Branch(  "PVx"				  , &PVx_,                              "PVx/D");
  bsTree_->Branch(  "PVy"				  , &PVy_,                              "PVy/D");
  bsTree_->Branch(  "PVz"				  , &PVz_,                              "PVz/D");
  bsTree_->Branch(  "PVerrx"			  , &PVerrx_,                           "PVerrx/D");
  bsTree_->Branch(  "PVerry"			  , &PVerry_,                           "PVerry/D");
  bsTree_->Branch(  "PVerrz"			  , &PVerrz_,                           "PVerrz/D");
  bsTree_->Branch(  "deltaRkaon1Jpsi"   , &deltaRkaon1Jpsi_,    "deltaRkaon1Jspi/D");
  bsTree_->Branch(  "deltaRkaon2Jpsi"   , &deltaRkaon2Jpsi_,    "deltaRkaon2Jspi/D");
  bsTree_->Branch( "deltaRBCandPion", &deltaRBCandPion_, "deltaRBCandPion/D");
  bsTree_->Branch(  "MuMuDCA", &MuMuDCA_,  "MuMuDCA/D");
  bsTree_->Branch("KaonsDCA", &KaonsDCA_,"KaonsDCA/D");
  bsTree_->Branch(  "JpsiMuMuDCA_beffit"			  , &JpsiMuMuDCA_beffit_,                           "JpsiMuMuDCA_beffit/D");
  bsTree_->Branch(  "MuMuDistance"			  , &MuMuDistance_,                           "MuMuDistance/D");
  bsTree_->Branch(  "MuMuDistanceSigma"			  , &MuMuDistanceSigma_,                           "MuMuDistanceSigma/D");
  bsTree_->Branch(  "MuDr1"			  , &MuDr1_,                           "MuDr1/D");
  bsTree_->Branch(  "MuDr2"			  , &MuDr2_,                           "MuDr2/D");
  bsTree_->Branch(  "BdMuMuDCA"			  , &BdMuMuDCA_,                           "BdMuMuDCA/D");
  bsTree_->Branch(  "BdMuMuDistance"              , &BdMuMuDistance_,                           "BdMuMuDistance/D");
  bsTree_->Branch(  "BdMuMuDistanceSigma"         , &BdMuMuDistanceSigma_,                           "BdMuMuDistanceSigma/D");
  bsTree_->Branch(  "BdMuDr1"			  , &BdMuDr1_,                           "BdMuDr1/D");
  bsTree_->Branch(  "BdMuDr2"			  , &BdMuDr2_,                           "BdMuDr2/D");

  bsTree_->Branch( "PVx_refit"   , &PVx_refit_   ,         "PVx_refit/D"   );
  bsTree_->Branch( "PVy_refit"   ,	 &PVy_refit_   ,	 "PVy_refit/D"   );
  bsTree_->Branch( "PVz_refit"   ,	 &PVz_refit_   ,	 "PVz_refit/D"   );

  bsTree_->Branch( "PVerrx_refit",	 &PVerrx_refit_,	 "PVerrx_refit/D");
  bsTree_->Branch( "PVerry_refit",	 &PVerry_refit_,	 "PVerry_refit/D");
  bsTree_->Branch( "PVerrz_refit",	 &PVerrz_refit_,	 "PVerrz_refit/D");


  bsTree_->Branch( "PVx_refit_cosTheta"   , &PVx_refit_cosTheta_   ,         "PVx_refit_cosTheta/D"   );
  bsTree_->Branch( "PVy_refit_cosTheta"   ,	 &PVy_refit_cosTheta_   ,	 "PVy_refit_cosTheta/D"   );
  bsTree_->Branch( "PVz_refit_cosTheta"   ,	 &PVz_refit_cosTheta_   ,	 "PVz_refit_cosTheta/D"   );

  bsTree_->Branch( "PVx_refit_closestZ"   , &PVx_refit_closestZ_   ,         "PVx_refit_closestZ/D"   );
  bsTree_->Branch( "PVy_refit_closestZ"   ,	 &PVy_refit_closestZ_   ,	 "PVy_refit_closestZ/D"   );
  bsTree_->Branch( "PVz_refit_closestZ"   ,	 &PVz_refit_closestZ_   ,	 "PVz_refit_closestZ/D"   );

  bsTree_->Branch( "BsSVx"   ,	 &BsSVx_  , "BsSVx/D"  );
  bsTree_->Branch( "BsSVy"   ,	 &BsSVy_  ,"BsSVy/D"  );
  bsTree_->Branch( "BsSVz"   ,	 &BsSVz_  ,"BsSVz/D"  );

  bsTree_->Branch( "JpsiSVx"   ,	 &JpsiSVx_  , "JpsiSVx/D"  );
  bsTree_->Branch( "JpsiSVy"   ,	 &JpsiSVy_  ,"JpsiSVy/D"  );
  bsTree_->Branch( "JpsiSVz"   ,	 &JpsiSVz_  ,"JpsiSVz/D"  );

  bsTree_->Branch( "isPV"   , &isPV_   ,         "isPV/I"   );
  bsTree_->Branch( "isBS"   , &isBS_   ,         "isBS/I"   );

  bsTree_->Branch( "runNumber"   , &runNumber_   ,         "runNumber/I"   );
  bsTree_->Branch( "eventNumber"   , &eventNumber_   ,       "eventNumber/i"   );
  bsTree_->Branch( "lumiSection"   , &lumiSection_   ,       "lumiSection/I"   );

  bsTree_->Branch( "PUinteraction"   , &PUinteraction_   ,         "PUinteraction/I"   );


  // bsTree_->Branch( "MuRecoPt2"                   ,&MuRecoPt2_ ,                      "MuRecoPt2[15]/D");
  // bsTree_->Branch( "MuRecoChg2"                   ,&MuRecoChg2_ ,                      "MuRecoChg2[15]/I");
  // bsTree_->Branch( "MuRecoEta2"                   ,&MuRecoEta2_ ,                      "MuRecoEta2[15]/D");
  // bsTree_->Branch( "MuRecoPhi2"                   ,&MuRecoPhi2_ ,                      "MuRecoPhi2[15]/D");






  bsTree_->Branch( "MuonType", &MuonType_, "MuonType/I");

  bsTree_->Branch(  "JpsiVtxProb"			  , &JpsiVtxProb_,                      "JpsiVtxProb/D");
  bsTree_->Branch(  "BdJpsiVtxProb"			  , &BdJpsiVtxProb_,                      "BdJpsiVtxProb/D");
  bsTree_->Branch(  "JpsiM_alone"			  , &JpsiM_alone_,                      "JpsiM_alone/D");
  bsTree_->Branch(  "JpsiPhi_alone"		  , &JpsiPhi_alone_,                    "JpsiPhi_alone/D");
  bsTree_->Branch(  "JpsiEta_alone"		  , &JpsiEta_alone_,                    "JpsiEta_alone/D");
  bsTree_->Branch(  "JpsiPt_alone"		  , &JpsiPt_alone_,                     "JpsiPt_alone/D");
  bsTree_->Branch(  "JpsiMu1Pt_alone"		  , &JpsiMu1Pt_alone_,                  "JpsiMu1Pt_alone/D");
  bsTree_->Branch(  "JpsiMu2Pt_alone"		  , &JpsiMu2Pt_alone_,                  "JpsiMu2Pt_alone/D");

  bsTree_->Branch(  "JpsiMu1Phi_alone"		  , &JpsiMu1Phi_alone_,                 "JpsiMu1Phi_alone/D");
  bsTree_->Branch(  "JpsiMu2Phi_alone"		  , &JpsiMu2Phi_alone_,                 "JpsiMu2Phi_alone/D");
  bsTree_->Branch(  "JpsiMu1Eta_alone"		  , &JpsiMu1Eta_alone_,                 "JpsiMu1Eta_alone/D");
  bsTree_->Branch(  "JpsiMu2Eta_alone"		  , &JpsiMu2Eta_alone_,                 "JpsiMu2Eta_alone/D");
  bsTree_->Branch(  "MuonCat1"		  , &MuonCat1_,               "MuonCat1/I");
  bsTree_->Branch(  "MuonCat2"		  , &MuonCat2_,               "MuonCat2/I");
  bsTree_->Branch( "MuonPairDR", &MuonPairDR_, "MuonPairDR/D" );

  bsTree_->Branch(  "JpsiMuonCat1"		  , &JpsiMuonCat1_ , "JpsiMuonCat1/I"  );
  bsTree_->Branch(  "JpsiMuonCat2"		  , &JpsiMuonCat2_ , "JpsiMuonCat2/I"  );

  /// Soft muon cuts ////

  bsTree_->Branch( "Mu1TrkBSDxy", &Mu1TrkBSDxy_ , "Mu1TrkBSDxy/D" );
  bsTree_->Branch( "Mu1TrkBSDz", &Mu1TrkBSDz_ ,"Mu1TrkBSDz/D");
  bsTree_->Branch( "Mu1PixelHits", &Mu1PixelHits_ , "Mu1PixelHits/I" );
  bsTree_->Branch( "Mu1TrackerHits", &Mu1TrackerHits_ , "Mu1TrackerHits/I" );
  bsTree_->Branch( "Mu1isGood", &Mu1isGood_ , "Mu1isGood/I" );
  bsTree_->Branch( "Mu1InnerTrkHighQuality", &Mu1InnerTrkHighQuality_ , "Mu1InnerTrkHighQuality/I" );


  bsTree_->Branch( "Mu2TrkBSDxy", &Mu2TrkBSDxy_,"Mu2TrkBSDxy/D" );
  bsTree_->Branch( "Mu2TrkBSDz", &Mu2TrkBSDz_ ,"Mu2TrkBSDxy/D");
  bsTree_->Branch( "Mu2PixelHits", &Mu2PixelHits_ , "Mu2PixelHits/I" );
  bsTree_->Branch( "Mu2TrackerHits", &Mu2TrackerHits_ , "Mu2TrackerHits/I" );
  bsTree_->Branch( "Mu2isGood", &Mu2isGood_ , "Mu2isGood/I" );
  bsTree_->Branch( "Mu2InnerTrkHighQuality", &Mu2InnerTrkHighQuality_ , "Mu2InnerTrkHighQuality/I" );


  bsTree_->Branch( "BpMu1TrkBSDxy", &BpMu1TrkBSDxy_ , "BpMu1TrkBSDxy/D" );
  bsTree_->Branch( "BpMu1TrkBSDz", &BpMu1TrkBSDz_ ,"BpMu1TrkBSDz/D");
  bsTree_->Branch( "BpMu1PixelHits", &BpMu1PixelHits_ , "BpMu1PixelHits/I" );
  bsTree_->Branch( "BpMu1TrackerHits", &BpMu1TrackerHits_ , "BpMu1TrackerHits/I" );
  bsTree_->Branch( "BpMu1isGood", &BpMu1isGood_ , "BpMu1isGood/I" );
  bsTree_->Branch( "BpMu1InnerTrkHighQuality", &BpMu1InnerTrkHighQuality_ , "BpMu1InnerTrkHighQuality/I" );


  bsTree_->Branch( "BpMu2TrkBSDxy", &BpMu2TrkBSDxy_,"BpMu2TrkBSDxy/D" );
  bsTree_->Branch( "BpMu2TrkBSDz", &BpMu2TrkBSDz_ ,"BpMu2TrkBSDxy/D");
  bsTree_->Branch( "BpMu2PixelHits", &BpMu2PixelHits_ , "BpMu2PixelHits/I" );
  bsTree_->Branch( "BpMu2TrackerHits", &BpMu2TrackerHits_ , "BpMu2TrackerHits/I" );
  bsTree_->Branch( "BpMu2isGood", &BpMu2isGood_ , "BpMu2isGood/I" );
  bsTree_->Branch( "BpMu2InnerTrkHighQuality", &BpMu2InnerTrkHighQuality_ , "BpMu2InnerTrkHighQuality/I" );



  bsTree_->Branch(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_,               "JpsiMuon1Cat_alone/I");
  bsTree_->Branch(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_,               "JpsiMuon2Cat_alone/I");
  bsTree_->Branch(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_,               "BdJpsiMuon1Cat_alone/I");
  bsTree_->Branch(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_,               "BdJpsiMuon2Cat_alone/I");
  bsTree_->Branch(  "JpsiMu1d0_alone"               , &JpsiMu1d0_alone_,                  "JpsiMu1d0_alone/D");
  bsTree_->Branch(  "JpsiMu2d0_alone"               , &JpsiMu2d0_alone_,                  "JpsiMu2d0_alone/D");
  bsTree_->Branch(  "JpsiMu1dz_alone"               , &JpsiMu1dz_alone_,                  "JpsiMu1dz_alone/D");
  bsTree_->Branch(  "JpsiMu2dz_alone"               , &JpsiMu2dz_alone_,                  "JpsiMu2dz_alone/D");
  bsTree_->Branch(  "JpsiMu1chi2_alone"             , &JpsiMu1chi2_alone_,                  "JpsiMu1chi2_alone/D");
  bsTree_->Branch(  "JpsiMu2chi2_alone"             , &JpsiMu2chi2_alone_,                  "JpsiMu2chi2_alone/D");
  bsTree_->Branch(  "JpsiMu1ndof_alone"             , &JpsiMu1ndof_alone_,                  "JpsiMu1ndof_alone/I");
  bsTree_->Branch(  "JpsiMu2ndof_alone"             , &JpsiMu2ndof_alone_,                  "JpsiMu2ndof_alone/I");
  bsTree_->Branch(  "JpsiMu1nHits_alone"            , &JpsiMu1nHits_alone_,                  "JpsiMu1nHits_alone/I");
  bsTree_->Branch(  "JpsiMu2nHits_alone"            , &JpsiMu2nHits_alone_,                  "JpsiMu2nHits_alone/I");
  bsTree_->Branch(  "JpsiMu1nPixHits_alone"         , &JpsiMu1nPixHits_alone_,               "JpsiMu1nPixHits_alone/I");
  bsTree_->Branch(  "JpsiMu2nPixHits_alone"         , &JpsiMu2nPixHits_alone_,               "JpsiMu2nPixHits_alone/I");

  bsTree_->Branch(  "K1Pt_beffit"			  , &K1Pt_beffit_,                       "K1Pt_beffit/D");
  bsTree_->Branch(  "K1Pz_beffit"			  , &K1Pz_beffit_,                       "K1Pz_beffit/D");
  bsTree_->Branch(  "K1Eta_beffit"			  , &K1Eta_beffit_,                      "K1Eta_beffit/D");
  bsTree_->Branch(  "K1Phi_beffit"			  , &K1Phi_beffit_,                      "K1Phi_beffit/D");

  bsTree_->Branch(  "K2Pt_beffit"			  , &K2Pt_beffit_,                       "K2Pt_beffit/D");
  bsTree_->Branch(  "K2Pz_beffit"			  , &K2Pz_beffit_,                       "K2Pz_beffit/D");
  bsTree_->Branch(  "K2Eta_beffit"			  , &K2Eta_beffit_,                      "K2Eta_beffit/D");
  bsTree_->Branch(  "K2Phi_beffit"			  , &K2Phi_beffit_,                      "K2Phi_beffit/D");

  bsTree_->Branch(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_,                       "Mu1Pt_beffit/D");
  bsTree_->Branch(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_,                       "Mu1Pz_beffit/D");
  bsTree_->Branch(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_,                      "Mu1Eta_beffit/D");
  bsTree_->Branch(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_,                      "Mu1Phi_beffit/D");
  bsTree_->Branch(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_,                       "Mu2Pt_beffit/D");
  bsTree_->Branch(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_,                       "Mu2Pz_beffit/D");
  bsTree_->Branch(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_,                      "Mu2Eta_beffit/D");
  bsTree_->Branch(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_,                      "Mu2Phi_beffit/D");

  bsTree_->Branch(  "BdMu1Pt_beffit"			  , &BdMu1Pt_beffit_,                       "BdMu1Pt_beffit/D");
  bsTree_->Branch(  "BdMu1Pz_beffit"			  , &BdMu1Pz_beffit_,                       "BdMu1Pz_beffit/D");
  bsTree_->Branch(  "BdMu1Eta_beffit"			  , &BdMu1Eta_beffit_,                      "BdMu1Eta_beffit/D");
  bsTree_->Branch(  "BdMu1Phi_beffit"			  , &BdMu1Phi_beffit_,                      "BdMu1Phi_beffit/D");
  bsTree_->Branch(  "BdMu2Pt_beffit"			  , &BdMu2Pt_beffit_,                       "BdMu2Pt_beffit/D");
  bsTree_->Branch(  "BdMu2Pz_beffit"			  , &BdMu2Pz_beffit_,                       "BdMu2Pz_beffit/D");
  bsTree_->Branch(  "BdMu2Eta_beffit"			  , &BdMu2Eta_beffit_,                      "BdMu2Eta_beffit/D");
  bsTree_->Branch(  "BdMu2Phi_beffit"			  , &BdMu2Phi_beffit_,                      "BdMu2Phi_beffit/D");

  bsTree_->Branch(  "BsFitChi2"			  , &BsFitChi2_,                        "BsFitChi2/D");
  bsTree_->Branch(  "BsFitNdof"			  , &BsFitNdof_,                        "BsFitNdof/I");
  bsTree_->Branch(  "BsFitVtxProb"		  , &BsFitVtxProb_,                     "BsFitVtxProb/D");

  bsTree_->Branch(  "BsFitM"			  , &BsFitM_,                           "BsFitM/D");
  bsTree_->Branch(  "K1Pt_fit"			  , &K1Pt_fit_,                           "K1Pt_fit/D");
  bsTree_->Branch(  "K2Pt_fit"			  , &K2Pt_fit_,                           "K2Pt_fit/D");
  bsTree_->Branch(  "PhiM_fit"			  , &PhiM_fit_,                           "PhiM_fit/D");
  bsTree_->Branch(  "BsFitEta"			  , &BsFitEta_,                         "BsFitEta/D");
  bsTree_->Branch(  "BsFitPt"			  , &BsFitPt_,                          "BsFitPt/D");
  bsTree_->Branch(  "BsFitPz"			  , &BsFitPz_,                          "BsFitPz/D");
  bsTree_->Branch(  "BsFitPhi"			  , &BsFitPhi_,                         "BsFitPhi/D");
  bsTree_->Branch(  "BsFitVtx_x"			  , &BsFitVtx_x_,                       "BsFitVtx_x/D");
  bsTree_->Branch(  "BsFitVtx_y"			  , &BsFitVtx_y_,                       "BsFitVtx_y/D");
  bsTree_->Branch(  "BsFitVtx_z"			  , &BsFitVtx_z_,                       "BsFitVtx_z/D");
  bsTree_->Branch(  "BsM_nofit"			  , &BsM_nofit_,                        "BsM_nofit/D");
  bsTree_->Branch(  "BsPhi_nofit"			  , &BsPhi_nofit_,                      "BsPhi_nofit/D");
  bsTree_->Branch(  "BsEta_nofit"			  , &BsEta_nofit_,                      "BsEta_nofit/D");
  bsTree_->Branch(  "BsPt_nofit"			  , &BsPt_nofit_,                       "BsPt_nofit/D");
  bsTree_->Branch(  "BsPz_nofit"			  , &BsPz_nofit_,                       "BsPz_nofit/D");
  bsTree_->Branch(  "JpsiM_nofit"			  , &JpsiM_nofit_,                      "JpsiM_nofit/D");
  bsTree_->Branch(  "JpsiPhi_nofit"		  , &JpsiPhi_nofit_,                    "JpsiPhi_nofit/D");
  bsTree_->Branch(  "JpsiEta_nofit"		  , &JpsiEta_nofit_,                    "JpsiEta_nofit/D");
  bsTree_->Branch(  "JpsiPt_nofit"		  , &JpsiPt_nofit_,                     "JpsiPt_nofit/D");
  bsTree_->Branch(  "JpsiPz_nofit"		  , &JpsiPz_nofit_,                     "JpsiPz_nofit/D");
  bsTree_->Branch(  "BdJpsiM_nofit"		  , &BdJpsiM_nofit_,                      "BdJpsiM_nofit/D");
  bsTree_->Branch(  "BdJpsiPhi_nofit"		  , &BdJpsiPhi_nofit_,                    "BdJpsiPhi_nofit/D");
  bsTree_->Branch(  "BdJpsiEta_nofit"		  , &BdJpsiEta_nofit_,                    "BdJpsiEta_nofit/D");
  bsTree_->Branch(  "BsSoftMuon1"		  , &BsSoftMuon1_,                    "BsSoftMuon1/I");
  bsTree_->Branch(  "BsSoftMuon2"		  , &BsSoftMuon2_,                    "BsSoftMuon2/I");

  bsTree_->Branch(  "BdJpsiPt_nofit"		  , &BdJpsiPt_nofit_,                     "BdJpsiPt_nofit/D");
  bsTree_->Branch(  "BdJpsiPz_nofit"		  , &BdJpsiPz_nofit_,                     "BdJpsiPz_nofit/D");
  bsTree_->Branch(  "PhiM_nofit"			  , &PhiM_nofit_,                       "PhiM_nofit/D");
  bsTree_->Branch(  "PhiPhi_nofit"		  , &PhiPhi_nofit_,                     "PhiPhi_nofit/D");
  bsTree_->Branch(  "PhiEta_nofit"		  , &PhiEta_nofit_,                     "PhiEta_nofit/D");
  bsTree_->Branch(  "PhiPt_nofit"			  , &PhiPt_nofit_,                      "PhiPt_nofit/D");
  bsTree_->Branch(  "PhiPz_nofit"			  , &PhiPz_nofit_,                      "PhiPz_nofit/D");
  bsTree_->Branch(  "K1Pt_nofit"			  , &K1Pt_nofit_,                       "K1Pt_nofit/D");
  bsTree_->Branch(  "K1Pz_nofit"			  , &K1Pz_nofit_,                       "K1Pz_nofit/D");

  bsTree_->Branch(  "K1HighPurityTrack"			  , &K1HighPurityTrack_  ,"K1HighPurityTrack/D");
  bsTree_->Branch(  "K2HighPurityTrack"			  , &K2HighPurityTrack_  ,"K2HighPurityTrack/D");

  bsTree_->Branch(  "K1Eta_nofit"			  , &K1Eta_nofit_,                      "K1Eta_nofit/D");
  bsTree_->Branch(  "K1Phi_nofit"			  , &K1Phi_nofit_,                      "K1Phi_nofit/D");
  bsTree_->Branch(  "K1Key_nofit"			  , &K1Key_nofit_,                       "K1Key_nofit/I");
  bsTree_->Branch(  "K2Eta_nofit"			  , &K2Eta_nofit_,                      "K2Eta_nofit/D");
  bsTree_->Branch(  "K2Pt_nofit"			  , &K2Pt_nofit_,                       "K2Pt_nofit/D");
  bsTree_->Branch(  "K2Pz_nofit"			  , &K2Pz_nofit_,                       "K2Pz_nofit/D");
  bsTree_->Branch(  "K2Phi_nofit"			  , &K2Phi_nofit_,                      "K2Phi_nofit/D");
  bsTree_->Branch(  "K2Key_nofit"			  , &K2Key_nofit_,                       "K2Key_nofit/I");
  bsTree_->Branch(  "K1Chi2"			  , &K1Chi2_,                           "K1Chi2/D");
  bsTree_->Branch(  "K1nHits"			  , &K1nHits_,                          "K1nHits/I");
  bsTree_->Branch(  "K2Chi2"			  , &K2Chi2_,                           "K2Chi2/D");
  bsTree_->Branch(  "K2nHits"			  , &K2nHits_,                          "K2nHits/I");
  bsTree_->Branch(  "K1pixH"			  , &K1pixH_,                           "K1pixH/I");
  bsTree_->Branch(  "K1trkH"			  , &K1trkH_,                           "K1trkH/I");
  bsTree_->Branch(  "K2pixH"			  , &K2pixH_,                           "K2pixH/I");
  bsTree_->Branch(  "K2trkH"			  , &K2trkH_,                           "K2trkH/I");
  bsTree_->Branch(  "Mu1Chi2"			  , &Mu1Chi2_,                          "Mu1Chi2/D");
  bsTree_->Branch(  "Mu1nHits"			  , &Mu1nHits_,                         "Mu1nHits/I");
  bsTree_->Branch(  "Mu2Chi2"			  , &Mu2Chi2_,                          "Mu2Chi2/D");
  bsTree_->Branch(  "Mu2nHits"			  , &Mu2nHits_,                         "Mu2nHits/I");
  bsTree_->Branch(  "Mu1pixH"			  , &Mu1pixH_,                          "Mu1pixH/I");
  bsTree_->Branch(  "Mu1trkH"			  , &Mu1trkH_,                          "Mu1trkH/I");
  bsTree_->Branch(  "Mu2pixH"			  , &Mu2pixH_,                          "Mu2pixH/I");
  bsTree_->Branch(  "Mu2trkH"			  , &Mu2trkH_,                          "Mu2trkH/I");

  bsTree_->Branch(  "Mu1d0"               , &Mu1d0_,                  "Mu1d0/D");
  bsTree_->Branch(  "Mu2d0"               , &Mu2d0_,                  "Mu2d0/D");
  bsTree_->Branch(  "Mu1dz"               , &Mu1dz_,                  "Mu1dz/D");
  bsTree_->Branch(  "Mu2dz"               , &Mu2dz_,                  "Mu2dz/D");

  bsTree_->Branch(  "costheta"			  , &costheta_,                         "costheta/D");
  bsTree_->Branch(  "phi"				  , &phi_,                              "phi/D");
  bsTree_->Branch(  "cospsi"			  , &cospsi_,                           "cospsi/D");
  bsTree_->Branch(  "Bdcostheta"			  , &Bdcostheta_,                         "Bdcostheta/D");
  bsTree_->Branch(  "Bdphi"		          , &Bdphi_,                              "Bdphi/D");
  bsTree_->Branch(  "Bdcospsi"			  , &Bdcospsi_,                           "Bdcospsi/D");
  bsTree_->Branch(  "BdcosthetaMC"		  , &BdcosthetaMC_,                         "BdcosthetaMC/D");
  bsTree_->Branch(  "BdphiMC"		          , &BdphiMC_,                              "BdphiMC/D");
  bsTree_->Branch(  "BdcospsiMC"			  , &BdcospsiMC_,                           "BdcospsiMC/D");

  bsTree_->Branch(  "AngleBsDecayLength"		  , &AngleBsDecayLength_,               "AngleBsDecayLength/D");
  bsTree_->Branch(  "CosDeltaAlpha"		  , &CosDeltaAlpha_,               "CosDeltaAlpha/D");
  bsTree_->Branch(  "BdCosDeltaAlpha"		  , &BdCosDeltaAlpha_,               "BdCosDeltaAlpha/D");
  bsTree_->Branch(  "isMatched"			  , &isMatched_,                        "isMatched/I");
  bsTree_->Branch(  "isMatchedBd"			  , &isMatchedBd_,                      "isMatchedBd/I");

  bsTree_->Branch(  "isMatchedJpsi"			  , &isMatchedJpsi_,                        "isMatchedJpsi/I");
  bsTree_->Branch(  "JpsiMu1mcId"			  , &JpsiMu1mcId_,                        "JpsiMu1mcId/I");
  bsTree_->Branch(  "JpsiMu2mcId"			  , &JpsiMu2mcId_,                        "JpsiMu2mcId/I");
  bsTree_->Branch(  "JpsiMu1momId"			  , &JpsiMu1momId_ , "JpsiMu1momId/I");
  bsTree_->Branch(  "JpsiMu2momId"			  , &JpsiMu2momId_,  "JpsiMu2momId/I");


  bsTree_->Branch(  "BsLxy"			  , &BsLxy_,                            "BsLxy/D");
  bsTree_->Branch(  "BsCt"                          , &BsCt_,                             "BsCt/D");
  bsTree_->Branch(  "BsCtErr"                       , &BsCtErr_,                             "BsCtErr/D");
  bsTree_->Branch(  "BsCt3D"                        , &BsCt3D_,                             "BsCt3D/D");
  bsTree_->Branch(  "BsCt2D"                        , &BsCt2D_,                             "BsCt2D/D");
  bsTree_->Branch(  "BsCt2DBS"                        , &BsCt2DBS_,                             "BsCt2DBS/D");
  bsTree_->Branch(  "BdCt2DBS"                        , &BdCt2DBS_,                             "BdCt2DBS/D");
  bsTree_->Branch(  "BdCt2DMC"                        , &BdCt2DMC_,                             "BdCt2DMC/D");
  bsTree_->Branch(  "BdCt3DMC"                        , &BdCt3DMC_,                             "BdCt3DMC/D");
  bsTree_->Branch(  "BsCtMPV"                       , &BsCtMPV_,                             "BsCtMPV/D");
  bsTree_->Branch(  "BsCtErr3D"                     , &BsCtErr3D_,                             "BsCtErr3D/D");
  bsTree_->Branch(  "BsCtErr2D"                     , &BsCtErr2D_,                             "BsCtErr2D/D");
  bsTree_->Branch(  "BsCtErr2DBS"                     , &BsCtErr2DBS_,                             "BsCtErr2DBS/D");
  bsTree_->Branch(  "BsCtErr2DClosestZ"               , &BsCtErr2DClosestZ_,              "BsCtErr2DClosestZ/D");




  bsTree_->Branch(  "BdCtErr2DBS"                     , &BdCtErr2DBS_,                             "BdCtErr2DBS/D");
  bsTree_->Branch(  "BsCtErrMPV"                    , &BsCtErrMPV_,                             "BsCtErrMPV/D");
  bsTree_->Branch(  "BsCt3Drefit"                   , &BsCt3Drefit_,                             "BsCt3Drefit/D");
  bsTree_->Branch(  "BsCt2Drefit"                   , &BsCt2Drefit_,                             "BsCt2Drefit/D");
  bsTree_->Branch(  "BsCtMPVrefit"                  , &BsCtMPVrefit_,                             "BsCtMPVrefit/D");
  bsTree_->Branch(  "BsCtErr2D2"                    , &BsCtErr2D2_,                             "BsCtErr2D2/D");
  bsTree_->Branch(  "BsCtErrMPV"                    , &BsCtErrMPV_,                             "BsCtErrMPV/D");
  bsTree_->Branch(  "BsCtErr3Drefit"                , &BsCtErr3Drefit_,                            "BsCtErr3Drefit/D");
  bsTree_->Branch(  "BsCtErr2Drefit"                , &BsCtErr2Drefit_,                            "BsCtErr2Drefit/D");
  bsTree_->Branch(  "BsCtErrMPVrefit"               , &BsCtErrMPVrefit_,                             "BsCtErrMPVrefit/D");

  bsTree_->Branch(  "BsLxyErr"			  , &BsLxyErr_,                            "BsLxyErr/D");
  bsTree_->Branch(  "JpsiNumberOfCandidates"          , &JpsiNumberOfCandidates_,            "JpsiNumberOfCandidates/I");
  bsTree_->Branch(  "JpsiGenNumberOfCandidates"        , &JpsiGenNumberOfCandidates_, "JpsiGenNumberOfCandidates/I");
  bsTree_->Branch(  "PhiNumberOfCandidatesBeforeFit"          , &PhiNumberOfCandidatesBeforeFit_,            "PhiNumberOfCandidatesBeforeFit/I");
  bsTree_->Branch(  "BsNumberOfCandidatesBeforeFit"          , &BsNumberOfCandidatesBeforeFit_,            "BsNumberOfCandidatesBeforeFit/I");
  bsTree_->Branch(  "BsNumberOfCandidatesAfterFit"          , &BsNumberOfCandidatesAfterFit_,            "BsNumberOfCandidatesAfterFit/I");
  bsTree_->Branch(  "BsNumberOfCandidatesAfterBestFit"          , &BsNumberOfCandidatesAfterBestFit_,            "BsNumberOfCandidatesAfterBestFit/I");

  bsTree_->Branch(  "BsErrX"			  , &BsErrX_,                           "BsErrX/D");
  bsTree_->Branch(  "BsErrY"			  , &BsErrY_,                           "BsErrY/D");
  bsTree_->Branch(  "BsErrXY"			  , &BsErrXY_,                          "BsErrXY/D");

  bsTree_->Branch(  "K1trkLay"			  , &K1trkLay_,                         "K1trkLay/I");
  bsTree_->Branch(  "K1muDTh"			  , &K1muDTh_,                          "K1muDTh/I");
  bsTree_->Branch(  "K1muCSCh"			  , &K1muCSCh_,                         "K1muCSCh/I");
  bsTree_->Branch(  "K1muRPCh"			  , &K1muRPCh_,                         "K1muRPCh/I");
  bsTree_->Branch(  "K2trkLay"			  , &K2trkLay_,                         "K2trkLay/I");
  bsTree_->Branch(  "K2muDTh"			  , &K2muDTh_,                          "K2muDTh/I");
  bsTree_->Branch(  "K2muCSCh"			  , &K2muCSCh_,                         "K2muCSCh/I");
  bsTree_->Branch(  "K2muRPCh"			  , &K2muRPCh_,                         "K2muRPCh/I");
  bsTree_->Branch(  "Mu1trkLay"			  , &Mu1trkLay_,                        "Mu1trkLay/I");
  bsTree_->Branch(  "Mu1muDTh"			  , &Mu1muDTh_,                         "Mu1muDTh/I");
  bsTree_->Branch(  "Mu1muCSCh"			  , &Mu1muCSCh_,                        "Mu1muCSCh/I");
  bsTree_->Branch(  "Mu1muRPCh"			  , &Mu1muRPCh_,                        "Mu1muRPCh/I");
  bsTree_->Branch(  "Mu2trkLay"			  , &Mu2trkLay_,                        "Mu2trkLay/I");
  bsTree_->Branch(  "Mu2muDTh"			  , &Mu2muDTh_,                         "Mu2muDTh/I");
  bsTree_->Branch(  "Mu2muCSCh"			  , &Mu2muCSCh_,                        "Mu2muCSCh/I");
  bsTree_->Branch(  "Mu2muRPCh"			  , &Mu2muRPCh_,                        "Mu2muRPCh/I");
  bsTree_->Branch(  "K1mcId"			  , &K1mcId_,                           "K1mcId/I");
  bsTree_->Branch(  "K1momId"			  , &K1momId_,                          "K1momId/I");
  bsTree_->Branch(  "K1gmomId"			  , &K1gmomId_,                         "K1gmomId/I");
  bsTree_->Branch(  "K2mcId"			  , &K2mcId_,                           "K2mcId/I");
  bsTree_->Branch(  "K2momId"			  , &K2momId_,                          "K2momId/I");
  bsTree_->Branch(  "K2gmomId"			  , &K2gmomId_,                         "K2gmomId/I");
  bsTree_->Branch(  "Mu1mcId"			  , &Mu1mcId_,                          "Mu1mcId/I");
  bsTree_->Branch(  "Mu1momId"			  , &Mu1momId_,                         "Mu1momId/I");
  bsTree_->Branch(  "Mu1gmomId"			  , &Mu1gmomId_,                        "Mu1gmomId/I");
  bsTree_->Branch(  "Mu2mcId"			  , &Mu2mcId_,                          "Mu2mcId/I");
  bsTree_->Branch(  "Mu2momId"			  , &Mu2momId_,                         "Mu2momId/I");
  bsTree_->Branch(  "Mu2gmomId"			  , &Mu2gmomId_,                        "Mu2gmomId/I");
  bsTree_->Branch(  "Mu1GlobalMuonPromptTight"      , &Mu1GlobalMuonPromptTight_,         "Mu1GlobalMuonPromptTight/I");
  bsTree_->Branch(  "Mu2GlobalMuonPromptTight"      , &Mu2GlobalMuonPromptTight_,         "Mu2GlobalMuonPromptTight/I");
  bsTree_->Branch(  "Mu1TrackerMuonArbitrated"      , &Mu1TrackerMuonArbitrated_,         "Mu1TrackerMuonArbitrated/I");
  bsTree_->Branch(  "Mu1TMLastStationTight"         , &Mu1TMLastStationTight_,            "Mu1TMLastStationTight/I");
  bsTree_->Branch(  "Mu1TMOneStationTight"          , &Mu1TMOneStationTight_,             "Mu1TMOneStationTight/I");
  bsTree_->Branch(  "Mu1TMLastStationOptimizedLowPtTight", &Mu1TMLastStationOptimizedLowPtTight_, "Mu1TMLastStationOptimizedLowPtTight/I");
  bsTree_->Branch(  "Mu1TMLastStationAngTight"      , &Mu1TMLastStationAngTight_,         "Mu1TMLastStationAngTight/I");
  bsTree_->Branch(  "Mu1TMOneStationAngTight"       , &Mu1TMOneStationAngTight_,          "Mu1TMOneStationAngTight/I");
  bsTree_->Branch(  "Mu1TMLastStationOptimizedBarrelLowPtTight", &Mu1TMLastStationOptimizedBarrelLowPtTight_, "Mu1TMLastStationOptimizedBarrelLowPtTight/I");
  bsTree_->Branch(  "Mu2TrackerMuonArbitrated"      , &Mu2TrackerMuonArbitrated_,         "Mu2TrackerMuonArbitrated/I");
  bsTree_->Branch(  "Mu2TMLastStationTight"         , &Mu2TMLastStationTight_,            "Mu2TMLastStationTight/I");
  bsTree_->Branch(  "Mu2TMOneStationTight"          , &Mu2TMOneStationTight_,             "Mu2TMOneStationTight/I");
  bsTree_->Branch(  "Mu2TMLastStationOptimizedLowPtTight", &Mu2TMLastStationOptimizedLowPtTight_, "Mu2TMLastStationOptimizedLowPtTight/I");
  bsTree_->Branch(  "Mu2TMLastStationAngTight"      , &Mu2TMLastStationAngTight_,         "Mu2TMLastStationAngTight/I");
  bsTree_->Branch(  "Mu2TMOneStationAngTight"       , &Mu2TMOneStationAngTight_,          "Mu2TMOneStationAngTight/I");
  bsTree_->Branch(  "Mu2TMLastStationOptimizedBarrelLowPtTight", &Mu2TMLastStationOptimizedBarrelLowPtTight_, "Mu2TMLastStationOptimizedBarrelLowPtTight/I");

  bsTree_->Branch(  "BsCt2DPVClosestZ"                     , &BsCt2DPVClosestZ_,                             "BsCt2DPVClosestZ/D");
  bsTree_->Branch(  "BsCt3DPVClosestZ"                     , &BsCt3DPVClosestZ_,                             "BsCt3DPVClosestZ/D");
  bsTree_->Branch(  "BsCtErr2DCostheta"                     , &BsCtErr2DCostheta_,                             "BsCtErr2DCostheta/D");
  bsTree_->Branch(  "BsCtErr2DCosthetaOld"                     , &BsCtErr2DCosthetaOld_,                             "BsCtErr2DCosthetaOld/D");

  bsTree_->Branch(  "BsCtErr2DBSOld"                     , &BsCtErr2DBSOld_,                             "BsCtErr2DBSOld/D");
  bsTree_->Branch("BsCt2DBSOld", &BsCt2DBSOld_, "BsCt2DBSOld/D");

  bsTree_->Branch(  "BsCtErr2DOld"                     , &BsCtErr2DOld_,                             "BsCtErr2DOld/D");
  bsTree_->Branch("BsCt2DOld", &BsCt2DOld_, "BsCt2DOld/D");

 ///
  bsTree_->Branch("BsCtErr2DBS_JpsiVtx", &BsCtErr2DBS_JpsiVtx_, "BsCtErr2DBS_JpsiVtx/D");
  bsTree_->Branch("BsCtErr2D_JpsiVtx", &BsCtErr2D_JpsiVtx_, "BsCtErr2D_JpsiVtx/D");
  bsTree_->Branch("BsCtErr2DClosestZ_JpsiVtx", &BsCtErr2DClosestZ_JpsiVtx_, "BsCtErr2DClosestZ_JpsiVtx/D");
  bsTree_->Branch("BsCtErr2DCostheta_JpsiVtx", &BsCtErr2DCostheta_JpsiVtx_, "BsCtErr2DCostheta_JpsiVtx/D");
 ///

  bsTree_->Branch("BsCtErr2DClosestZOld", &BsCtErr2DClosestZOld_ ,"BsCtErr2DClosestZOld/D");
  bsTree_->Branch("BsCt2DPVClosestZOld", &BsCt2DPVClosestZOld_ ,"BsCt2DPVClosestZOld/D");

  bsTree_->Branch("BsCt3DPVCosTheta", &BsCt3DPVCosTheta_, "BsCt3DPVCosTheta/D");
  bsTree_->Branch("BsCt2DPVCosTheta", &BsCt2DPVCosTheta_, "BsCt2DPVCosTheta/D");
  bsTree_->Branch("BsCt2DPVCosThetaOld", &BsCt2DPVCosThetaOld_, "BsCt2DPVCosThetaOld/D");
  bsTree_->Branch(  "BsDist3d"			  , &BsDist3d_,                         "BsDist3d/D");
  bsTree_->Branch(  "BsDist3dErr"			  , &BsDist3dErr_,                      "BsDist3dErr/D");
  bsTree_->Branch(  "BsTime3d"			  , &BsTime3d_,                         "BsTime3d/D");
  bsTree_->Branch(  "BsTime3dErr"			  , &BsTime3dErr_,                      "BsTime3dErr/D");
  bsTree_->Branch(  "BsDist2d"			  , &BsDist2d_,                         "BsDist2d/D");
  bsTree_->Branch(  "BsDist2dErr"			  , &BsDist2dErr_,                      "BsDist2dErr/D");
  bsTree_->Branch(  "BsTime2d"			  , &BsTime2d_,                         "BsTime2d/D");
  bsTree_->Branch(  "BsTime2dErr"			  , &BsTime2dErr_,                      "BsTime2dErr/D");
  bsTree_->Branch(  "dedxTrk"			  , &dedxTrk_,                          "dedxTrk/D");
  bsTree_->Branch(  "errdedxTrk"			  , &errdedxTrk_,                       "errdedxTrk/D");
  bsTree_->Branch(  "numdedxTrk"			  , &numdedxTrk_,                       "numdedxTrk/I");
  bsTree_->Branch(  "iPassedCutIdent"		  , &iPassedCutIdent_,                  "iPassedCutIdent/I");
  bsTree_->Branch(  "iPassedCutIdentBd"		  , &iPassedCutIdentBd_,                "iPassedCutIdentBd/I");
  bsTree_->Branch(  "BdTrack1Charge"		  , &BdTrack1Charge_,                "BdTrack1Charge/I");
  bsTree_->Branch(  "K1Fit_par"			  , K1Fit_par_,                         "K1Fit_par[7]/D");
  bsTree_->Branch(  "K2Fit_par"			  , K2Fit_par_,                         "K2Fit_par[7]/D");
  bsTree_->Branch(  "K1Fit_sigX"			  , &K1Fit_sigX_,                       "K1Fit_sigX/D");
  bsTree_->Branch(  "K1Fit_sigY"			  , &K1Fit_sigY_,                       "K1Fit_sigY/D");
  bsTree_->Branch(  "K1Fit_sigZ"			  , &K1Fit_sigZ_,                       "K1Fit_sigZ/D");
  bsTree_->Branch(  "K2Fit_sigX"			  , &K2Fit_sigX_,                       "K2Fit_sigX/D");
  bsTree_->Branch(  "K2Fit_sigY"			  , &K2Fit_sigY_,                       "K2Fit_sigY/D");
  bsTree_->Branch(  "K2Fit_sigZ"			  , &K2Fit_sigZ_,                       "K2Fit_sigZ/D");
  bsTree_->Branch(  "K1Fit_sigPX"			  , &K1Fit_sigPX_,                      "K1Fit_sigPX/D");
  bsTree_->Branch(  "K1Fit_sigPY"			  , &K1Fit_sigPY_,                      "K1Fit_sigPY/D");
  bsTree_->Branch(  "K1Fit_sigPZ"			  , &K1Fit_sigPZ_,                      "K1Fit_sigPZ/D");
  bsTree_->Branch(  "K2Fit_sigPX"			  , &K2Fit_sigPX_,                      "K2Fit_sigPX/D");
  bsTree_->Branch(  "K2Fit_sigPY"			  , &K2Fit_sigPY_,                      "K2Fit_sigPY/D");
  bsTree_->Branch(  "K2Fit_sigPZ"			  , &K2Fit_sigPZ_,                      "K2Fit_sigPZ/D");

  bsTree_->Branch(  "GenNumberOfBdecays"		  , &GenNumberOfBdecays_,               "GenNumberOfBdecays/I");
  bsTree_->Branch(  "BmesonsId"			  , BmesonsId_,                         "BmesonsId[10]/I");
  bsTree_->Branch(  "BDauIdMC"			  , BDauIdMC_,                          "BDauIdMC[10][15]/I");
  bsTree_->Branch(  "BDauDauIdMC"			  , BDauDauIdMC_,                     	"BDauDauIdMC[10][15][10]/I");
  bsTree_->Branch(  "GenNumberOfDaughters"	  , GenNumberOfDaughters_,              "GenNumberOfDaughters[10]/I");
  bsTree_->Branch(  "GenNumberOfDaughtersDaughters" , GenNumberOfDaughtersDaughters_,     "GenNumberOfDaughtersDaughters[10][15]/I");
  bsTree_->Branch(  "BDauMMC"			  , BDauMMC_,                           "BDauMMC[10][15]/D");
  bsTree_->Branch(  "BDauPtMC"			  , BDauPtMC_,                          "BDauPtMC[10][15]/D");
  bsTree_->Branch(  "BDauPzMC"			  , BDauPzMC_,                          "BDauPzMC[10][15]/D");
  bsTree_->Branch(  "BDauEtaMC"			  , BDauEtaMC_,                         "BDauEtaMC[10][15]/D");
  bsTree_->Branch(  "BDauPhiMC"			  , BDauPhiMC_,                         "BDauPhiMC[10][15]/D");
  bsTree_->Branch(  "BDauDauMMC"			  , BDauDauMMC_,                      	"BDauDauMMC[10][15][10]/D");
  bsTree_->Branch(  "BDauDauPtMC"			  , BDauDauPtMC_,                     	"BDauDauPtMC[10][15][10]/D");
  bsTree_->Branch(  "BDauDauPzMC"			  , BDauDauPzMC_,                     	"BDauDauPzMC[10][15][10]/D");
  bsTree_->Branch(  "BDauDauEtaMC"		  , BDauDauEtaMC_,                    	"BDauDauEtaMC[10][15][10]/D");
  bsTree_->Branch(  "BDauDauPhiMC"		  , BDauDauPhiMC_,                    	"BDauDauPhiMC[10][15][10]/D");
  bsTree_->Branch(  "BMMC"			  , BMMC_,                              "BMMC[10]/D");
  bsTree_->Branch(  "BsCt3DMC"			  , &BsCt3DMC_,                              "BsCt3DMC/D");
  bsTree_->Branch(  "BsCt2DMC"			  , &BsCt2DMC_,                              "BsCt2DMC/D");
  bsTree_->Branch(  "BpCt2DMC"			  , &BpCt2DMC_,                              "BpCt2DMC/D");
  bsTree_->Branch(  "BsIniFlavour"			  , &BsIniFlavour_,                              "BsIniFlavour/I");
  bsTree_->Branch(  "BsEndFlavour"			  , &BsEndFlavour_,                              "BsEndFlavour/I");
  bsTree_->Branch(  "BdIniFlavour"                        , &BdIniFlavour_,                              "BdIniFlavour/I");
  bsTree_->Branch(  "BdEndFlavour"			  , &BdEndFlavour_,                              "BdEndFlavour/I");
  bsTree_->Branch(  "ChannelID"			  , &ChannelID_,                              "ChannelID/I");
  bsTree_->Branch(  "BdChannelID"			  , &BdChannelID_,                              "BdChannelID/I");
  bsTree_->Branch(  "BPtMC"			  , BPtMC_,                             "BPtMC[10]/D");
  bsTree_->Branch(  "BPxMC"			  , BPxMC_,                             "BPxMC[10]/D");
  bsTree_->Branch(  "BPyMC"			  , BPyMC_,                             "BPyMC[10]/D");
  bsTree_->Branch(  "BPzMC"			  , BPzMC_,                             "BPzMC[10]/D");
  bsTree_->Branch(  "BEtaMC"			  , BEtaMC_,                            "BEtaMC[10]/D");
  bsTree_->Branch(  "BPhiMC"			  , BPhiMC_,                            "BPhiMC[10]/D");

  bsTree_->Branch(  "costhetaMC"			  , costhetaMC_,                        "costhetaMC[10]/D");
  bsTree_->Branch(  "phiMC"			  , phiMC_,                             "phiMC[10]/D");
  bsTree_->Branch(  "cospsiMC"			  , cospsiMC_,                          "cospsiMC[10]/D");
  bsTree_->Branch(  "BscosthetaMC"			  , &BscosthetaMC_,                        "BscosthetaMC/D");
  bsTree_->Branch(  "BsphiMC"			  , &BsphiMC_,                             "BsphiMC/D");
  bsTree_->Branch(  "BscospsiMC"			  , &BscospsiMC_,                          "BscospsiMC/D");

  bsTree_->Branch(  "BVtxMC_x" ,   BVtxMC_x_  , "BVtxMC_x[10]/D" );
  bsTree_->Branch(  "BVtxMC_y" ,	 BVtxMC_y_  , "BVtxMC_y[10]/D" );
  bsTree_->Branch(  "BVtxMC_z" ,	 BVtxMC_z_  , "BVtxMC_z[10]/D" );
  bsTree_->Branch(  "BSVtxMC_x",	 BSVtxMC_x_ , "BSVtxMC_x[10]/D");
  bsTree_->Branch(  "BSVtxMC_y",	 BSVtxMC_y_ , "BSVtxMC_y[10]/D");
  bsTree_->Branch(  "BSVtxMC_z",	 BSVtxMC_z_ , "BSVtxMC_z[10]/D");
  bsTree_->Branch(  "BLxy_MC"  ,	 BLxy_MC_   , "BLxy_MC[10]/D"  );
  bsTree_->Branch(  "BCt_MC"   ,	 BCt_MC_    , "BCt_MC[10]/D"   );
  bsTree_->Branch(  "BCt_MC2D"   ,	 BCt_MC2D_    , "BCt_MC2D[10]/D"   );
  bsTree_->Branch(  "BCt_MC3D"   ,	 BCt_MC3D_    , "BCt_MC3D[10]/D"   );


 bsTree_->Branch(  "genBpPV_z"			  , &genBpPV_z_  ,"genBpPV_z/D");
  bsTree_->Branch(  "genBpPV_y"			  , &genBpPV_y_  ,"genBpPV_y/D");
  bsTree_->Branch(  "genBpPV_x"			  , &genBpPV_x_  ,"genBpPV_x/D");

  bsTree_->Branch(  "genBsPV_z"			  , &genBsPV_z_  ,"genBsPV_z/D");
  bsTree_->Branch(  "genBsPV_y"			  , &genBsPV_y_  ,"genBsPV_y/D");
  bsTree_->Branch(  "genBsPV_x"			  , &genBsPV_x_  ,"genBsPV_x/D");

  bsTree_->Branch(  "genBsVtx_z"			  , &genBsVtx_z_,                       "genBsVtx_z/D");
  bsTree_->Branch(  "genBsVtx_y"			  , &genBsVtx_y_,                       "genBsVtx_y/D");
  bsTree_->Branch(  "genBsVtx_x"			  , &genBsVtx_x_,                       "genBsVtx_x/D");
  bsTree_->Branch(  "genBsSVtx_z"			  , &genBsSVtx_z_,                      "genBsSVtx_z/D");
  bsTree_->Branch(  "genBsSVtx_y" 		  , &genBsSVtx_y_,                      "genBsSVtx_y/D");
  bsTree_->Branch(  "genBsSVtx_x"			  , &genBsSVtx_x_,                      "genBsSVtx_x/D");
  bsTree_->Branch(  "isGenJpsiEvent"		  , &isGenJpsiEvent_,                   "isGenJpsiEvent/I");
  bsTree_->Branch(  "BdFitChi2_Hyp1"		  , &BdFitChi2_Hyp1_,                   "BdFitChi2_Hyp1/D");
  bsTree_->Branch(  "BdFitNdof_Hyp1"		  , &BdFitNdof_Hyp1_,                   "BdFitNdof_Hyp1/I");
  bsTree_->Branch(  "BdFitVtxProb_Hyp1"		  , &BdFitVtxProb_Hyp1_,                "BdFitVtxProb_Hyp1/D");
  bsTree_->Branch(  "BdFitVtxProb"		  , &BdFitVtxProb_,                "BdFitVtxProb/D");


  bsTree_->Branch(  "BdFitM_Hyp1"			  , &BdFitM_Hyp1_,                      "BdFitM_Hyp1/D");
  bsTree_->Branch(  "BdFitEta_Hyp1"		  , &BdFitEta_Hyp1_,                    "BdFitEta_Hyp1/D");
  bsTree_->Branch(  "BdFitPt_Hyp1"		  , &BdFitPt_Hyp1_,                     "BdFitPt_Hyp1/D");
  bsTree_->Branch(  "BdFitPz_Hyp1"		  , &BdFitPz_Hyp1_,                     "BdFitPz_Hyp1/D");
  bsTree_->Branch(  "BdFitPhi_Hyp1"		  , &BdFitPhi_Hyp1_,                    "BdFitPhi_Hyp1/D");
  bsTree_->Branch(  "BdFitVtx_x_Hyp1"		  , &BdFitVtx_x_Hyp1_,                  "BdFitVtx_x_Hyp1/D");
  bsTree_->Branch(  "BdFitVtx_y_Hyp1"		  , &BdFitVtx_y_Hyp1_,                  "BdFitVtx_y_Hyp1/D");
  bsTree_->Branch(  "BdFitVtx_z_Hyp1"		  , &BdFitVtx_z_Hyp1_,                  "BdFitVtx_z_Hyp1/D");
  bsTree_->Branch(  "BdM_nofit"		  , &BdM_nofit_,                   "BdM_nofit/D");
  bsTree_->Branch(  "BdPhi_nofit"		  , &BdPhi_nofit_,                 "BdPhi_nofit/D");
  bsTree_->Branch(  "BdEta_nofit"		  , &BdEta_nofit_,                 "BdEta_nofit/D");
  bsTree_->Branch(  "BdPt_nofit"		  , &BdPt_nofit_,                  "BdPt_nofit/D");
  bsTree_->Branch(  "BdPz_nofit"		  , &BdPz_nofit_,                  "BdPz_nofit/D");
  bsTree_->Branch(  "KstarMass_nofit_Hyp1"	  , &KstarMass_nofit_Hyp1_,             "KstarMass_nofit_Hyp1/D");
  bsTree_->Branch(  "KstarMass_nofit_Hyp2"	  , &KstarMass_nofit_Hyp2_,             "KstarMass_nofit_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_par_Hyp1"		  , BdK1_kpi_par_Hyp1_,                 "BdK1_kpi_par_Hyp1[7]/D");
  bsTree_->Branch(  "BdK2_kpi_par_Hyp1"		  , BdK2_kpi_par_Hyp1_,                 "BdK2_kpi_par_Hyp1[7]/D");
  bsTree_->Branch(  "BdK1_kpi_sigX_Hyp1"		  , &BdK1_kpi_sigX_Hyp1_,               "BdK1_kpi_sigX_Hyp1/D");
  bsTree_->Branch(  "BdK1_kpi_sigY_Hyp1"		  , &BdK1_kpi_sigY_Hyp1_,               "BdK1_kpi_sigY_Hyp1/D");
  bsTree_->Branch(  "BdK1_kpi_sigZ_Hyp1"		  , &BdK1_kpi_sigZ_Hyp1_,               "BdK1_kpi_sigZ_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigX_Hyp1"		  , &BdK2_kpi_sigX_Hyp1_,               "BdK2_kpi_sigX_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigY_Hyp1"		  , &BdK2_kpi_sigY_Hyp1_,               "BdK2_kpi_sigY_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigZ_Hyp1"		  , &BdK2_kpi_sigZ_Hyp1_,               "BdK2_kpi_sigZ_Hyp1/D");
  bsTree_->Branch(  "BdK1_kpi_sigPX_Hyp1"		  , &BdK1_kpi_sigPX_Hyp1_,              "BdK1_kpi_sigPX_Hyp1/D");
  bsTree_->Branch(  "BdK1_kpi_sigPY_Hyp1"		  , &BdK1_kpi_sigPY_Hyp1_,              "BdK1_kpi_sigPY_Hyp1/D");
  bsTree_->Branch(  "BdK1_kpi_sigPZ_Hyp1"		  , &BdK1_kpi_sigPZ_Hyp1_,              "BdK1_kpi_sigPZ_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigPX_Hyp1"		  , &BdK2_kpi_sigPX_Hyp1_,              "BdK2_kpi_sigPX_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigPY_Hyp1"		  , &BdK2_kpi_sigPY_Hyp1_,              "BdK2_kpi_sigPY_Hyp1/D");
  bsTree_->Branch(  "BdK2_kpi_sigPZ_Hyp1"		  , &BdK2_kpi_sigPZ_Hyp1_,              "BdK2_kpi_sigPZ_Hyp1/D");
  bsTree_->Branch(  "BdFitChi2_Hyp2"		  , &BdFitChi2_Hyp2_,                   "BdFitChi2_Hyp2/D");
  bsTree_->Branch(  "BdFitNdof_Hyp2"		  , &BdFitNdof_Hyp2_,                   "BdFitNdof_Hyp2/I");
  bsTree_->Branch(  "BdFitVtxProb_Hyp2"		  , &BdFitVtxProb_Hyp2_,                "BdFitVtxProb_Hyp2/D");
  bsTree_->Branch(  "BdFitM_Hyp2"			  , &BdFitM_Hyp2_,                      "BdFitM_Hyp2/D");
  bsTree_->Branch(  "BdFitEta_Hyp2"		  , &BdFitEta_Hyp2_,                    "BdFitEta_Hyp2/D");
  bsTree_->Branch(  "BdFitPt_Hyp2"		  , &BdFitPt_Hyp2_,                     "BdFitPt_Hyp2/D");
  bsTree_->Branch(  "BdFitPz_Hyp2"		  , &BdFitPz_Hyp2_,                     "BdFitPz_Hyp2/D");
  bsTree_->Branch(  "BdFitPhi_Hyp2"		  , &BdFitPhi_Hyp2_,                    "BdFitPhi_Hyp2/D");
  bsTree_->Branch(  "BdFitVtx_x_Hyp2"		  , &BdFitVtx_x_Hyp2_,                  "BdFitVtx_x_Hyp2/D");
  bsTree_->Branch(  "BdFitVtx_y_Hyp2"		  , &BdFitVtx_y_Hyp2_,                  "BdFitVtx_y_Hyp2/D");
  bsTree_->Branch(  "BdFitVtx_z_Hyp2"		  , &BdFitVtx_z_Hyp2_,                  "BdFitVtx_z_Hyp2/D");
  bsTree_->Branch(  "BdNumberOfCandidates"          , &BdNumberOfCandidates_,            "BdNumberOfCandidates/I");
  bsTree_->Branch(  "BdK1_kpi_par_Hyp2"		  , BdK1_kpi_par_Hyp2_,                 "BdK1_kpi_par_Hyp2[7]/D");
  bsTree_->Branch(  "BdK2_kpi_par_Hyp2"		  , BdK2_kpi_par_Hyp2_,                 "BdK2_kpi_par_Hyp2[7]/D");
  bsTree_->Branch(  "BdK1_kpi_sigX_Hyp2"		  , &BdK1_kpi_sigX_Hyp2_,               "BdK1_kpi_sigX_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_sigY_Hyp2"		  , &BdK1_kpi_sigY_Hyp2_,               "BdK1_kpi_sigY_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_sigZ_Hyp2"		  , &BdK1_kpi_sigZ_Hyp2_,               "BdK1_kpi_sigZ_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigX_Hyp2"		  , &BdK2_kpi_sigX_Hyp2_,               "BdK2_kpi_sigX_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigY_Hyp2"		  , &BdK2_kpi_sigY_Hyp2_,               "BdK2_kpi_sigY_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigZ_Hyp2"		  , &BdK2_kpi_sigZ_Hyp2_,               "BdK2_kpi_sigZ_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_sigPX_Hyp2"		  , &BdK1_kpi_sigPX_Hyp2_,              "BdK1_kpi_sigPX_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_sigPY_Hyp2"		  , &BdK1_kpi_sigPY_Hyp2_,              "BdK1_kpi_sigPY_Hyp2/D");
  bsTree_->Branch(  "BdK1_kpi_sigPZ_Hyp2"		  , &BdK1_kpi_sigPZ_Hyp2_,              "BdK1_kpi_sigPZ_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigPX_Hyp2"		  , &BdK2_kpi_sigPX_Hyp2_,              "BdK2_kpi_sigPX_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigPY_Hyp2"		  , &BdK2_kpi_sigPY_Hyp2_,              "BdK2_kpi_sigPY_Hyp2/D");
  bsTree_->Branch(  "BdK2_kpi_sigPZ_Hyp2"		  , &BdK2_kpi_sigPZ_Hyp2_,              "BdK2_kpi_sigPZ_Hyp2/D");
  bsTree_->Branch(  "BdK1Pt_nofit" 		  , &BdK1Pt_nofit_,                     "BdK1Pt_nofit/D");
  bsTree_->Branch(  "BdK1Pz_nofit" 		  , &BdK1Pz_nofit_,                     "BdK1Pz_nofit/D");
  bsTree_->Branch(  "BdK1Eta_nofit" 		  , &BdK1Eta_nofit_,                    "BdK1Eta_nofit/D");
  bsTree_->Branch(  "BdK1Phi_nofit" 		  , &BdK1Phi_nofit_,                    "BdK1Phi_nofit/D");
  bsTree_->Branch(  "BdK1Key_nofit" 		  , &BdK1Key_nofit_,                     "BdK1Key_nofit/I");
  bsTree_->Branch(  "BdK2Pt_nofit" 		  , &BdK2Pt_nofit_,                     "BdK2Pt_nofit/D");
  bsTree_->Branch(  "BdK2Pz_nofit" 		  , &BdK2Pz_nofit_,                     "BdK2Pz_nofit/D");
  bsTree_->Branch(  "BdK2Eta_nofit" 		  , &BdK2Eta_nofit_,                    "BdK2Eta_nofit/D");
  bsTree_->Branch(  "BdK2Phi_nofit" 		  , &BdK2Phi_nofit_,                    "BdK2Phi_nofit/D");
  bsTree_->Branch(  "BdK2Key_nofit" 		  , &BdK2Key_nofit_,                     "BdK2Key_nofit/I");

  bsTree_->Branch(  "BdPVx_refit" , &BdPVx_refit_ ,  "BdPVx_refit/D" );
  bsTree_->Branch(  "BdPVy_refit" ,&BdPVy_refit_ ,    "BdPVy_refit/D");
  bsTree_->Branch(  "BdPVz_refit" , &BdPVz_refit_ ,  "BdPVz_refit/D" );
  bsTree_->Branch(  "BdPVerrx_refit" , &BdPVerrx_refit_ , "BdPVerrx_refit/D" );
  bsTree_->Branch(  "BdPVerry_refit" ,&BdPVerry_refit_ ,"BdPVerry_refit/D" );
  bsTree_->Branch(  "BdPVerrz_refit" , &BdPVerrz_refit_ , "BdPVerrz_refit/D" );

  bsTree_->Branch(  "BdLxy"			  , &BdLxy_,                            "BdLxy/D");
  bsTree_->Branch(  "BdLxyErr"			  , &BdLxyErr_,                            "BdLxyErr/D");
  bsTree_->Branch(  "BdErrX"			  , &BdErrX_,                           "BdErrX/D");
  bsTree_->Branch(  "BdErrY"			  , &BdErrY_,                           "BdErrY/D");
  bsTree_->Branch(  "BdErrXY"			  , &BdErrXY_,                          "BdErrXY/D");
  bsTree_->Branch(  "BdCt"			  , &BdCt_,                             "BdCt/D");
  bsTree_->Branch(  "BdCtErr"			  , &BdCtErr_,                           "BdCtErr/D");
  bsTree_->Branch(  "BdDist3d"			  , &BdDist3d_,                         "BdDist3d/D");
  bsTree_->Branch(  "BdDist3dErr"			  , &BdDist3dErr_,                      "BdDist3dErr/D");
  bsTree_->Branch(  "BdTime3d"			  , &BdTime3d_,                         "BdTime3d/D");
  bsTree_->Branch(  "BdTime3dErr"			  , &BdTime3dErr_,                      "BdTime3dErr/D");
  bsTree_->Branch(  "BdDist2d"			  , &BdDist2d_,                         "BdDist2d/D");
  bsTree_->Branch(  "BdDist2dErr"			  , &BdDist2dErr_,                      "BdDist2dErr/D");
  bsTree_->Branch(  "BdTime2d"			  , &BdTime2d_,                         "BdTime2d/D");
  bsTree_->Branch(  "BdTime2dErr"                   , &BdTime2dErr_,                      "BdTime2dErr/D");
  bsTree_->Branch("BdCt2DPVCosTheta", &BdCt2DPVCosTheta_ ,"BdCt2DPVCosTheta/D");
  bsTree_->Branch("Bdt2DPVCosTheta", &Bdt2DPVCosTheta_ ,"Bdt2DPVCosTheta/D");
  bsTree_->Branch("BdCtErr2DCostheta", &BdCtErr2DCostheta_ ,"BdCtErr2DCostheta/D");
  bsTree_->Branch("BdtErr2DCostheta", &BdtErr2DCostheta_ ,"BdtErr2DCostheta/D");

  bsTree_->Branch(  "BdSoftMuon1"        , &BdSoftMuon1_   ,"BdSoftMuon1/I" );
  bsTree_->Branch(  "BdSoftMuon2"        , &BdSoftMuon2_  ,"BdSoftMuon2/I");

  bsTree_->Branch(  "BdMCKstarPion"        , &BdMCKstarPion_   ,"BdMCKstarPion/I" );
  bsTree_->Branch(  "BdMCKstarKaon"        , &BdMCKstarKaon_  ,"BdMCKstarKaon/I");

  bsTree_->Branch(  "BdK1mcId"     , &BdK1mcId_     ,  "BdK1mcId/I"    );
  bsTree_->Branch(  "BdK1momId"	  ,	 &BdK1momId_    ,  "BdK1momId/I"   );
  bsTree_->Branch(  "BdK1gmomId"	  ,	 &BdK1gmomId_   ,  "BdK1gmomId/I"  );
  bsTree_->Branch(  "BdK2mcId"	  ,	 &BdK2mcId_     ,  "BdK2mcId/I"    );
  bsTree_->Branch(  "BdK2momId"	  ,	 &BdK2momId_    ,  "BdK2momId/I"   );
  bsTree_->Branch(  "BdK2gmomId"	  ,	 &BdK2gmomId_   ,  "BdK2gmomId/I"  );
  bsTree_->Branch(  "BdMu1mcId"	  ,	 &BdMu1mcId_    ,  "BdMu1mcId/I"   );
  bsTree_->Branch(  "BdMu1momId"	  ,	 &BdMu1momId_   ,  "BdMu1momId/I"  );
  bsTree_->Branch(  "BdMu1gmomId"  ,	 &BdMu1gmomId_  ,  "BdMu1gmomId/I" );
  bsTree_->Branch(  "BdMu2mcId"	  ,	 &BdMu2mcId_    ,  "BdMu2mcId/I"   );
  bsTree_->Branch(  "BdMu2momId"	  ,	 &BdMu2momId_   ,  "BdMu2momId/I"  );
  bsTree_->Branch(  "BdMu2gmomId"  ,	 &BdMu2gmomId_  ,  "BdMu2gmomId/I" );

}

BsToJpsiPhiRootTree::~BsToJpsiPhiRootTree()
{}

void BsToJpsiPhiRootTree::writeFile()
{
  bsFile_->Write();
  bsFile_->Close();

}

void BsToJpsiPhiRootTree::resetEntries()
{

//   PVTrkCharge_->clear();
//   PVTrkPt_->clear();
//   PVTrkEta_->clear();
//   PVTrkPhi_->clear();
//   BpPVTrkCharge_->clear();
//   BpPVTrkPt_->clear();
//   BpPVTrkEta_->clear();
//   BpPVTrkPhi_->clear();
//   BJetTrkCharge_->clear();
//   BJetTrkPt_->clear();
//   BpBJetTrkCharge_->clear();
//   BpBJetTrkPt_->clear();

  TrackMultiplicityBd_ = -9999999;
  TrackMultiplicityBp_ = -9999999;
  TrackMultiplicity_ = -9999999;
  MuonMultiplicity_ = -9999999;
  ElectronMultiplicity_ = -9999999;

  BJetParton_ = -9999999;
  BJetEta_= -9999999;
  BJetPhi_ = -9999999;
  BJetPt_ = -9999999;
  JetBTagProb_ = -9999999;

  BpBJetParton_ = -9999999;
  BpBJetEta_= -9999999;
  BpBJetPhi_ = -9999999;
  BpBJetPt_ = -9999999;
  BpJetBTagProb_ = -9999999;


  BplusCharge_=-9999999;
  BplusMu1Eta_=-9999999;
  BplusMu2Eta_=-9999999;
  BpSoftMuon1_=-9999999;
  BpSoftMuon2_=-9999999;
  JpsiPtbplus_fit_=-9999999;
  BdMuonCat1_=-9999999;
  BdMuonCat2_=-9999999;
  BplusM_fit_=-9999999;
  BplusVtxProb_=-9999999;
  BplusChi2_=-9999999;
  BplusPt_=-9999999;
  BplusPtot_=-9999999;
  KplusPt_=-9999999;
  KplusHighPurityTrack_=-9999999;
  KplusPtot_=-9999999;
  BplusMu1Pt_=-9999999;
  BplusMu2Pt_=-9999999;
  BpJpsiVtxProb_=-9999999;
  BpCosDeltaAlpha_=-9999999;
  BpMuMuDCA_=-9999999;
  BpMuMuDistance_=-9999999;
  BpMuMuDistanceSigma_=-9999999;
  BpMuDr1_=-9999999;
  BpMuDr2_=-9999999;
  BpMuonCat1_=-9999999;
  BpMuonCat2_=-9999999;
  BplusMu1Ptot_=-9999999;
  BplusMu2Ptot_=-9999999;
  BplusEta_=-9999999;
  BplusPhi_=-9999999;
  JpsiMass_bplus_=-9999999;
  JpsiPt_bplus_=-9999999;
  BplusPVindex_=-9999999;
  BpCt2DPVCosTheta_=-9999999;
  BpCt3DPVCosTheta_=-9999999;
  BpCt2D_=-9999999;
  BpCt3D_=-9999999;
  BpCt2DBS_=-9999999;
  BpCt3DPVClosestZ_=-9999999;
  BpCt2DPVClosestZ_=-9999999;
  Mu1SoftID_=-9999999;
  Mu2SoftID_=-9999999;
  triggerbit_HLTmu4TkTkDis_ = -9999999;
  triggerbit_HLTmu4TkDis_ = -9999999;
  triggerbit_HLTDimuon0JpsiMuon_ = -9999999;
  //Dimuon0JpsiMuon
  //triggerbit_HLTmu3Tk_ = -9999999;
  //triggerbit_HLTmu5_ = -9999999;
  BpPx_=-9999999;
  BpPy_=-9999999;
  BpPz_=-9999999;

  BpSVx_=-9999999;
  BpSVy_=-9999999;
  BpSVz_=-9999999;

  BpPVx_refit_=-9999999;
  BpPVy_refit_=-9999999;
  BpPVz_refit_=-9999999;

  BpPVx_refit_closestZ_=-9999999;
  BpPVy_refit_closestZ_=-9999999;
  BpPVz_refit_closestZ_=-9999999;

  BpPVx_refit_cosTheta_=-9999999;
  BpPVy_refit_cosTheta_=-9999999;
  BpPVz_refit_cosTheta_=-9999999;

  BpPVCL_ClosestZ_=-9999999;
  BpPVCL_CosTheta_=-9999999;
  BpPVCL_=-9999999;
  BpJpsiVtxCL_=-9999999;
  BpVtxCL_=-9999999;
  BpdRKaonJpsi_ = -9999999;
  BpCtErr2DCostheta_=-9999999;
  BpCtErr2DBS_=-9999999;
  BpCtErr2D_=-9999999;
  BpCtErr2DClosestZ_=-9999999;
  
  IP3DKandJpsiVtx_=-9999999;
  IP3DKandJpsiVtxErr_=-9999999;
  BpmatchDoubleMu01_=-9999999;
  BpmatchDoubleMu02_=-9999999;
  BpmatchDoubleMu41_=-9999999;
  BpmatchDoubleMu42_=-9999999;
  BpmatchDoubleMu41_v12_=-9999999;
  BpmatchDoubleMu42_v12_=-9999999;
  BpmatchDoubleMu01DiMuon0_=-9999999;
  BpmatchDoubleMu02DiMuon0_=-9999999;
  BplusKmcId_=-9999999;
  BplusKmomId_=-9999999;
  BplusMu1mcId_=-9999999;
  BplusMu1momId_=-9999999;
  BplusMu1gmomId_=-9999999;
  BplusMu2mcId_=-9999999;
  BplusMu2momId_=-9999999;
  BplusMu2gmomId_=-9999999;
  isMatchedBplus_=-9999999;


  JpsiMu2mcId_=-9999999;
  JpsiMu1mcId_=-9999999;
  JpsiMu2momId_=-9999999;
  JpsiMu1momId_=-9999999;
  isMatchedJpsi_=-9999999;

  BplusDecayChannel_=-9999999;

  runNumber_ =  -9999999;
  eventNumber_ =  -9999999;
  lumiSection_ =  -9999999;

  BcP_ =  -9999999;
  BcIP3D_ =  -9999999;
  BcCosAlpha_ =  -9999999;
  BcCt_ =  -9999999;
  BcM_ =  -9999999;
  BcMC_ =  -9999999;
  BcAng_ =  -9999999;
  BcAngPi_ =  -9999999;
  BcProb_ =  -9999999;
  BctrackPt_ = -9999999;

  PUinteraction_ = -9999999;

  ihaveajpsi_=  -9999999;
  BsCowboy_=  -9999999;
  BdCowboy_=  -9999999;
  BsPhiVtxProb_=  -9999999;
  BsMu1QualityG_=  -9999999;
  BsMu2QualityG_=  -9999999;
  BsMu1QualityT_=  -9999999;
  BsMu2QualityT_=  -9999999;
  BdMu1QualityG_=  -9999999;
  BdMu2QualityG_=  -9999999;
  BdMu1QualityT_=  -9999999;
  BdMu2QualityT_=  -9999999;


  NVertices_ = -9999999;
 

 JpsiMuonMatch41_= -9999999;
 JpsiMuonMatch42_= -9999999;

 JpsiCosDeltaAlpha_ = -9999999;
 JpsiLxySigma_ = -9999999;
 JpsiLxy_ = -9999999;
 JpsiGenLxy_ = -9999999;
 JpsiGenLxyOld_ = -9999999;
 JpsiGenLxyOverPt_= -9999999;
 JpsiLxyOverPt_= -9999999;
 JpsiGenPt_ = -9999999;

 JpsiCosAlphaMC_ = -9999999;
 JpsiMu1EtaMC_= -9999999;
 JpsiMu2EtaMC_= -9999999;

 JpsiMu2PtMC_= -9999999;
 JpsiMu1PtMC_= -9999999;

 JpsiMu2PhiMC_= -9999999;
 JpsiMu1PhiMC_= -9999999;

 JpsiGenPVz_ = -9999999;
 JpsiGenPVy_ = -9999999;
 JpsiGenPVx_ = -9999999;

 deltaRkaon1Jpsi_ = -9999999;
 deltaRkaon2Jpsi_ = -9999999;
 deltaRBCandPion_ = -9999999;
 JpsiTrigVtxProb_ = -9999999;
 JpsiMu1TrackerMuonArbitrated_  = -9999999;
 JpsiMu2TrackerMuonArbitrated_  = -9999999;

  for(int i = 0; i<30; i++){
    SVZpos_[i] = -9999999;
    PVZpos_[i] = -9999999;
    PVAbsPt_[i] = -9999999;
    NTracksInPV_[i] = -9999999;
  }

  BsLxy3DMC_ = -9999999;
  BsLxy2DMC_ = -9999999;
  BsPMC_ = -9999999;


  BsMassMC_ = -9999999;
  BsCt2DMC_TrueMass_ = -9999999;
  BsCosAlphaMC_  = -9999999;
  BsJpsiMu2PtMC_  = -9999999;
  BsJpsiMu1PtMC_  = -9999999;
  BsJpsiMu1EtaMC_  = -9999999;
  BsJpsiMu2EtaMC_  = -9999999;
  BsJpsiMu1PhiMC_  = -9999999;
  BsJpsiMu2PhiMC_  = -9999999;
  BsJpsiPtMC_  = -9999999;
  BsJpsiMassMC_  = -9999999;
  BsPhiPtMC_  = -9999999;
  BsPhiMassMC_  = -9999999;
  BsPhiK1PtMC_  = -9999999;
  BsPhiK2PtMC_  = -9999999;

  BSx_ = -9999999;
  BSy_ = -9999999;
  BSz_ = -9999999;
  BSdx_ = -9999999;
  BSdy_ = -9999999;
  BSdz_ = -9999999;
  BSsigmaZ_ = -9999999;
  BSdsigmaZ_ = -9999999;

  BSdxdz_ = -9999999;
  BSdydz_ = -9999999;




  PVx_ = -9999999;
  PVy_ = -9999999;
  PVz_ = -9999999;
  PVerrx_ = -9999999;
  PVerry_ = -9999999;
  PVerrz_ = -9999999;


  
  JpsiMuMuDCA_beffit_ =  -9999999;
  MuMuDCA_ = -9999999;
  KaonsDCA_= -9999999;
  MuMuDistance_ = -9999999;
  MuMuDistanceSigma_ = -9999999;
  MuDr1_ = -9999999;
  MuDr2_ = -9999999;
  BdMuMuDCA_ = -9999999;
  BdMuMuDistance_ = -9999999;
  BdMuMuDistanceSigma_ = -9999999;
  BdMuDr1_ = -9999999;
  BdMuDr2_ = -9999999;

  PVx_refit_= -9999999;
  PVy_refit_= -9999999;
  PVz_refit_= -9999999;
  PVerrx_refit_= -9999999;
  PVerry_refit_= -9999999;
  PVerrz_refit_= -9999999;

  PVx_refit_closestZ_ = -9999999;
  PVy_refit_closestZ_ = -9999999;
  PVz_refit_closestZ_ = -9999999;

  PVx_refit_cosTheta_ = -9999999;
  PVy_refit_cosTheta_ = -9999999;
  PVz_refit_cosTheta_ = -9999999;

	BsSVx_ = -9999999;
	BsSVy_ = -9999999;
	BsSVz_ = -9999999;

	JpsiSVx_ = -9999999;
	JpsiSVy_ = -9999999;
	JpsiSVz_ = -9999999;

  

  BdMu1Pt_beffit_= -9999999;
  BdMu1Pz_beffit_= -9999999;
  BdMu1Eta_beffit_= -9999999;
  BdMu1Phi_beffit_= -9999999;
  BdMu2Pt_beffit_= -9999999;
  BdMu2Pz_beffit_= -9999999;
  BdMu2Eta_beffit_= -9999999;
  BdMu2Phi_beffit_= -9999999;
  BdJpsiM_nofit_= -9999999;
  BdJpsiEta_nofit_= -9999999;
  BdJpsiPhi_nofit_= -9999999;
  BdJpsiPt_nofit_= -9999999;
  BdJpsiPz_nofit_= -9999999;
  BsSoftMuon1_ = -9999999;
  BsSoftMuon2_ = -9999999;

  JpsiVtxProb_ = -9999999;
  BdJpsiVtxProb_ = -9999999;
  JpsiM_alone_ = -9999999;
  JpsiPhi_alone_ = -9999999;
  JpsiEta_alone_ = -9999999;
  JpsiPt_alone_ = -9999999;
  JpsiMu1Pt_alone_ = -9999999;
  JpsiMu2Pt_alone_ = -9999999;
  JpsiMu1Phi_alone_ = -9999999;
  JpsiMu2Phi_alone_ = -9999999;
  JpsiMu1Eta_alone_ = -9999999;
  JpsiMu2Eta_alone_ = -9999999;
  MuonCat1_ = -9999999;
  MuonCat2_ = -9999999;
  JpsiMuonCat1_ = -9999999;
  JpsiMuonCat2_ = -9999999;

  JpsiMuon1Cat_alone_ = -9999999;
  JpsiMuon2Cat_alone_ = -9999999;
  BdJpsiMuon1Cat_alone_ = -9999999;
  BdJpsiMuon2Cat_alone_ = -9999999;

  Mu1TrkBSDxy_= -9999999;
  Mu1TrkBSDz_= -9999999;
  Mu1PixelHits_ = -9999999;
  Mu1TrackerHits_ = -9999999;
  Mu1isGood_ = -9999999;
  Mu1InnerTrkHighQuality_ = -9999999;

  Mu2TrkBSDxy_= -9999999;
  Mu2TrkBSDz_= -9999999;
  Mu2PixelHits_ = -9999999;
  Mu2TrackerHits_ = -9999999;
  Mu2isGood_ = -9999999;
  Mu2InnerTrkHighQuality_ = -9999999;

  BpMu1TrkBSDxy_= -9999999;
  BpMu1TrkBSDz_= -9999999;
  BpMu1PixelHits_ = -9999999;
  BpMu1TrackerHits_ = -9999999;
  BpMu1isGood_ = -9999999;
  BpMu1InnerTrkHighQuality_ = -9999999;

  BpMu2TrkBSDxy_= -9999999;
  BpMu2TrkBSDz_= -9999999;
  BpMu2PixelHits_ = -9999999;
  BpMu2TrackerHits_ = -9999999;
  BpMu2isGood_ = -9999999;
  BpMu2InnerTrkHighQuality_ = -9999999;



  BpMuonPairDR_= -9999999;
  MuonPairDR_= -9999999;
  JpsiMu1d0_alone_ = -9999999;
  JpsiMu2d0_alone_ = -9999999;
  JpsiMu1dz_alone_ = -9999999;
  JpsiMu2dz_alone_ = -9999999;
  JpsiMu1chi2_alone_ = -9999999;
  JpsiMu2chi2_alone_ = -9999999;
  JpsiMu1ndof_alone_ = -9999999;
  JpsiMu2ndof_alone_ = -9999999;
  JpsiMu1nHits_alone_ = -9999999;
  JpsiMu2nHits_alone_ = -9999999;
  JpsiMu1nPixHits_alone_ = -9999999;
  JpsiMu2nPixHits_alone_ = -9999999;

  K1Pt_beffit_ = -9999999;
  K1Pz_beffit_ = -9999999;
  K1Eta_beffit_ = -9999999;
  K1Phi_beffit_ = -9999999;
  K2Pt_beffit_ = -9999999;
  K2Pz_beffit_ = -9999999;
  K2Eta_beffit_ = -9999999;
  K2Phi_beffit_ = -9999999;

  Mu1Pt_beffit_ = -9999999;
  Mu1Pz_beffit_ = -9999999;
  Mu1Eta_beffit_ = -9999999;
  Mu1Phi_beffit_ = -9999999;
  Mu2Pt_beffit_ = -9999999;
  Mu2Pz_beffit_ = -9999999;
  Mu2Eta_beffit_ = -9999999;
  Mu2Phi_beffit_ = -9999999;

  BsFitChi2_ = -9999999;
  BsFitNdof_ = -9999999;
  BsFitVtxProb_ = -9999999;
  JpsiNumberOfCandidates_ = 0;
  JpsiGenNumberOfCandidates_=0;
  PhiNumberOfCandidatesBeforeFit_ =  0;
  BsNumberOfCandidatesBeforeFit_ =  0;
  BsNumberOfCandidatesAfterFit_ =  0;
  BsNumberOfCandidatesAfterBestFit_ =  0;
  BsFitM_ = -9999999;
  K1Pt_fit_ = -9999999;
  K2Pt_fit_ = -9999999;
  PhiM_fit_ = -9999999;
  BsFitEta_ = -9999999;
  BsFitPt_ = -9999999;
  BsFitPz_ = -9999999;
  BsFitPhi_ = -9999999;
  BsFitVtx_x_ = -9999999;
  BsFitVtx_y_ = -9999999;
  BsFitVtx_z_ = -9999999;
  BsM_nofit_ = -9999999;
  BsPhi_nofit_ = -9999999;
  BsEta_nofit_ = -9999999;
  BsPt_nofit_ = -9999999;
  BsPz_nofit_ = -9999999;
  JpsiM_nofit_ = -9999999;
  JpsiPhi_nofit_ = -9999999;
  JpsiEta_nofit_ = -9999999;
  JpsiPt_nofit_ = -9999999;
  JpsiPz_nofit_ = -9999999;
  PhiM_nofit_ = -9999999;
  PhiPhi_nofit_ = -9999999;
  PhiEta_nofit_ = -9999999;
  PhiPt_nofit_ = -9999999;
  PhiPz_nofit_ = -9999999;
  K1Pt_nofit_ = -9999999;
  K1Pz_nofit_ = -9999999;
  K1HighPurityTrack_ = -9999999;
  K2HighPurityTrack_ = -9999999;

  K1Eta_nofit_ = -9999999;
  K1Phi_nofit_ = -9999999;
  K1Key_nofit_ = -9999999;
  K2Eta_nofit_ = -9999999;
  K2Pt_nofit_ = -9999999;
  K2Pz_nofit_ = -9999999;
  K2Phi_nofit_ = -9999999;
  K2Key_nofit_ = -9999999;
  K1Chi2_ = -9999999;
  K1nHits_ = -9999999;
  K2Chi2_ = -9999999;
  K2nHits_ = -9999999;
  K1pixH_ = -9999999;
  K1trkH_ = -9999999;
  K2pixH_ = -9999999;
  K2trkH_ = -9999999;
  Mu1Chi2_ = -9999999;
  Mu1nHits_ = -9999999;
  Mu2Chi2_ = -9999999;
  Mu2nHits_ = -9999999;
  Mu1pixH_ = -9999999;
  Mu1trkH_ = -9999999;
  Mu2pixH_ = -9999999;
  Mu2trkH_ = -9999999;
  costheta_ = -9999999;
  phi_ = -9999999;
  cospsi_ = -9999999;
  Bdcostheta_ = -9999999;
  Bdphi_ = -9999999;
  Bdcospsi_ = -9999999;
  BdcosthetaMC_ = -9999999;
  BdphiMC_ = -9999999;
  BdcospsiMC_ = -9999999;
  AngleBsDecayLength_ = -9999999;
  CosDeltaAlpha_ = -9999999;
  BdCosDeltaAlpha_ = -9999999;

  isPV_ = -9999999;
  isBS_ = -9999999;

  isMatched_ = -9999999;
  isMatchedBd_ = -9999999;
  BsLxy_ = -9999999;
  BsLxyErr_ = -9999999;
  BsErrX_ = -9999999;
  BsErrY_ = -9999999;
  BsErrXY_ = -9999999;
  BsCt_ = -9999999;
  BsCtErr_ = -9999999;
  BsCt3D_ = -9999999;
  BsCt2D_ = -9999999;
  BsCt2DBS_ = -9999999;
  BdCt2DBS_ = -9999999;
  BdCt2DMC_ = -9999999;
  BdCt3DMC_ = -9999999;
  BsCtMPV_ = -9999999;
  BsCt3Drefit_ = -9999999;
  BsCt2Drefit_ = -9999999;
  BsCtMPVrefit_ = -9999999;
  BsCtErr3D_ = -9999999;
  BsCtErr2D_ = -9999999;
  BsCtErr2DBS_ = -9999999;
  BsCtErr2DClosestZ_= -9999999;
  BdCtErr2DBS_ = -9999999;
  BsCtErr2D2_ = -9999999;
  BsCtErrMPV_ = -9999999;
  BsCtErr3Drefit_ = -9999999;
  BsCtErr2Drefit_ = -9999999;
  BsCtErrMPVrefit_ = -9999999;

  BsCtErr2DBSOld_ = -9999999;
  BsCt2DBSOld_= -9999999;

  BsCtErr2DBS_JpsiVtx_ = -9999999;
  BsCtErr2D_JpsiVtx_ = -9999999;
  BsCtErr2DClosestZ_JpsiVtx_ = -9999999;
  BsCtErr2DCostheta_JpsiVtx_ = -9999999;



  BsCtErr2DClosestZOld_ = -9999999;
  BsCt2DPVClosestZOld_= -9999999;

  BsCtErr2DOld_= -9999999;
  BsCt2DOld_= -9999999;

  K1trkLay_ = -9999999;
  K1muDTh_ = -9999999;
  K1muCSCh_ = -9999999;
  K1muRPCh_ = -9999999;
  K2trkLay_ = -9999999;
  K2muDTh_ = -9999999;
  K2muCSCh_ = -9999999;
  K2muRPCh_ = -9999999;
  Mu1trkLay_ = -9999999;
  Mu1muDTh_ = -9999999;
  Mu1muCSCh_ = -9999999;
  Mu1muRPCh_ = -9999999;
  Mu2trkLay_ = -9999999;
  Mu2muDTh_ = -9999999;
  Mu2muCSCh_ = -9999999;
  Mu2muRPCh_ = -9999999;
  K1mcId_ = -9999999;
  K1momId_ = -9999999;
  K1gmomId_ = -9999999;
  K2mcId_ = -9999999;
  K2momId_ = -9999999;
  K2gmomId_ = -9999999;
  Mu1mcId_ = -9999999;
  Mu1momId_ = -9999999;
  Mu1gmomId_ = -9999999;
  Mu2mcId_ = -9999999;
  Mu2momId_ = -9999999;
  Mu2gmomId_ = -9999999;

  Mu1d0_ = -9999999;
  Mu2d0_ = -9999999;
  Mu1dz_ = -9999999;
  Mu2dz_ = -9999999;
  Mu1GlobalMuonPromptTight_=-9999999;
  Mu2GlobalMuonPromptTight_=-9999999;
  Mu1TrackerMuonArbitrated_=-9999999;
  Mu1TMLastStationTight_=-9999999;
  Mu1TMOneStationTight_=-9999999;
  Mu1TMLastStationOptimizedLowPtTight_=-9999999;
  Mu1TMLastStationAngTight_=-9999999;
  Mu1TMOneStationAngTight_=-9999999;
  Mu1TMLastStationOptimizedBarrelLowPtTight_=-9999999;
  Mu2TrackerMuonArbitrated_=-9999999;
  Mu2TMLastStationTight_=-9999999;
  Mu2TMOneStationTight_=-9999999;
  Mu2TMLastStationOptimizedLowPtTight_=-9999999;
  Mu2TMLastStationAngTight_=-9999999;
  Mu2TMOneStationAngTight_=-9999999;
  Mu2TMLastStationOptimizedBarrelLowPtTight_=-9999999;

  BsCtErr2DCosthetaOld_ = -9999999;
  BsCtErr2DCostheta_ = -9999999;
  BsCt3DPVCosTheta_ = -9999999;
  BsCt2DPVCosTheta_ =-9999999;
  BsCt2DPVCosThetaOld_ =-9999999;
  BsCt2DPVClosestZ_=-9999999;
  BsCt3DPVClosestZ_=-9999999;
 
  BsDist3d_ = -9999999;
  BsDist3dErr_ = -9999999;
  BsTime3d_ = -9999999;
  BsTime3dErr_ = -9999999;
  BsDist2d_ = -9999999;
  BsDist2dErr_ = -9999999;
  BsTime2d_ = -9999999;
  BsTime2dErr_ = -9999999;
  dedxTrk_ = -9999999;
  errdedxTrk_ = -9999999;
  numdedxTrk_ = -9999999;
  iPassedCutIdent_ = -9999999;
  iPassedCutIdentBd_ = -9999999;
  BdTrack1Charge_ = -9999999;

  K1Fit_sigX_ = -9999999;
  K1Fit_sigY_ = -9999999;
  K1Fit_sigZ_ = -9999999;
  K2Fit_sigX_ = -9999999;
  K2Fit_sigY_ = -9999999;
  K2Fit_sigZ_ = -9999999;
  K1Fit_sigPX_ = -9999999;
  K1Fit_sigPY_ = -9999999;
  K1Fit_sigPZ_ = -9999999;
  K2Fit_sigPX_ = -9999999;
  K2Fit_sigPY_ = -9999999;
  K2Fit_sigPZ_ = -9999999;

  GenNumberOfBdecays_ = -9999999;

  genBpPV_z_ = -9999999;
  genBpPV_y_ = -9999999;
  genBpPV_x_ = -9999999;  

  genBsPV_z_ = -9999999;
  genBsPV_y_ = -9999999;
  genBsPV_x_ = -9999999;  

  genBsVtx_z_ = -9999999;
  genBsVtx_y_ = -9999999;
  genBsVtx_x_ = -9999999;
  genBsSVtx_z_ = -9999999;
  genBsSVtx_y_ = -9999999;
  genBsSVtx_x_ = -9999999;
  isGenJpsiEvent_ = -9999999;
  BdFitChi2_Hyp1_ = -9999999;
  BdFitNdof_Hyp1_ = -9999999;
  BdFitVtxProb_Hyp1_ = -9999999;
  BdFitVtxProb_ = -9999999;
  BdFitM_Hyp1_ = -9999999;
  BdFitEta_Hyp1_ = -9999999;
  BdFitPt_Hyp1_ = -9999999;
  BdFitPz_Hyp1_ = -9999999;
  BdFitPhi_Hyp1_ = -9999999;
  BdFitVtx_x_Hyp1_ = -9999999;
  BdFitVtx_y_Hyp1_ = -9999999;
  BdFitVtx_z_Hyp1_ = -9999999;
  BdM_nofit_ = -9999999;
  BdPhi_nofit_ = -9999999;
  BdEta_nofit_ = -9999999;
  BdPt_nofit_ = -9999999;
  BdPz_nofit_ = -9999999;
  KstarMass_nofit_Hyp1_ = -9999999;
  KstarMass_nofit_Hyp2_ = -9999999;
  BdK1_kpi_sigX_Hyp1_ = -9999999;
  BdK1_kpi_sigY_Hyp1_ = -9999999;
  BdK1_kpi_sigZ_Hyp1_ = -9999999;
  BdK2_kpi_sigX_Hyp1_ = -9999999;
  BdK2_kpi_sigY_Hyp1_ = -9999999;
  BdK2_kpi_sigZ_Hyp1_ = -9999999;
  BdK1_kpi_sigPX_Hyp1_ = -9999999;
  BdK1_kpi_sigPY_Hyp1_ = -9999999;
  BdK1_kpi_sigPZ_Hyp1_ = -9999999;
  BdK2_kpi_sigPX_Hyp1_ = -9999999;
  BdK2_kpi_sigPY_Hyp1_ = -9999999;
  BdK2_kpi_sigPZ_Hyp1_ = -9999999;
  BdFitChi2_Hyp2_ = -9999999;
  BdFitNdof_Hyp2_ = -9999999;
  BdFitVtxProb_Hyp2_ = -9999999;
  BdFitM_Hyp2_ = -9999999;
  BdFitEta_Hyp2_ = -9999999;
  BdFitPt_Hyp2_ = -9999999;
  BdFitPz_Hyp2_ = -9999999;
  BdFitPhi_Hyp2_ = -9999999;
  BdFitVtx_x_Hyp2_ = -9999999;
  BdFitVtx_y_Hyp2_ = -9999999;
  BdFitVtx_z_Hyp2_ = -9999999;
  BdNumberOfCandidates_ =  0;

  BdPVx_refit_    = -9999999;
  BdPVy_refit_    = -9999999;
  BdPVz_refit_    = -9999999;

  BdPVerrx_refit_ = -9999999;
  BdPVerry_refit_ = -9999999;
  BdPVerrz_refit_ = -9999999;

  BdK1_kpi_sigX_Hyp2_ = -9999999;
  BdK1_kpi_sigY_Hyp2_ = -9999999;
  BdK1_kpi_sigZ_Hyp2_ = -9999999;
  BdK2_kpi_sigX_Hyp2_ = -9999999;
  BdK2_kpi_sigY_Hyp2_ = -9999999;
  BdK2_kpi_sigZ_Hyp2_ = -9999999;
  BdK1_kpi_sigPX_Hyp2_ = -9999999;
  BdK1_kpi_sigPY_Hyp2_ = -9999999;
  BdK1_kpi_sigPZ_Hyp2_ = -9999999;
  BdK2_kpi_sigPX_Hyp2_ = -9999999;
  BdK2_kpi_sigPY_Hyp2_ = -9999999;
  BdK2_kpi_sigPZ_Hyp2_ = -9999999;
  BdK1Pt_nofit_ = -9999999;
  BdK1Pz_nofit_ = -9999999;
  BdK1Eta_nofit_ = -9999999;
  BdK1Phi_nofit_ = -9999999;
  BdK1Key_nofit_ = -9999999;
  BdK2Pt_nofit_ = -9999999;
  BdK2Pz_nofit_ = -9999999;
  BdK2Eta_nofit_ = -9999999;
  BdK2Phi_nofit_ = -9999999;
  BdK2Key_nofit_ = -9999999;
  BdLxy_ = -9999999;
  BdLxyErr_ = -9999999;
  BdErrX_ = -9999999;
  BdErrY_ = -9999999;
  BdErrXY_ = -9999999;
  BdCt_ = -9999999;
  BdCtErr_ = -9999999;
  BdDist3d_ = -9999999;
  BdDist3dErr_ = -9999999;
  BdTime3d_ = -9999999;
  BdTime3dErr_ = -9999999;
  BdDist2d_ = -9999999;
  BdDist2dErr_ = -9999999;
  BdTime2d_ = -9999999;
  BdTime2dErr_ = -9999999;
  BdCt2DPVCosTheta_ = -9999999;
  Bdt2DPVCosTheta_  = -9999999;
  BdCtErr2DCostheta_ = -9999999;
  BdtErr2DCostheta_ = -9999999;

  BdSoftMuon1_ =  -9999999;
  BdSoftMuon2_ =  -9999999;
  BdMCKstarKaon_ = -9999999;
  BdMCKstarPion_ = -9999999;  
  BdK1mcId_    = -9999999;
  BdK1momId_   = -9999999;
  BdK1gmomId_  = -9999999;
  BdK2mcId_    = -9999999;
  BdK2momId_   = -9999999;
  BdK2gmomId_  = -9999999;
  BdMu1mcId_   = -9999999;
  BdMu1momId_  = -9999999;
  BdMu1gmomId_ = -9999999;
  BdMu2mcId_   = -9999999;
  BdMu2momId_  = -9999999;
  BdMu2gmomId_ = -9999999;

  BsCt3DMC_ = -9999999;
  BsCt2DMC_ = -9999999;
  BpCt2DMC_ = -9999999;
  BsIniFlavour_ = -9999999;
  BsEndFlavour_ = -9999999;
  BdIniFlavour_ = -9999999;

  BdEndFlavour_ = -9999999;
  ChannelID_ = -9999999;
  BdChannelID_ = -9999999;
  BscosthetaMC_ = -9999999;
  BsphiMC_ = -9999999;
  BscospsiMC_ = -9999999;

  for(int i=0; i<7; i++){
    K1Fit_par_[i] = -9999999;
    K2Fit_par_[i] = -9999999;

    BdK1_kpi_par_Hyp1_[i] = -9999999;
    BdK2_kpi_par_Hyp1_[i] = -9999999;
    BdK1_kpi_par_Hyp2_[i] = -9999999;
    BdK2_kpi_par_Hyp2_[i] = -9999999;
  }




  for(int i=0; i<10; i++){
    BmesonsId_[i]  =  -9999999;

    costhetaMC_[i] = -9999999;
    phiMC_[i] = -9999999;
    cospsiMC_[i] = -9999999;

    BMMC_[i] =  -9999999;
    BPtMC_[i] =  -9999999;
    BPxMC_[i] =  -9999999;
    BPyMC_[i] =  -9999999;
    BPzMC_[i] =  -9999999;
    BEtaMC_[i] =  -9999999;
    BPhiMC_[i] =  -9999999;

    BVtxMC_x_[i] =  -9999999;
    BVtxMC_y_[i]  =  -9999999;
    BVtxMC_z_[i]  =  -9999999;
    BSVtxMC_x_[i] =  -9999999;
    BSVtxMC_y_[i] =  -9999999;
    BSVtxMC_z_[i] =  -9999999;
    BLxy_MC_[i]   =  -9999999;
    BCt_MC_[i]    =  -9999999;
    BCt_MC2D_[i]    =  -9999999;
    BCt_MC3D_[i]    =  -9999999;

    GenNumberOfDaughters_[i] =  -9999999;

    for(int j=0;j<15;j++){
      BDauIdMC_[i][j]= -9999999;

      BDauMMC_[i][j]= -9999999;
      BDauPtMC_[i][j]= -9999999;
      BDauPzMC_[i][j]= -9999999;
      BDauEtaMC_[i][j]= -9999999;
      BDauPhiMC_[i][j]= -9999999;

      GenNumberOfDaughtersDaughters_[i][j] =  -9999999;

      for(int k=0; k<10; k++){
        BDauDauIdMC_[i][j][k]= -9999999;

        BDauDauMMC_[i][j][k]= -9999999;
        BDauDauPtMC_[i][j][k]= -9999999;
        BDauDauPzMC_[i][j][k]= -9999999;
        BDauDauEtaMC_[i][j][k]= -9999999;
        BDauDauPhiMC_[i][j][k]= -9999999;
      }
    }
  }

}

void BsToJpsiPhiRootTree::getDeDx(const double f1, const double f2, const int f3)
{
  dedxTrk_ = f1;
  errdedxTrk_ = f2;
  numdedxTrk_ = f3;
}


void BsToJpsiPhiRootTree::getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff,
                                 const double gg, const double hh, const double ii)
{
  BSx_ = aa;
  BSy_ = bb;
  BSz_ = cc;
  PVx_ = dd;
  PVy_ = ee;
  PVz_ = ff;
  PVerrx_ = gg;
  PVerry_ = hh;
  PVerrz_ = ii;
}


void BsToJpsiPhiRootTree::getAngles(const double aa, const double bb, const double cc, const double dd)
{
  costheta_ = aa;
  phi_ = bb;
  cospsi_ = cc;
  AngleBsDecayLength_ = dd;
}


void BsToJpsiPhiRootTree::fill()
{
  bsTree_->Fill();
}


void BsToJpsiPhiRootTree::readTree(const std::string filename)
{

  // open root file
  bsFile_ = new TFile (filename.c_str(), "READ" );

  // create tree structure
  bsTree_ =  (TTree*) bsFile_->Get("BsTree");

  setBranchAddresses();
}

void BsToJpsiPhiRootTree::readTree(std::vector<std::string> filenames){
  TChain * myChain = new TChain("BsTree");
  //  for(int i=0;i<filenames.size();++i) {
    //    std::string filenome = filenames[i] ;
  //    myChain->Add(filenome,-1);
  //  }

  for(std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++){
    myChain->Add( (*it).c_str());
  }

  bsTree_ = myChain;
  setBranchAddresses();
}

void BsToJpsiPhiRootTree::setBranchAddresses(){

  bsTree_->SetBranchAddress( "BplusCharge", &BplusCharge_ );
  bsTree_->SetBranchAddress( "BplusMu1Eta", &BplusMu1Eta_ );
  bsTree_->SetBranchAddress( "BplusMu2Eta", &BplusMu2Eta_ );
  bsTree_->SetBranchAddress("BpSoftMuon1",&BpSoftMuon1_);
  bsTree_->SetBranchAddress("BpSoftMuon2",&BpSoftMuon2_);
  bsTree_->SetBranchAddress( "JpsiPtbplus_fit", &JpsiPtbplus_fit_ );

  bsTree_->SetBranchAddress( "BplusPt", &BplusPt_ );
  bsTree_->SetBranchAddress( "BplusPtot", &BplusPtot_ );
  bsTree_->SetBranchAddress( "KplusPt", &KplusPt_ );
  bsTree_->SetBranchAddress( "KplusPtot", &KplusPtot_ );
  bsTree_->SetBranchAddress( "KplusHighPurityTrack", &KplusHighPurityTrack_ );
  bsTree_->SetBranchAddress( "BplusMu1Pt", &BplusMu1Pt_ );
  bsTree_->SetBranchAddress( "BplusMu2Pt", &BplusMu2Pt_ );
  bsTree_->SetBranchAddress( "BpJpsiVtxProb", &BpJpsiVtxProb_ );
  bsTree_->SetBranchAddress( "BpCosDeltaAlpha", &BpCosDeltaAlpha_ );
  bsTree_->SetBranchAddress( "BpMuMuDCA", &BpMuMuDCA_ );
  bsTree_->SetBranchAddress( "BpMuMuDistance", &BpMuMuDistance_ );
  bsTree_->SetBranchAddress( "BpMuMuDistanceSigma", &BpMuMuDistanceSigma_ );
  bsTree_->SetBranchAddress( "BpMuDr1", &BpMuDr1_ );
  bsTree_->SetBranchAddress( "BpMuDr2", &BpMuDr2_ );
  bsTree_->SetBranchAddress( "BpMuonCat1", &BpMuonCat1_ );
  bsTree_->SetBranchAddress( "BpMuonCat2", &BpMuonCat2_ );
  bsTree_->SetBranchAddress( "BplusMu1Ptot", &BplusMu1Ptot_ );
  bsTree_->SetBranchAddress( "BplusMu2Ptot", &BplusMu2Ptot_ );
  bsTree_->SetBranchAddress( "BplusEta", &BplusEta_ );
  bsTree_->SetBranchAddress( "BplusPhi", &BplusPhi_ );
  bsTree_->SetBranchAddress( "JpsiMass_bplus", &JpsiMass_bplus_ );
  bsTree_->SetBranchAddress( "JpsiPt_bplus", &JpsiPt_bplus_ );
  bsTree_->SetBranchAddress( "BplusPVindex", &BplusPVindex_ );
  bsTree_->SetBranchAddress("BpCt2D",&BpCt2D_);
  bsTree_->SetBranchAddress("BpCt3D",&BpCt3D_);
  bsTree_->SetBranchAddress("BpCt2DBS",&BpCt2DBS_);
  bsTree_->SetBranchAddress( "BpCt3DPVCosTheta", &BpCt3DPVCosTheta_ );
  bsTree_->SetBranchAddress( "BpCt2DPVCosTheta", &BpCt2DPVCosTheta_ );
  bsTree_->SetBranchAddress( "BpCtErr2DCostheta", &BpCtErr2DCostheta_ );
  bsTree_->SetBranchAddress("BpCtErr2DBS",&BpCtErr2DBS_);
  bsTree_->SetBranchAddress("BpCt2DPVClosestZ",&BpCt2DPVClosestZ_);
  bsTree_->SetBranchAddress("BpCt3DPVClosestZ",&BpCt3DPVClosestZ_);
  bsTree_->SetBranchAddress("BpCtErr2D",&BpCtErr2D_);
  bsTree_->SetBranchAddress("BpCtErr2DClosestZ",&BpCtErr2DClosestZ_);
 

  bsTree_->SetBranchAddress("BpPx",&BpPx_);
  bsTree_->SetBranchAddress("BpPy",&BpPy_);
  bsTree_->SetBranchAddress("BpPz",&BpPz_);

  bsTree_->SetBranchAddress("BpSVx",&BpSVx_);
  bsTree_->SetBranchAddress("BpSVy",&BpSVy_);
  bsTree_->SetBranchAddress("BpSVz",&BpSVz_);

  bsTree_->SetBranchAddress("BpPVx_refit",&BpPVx_refit_);
  bsTree_->SetBranchAddress("BpPVy_refit",&BpPVy_refit_);
  bsTree_->SetBranchAddress("BpPVz_refit",&BpPVz_refit_);

  bsTree_->SetBranchAddress("BpPVx_refit_closestZ",&BpPVx_refit_closestZ_);
  bsTree_->SetBranchAddress("BpPVy_refit_closestZ",&BpPVy_refit_closestZ_);
  bsTree_->SetBranchAddress("BpPVz_refit_closestZ",&BpPVz_refit_closestZ_);

  bsTree_->SetBranchAddress("BpPVx_refit_cosTheta",&BpPVx_refit_cosTheta_);
  bsTree_->SetBranchAddress("BpPVy_refit_cosTheta",&BpPVy_refit_cosTheta_);
  bsTree_->SetBranchAddress("BpPVz_refit_cosTheta",&BpPVz_refit_cosTheta_);

  bsTree_->SetBranchAddress("BpPVCL_ClosestZ",&BpPVCL_ClosestZ_);
  bsTree_->SetBranchAddress("BpPVCL_CosTheta",&BpPVCL_CosTheta_);
  bsTree_->SetBranchAddress("BpPVCL",&BpPVCL_);
  bsTree_->SetBranchAddress("BpJpsiVtxCL",&BpJpsiVtxCL_);
  bsTree_->SetBranchAddress("BpVtxCL",&BpVtxCL_);
  bsTree_->SetBranchAddress("BpdRKaonJpsi",&BpdRKaonJpsi_);

  bsTree_->SetBranchAddress( "IP3DKandJpsiVtx", &IP3DKandJpsiVtx_ );
  bsTree_->SetBranchAddress( "IP3DKandJpsiVtxErr", &IP3DKandJpsiVtxErr_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu01", &BpmatchDoubleMu01_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu02", &BpmatchDoubleMu02_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu41", &BpmatchDoubleMu41_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu42", &BpmatchDoubleMu42_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu41_v12", &BpmatchDoubleMu41_v12_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu42_v12", &BpmatchDoubleMu42_v12_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu01DiMuon0", &BpmatchDoubleMu01DiMuon0_ );
  bsTree_->SetBranchAddress( "BpmatchDoubleMu02DiMuon0", &BpmatchDoubleMu02DiMuon0_ );
  bsTree_->SetBranchAddress( "BplusKmcId", &BplusKmcId_ );
  bsTree_->SetBranchAddress( "BplusKmomId", &BplusKmomId_ );
  bsTree_->SetBranchAddress( "BplusMu1mcId", &BplusMu1mcId_ );
  bsTree_->SetBranchAddress( "BplusMu1momId", &BplusMu1momId_ );
  bsTree_->SetBranchAddress( "BplusMu1gmomId", &BplusMu1gmomId_ );
  bsTree_->SetBranchAddress( "BplusMu2mcId", &BplusMu2mcId_ );
  bsTree_->SetBranchAddress( "BplusMu2momId", &BplusMu2momId_ );
  bsTree_->SetBranchAddress( "BplusMu2gmomId", &BplusMu2gmomId_ );
  bsTree_->SetBranchAddress( "isMatchedBplus", &isMatchedBplus_ );

  bsTree_->SetBranchAddress( "BplusDecayChannel", &BplusDecayChannel_ );

  bsTree_->SetBranchAddress(  "isMatchedJpsi"		  , &isMatchedJpsi_);
  bsTree_->SetBranchAddress(  "JpsiMu1mcId"			  , &JpsiMu1mcId_);
  bsTree_->SetBranchAddress(  "JpsiMu2mcId"			  , &JpsiMu2mcId_);
  bsTree_->SetBranchAddress(  "JpsiMu1momId"			  , &JpsiMu1momId_);
  bsTree_->SetBranchAddress(  "JpsiMu2momId"			  , &JpsiMu2momId_);

  bsTree_->SetBranchAddress(  "BcP"             , &BcP_);
  bsTree_->SetBranchAddress(  "BcCosAlpha"             , &BcCosAlpha_);
  bsTree_->SetBranchAddress(  "BcIP3D"             , &BcIP3D_);
  bsTree_->SetBranchAddress(  "BcCt"             , &BcCt_);
  bsTree_->SetBranchAddress(  "BcM"             , &BcM_);
  bsTree_->SetBranchAddress(  "BcMC"             , &BcMC_);
  bsTree_->SetBranchAddress(  "BcAng"             , &BcAng_);
  bsTree_->SetBranchAddress(  "BcAngPi"             , &BcAngPi_);
  bsTree_->SetBranchAddress(  "BcProb"             , &BcProb_);
  bsTree_->SetBranchAddress(  "BctrackPt"             , &BctrackPt_);

  bsTree_->SetBranchAddress("BsLxy3DMC", &BsLxy3DMC_);
  bsTree_->SetBranchAddress("BsLxy2DMC", &BsLxy2DMC_);
  bsTree_->SetBranchAddress("BsPMC", &BsPMC_);
  bsTree_->SetBranchAddress(  "PVZpos"            , &PVZpos_ );
  bsTree_->SetBranchAddress(  "NTracksInPV"       , &NTracksInPV_  );
  bsTree_->SetBranchAddress( "PVAbsPt"           , &PVAbsPt_);
  bsTree_->SetBranchAddress(  "SVZpos"            , &SVZpos_ );
  bsTree_->SetBranchAddress("BsPtMC", &BsPtMC_);
  bsTree_->SetBranchAddress("FirstBsMCZpos", &FirstBsMCZpos_);

  bsTree_->SetBranchAddress(  "runNumber"             , &runNumber_  );
  bsTree_->SetBranchAddress(  "eventNumber"             , &eventNumber_  );
  bsTree_->SetBranchAddress(  "lumiSection"             , &lumiSection_  );

  bsTree_->SetBranchAddress(  "PUinteraction"             , &PUinteraction_  );
  //----------------------------------------------------------------------------------------------------------------------muonsoftid
  bsTree_->SetBranchAddress( "Mu1SoftID", &Mu1SoftID_);
  bsTree_->SetBranchAddress( "Mu2SoftID", &Mu2SoftID_);
  bsTree_->SetBranchAddress( "deltaRkaon1Jpsi", &deltaRkaon1Jpsi_);
  bsTree_->SetBranchAddress( "deltaRkaon2Jpsi", &deltaRkaon2Jpsi_);
  bsTree_->SetBranchAddress( "deltaRBCandJpsi", &deltaRBCandPion_);
  //---------------------------------------------------------------------------------------------------trigger
  bsTree_->SetBranchAddress(  "triggerbit_HLTmu4TkDis"             , &triggerbit_HLTmu4TkTkDis_);
  bsTree_->SetBranchAddress(  "triggerbit_HLTmu4TkTkDis"		  , &triggerbit_HLTmu4TkTkDis_  );
  bsTree_->SetBranchAddress(  "triggerbit_HLTDimuon0JpsiMuon"                  , &triggerbit_HLTDimuon0JpsiMuon_  );
  //Dimuon0JpsiMuon
  bsTree_->SetBranchAddress(  "ihaveajpsi"             , &ihaveajpsi_);
  bsTree_->SetBranchAddress(  "BsCowboy"             , &BsCowboy_);
  bsTree_->SetBranchAddress(  "BdCowboy"             , &BdCowboy_);
  bsTree_->SetBranchAddress(  "BsPhiVtxProb"             , &BsPhiVtxProb_);
  bsTree_->SetBranchAddress(  "BsMu1QualityG"             , &BsMu1QualityG_);
  bsTree_->SetBranchAddress(  "BsMu2QualityG"             , &BsMu2QualityG_);
  bsTree_->SetBranchAddress(  "BsMu1QualityT"             , &BsMu1QualityT_);
  bsTree_->SetBranchAddress(  "BsMu2QualityT"             , &BsMu2QualityT_);
  bsTree_->SetBranchAddress(  "BdMu1QualityG"             , &BdMu1QualityG_);
  bsTree_->SetBranchAddress(  "BdMu2QualityG"             , &BdMu2QualityG_);
  bsTree_->SetBranchAddress(  "BdMu1QualityT"             , &BdMu1QualityT_);
  bsTree_->SetBranchAddress(  "BdMu2QualityT"             , &BdMu2QualityT_);

  bsTree_->SetBranchAddress(  "NVertices"             , &NVertices_  );


 bsTree_->SetBranchAddress(  "JpsiMuonMatch41"         , &JpsiMuonMatch41_);
 bsTree_->SetBranchAddress(  "JpsiMuonMatch42"         , &JpsiMuonMatch42_);

 bsTree_->SetBranchAddress(  "JpsiCosDeltaAlpha"         , &JpsiCosDeltaAlpha_);
 bsTree_->SetBranchAddress(  "JpsiLxySigma"         , &JpsiLxySigma_);
 bsTree_->SetBranchAddress(  "JpsiLxy"         , &JpsiLxy_);


 bsTree_->SetBranchAddress(  "JpsiGenPt"         , &JpsiGenPt_);
 bsTree_->SetBranchAddress(  "JpsiGenLxy"         , &JpsiGenLxy_);
 bsTree_->SetBranchAddress(  "JpsiGenLxyOld"         , &JpsiGenLxyOld_);
 bsTree_->SetBranchAddress(  "JpsiGenLxyOverPt"         , &JpsiGenLxyOverPt_);
 bsTree_->SetBranchAddress(  "JpsiLxyOverPt"         , &JpsiLxyOverPt_);

 bsTree_->SetBranchAddress(  "JpsiMu1EtaMC"         , &JpsiMu1EtaMC_);
 bsTree_->SetBranchAddress(  "JpsiMu2EtaMC"         , &JpsiMu2EtaMC_);

 bsTree_->SetBranchAddress(  "JpsiMu1PtMC"         , &JpsiMu1PtMC_);
 bsTree_->SetBranchAddress(  "JpsiMu2PtMC"         , &JpsiMu2PtMC_);

 bsTree_->SetBranchAddress(  "JpsiMu1PhiMC"         , &JpsiMu1PhiMC_);
 bsTree_->SetBranchAddress(  "JpsiMu2PhiMC"         , &JpsiMu2PhiMC_);
 bsTree_->SetBranchAddress(  "JpsiCosAlphaMC"         , &JpsiCosAlphaMC_);



 bsTree_->SetBranchAddress(  "JpsiGenPVy"         , &JpsiGenPVy_);
 bsTree_->SetBranchAddress(  "JpsiGenPVx"         , &JpsiGenPVx_);
 bsTree_->SetBranchAddress(  "JpsiGenPVz"         , &JpsiGenPVz_);

 bsTree_->SetBranchAddress(  "JpsiTrigVtxProb"         , &JpsiTrigVtxProb_); 

 bsTree_->SetBranchAddress(  "JpsiMu1TrackerMuonArbitrated"         , &JpsiMu1TrackerMuonArbitrated_ ); 
 bsTree_->SetBranchAddress(  "JpsiMu2TrackerMuonArbitrated"         , &JpsiMu2TrackerMuonArbitrated_ ); 


  bsTree_->SetBranchAddress(  "BSx"				  , &BSx_  );
  bsTree_->SetBranchAddress(  "BSy"				  , &BSy_  );
  bsTree_->SetBranchAddress(  "BSz"				  , &BSz_  );
  bsTree_->SetBranchAddress(  "BSdx"                                , &BSdx_  );
  bsTree_->SetBranchAddress(  "BSdy"                                , &BSdy_  );
  bsTree_->SetBranchAddress(  "BSdz"                                , &BSdz_  );

  bsTree_->SetBranchAddress(  "BSdxdz"                           , &BSdxdz_);
  bsTree_->SetBranchAddress(  "BSdydz"                           , &BSdydz_);


  bsTree_->SetBranchAddress(  "BSsigmaZ"                            , &BSsigmaZ_  );
  bsTree_->SetBranchAddress(  "BSdsigmaZ"                           , &BSdsigmaZ_  );
  bsTree_->SetBranchAddress(  "PVx"				  , &PVx_  );
  bsTree_->SetBranchAddress(  "PVy"				  , &PVy_  );
  bsTree_->SetBranchAddress(  "PVz"				  , &PVz_  );
  bsTree_->SetBranchAddress(  "PVerrx"			  , &PVerrx_  );
  bsTree_->SetBranchAddress(  "PVerry"			  , &PVerry_  );
  bsTree_->SetBranchAddress(  "PVerrz"			  , &PVerrz_  );
  bsTree_->SetBranchAddress(  "JpsiMuMuDCA_beffit"			  , &JpsiMuMuDCA_beffit_);
  bsTree_->SetBranchAddress(  "MuMuDCA"			  , &MuMuDCA_  );
  bsTree_->SetBranchAddress(  "KaonsDCA"                   , &KaonsDCA_  );
  bsTree_->SetBranchAddress(  "MuMuDistance"			  , &MuMuDistance_  );
  bsTree_->SetBranchAddress(  "MuMuDistanceSigma"			  , &MuMuDistanceSigma_  );
  bsTree_->SetBranchAddress(  "MuDr1"			  , &MuDr1_  );
  bsTree_->SetBranchAddress(  "MuDr2"			  , &MuDr2_  );
  bsTree_->SetBranchAddress(  "BdMuMuDCA"			  , &BdMuMuDCA_  );
  bsTree_->SetBranchAddress(  "BdMuMuDistance"			  , &BdMuMuDistance_  );
  bsTree_->SetBranchAddress(  "BdMuMuDistanceSigma"			  , &BdMuMuDistanceSigma_  );
  bsTree_->SetBranchAddress(  "BdMuDr1"			  , &BdMuDr1_  );
  bsTree_->SetBranchAddress(  "BdMuDr2"			  , &BdMuDr2_  );

  bsTree_->SetBranchAddress(  "isPV"				  , &isPV_  );
  bsTree_->SetBranchAddress(  "isBS"				  , &isBS_  );

  bsTree_->SetBranchAddress( "PVx_refit"   , &PVx_refit_      );
  bsTree_->SetBranchAddress( "PVy_refit"   ,	 &PVy_refit_        );
  bsTree_->SetBranchAddress( "PVz_refit"   ,	 &PVz_refit_        );
  bsTree_->SetBranchAddress( "PVerrx_refit",	 &PVerrx_refit_     );
  bsTree_->SetBranchAddress( "PVerry_refit",	 &PVerry_refit_     );
  bsTree_->SetBranchAddress( "PVerrz_refit",	 &PVerrz_refit_     );


  bsTree_->SetBranchAddress( "PVx_refit_cosTheta"   , &PVx_refit_cosTheta_   );
  bsTree_->SetBranchAddress( "PVy_refit_cosTheta"   ,	 &PVy_refit_cosTheta_  );
  bsTree_->SetBranchAddress( "PVz_refit_cosTheta"   ,	 &PVz_refit_cosTheta_  );

  bsTree_->SetBranchAddress( "PVx_refit_closestZ"   , &PVx_refit_closestZ_   );
  bsTree_->SetBranchAddress( "PVy_refit_closestZ"   ,	 &PVy_refit_closestZ_  );
  bsTree_->SetBranchAddress( "PVz_refit_closestZ"   ,	 &PVz_refit_closestZ_  );

  bsTree_->SetBranchAddress( "BsSVx"   ,	 &BsSVx_  );
  bsTree_->SetBranchAddress( "BsSVy"   ,	 &BsSVy_  );
  bsTree_->SetBranchAddress( "BsSVz"   ,	 &BsSVz_  );


  bsTree_->SetBranchAddress(  "JpsiSVx"   ,	 &JpsiSVx_    );
  bsTree_->SetBranchAddress( "JpsiSVy"   ,	 &JpsiSVy_   );
  bsTree_->SetBranchAddress(  "JpsiSVz"   ,	 &JpsiSVz_    );

  bsTree_->Branch(  "triggerbit_HLTmu4TkTkDis"             , &triggerbit_HLTmu4TkTkDis_,                "triggerbit_HLTmu4TkTkDis/I");
  bsTree_->Branch(  "triggerbit_HLTmu4TkDis"		  , &triggerbit_HLTmu4TkDis_,                "triggerbit_HLTmu4TkDis/I");
 bsTree_->Branch(  "triggerbit_HLTDimuon0JpsiMuon"              , &triggerbit_HLTDimuon0JpsiMuon_,                "triggerbit_HLTDimuon0JpsiMuon/I");
 //bsTree_->Branch(  "triggerbit_HLTmu7"		  , &triggerbit_HLTmu7_,                "triggerbit_HLTmu7/I");
  //bsTree_->Branch(  "triggerbit_HLTdoubleIsoMu3"	  , &triggerbit_HLTdoubleIsoMu3_,       "triggerbit_HLTdoubleIsoMu3/I");


  bsTree_->SetBranchAddress( "BdMu1Pt_beffit", &BdMu1Pt_beffit_);
  bsTree_->SetBranchAddress( "BdMu1Pz_beffit", &BdMu1Pz_beffit_);
  bsTree_->SetBranchAddress( "BdMu1Eta_beffit", &BdMu1Eta_beffit_);
  bsTree_->SetBranchAddress( "BdMu1Phi_beffit", &BdMu1Phi_beffit_);
  bsTree_->SetBranchAddress( "BdMu2Pt_beffit", &BdMu2Pt_beffit_);
  bsTree_->SetBranchAddress( "BdMu2Pz_beffit", &BdMu2Pz_beffit_);
  bsTree_->SetBranchAddress( "BdMu2Eta_beffit", &BdMu2Eta_beffit_);
  bsTree_->SetBranchAddress( "BdMu2Phi_beffit", &BdMu2Phi_beffit_);
  bsTree_->SetBranchAddress( "BdJpsiM_nofit", &BdJpsiM_nofit_);
  bsTree_->SetBranchAddress( "BdJpsiEta_nofit", &BdJpsiEta_nofit_);
  bsTree_->SetBranchAddress( "BdJpsiPhi_nofit", &BdJpsiPhi_nofit_);
  bsTree_->SetBranchAddress( "BdJpsiPt_nofit", &BdJpsiPt_nofit_);
  bsTree_->SetBranchAddress( "BdJpsiPz_nofit", &BdJpsiPz_nofit_);
  bsTree_->SetBranchAddress(  "BsSoftMuon1"		  , &BsSoftMuon1_);
  bsTree_->SetBranchAddress(  "BsSoftMuon2"		  , &BsSoftMuon2_);

  bsTree_->SetBranchAddress( "MuonType", &MuonType_ );

  bsTree_->SetBranchAddress(  "JpsiVtxProb"			  , &JpsiVtxProb_  );
  bsTree_->SetBranchAddress(  "BdJpsiVtxProb"			  , &BdJpsiVtxProb_  );
  bsTree_->SetBranchAddress(  "JpsiM_alone"			  , &JpsiM_alone_  );
  bsTree_->SetBranchAddress(  "JpsiPhi_alone"		  , &JpsiPhi_alone_  );
  bsTree_->SetBranchAddress(  "JpsiEta_alone"		  , &JpsiEta_alone_  );
  bsTree_->SetBranchAddress(  "JpsiPt_alone"		  , &JpsiPt_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1Pt_alone"		  , &JpsiMu1Pt_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2Pt_alone"		  , &JpsiMu2Pt_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1Phi_alone"		  , &JpsiMu1Phi_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2Phi_alone"		  , &JpsiMu2Phi_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1Eta_alone"		  , &JpsiMu1Eta_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2Eta_alone"		  , &JpsiMu2Eta_alone_  );
  bsTree_->SetBranchAddress(  "MuonCat1"		  , &MuonCat1_  );
  bsTree_->SetBranchAddress(  "MuonCat2"		  , &MuonCat2_  );


  bsTree_->SetBranchAddress( "MuonPairDR", &MuonPairDR_ );
  bsTree_->SetBranchAddress( "BpMuonPairDR", &BpMuonPairDR_ );
  bsTree_->SetBranchAddress(  "JpsiMuonCat1"		  , &JpsiMuonCat1_  );
  bsTree_->SetBranchAddress(  "JpsiMuonCat2"		  , &JpsiMuonCat2_  );

	// soft muon variables

  bsTree_->SetBranchAddress( "Mu1TrkBSDxy", &Mu1TrkBSDxy_ );
  bsTree_->SetBranchAddress( "Mu1TrkBSDz", &Mu1TrkBSDz_ );
  bsTree_->SetBranchAddress( "Mu1PixelHits", &Mu1PixelHits_  );
  bsTree_->SetBranchAddress( "Mu1TrackerHits", &Mu1TrackerHits_  );
  bsTree_->SetBranchAddress( "Mu1isGood", &Mu1isGood_  );
  bsTree_->SetBranchAddress( "Mu1InnerTrkHighQuality", &Mu1InnerTrkHighQuality_ );

  bsTree_->SetBranchAddress( "Mu2TrkBSDxy", &Mu2TrkBSDxy_ );
  bsTree_->SetBranchAddress( "Mu2TrkBSDz", &Mu2TrkBSDz_ );
  bsTree_->SetBranchAddress( "Mu2PixelHits", &Mu2PixelHits_  );
  bsTree_->SetBranchAddress( "Mu2TrackerHits", &Mu2TrackerHits_ );
  bsTree_->SetBranchAddress( "Mu2isGood", &Mu2isGood_ );
  bsTree_->SetBranchAddress( "Mu2InnerTrkHighQuality", &Mu2InnerTrkHighQuality_ );


  bsTree_->SetBranchAddress( "BpMu1TrkBSDxy", &BpMu1TrkBSDxy_ );
  bsTree_->SetBranchAddress( "BpMu1TrkBSDz", &BpMu1TrkBSDz_ );
  bsTree_->SetBranchAddress( "BpMu1PixelHits", &BpMu1PixelHits_  );
  bsTree_->SetBranchAddress( "BpMu1TrackerHits", &BpMu1TrackerHits_  );
  bsTree_->SetBranchAddress( "BpMu1isGood", &BpMu1isGood_  );
  bsTree_->SetBranchAddress( "BpMu1InnerTrkHighQuality", &BpMu1InnerTrkHighQuality_ );

  bsTree_->SetBranchAddress( "BpMu2TrkBSDxy", &BpMu2TrkBSDxy_ );
  bsTree_->SetBranchAddress( "BpMu2TrkBSDz", &BpMu2TrkBSDz_ );
  bsTree_->SetBranchAddress( "BpMu2PixelHits", &BpMu2PixelHits_  );
  bsTree_->SetBranchAddress( "BpMu2TrackerHits", &BpMu2TrackerHits_ );
  bsTree_->SetBranchAddress( "BpMu2isGood", &BpMu2isGood_ );
  bsTree_->SetBranchAddress( "BpMu2InnerTrkHighQuality", &BpMu2InnerTrkHighQuality_ );


  bsTree_->SetBranchAddress(  "JpsiMuon1Cat_alone"		  , &JpsiMuon1Cat_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMuon2Cat_alone"		  , &JpsiMuon2Cat_alone_  );
  bsTree_->SetBranchAddress(  "BdJpsiMuon1Cat_alone"		  , &BdJpsiMuon1Cat_alone_  );
  bsTree_->SetBranchAddress(  "BdJpsiMuon2Cat_alone"		  , &BdJpsiMuon2Cat_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1d0_alone"             , &JpsiMu1d0_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2d0_alone"             , &JpsiMu2d0_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1dz_alone"             , &JpsiMu1dz_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2dz_alone"             , &JpsiMu2dz_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1chi2_alone"           , &JpsiMu1chi2_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2chi2_alone"           , &JpsiMu2chi2_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1ndof_alone"           , &JpsiMu1ndof_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2ndof_alone"           , &JpsiMu2ndof_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1nHits_alone"                  , &JpsiMu1nHits_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2nHits_alone"                  , &JpsiMu2nHits_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu1nPixHits_alone"               , &JpsiMu1nPixHits_alone_  );
  bsTree_->SetBranchAddress(  "JpsiMu2nPixHits_alone"               , &JpsiMu2nPixHits_alone_  );

  bsTree_->SetBranchAddress(  "BsFitChi2"			  , &BsFitChi2_  );
  bsTree_->SetBranchAddress(  "BsFitNdof"			  , &BsFitNdof_  );
  bsTree_->SetBranchAddress(  "BsFitVtxProb"		  , &BsFitVtxProb_  );
  bsTree_->SetBranchAddress(  "JpsiNumberOfCandidates"        , &JpsiNumberOfCandidates_);
  bsTree_->SetBranchAddress(  "JpsiGenNumberOfCandidates"        , &JpsiGenNumberOfCandidates_);
  bsTree_->SetBranchAddress(  "PhiNumberOfCandidatesBeforeFit"        , &PhiNumberOfCandidatesBeforeFit_);
  bsTree_->SetBranchAddress(  "BsNumberOfCandidatesBeforeFit"        , &BsNumberOfCandidatesBeforeFit_);
  bsTree_->SetBranchAddress(  "BsNumberOfCandidatesAfterFit"        , &BsNumberOfCandidatesAfterFit_);
  bsTree_->SetBranchAddress(  "BsNumberOfCandidatesAfterBestFit"        , &BsNumberOfCandidatesAfterBestFit_);

  bsTree_->SetBranchAddress(  "K1Pt_beffit"			  , &K1Pt_beffit_  );
  bsTree_->SetBranchAddress(  "K1Pz_beffit"			  , &K1Pz_beffit_  );
  bsTree_->SetBranchAddress(  "K1Eta_beffit"			  , &K1Eta_beffit_  );
  bsTree_->SetBranchAddress(  "K1Phi_beffit"			  , &K1Phi_beffit_  );
  bsTree_->SetBranchAddress(  "K2Pt_beffit"			  , &K2Pt_beffit_  );
  bsTree_->SetBranchAddress(  "K2Pz_beffit"			  , &K2Pz_beffit_  );
  bsTree_->SetBranchAddress(  "K2Eta_beffit"			  , &K2Eta_beffit_  );
  bsTree_->SetBranchAddress(  "K2Phi_beffit"			  , &K2Phi_beffit_  );

  bsTree_->SetBranchAddress(  "Mu1Pt_beffit"			  , &Mu1Pt_beffit_  );
  bsTree_->SetBranchAddress(  "Mu1Pz_beffit"			  , &Mu1Pz_beffit_  );
  bsTree_->SetBranchAddress(  "Mu1Eta_beffit"			  , &Mu1Eta_beffit_  );
  bsTree_->SetBranchAddress(  "Mu1Phi_beffit"			  , &Mu1Phi_beffit_  );
  bsTree_->SetBranchAddress(  "Mu2Pt_beffit"			  , &Mu2Pt_beffit_  );
  bsTree_->SetBranchAddress(  "Mu2Pz_beffit"			  , &Mu2Pz_beffit_  );
  bsTree_->SetBranchAddress(  "Mu2Eta_beffit"			  , &Mu2Eta_beffit_  );
  bsTree_->SetBranchAddress(  "Mu2Phi_beffit"			  , &Mu2Phi_beffit_  );

  bsTree_->SetBranchAddress(  "Mu1d0"             , &Mu1d0_  );
  bsTree_->SetBranchAddress(  "Mu2d0"             , &Mu2d0_  );
  bsTree_->SetBranchAddress(  "Mu1dz"             , &Mu1dz_  );
  bsTree_->SetBranchAddress(  "Mu2dz"             , &Mu2dz_  );

  bsTree_->SetBranchAddress(  "BsFitM"			  , &BsFitM_  );
  bsTree_->SetBranchAddress(  "K1Pt_fit"			  , &K1Pt_fit_  );
  bsTree_->SetBranchAddress(  "K2Pt_fit"			  , &K2Pt_fit_  );
  bsTree_->SetBranchAddress(  "PhiM_fit"			  , &PhiM_fit_  );
  bsTree_->SetBranchAddress(  "BsFitEta"			  , &BsFitEta_  );
  bsTree_->SetBranchAddress(  "BsFitPt"			  , &BsFitPt_  );
  bsTree_->SetBranchAddress(  "BsFitPz"			  , &BsFitPz_  );
  bsTree_->SetBranchAddress(  "BsFitPhi"			  , &BsFitPhi_  );
  bsTree_->SetBranchAddress(  "BsFitVtx_x"			  , &BsFitVtx_x_  );
  bsTree_->SetBranchAddress(  "BsFitVtx_y"			  , &BsFitVtx_y_  );
  bsTree_->SetBranchAddress(  "BsFitVtx_z"			  , &BsFitVtx_z_  );
  bsTree_->SetBranchAddress(  "BsM_nofit"			  , &BsM_nofit_  );
  bsTree_->SetBranchAddress(  "BsPhi_nofit"			  , &BsPhi_nofit_  );
  bsTree_->SetBranchAddress(  "BsEta_nofit"			  , &BsEta_nofit_  );
  bsTree_->SetBranchAddress(  "BsPt_nofit"			  , &BsPt_nofit_  );
  bsTree_->SetBranchAddress(  "BsPz_nofit"			  , &BsPz_nofit_  );
  bsTree_->SetBranchAddress(  "JpsiM_nofit"			  , &JpsiM_nofit_  );
  bsTree_->SetBranchAddress(  "JpsiPhi_nofit"		  , &JpsiPhi_nofit_  );
  bsTree_->SetBranchAddress(  "JpsiEta_nofit"		  , &JpsiEta_nofit_  );
  bsTree_->SetBranchAddress(  "JpsiPt_nofit"		  , &JpsiPt_nofit_  );
  bsTree_->SetBranchAddress(  "JpsiPz_nofit"		  , &JpsiPz_nofit_  );
  bsTree_->SetBranchAddress(  "PhiM_nofit"			  , &PhiM_nofit_  );
  bsTree_->SetBranchAddress(  "PhiPhi_nofit"		  , &PhiPhi_nofit_  );
  bsTree_->SetBranchAddress(  "PhiEta_nofit"		  , &PhiEta_nofit_  );
  bsTree_->SetBranchAddress(  "PhiPt_nofit"			  , &PhiPt_nofit_  );
  bsTree_->SetBranchAddress(  "PhiPz_nofit"			  , &PhiPz_nofit_  );
  bsTree_->SetBranchAddress(  "K1Pt_nofit"			  , &K1Pt_nofit_  );
  bsTree_->SetBranchAddress(  "K1HighPurityTrack"			  , &K1HighPurityTrack_  );
  bsTree_->SetBranchAddress(  "K2HighPurityTrack"			  , &K2HighPurityTrack_  );

  bsTree_->SetBranchAddress(  "K1Pz_nofit"			  , &K1Pz_nofit_  );
  bsTree_->SetBranchAddress(  "K1Eta_nofit"			  , &K1Eta_nofit_  );
  bsTree_->SetBranchAddress(  "K1Phi_nofit"			  , &K1Phi_nofit_  );
  bsTree_->SetBranchAddress(  "K1Key_nofit"			  , &K1Key_nofit_  );
  bsTree_->SetBranchAddress(  "K2Eta_nofit"			  , &K2Eta_nofit_  );
  bsTree_->SetBranchAddress(  "K2Pt_nofit"			  , &K2Pt_nofit_  );
  bsTree_->SetBranchAddress(  "K2Pz_nofit"			  , &K2Pz_nofit_  );
  bsTree_->SetBranchAddress(  "K2Phi_nofit"			  , &K2Phi_nofit_  );
  bsTree_->SetBranchAddress(  "K2Key_nofit"			  , &K2Key_nofit_  );
  bsTree_->SetBranchAddress(  "K1Chi2"			  , &K1Chi2_  );
  bsTree_->SetBranchAddress(  "K1nHits"			  , &K1nHits_  );
  bsTree_->SetBranchAddress(  "K2Chi2"			  , &K2Chi2_  );
  bsTree_->SetBranchAddress(  "K2nHits"			  , &K2nHits_  );
  bsTree_->SetBranchAddress(  "K1pixH"			  , &K1pixH_  );
  bsTree_->SetBranchAddress(  "K1trkH"			  , &K1trkH_  );
  bsTree_->SetBranchAddress(  "K2pixH"			  , &K2pixH_  );
  bsTree_->SetBranchAddress(  "K2trkH"			  , &K2trkH_  );
  bsTree_->SetBranchAddress(  "Mu1Chi2"			  , &Mu1Chi2_  );
  bsTree_->SetBranchAddress(  "Mu1nHits"			  , &Mu1nHits_  );
  bsTree_->SetBranchAddress(  "Mu2Chi2"			  , &Mu2Chi2_  );
  bsTree_->SetBranchAddress(  "Mu2nHits"			  , &Mu2nHits_  );
  bsTree_->SetBranchAddress(  "Mu1pixH"			  , &Mu1pixH_  );
  bsTree_->SetBranchAddress(  "Mu1trkH"			  , &Mu1trkH_  );
  bsTree_->SetBranchAddress(  "Mu2pixH"			  , &Mu2pixH_  );
  bsTree_->SetBranchAddress(  "Mu2trkH"			  , &Mu2trkH_  );
  bsTree_->SetBranchAddress(  "costheta"			  , &costheta_  );
  bsTree_->SetBranchAddress(  "cospsi"			  , &cospsi_  );
  bsTree_->SetBranchAddress(  "phi"				  , &phi_  );
  bsTree_->SetBranchAddress(  "Bdcospsi"			  , &Bdcospsi_  );
  bsTree_->SetBranchAddress(  "Bdcostheta"			  , &Bdcostheta_  );
  bsTree_->SetBranchAddress(  "Bdphi"				  , &Bdphi_  );
  bsTree_->SetBranchAddress(  "Bdcospsi"			  , &Bdcospsi_  );
  bsTree_->SetBranchAddress(  "BdcospsiMC"			  , &BdcospsiMC_  );
  bsTree_->SetBranchAddress(  "BdcosthetaMC"			  , &BdcosthetaMC_  );
  bsTree_->SetBranchAddress(  "BdphiMC"				  , &BdphiMC_  );
  bsTree_->SetBranchAddress(  "AngleBsDecayLength"		  , &AngleBsDecayLength_  );
  bsTree_->SetBranchAddress(  "CosDeltaAlpha"		  , &CosDeltaAlpha_  );
  bsTree_->SetBranchAddress(  "BdCosDeltaAlpha"		  , &BdCosDeltaAlpha_  );
  bsTree_->SetBranchAddress(  "isMatched"			  , &isMatched_  );
  bsTree_->SetBranchAddress(  "isMatchedBd"			  , &isMatchedBd_  );
  bsTree_->SetBranchAddress(  "BsLxy"			  , &BsLxy_  );
  bsTree_->SetBranchAddress(  "BsLxyErr"			  , &BsLxyErr_  );
  bsTree_->SetBranchAddress(  "BsErrX"			  , &BsErrX_  );
  bsTree_->SetBranchAddress(  "BsErrY"			  , &BsErrY_  );
  bsTree_->SetBranchAddress(  "BsErrXY"			  , &BsErrXY_  );
  bsTree_->SetBranchAddress(  "BsCt"			  , &BsCt_  );
  bsTree_->SetBranchAddress(  "BsCtErr"			  , &BsCtErr_  );
  bsTree_->SetBranchAddress(  "BsCt3D"			  , &BsCt3D_  );
  bsTree_->SetBranchAddress(  "BsCt2D"			  , &BsCt2D_  );
  bsTree_->SetBranchAddress(  "BsCt2DBS"			  , &BsCt2DBS_  );
  bsTree_->SetBranchAddress(  "BdCt2DBS"			  , &BdCt2DBS_  );
  bsTree_->SetBranchAddress(  "BdCt2DMC"			  , &BdCt2DMC_  );
  bsTree_->SetBranchAddress(  "BdCt3DMC"			  , &BdCt3DMC_  );
  bsTree_->SetBranchAddress(  "BsCtMPV"			  , &BsCtMPV_  );
  bsTree_->SetBranchAddress(  "BsCt3Drefit"		  , &BsCt3Drefit_  );
  bsTree_->SetBranchAddress(  "BsCt2Drefit"		  , &BsCt2Drefit_  );
  bsTree_->SetBranchAddress(  "BsCtMPVrefit"		  , &BsCtMPVrefit_  );
  bsTree_->SetBranchAddress(  "BsCtErr3D"			  , &BsCtErr3D_  );
  bsTree_->SetBranchAddress(  "BsCtErr2D"			  , &BsCtErr2D_  );
  bsTree_->SetBranchAddress(  "BsCtErr2DBS"		, &BsCtErr2DBS_  );
  bsTree_->SetBranchAddress(  "BsCtErr2DClosestZ" , &BsCtErr2DClosestZ_  );
  bsTree_->SetBranchAddress( "BsCtErr2DBSOld"     , &BsCtErr2DBSOld_);
  bsTree_->SetBranchAddress("BsCt2DBSOld", &BsCt2DBSOld_);
  bsTree_->SetBranchAddress("BsCtErr2DClosestZOld", &BsCtErr2DClosestZOld_);
  bsTree_->SetBranchAddress("BsCt2DPVClosestZOld", &BsCt2DPVClosestZOld_);

  bsTree_->SetBranchAddress( "BsCtErr2DOld"  , &BsCtErr2DOld_);
  bsTree_->SetBranchAddress("BsCt2DOld", &BsCt2DOld_);

  bsTree_->SetBranchAddress("BsCtErr2DBS_JpsiVtx", &BsCtErr2DBS_JpsiVtx_);
  bsTree_->SetBranchAddress("BsCtErr2D_JpsiVtx", &BsCtErr2D_JpsiVtx_);
  bsTree_->SetBranchAddress("BsCtErr2DClosestZ_JpsiVtx", &BsCtErr2DClosestZ_JpsiVtx_);
  bsTree_->SetBranchAddress("BsCtErr2DCostheta_JpsiVtx", &BsCtErr2DCostheta_JpsiVtx_);


  bsTree_->SetBranchAddress(  "BdCtErr2DBS"			  , &BdCtErr2DBS_  );
  bsTree_->SetBranchAddress(  "BsCtErr2D2"		  , &BsCtErr2D2_  );
  bsTree_->SetBranchAddress(  "BsCtErrMPV"		  , &BsCtErrMPV_  );
  bsTree_->SetBranchAddress(  "BsCtErr3Drefit"		  , &BsCtErr3Drefit_  );
  bsTree_->SetBranchAddress(  "BsCtErr2Drefit"		  , &BsCtErr2Drefit_  );
  bsTree_->SetBranchAddress(  "BsCtErrMPVrefit"		  , &BsCtErrMPVrefit_  );
  bsTree_->SetBranchAddress(  "K1trkLay"			  , &K1trkLay_  );
  bsTree_->SetBranchAddress(  "K1muDTh"			  , &K1muDTh_  );
  bsTree_->SetBranchAddress(  "K1muCSCh"			  , &K1muCSCh_  );
  bsTree_->SetBranchAddress(  "K1muRPCh"			  , &K1muRPCh_  );
  bsTree_->SetBranchAddress(  "K2trkLay"			  , &K2trkLay_  );
  bsTree_->SetBranchAddress(  "K2muDTh"			  , &K2muDTh_  );
  bsTree_->SetBranchAddress(  "K2muCSCh"			  , &K2muCSCh_  );
  bsTree_->SetBranchAddress(  "K2muRPCh"			  , &K2muRPCh_  );
  bsTree_->SetBranchAddress(  "Mu1trkLay"			  , &Mu1trkLay_  );
  bsTree_->SetBranchAddress(  "Mu1muDTh"			  , &Mu1muDTh_  );
  bsTree_->SetBranchAddress(  "Mu1muCSCh"			  , &Mu1muCSCh_  );
  bsTree_->SetBranchAddress(  "Mu1muRPCh"			  , &Mu1muRPCh_  );
  bsTree_->SetBranchAddress(  "Mu2trkLay"			  , &Mu2trkLay_  );
  bsTree_->SetBranchAddress(  "Mu2muDTh"			  , &Mu2muDTh_  );
  bsTree_->SetBranchAddress(  "Mu2muCSCh"			  , &Mu2muCSCh_  );
  bsTree_->SetBranchAddress(  "Mu2muRPCh"			  , &Mu2muRPCh_  );
  bsTree_->SetBranchAddress(  "K1mcId"			  , &K1mcId_  );
  bsTree_->SetBranchAddress(  "K1momId"			  , &K1momId_  );
  bsTree_->SetBranchAddress(  "K1gmomId"			  , &K1gmomId_  );
  bsTree_->SetBranchAddress(  "K2mcId"			  , &K2mcId_  );
  bsTree_->SetBranchAddress(  "K2momId"			  , &K2momId_  );
  bsTree_->SetBranchAddress(  "K2gmomId"			  , &K2gmomId_  );
  bsTree_->SetBranchAddress(  "Mu1mcId"			  , &Mu1mcId_  );
  bsTree_->SetBranchAddress(  "Mu1momId"			  , &Mu1momId_  );
  bsTree_->SetBranchAddress(  "Mu1gmomId"			  , &Mu1gmomId_  );
  bsTree_->SetBranchAddress(  "Mu2mcId"			  , &Mu2mcId_  );
  bsTree_->SetBranchAddress(  "Mu2momId"			  , &Mu2momId_  );
  bsTree_->SetBranchAddress(  "Mu2gmomId"			  , &Mu2gmomId_  );
  bsTree_->SetBranchAddress(  "Mu1GlobalMuonPromptTight"    , &Mu1GlobalMuonPromptTight_   );
  bsTree_->SetBranchAddress(  "Mu2GlobalMuonPromptTight"    , &Mu2GlobalMuonPromptTight_  );
  bsTree_->SetBranchAddress(  "Mu1TrackerMuonArbitrated"    , &Mu1TrackerMuonArbitrated_  );
  bsTree_->SetBranchAddress(  "Mu1TMLastStationTight"       , &Mu1TMLastStationTight_  );
  bsTree_->SetBranchAddress(  "Mu1TMOneStationTight"         , &Mu1TMOneStationTight_  );
  bsTree_->SetBranchAddress(  "Mu1TMLastStationOptimizedLowPtTight", &Mu1TMLastStationOptimizedLowPtTight_  );
  bsTree_->SetBranchAddress(  "Mu1TMLastStationAngTight"    , &Mu1TMLastStationAngTight_  );
  bsTree_->SetBranchAddress(  "Mu1TMOneStationAngTight"     , &Mu1TMOneStationAngTight_  );
  bsTree_->SetBranchAddress(  "Mu1TMLastStationOptimizedBarrelLowPtTight", &Mu1TMLastStationOptimizedBarrelLowPtTight_  );
  bsTree_->SetBranchAddress(  "Mu2TrackerMuonArbitrated"    , &Mu2TrackerMuonArbitrated_  );
  bsTree_->SetBranchAddress(  "Mu2TMLastStationTight"       , &Mu2TMLastStationTight_  );
  bsTree_->SetBranchAddress(  "Mu2TMOneStationTight"         , &Mu2TMOneStationTight_  );
  bsTree_->SetBranchAddress(  "Mu2TMLastStationOptimizedLowPtTight", &Mu2TMLastStationOptimizedLowPtTight_  );
  bsTree_->SetBranchAddress(  "Mu2TMLastStationAngTight"    , &Mu2TMLastStationAngTight_  );
  bsTree_->SetBranchAddress(  "Mu2TMOneStationAngTight"     , &Mu2TMOneStationAngTight_  );
  bsTree_->SetBranchAddress(  "Mu2TMLastStationOptimizedBarrelLowPtTight", &Mu2TMLastStationOptimizedBarrelLowPtTight_  );

  bsTree_->SetBranchAddress(  "BsCtErr2DCosthetaOld"                   , &BsCtErr2DCosthetaOld_  );
  bsTree_->SetBranchAddress(  "BsCtErr2DCostheta"                   , &BsCtErr2DCostheta_  );
  bsTree_->SetBranchAddress("BsCt3DPVCosTheta", &BsCt3DPVCosTheta_ );
  bsTree_->SetBranchAddress(  "BsCt2DPVClosestZ", &BsCt2DPVClosestZ_);
  bsTree_->SetBranchAddress( "BsCt3DPVClosestZ" , &BsCt3DPVClosestZ_);
//triggerbit_HLTDimuon0JpsiMuon
  bsTree_->SetBranchAddress("BsCt2DPVCosThetaOld", &BsCt2DPVCosThetaOld_ );
  bsTree_->SetBranchAddress("BsCt2DPVCosTheta", &BsCt2DPVCosTheta_ );
  bsTree_->SetBranchAddress(  "BsDist3d"			  , &BsDist3d_  );
  bsTree_->SetBranchAddress(  "BsDist3dErr"			  , &BsDist3dErr_  );
  bsTree_->SetBranchAddress(  "BsTime3d"			  , &BsTime3d_  );
  bsTree_->SetBranchAddress(  "BsTime3dErr"			  , &BsTime3dErr_  );
  bsTree_->SetBranchAddress(  "BsDist2d"			  , &BsDist2d_  );
  bsTree_->SetBranchAddress(  "BsDist2dErr"			  , &BsDist2dErr_  );
  bsTree_->SetBranchAddress(  "BsTime2d"			  , &BsTime2d_  );
  bsTree_->SetBranchAddress(  "BsTime2dErr"			  , &BsTime2dErr_  );
  bsTree_->SetBranchAddress(  "dedxTrk"			  , &dedxTrk_  );
  bsTree_->SetBranchAddress(  "errdedxTrk"			  , &errdedxTrk_  );
  bsTree_->SetBranchAddress(  "numdedxTrk"			  , &numdedxTrk_  );
  bsTree_->SetBranchAddress(  "iPassedCutIdent"		  , &iPassedCutIdent_  );
  bsTree_->SetBranchAddress(  "iPassedCutIdentBd"		  , &iPassedCutIdentBd_  );
  bsTree_->SetBranchAddress(  "BdTrack1Charge"		  , &BdTrack1Charge_  );
  bsTree_->SetBranchAddress(  "K1Fit_par"			  , K1Fit_par_  );
  bsTree_->SetBranchAddress(  "K2Fit_par"			  , K2Fit_par_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigX"			  , &K1Fit_sigX_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigY"			  , &K1Fit_sigY_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigZ"			  , &K1Fit_sigZ_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigX"			  , &K2Fit_sigX_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigY"			  , &K2Fit_sigY_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigZ"			  , &K2Fit_sigZ_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigPX"			  , &K1Fit_sigPX_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigPY"			  , &K1Fit_sigPY_  );
  bsTree_->SetBranchAddress(  "K1Fit_sigPZ"			  , &K1Fit_sigPZ_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigPX"			  , &K2Fit_sigPX_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigPY"			  , &K2Fit_sigPY_  );
  bsTree_->SetBranchAddress(  "K2Fit_sigPZ"			  , &K2Fit_sigPZ_  );

  bsTree_->SetBranchAddress(  "GenNumberOfBdecays"		  , &GenNumberOfBdecays_  );
  bsTree_->SetBranchAddress(  "BmesonsId"			  , BmesonsId_  );
  bsTree_->SetBranchAddress(  "BDauIdMC"			  , BDauIdMC_  );
  bsTree_->SetBranchAddress(  "BDauDauIdMC"			  , BDauDauIdMC_  );
  bsTree_->SetBranchAddress(  "GenNumberOfDaughters"	  , GenNumberOfDaughters_  );
  bsTree_->SetBranchAddress(  "GenNumberOfDaughtersDaughters" , GenNumberOfDaughtersDaughters_  );
  bsTree_->SetBranchAddress(  "BDauMMC"			  , BDauMMC_  );
  bsTree_->SetBranchAddress(  "BDauPtMC"			  , BDauPtMC_  );
  bsTree_->SetBranchAddress(  "BDauPzMC"			  , BDauPzMC_  );
  bsTree_->SetBranchAddress(  "BDauEtaMC"			  , BDauEtaMC_  );
  bsTree_->SetBranchAddress(  "BDauPhiMC"			  , BDauPhiMC_  );
  bsTree_->SetBranchAddress(  "BDauDauMMC"			  , BDauDauMMC_  );
  bsTree_->SetBranchAddress(  "BDauDauPtMC"			  , BDauDauPtMC_  );
  bsTree_->SetBranchAddress(  "BDauDauPzMC"			  , BDauDauPzMC_  );
  bsTree_->SetBranchAddress(  "BDauDauEtaMC"		  , BDauDauEtaMC_  );
  bsTree_->SetBranchAddress(  "BDauDauPhiMC"		  , BDauDauPhiMC_  );
  bsTree_->SetBranchAddress(  "BMMC"			  , BMMC_  );
  bsTree_->SetBranchAddress(  "BsCt2DMC"			  , &BsCt2DMC_  );
  bsTree_->SetBranchAddress(  "BpCt2DMC"			  , &BpCt2DMC_  );
  bsTree_->SetBranchAddress(  "BsCt3DMC"			  , &BsCt3DMC_  );
  bsTree_->SetBranchAddress(  "BsIniFlavour"			  , &BsIniFlavour_  );
  bsTree_->SetBranchAddress(  "BsEndFlavour"			  , &BsEndFlavour_  );
  bsTree_->SetBranchAddress(  "BdIniFlavour"                      , &BdIniFlavour_  );
  bsTree_->SetBranchAddress(  "BdEndFlavour"			  , &BdEndFlavour_  );
  bsTree_->SetBranchAddress(  "ChannelID"			  , &ChannelID_  );
  bsTree_->SetBranchAddress(  "BdChannelID"			  , &BdChannelID_  );
  bsTree_->SetBranchAddress(  "BPtMC"			  , BPtMC_  );
  bsTree_->SetBranchAddress(  "BPxMC"			  , BPxMC_  );
  bsTree_->SetBranchAddress(  "BPyMC"			  , BPyMC_  );
  bsTree_->SetBranchAddress(  "BPzMC"			  , BPzMC_  );
  bsTree_->SetBranchAddress(  "BEtaMC"			  , BEtaMC_  );
  bsTree_->SetBranchAddress(  "BPhiMC"			  , BPhiMC_  );
  bsTree_->SetBranchAddress(  "costhetaMC"		  , costhetaMC_);
  bsTree_->SetBranchAddress(  "phiMC"			  , phiMC_);
  bsTree_->SetBranchAddress(  "cospsiMC"			  , cospsiMC_);
  bsTree_->SetBranchAddress(  "BscosthetaMC"		  , &BscosthetaMC_);
  bsTree_->SetBranchAddress(  "BsphiMC"			  , &BsphiMC_);
  bsTree_->SetBranchAddress(  "BscospsiMC"			  , &BscospsiMC_);

  bsTree_->SetBranchAddress(  "BVtxMC_x" , BVtxMC_x_);
  bsTree_->SetBranchAddress(  "BVtxMC_y" , BVtxMC_y_);
  bsTree_->SetBranchAddress(  "BVtxMC_z" , BVtxMC_z_);
  bsTree_->SetBranchAddress(  "BSVtxMC_x", BSVtxMC_x_);
  bsTree_->SetBranchAddress(  "BSVtxMC_y", BSVtxMC_y_);
  bsTree_->SetBranchAddress(  "BSVtxMC_z", BSVtxMC_z_);
  bsTree_->SetBranchAddress(  "BLxy_MC"  , BLxy_MC_);
  bsTree_->SetBranchAddress(  "BCt_MC"   , BCt_MC_);
  bsTree_->SetBranchAddress(  "BCt_MC2D"   , BCt_MC2D_);
  bsTree_->SetBranchAddress(  "BCt_MC3D"   , BCt_MC3D_);

  bsTree_->SetBranchAddress(  "genBsPV_z"			  , &genBsPV_z_  );
  bsTree_->SetBranchAddress(  "genBsPV_y"			  , &genBsPV_y_  );
  bsTree_->SetBranchAddress(  "genBsPV_x"			  , &genBsPV_x_  );

  bsTree_->SetBranchAddress(  "genBpPV_z"			  , &genBpPV_z_  );
  bsTree_->SetBranchAddress(  "genBpPV_y"			  , &genBpPV_y_  );
  bsTree_->SetBranchAddress(  "genBpPV_x"			  , &genBpPV_x_  );

  bsTree_->SetBranchAddress(  "genBsVtx_z"			  , &genBsVtx_z_  );
  bsTree_->SetBranchAddress(  "genBsVtx_y"			  , &genBsVtx_y_  );
  bsTree_->SetBranchAddress(  "genBsVtx_x"			  , &genBsVtx_x_  );
  bsTree_->SetBranchAddress(  "genBsSVtx_z"			  , &genBsSVtx_z_  );
  bsTree_->SetBranchAddress(  "genBsSVtx_y" 		  , &genBsSVtx_y_  );
  bsTree_->SetBranchAddress(  "genBsSVtx_x"			  , &genBsSVtx_x_  );

  bsTree_->SetBranchAddress(  "isGenJpsiEvent"		  , &isGenJpsiEvent_  );
  bsTree_->SetBranchAddress(  "BdFitChi2_Hyp1"		  , &BdFitChi2_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitNdof_Hyp1"		  , &BdFitNdof_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitVtxProb_Hyp1"		  , &BdFitVtxProb_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitVtxProb"		  , &BdFitVtxProb_  );
  bsTree_->SetBranchAddress(  "BdFitM_Hyp1"			  , &BdFitM_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitEta_Hyp1"		  , &BdFitEta_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitPt_Hyp1"		  , &BdFitPt_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitPz_Hyp1"		  , &BdFitPz_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitPhi_Hyp1"		  , &BdFitPhi_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_x_Hyp1"		  , &BdFitVtx_x_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_y_Hyp1"		  , &BdFitVtx_y_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_z_Hyp1"		  , &BdFitVtx_z_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdM_nofit"		  , &BdM_nofit_  );
  bsTree_->SetBranchAddress(  "BdPhi_nofit"		  , &BdPhi_nofit_  );
  bsTree_->SetBranchAddress(  "BdEta_nofit"		  , &BdEta_nofit_  );
  bsTree_->SetBranchAddress(  "BdPt_nofit"		  , &BdPt_nofit_  );
  bsTree_->SetBranchAddress(  "BdPz_nofit"		  , &BdPz_nofit_  );
  bsTree_->SetBranchAddress(  "KstarMass_nofit_Hyp1"	  , &KstarMass_nofit_Hyp1_  );
  bsTree_->SetBranchAddress(  "KstarMass_nofit_Hyp2"	  , &KstarMass_nofit_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_par_Hyp1"		  , BdK1_kpi_par_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_par_Hyp1"		  , BdK2_kpi_par_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigX_Hyp1"		  , &BdK1_kpi_sigX_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigY_Hyp1"		  , &BdK1_kpi_sigY_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigZ_Hyp1"		  , &BdK1_kpi_sigZ_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigX_Hyp1"		  , &BdK2_kpi_sigX_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigY_Hyp1"		  , &BdK2_kpi_sigY_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigZ_Hyp1"		  , &BdK2_kpi_sigZ_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPX_Hyp1"		  , &BdK1_kpi_sigPX_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPY_Hyp1"		  , &BdK1_kpi_sigPY_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPZ_Hyp1"		  , &BdK1_kpi_sigPZ_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPX_Hyp1"		  , &BdK2_kpi_sigPX_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPY_Hyp1"		  , &BdK2_kpi_sigPY_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPZ_Hyp1"		  , &BdK2_kpi_sigPZ_Hyp1_  );
  bsTree_->SetBranchAddress(  "BdFitChi2_Hyp2"		  , &BdFitChi2_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitNdof_Hyp2"		  , &BdFitNdof_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitVtxProb_Hyp2"		  , &BdFitVtxProb_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitM_Hyp2"			  , &BdFitM_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitEta_Hyp2"		  , &BdFitEta_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitPt_Hyp2"		  , &BdFitPt_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitPz_Hyp2"		  , &BdFitPz_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitPhi_Hyp2"		  , &BdFitPhi_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_x_Hyp2"		  , &BdFitVtx_x_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_y_Hyp2"		  , &BdFitVtx_y_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdFitVtx_z_Hyp2"		  , &BdFitVtx_z_Hyp2_  );

  bsTree_->SetBranchAddress(  "BdPVx_refit", &BdPVx_refit_   );
  bsTree_->SetBranchAddress(  "BdPVy_refit", &BdPVy_refit_   );
  bsTree_->SetBranchAddress(  "BdPVz_refit", &BdPVz_refit_   );
  bsTree_->SetBranchAddress(  "BdPVerrx_refit", &BdPVerrx_refit_);
  bsTree_->SetBranchAddress(  "BdPVerry_refit", &BdPVerry_refit_);
  bsTree_->SetBranchAddress(  "BdPVerrz_refit", &BdPVerrz_refit_);
  bsTree_->SetBranchAddress(  "BdNumberOfCandidates"        , &BdNumberOfCandidates_);
  bsTree_->SetBranchAddress(  "BdK1_kpi_par_Hyp2"		  , BdK1_kpi_par_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_par_Hyp2"		  , BdK2_kpi_par_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigX_Hyp2"		  , &BdK1_kpi_sigX_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigY_Hyp2"		  , &BdK1_kpi_sigY_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigZ_Hyp2"		  , &BdK1_kpi_sigZ_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigX_Hyp2"		  , &BdK2_kpi_sigX_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigY_Hyp2"		  , &BdK2_kpi_sigY_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigZ_Hyp2"		  , &BdK2_kpi_sigZ_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPX_Hyp2"		  , &BdK1_kpi_sigPX_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPY_Hyp2"		  , &BdK1_kpi_sigPY_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1_kpi_sigPZ_Hyp2"		  , &BdK1_kpi_sigPZ_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPX_Hyp2"		  , &BdK2_kpi_sigPX_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPY_Hyp2"		  , &BdK2_kpi_sigPY_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK2_kpi_sigPZ_Hyp2"		  , &BdK2_kpi_sigPZ_Hyp2_  );
  bsTree_->SetBranchAddress(  "BdK1Pt_nofit" 		  , &BdK1Pt_nofit_  );
  bsTree_->SetBranchAddress(  "BdK1Pz_nofit" 		  , &BdK1Pz_nofit_  );
  bsTree_->SetBranchAddress(  "BdK1Eta_nofit" 		  , &BdK1Eta_nofit_  );
  bsTree_->SetBranchAddress(  "BdK1Phi_nofit" 		  , &BdK1Phi_nofit_  );
  bsTree_->SetBranchAddress(  "BdK1Key_nofit" 		  , &BdK1Key_nofit_  );
  bsTree_->SetBranchAddress(  "BdK2Pt_nofit" 		  , &BdK2Pt_nofit_  );
  bsTree_->SetBranchAddress(  "BdK2Pz_nofit" 		  , &BdK2Pz_nofit_  );
  bsTree_->SetBranchAddress(  "BdK2Eta_nofit" 		  , &BdK2Eta_nofit_  );
  bsTree_->SetBranchAddress(  "BdK2Phi_nofit" 		  , &BdK2Phi_nofit_  );
  bsTree_->SetBranchAddress(  "BdK2Key_nofit" 		  , &BdK2Key_nofit_  );
  bsTree_->SetBranchAddress(  "BdLxy"			  , &BdLxy_  );
  bsTree_->SetBranchAddress(  "BdLxyErr"			  , &BdLxyErr_  );
  bsTree_->SetBranchAddress(  "BdErrX"			  , &BdErrX_  );
  bsTree_->SetBranchAddress(  "BdErrY"			  , &BdErrY_  );
  bsTree_->SetBranchAddress(  "BdErrXY"			  , &BdErrXY_  );
  bsTree_->SetBranchAddress(  "BdCt"			  , &BdCt_  );
  bsTree_->SetBranchAddress(  "BdCtErr"			  , &BdCtErr_  );
  bsTree_->SetBranchAddress(  "BdDist3d"			  , &BdDist3d_  );
  bsTree_->SetBranchAddress(  "BdDist3dErr"			  , &BdDist3dErr_  );
  bsTree_->SetBranchAddress(  "BdTime3d"			  , &BdTime3d_  );
  bsTree_->SetBranchAddress(  "BdTime3dErr"			  , &BdTime3dErr_  );
  bsTree_->SetBranchAddress(  "BdDist2d"			  , &BdDist2d_  );
  bsTree_->SetBranchAddress(  "BdDist2dErr"			  , &BdDist2dErr_  );
  bsTree_->SetBranchAddress(  "BdTime2d"			  , &BdTime2d_  );
  bsTree_->SetBranchAddress(  "BdTime2dErr"                   , &BdTime2dErr_  );
  bsTree_->SetBranchAddress("BdCt2DPVCosTheta", &BdCt2DPVCosTheta_ );
  bsTree_->SetBranchAddress("Bdt2DPVCosTheta", &Bdt2DPVCosTheta_ );
  bsTree_->SetBranchAddress("BdCtErr2DCostheta", &BdCtErr2DCostheta_ );
  bsTree_->SetBranchAddress("BdtErr2DCostheta", &BdtErr2DCostheta_ );

  bsTree_->SetBranchAddress(  "BdSoftMuon1"                   , &BdSoftMuon1_  );
  bsTree_->SetBranchAddress(  "BdSoftMuon2"                   , &BdSoftMuon2_  );

  bsTree_->SetBranchAddress(  "BdMCKstarKaon"                   , &BdMCKstarKaon_  );
  bsTree_->SetBranchAddress(  "BdMCKstarPion"                   , &BdMCKstarPion_  );
  bsTree_->SetBranchAddress(  "BdK1mcId"           , &BdK1mcId_     );
  bsTree_->SetBranchAddress(  "BdK1momId"	  ,	 &BdK1momId_    );
  bsTree_->SetBranchAddress(  "BdK1gmomId"	  ,	 &BdK1gmomId_   );
  bsTree_->SetBranchAddress(  "BdK2mcId"	          ,	 &BdK2mcId_     );
  bsTree_->SetBranchAddress(  "BdK2momId"	  ,	 &BdK2momId_    );
  bsTree_->SetBranchAddress(  "BdK2gmomId"	  ,	 &BdK2gmomId_   );
  bsTree_->SetBranchAddress(  "BdMu1mcId"	  ,	 &BdMu1mcId_    );
  bsTree_->SetBranchAddress(  "BdMu1momId"	  ,	 &BdMu1momId_   );
  bsTree_->SetBranchAddress(  "BdMu1gmomId"        ,	 &BdMu1gmomId_  );
  bsTree_->SetBranchAddress(  "BdMu2mcId"	  ,	 &BdMu2mcId_    );
  bsTree_->SetBranchAddress(  "BdMu2momId"	  ,	 &BdMu2momId_   );
  bsTree_->SetBranchAddress(  "BdMu2gmomId"        ,	 &BdMu2gmomId_  );


}
