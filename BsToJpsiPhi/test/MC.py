import FWCore.ParameterSet.Config as cms


process = cms.Process("MUMU")

#-- NUMBER OF EVENTS --#
process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

#-- SOURCE FILES --#
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/BcToBsPi_JpsiPhiPi_MuMuKKPi_JpsiPhiFilter_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/4E457795-DA2E-7340-9F04-9FD229FE5465.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_DGamma0_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FEE2F092-6C44-E811-BAAE-AC1F6B1AF002.root',
#'root://cms-xrd-global.cern.ch//store/data/Run2017C/Charmonium/MINIAOD/12Sep2017-v1/70000/B6F7B207-48A7-E711-9229-FA163EC5EF3C.root',		
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/3A0A3362-36AC-E711-BBED-FA163ED8E79E.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/7A22593E-4AAC-E711-84F7-02163E0176A5.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_N1_94X_mc2017_realistic_v14-v1/90000/D6A9C35A-2B42-E811-B248-008CFAC93F40.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_N1_94X_mc2017_realistic_v14-v1/70000/FE34FCB0-EE41-E811-B5AD-008CFAC93B94.root',

#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/14422CC4-3BAD-E711-85D2-FA163E3BCC72.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/3A0A3362-36AC-E711-BBED-FA163ED8E79E.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/4469C355-B2AC-E711-BE8E-02163E012D69.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/4CD344EC-3BAD-E711-BD95-FA163E3201B4.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/549DB69D-3CAD-E711-849F-0CC47A7C3422.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/5AC2A21E-64AC-E711-B29B-FA163E7625E2.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6E47188D-3BAD-E711-8B3F-0CC47A6C1874.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/7A22593E-4AAC-E711-84F7-02163E0176A5.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/80125A7C-7EAC-E711-8CDA-FA163ECC515A.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/8E5FA699-64AC-E711-9772-FA163EE771A0.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/949F3A31-6BAC-E711-B38C-FA163EF027B6.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A0126530-5EAC-E711-98A1-02163E0164E3.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A60BFA1C-3CAD-E711-AE31-008CFAC91B60.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/B8DB27C5-57AC-E711-9112-FA163E7B5756.root'


#trigger v9
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v3/20000/0484C052-468E-E311-BE47-E0CB4E29C4D0.root'															
									 #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v3/20000/083A7164-518E-E311-9E3D-E0CB4E1A1180.root'
									 #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToJPsiPhi_2K2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v1/0000/00B97FCB-39D9-E111-BD80-001A645C185E.root'
									#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v3/20000/02F84F4C-5F8E-E311-B85D-20CF305B0585.root'
									 #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BdToKstarJPsi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/0002ECCA-AD4A-E311-99EE-1CC1DE050110.root'

#triggers v9, v12 B+
										#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BuJpsiK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/00086026-C705-E411-B655-90E6BAE8CC0C.root'

#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BuJpsiK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/0CCCA02D-AE06-E411-90D5-52540017A56F.root'

##trigger v12
									#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v3/20000/1E45701D-288E-E311-B50D-0025907B4EDC.root'
								#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BuJpsiK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/00086026-C705-E411-B655-90E6BAE8CC0C.root'
									#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BuJpsiK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/04B650A2-8505-E411-A987-002590D0AFFE.root' # ongelmallinen file -481,65 	

                     		#'file:/afs/cern.ch/work/t/terhi/private/BsMCtest.root'
									#'file:/afs/cern.ch/work/t/terhi/private/BuMCtest.root'
									#'file:/tmp/terhi/Bstest.root'
									#'file:/tmp/terhi/LambdaB.root'
									#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BsToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v1/0000/005DE3B0-FDDC-E111-9812-00266CFFC198.root'

                            )
)
#Global tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = "92X_upgrade2017_realistic_v10"


#-- LOGGER --#
process.MessageLogger = cms.Service(
                                    "MessageLogger",
                                    #destinations = cms.untracked.vstring('BuToJPsiK_log.txt'),#qui si decide il file di output
                                    #destinations = cms.untracked.vstring('BsMCtest_log.txt'),#qui si decide il file di output
                                    #destinations = cms.untracked.vstring('BdToJPsiKstar_log.txt'),#qui si decide il file di output
                                    default = cms.untracked.PSet( reportEvery = cms.untracked.int32(1000) ) #qui riporta i messaggi ogni ()
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#-- GEOMETRY + B-FIELD --#
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


#-- GLOBAL TAG --#
#process.GlobalTag.globaltag = cms.string('START53_V19F::All')
#process.GlobalTag.globaltag = cms.string('92X_upgrade2017_realistic_v10')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12', '')

#-- PAT LAYER 0+1 --#
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *
#usePFIso( process )

#-- PAT OVERLAP mu/ele --#
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

#-- PAT MC MATCHING mu --#
#MUON MC-MATCHING VALUES FROM BsMuMu MUON-ID STUDIES
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.matched = cms.InputTag("prunedGenParticles")
process.muonMatch.maxDeltaR = cms.double(0.12)
process.muonMatch.maxDPtRel = cms.double(0.3)
process.muonMatch.checkCharge = cms.bool(True)
process.muonMatch.resolveAmbiguities = cms.bool(True)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT MC MATCHING ele --#
#ELECTRON MC-MATCHING VALUES FROM: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
# object        electron        photon  muon    tau to jet      jet
# maxDPtRel     0.5             1.0     0.5     3.0             3.0
# maxDeltaR     0.5             0.2     0.5     0.1             0.4
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi")
process.electronMatch.matched = cms.InputTag("prunedGenParticles")
process.electronMatch.maxDeltaR = cms.double(0.5)
process.electronMatch.maxDPtRel = cms.double(0.5)
process.electronMatch.checkCharge = cms.bool(True)
process.electronMatch.resolveAmbiguities = cms.bool(True)
process.electronMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT MC MATCHING jet --#
#JET MC-MATCHING VALUES FROM: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
# object        electron        photon  muon    tau to jet      jet
# maxDPtRel     0.5             1.0     0.5     3.0             3.0
# maxDeltaR     0.5             0.2     0.5     0.1             0.4
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
process.patJetPartonMatch.src     = cms.InputTag("slimmedJets")
process.patJetPartonMatch.matched = cms.InputTag("prunedGenParticles")
process.patJetPartonMatch.maxDeltaR = cms.double(0.4)
process.patJetPartonMatch.maxDPtRel = cms.double(3.0)
process.patJetPartonMatch.checkCharge = cms.bool(False)
process.patJetPartonMatch.resolveAmbiguities = cms.bool(True)
process.patJetPartonMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT TRACKS --#
from PhysicsTools.PatAlgos.tools.trackTools import *

#kaons
#process.allKTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
 #                                   src = cms.InputTag("isolatedTracks"),
 #                                   particleType = cms.string('K+')
#                                    )
#process.kTracks = cms.EDFilter("CandViewRefSelector",
#                               src = cms.InputTag("allKTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
 #                              cut = cms.string("pt > 0.5 & abs(eta) < 2.4")
  #                             )

#pions
#process.allPiTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
 #                                   src = cms.InputTag("isolatedTracks"),
 #                                   particleType = cms.string('pi+')
  #                                  )
#process.piTracks = cms.EDFilter("CandViewRefSelector",
 #                              src = cms.InputTag("allPiTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
  #                             cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
   #                            )

#-- SIM HITS --#
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
                                            #src = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
                                            #setStatus  = cms.int32(8),
                                            ##   particleTypes = cms.vstring(""),
                                            #filter = cms.vstring("pt > 0.0"),  # just for testing
                                            #genParticles = cms.InputTag("genParticles")
                                            #)

#-- PAT ELECTRONS --#
#process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
#process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
#process.mvaID = cms.Sequence(  process.mvaTrigV0 +  process.mvaTrigNoIPV0 + process.mvaNonTrigV0 + process.simpleEleIdSequence )

#process.patElectrons.useParticleFlow  = cms.bool(True)
#process.patElectrons.pfElectronSource = cms.InputTag("packedPFCandidates")
#process.patElectrons.embedTrack       = cms.bool(True)

#process.patElectrons.electronIDSources = cms.PSet(
    #MVA
 #   mvaTrigV0 = cms.InputTag("mvaTrigV0"),
 #   mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
 #   mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
    #non-MVA
    #eidVeto = cms.InputTag("eidVeto"),
 #   eidTight = cms.InputTag("eidTight"),
 #   eidLoose = cms.InputTag("eidLoose"),
 #   eidRobustTight = cms.InputTag("eidRobustTight"),
 #   eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
 #   eidRobustLoose = cms.InputTag("eidRobustLoose"),
 #   simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
 #   simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
 #   simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
 #   simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
 #   simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
 #   simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
 #   simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
 #   simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
 #   simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
 #   simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
 #   simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
 #   simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
 #   )

#-- PAT CONVERSION --#
#process.patConversions = cms.EDProducer("PATConversionProducer",
 #                                       electronSource = cms.InputTag("cleanPatElectrons")  ,
                                        #muonSource = cms.InputTag("cleanPatMuons")
                                        # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,
#)

#-- PAT JETS --#
import PhysicsTools.PatAlgos.tools.jetTools as jetTools

#process.patJets.addTagInfos = cms.bool(True)

#jetTools.switchJetCollection(process,
 #                            cms.InputTag('slimmedJetsAK8'),
 #                            doJTA              = True,
 #                            doBTagging         = True,
 #                            btagInfo           = cms.vstring('impactParameterTagInfos','secondaryVertexTagInfos','softPFMuonsTagInfos','softPFElectronsTagInfos'),
 #                            btagdiscriminators = cms.vstring('jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags'),
  #                           jetCorrLabel     = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),
                             #### data:                   jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']),
   #                          doType1MET       = False,
   #                          genJetCollection = cms.InputTag("slimmedGenJetsAK8"),
   #                          doJetID          = True,
   #                          jetIdLabel       = "ak8",
  #                           outputModules    = [],
#)

#-- ANALYZER TAGS AND PARAMETERS --#
process.bsVertexAnalysis = cms.EDAnalyzer("BsToJpsiPhiAnalysis",
                                          isMCstudy                     = cms.bool( True ),
                                          #TriggerTag                    = cms.InputTag("TriggerResults"),
                                          genParticlesLabel             = cms.InputTag("prunedGenParticles"),                                        
                                          MuonTag                       = cms.InputTag("slimmedMuons"),
                                          JetCollection                 = cms.InputTag("slimmedJets"),
                                          PUInfo                        = cms.InputTag("slimmedAddPileupInfo"),
                                          vertexBeamSpot                = cms.InputTag("offlineBeamSpot"),
                                          primaryvertex                 = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                          triggerresults                = cms.InputTag("TriggerResults",'',"HLT"),
                                          ElectronTag                   = cms.InputTag("slimmedElectrons"),
                                          track                         = cms.InputTag("packedPFCandidates"),
                                          isotrack                      = cms.InputTag("isolatedTracks"),
                                          StoreDeDxInfo                 = cms.bool( False ),
                                          JpsiMassWindowBeforeFit       = cms.double(0.31), #leave this selection looser than the trigger one for the efficiency calculation
                                          JpsiMassWindowAfterFit        = cms.double(0.150),
                                          MuonPtCut                     = cms.double(4),
                                          JpsiPtCut                     = cms.double(7),
                                          KaonTrackPtCut                = cms.double(0.7),
                                          BdKaonTrackPtCut              = cms.double(0.6),
                                          PhiMassWindowBeforeFit        = cms.double(0.03),
                                          PhiMassWindowAfterFit         = cms.double(0.02),
                                          BsLowerMassCutBeforeFit       = cms.double(4.5),
                                          BsUpperMassCutBeforeFit       = cms.double(6),
                                          BsLowerMassCutAfterFit        = cms.double(5),
                                          BsUpperMassCutAfterFit        = cms.double(6),
                                          KstarMassWindowBeforeFit      = cms.double(0.2),
                                          KstarMassWindowAfterFit       = cms.double(0.15),
                                          BdLowerMassCutBeforeFit       = cms.double(4.5),
                                          BdUpperMassCutBeforeFit       = cms.double(6),
                                          BdLowerMassCutAfterFit        = cms.double(4.9),
                                          BdUpperMassCutAfterFit        = cms.double(5.7),
                                          verbose                       = cms.bool( False ),
                                          TestVerbose                   = cms.bool( False ),
                                          BsPDGMass = cms.double(5.3699),
                                          BdPDGMass = cms.double(5.2794),
                                          BpPDGMass = cms.double(5.2790),
                                          #outputFile                   = cms.untracked.string("BuToJPsiK_AllOktest.root"),
                                          #outputFile                   = cms.untracked.string("BuToJPsiKCt_Dimuon8Trig.root"),
                                          outputFile                   = cms.untracked.string("BsJpsiPhi_dG0_data.root"),#BsJpsiPhi_dG0.root
                                          #outputFile                   = cms.untracked.string("BdToJPsiKstar.root"),
)

#-- PAT MUONS --#
#process.patMuons.useParticleFlow   = cms.bool(True)
process.patMuons.pfMuonSource      = cms.InputTag("packedPFCandidates")
process.patMuons.embedPFCandidate  = cms.bool(True)

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()

from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
addMCinfo(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR        = 0.1
process.muonMatchHLTL3.maxDPtRel        = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR  = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel  = 10.0
process.muonMatchHLTTrackMu.maxDeltaR   = 0.1
process.muonMatchHLTTrackMu.maxDPtRel   = 10.0

process.patMuons = cms.EDFilter("PATMuonSelector",
                                src = cms.InputTag("patMuonsWithTrigger"),
                                # cut = cms.string("p>2 && abs(eta)<2.4"),
                                cut = cms.string("p>0 && abs(eta)<1000"),
)

process.bphAnalyzer = cms.EDAnalyzer('TestBPHRecoDecay',
#    verbose = cms.untracked.string('t')
)


#-- WRAPPING UP --#
#process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)

#process.pat = cms.Path(process.patDefaultSequence)

process.ntup = cms.Path( process.bsVertexAnalysis )

#process.schedule = cms.Schedule(process.pat, process.ntup )
#process.p = mumu.Path(process.VertexAnalysis)
#######               DUMP Completo                         ########
####################################################################
#temp = process.dumpPython()
#outputfile = file("Complete.py",'w')
#outputfile.write(temp)
#outputfile.close()
