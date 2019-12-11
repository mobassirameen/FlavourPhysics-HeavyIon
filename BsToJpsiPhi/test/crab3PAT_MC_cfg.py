from WMCore.Configuration import Configuration
config = Configuration()

###### TARKISTA SEURAAVAT ASIAT #####
# 1. datasetpath on oikea
# 2. datasetpathissa maariteltya datasettia vastaa oikea JSON
# 3. data tallentuu haluttuun kansioon
# 4. output tiedoston nimi on oikea
# 5. MC.py tiedostossa maaritelty global tag on datasetin global tag
# 6. MC.py tiedostossa ei ole kommentoimattomia datatiedostoja seka eventtien maara on -1
# 7. tarkista ignorelocality & whitelist!!


config.section_('General')
config.General.requestName = 'BsToJPsiPhi_AllTrig_dG0'
#request name is the name of the folder where crab log is saved
config.General.workArea = 'Bs_Data_05_MAR'
config.General.transferOutputs = True
config.General.transferLogs = False


config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/m/mumuhamm/private/CMSSW_9_4_0/src/HeavyFlavorAnalysis/BsToJpsiPhi/test/MC.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['BsJpsiPhi_dG0_data.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

#config.Data.inputDataset = '/LambdaBToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM' 

#config.Data.inputDataset = '/BpToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
#config.Data.inputDataset = '/B0ToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'

#config.Data.inputDataset = '/BsToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' 

#config.Data.inputDataset = '/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/Summer12_DR53X-PU_RD2_START53_V19F-v3/AODSIM' # sample w mu pt + k pt cuts and dG > 0 
#config.Data.inputDataset = '/BsToJpsiPhiV2_BFilter_TuneZ2star_8TeV-pythia6-evtgen/Summer12DR53X-PU_RD2_20150709_START53_V19F-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/BsToJPsiPhi_2K2MuFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'  # sample w/o mu pt + k pt cuts and dG = 0
#config.Data.inputDataset = '/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17MiniAOD-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/MINIAODSIM'  # sample w mu pt + k pt cuts and dG = 0
#==================================================================================================For data next two lines
config.Data.lumiMask = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
#config.Data.runRange = '299370-300238' # '193093-194075'
#bs signal 1 2017
#config.Data.inputDataset = '/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIFall17MiniAODv2-PU2017_12Apr2018_N1_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#bd signal 1 2017
#config.Data.inputDataset = '/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIFall17MiniAODv2-PU2017_12Apr2018_N1_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#bd signal 2 2017
#config.Data.inputDataset = '/BsToJpsiPhi_BMuonFilter_DGamma0_TuneCP5_13TeV-pythia8-evtgen/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/Charmonium/Run2017C-12Sep2017-v1/MINIAOD'
#config.Data.inputDataset = '/store/data/Run2017E/Charmonium/MINIAOD/31Mar2018-v1/90000/FE953E8E-1838-E811-85F2-0025905C42A8.root'
config.Data.inputDataset = '/Charmonium/Run2017E-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset = '/BuJpsiK_TuneZ2star_8TeV_Pythia6/Summer12_DR53X-PU_RD2_START53_V19F-v1/AODSIM'	
config.Data.unitsPerJob     = 7000#5000
config.Data.totalUnits      = 70000000
config.Data.splitting       = 'EventAwareLumiBased'
config.JobType.maxMemoryMB = 8000
#config.Data.splitting = 'LumiBased'

#config.Data.splitting = 'Automatic'#'FileBased'
#config.Data.unitsPerJob = 1
#NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/user/mumuhamm/BsSecond'
config.Data.publication = False

############## 

config.Data.ignoreLocality = False

config.section_("Site")
#config.Site.whitelist = ['T2_FI_HIP']
#config.Site.blacklist = ['T2_FI_HIP']
config.Site.storageSite = 'T2_IN_TIFR'





