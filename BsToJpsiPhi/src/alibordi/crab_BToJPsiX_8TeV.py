from WMCore.Configuration import Configuration

config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'BuToJPsiX_8TeV'

config.section_('JobType')
config.JobType.pluginName = 'Analysis' #'PrivateMC'
config.JobType.psetName = 'btokstarmumu_Run2012.py'
config.JobType.outputFiles = ['BToKstarMuMu.root']

config.section_('Data')
config.Data.inputDataset = '/BpToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.unitsPerJob = 5
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFN = '/store/user/pchen/BToKstarMuMu/ntp/v3p2/BuToJPsiX_8TeV'

config.section_('User')
config.User.voGroup = 'twcms'

config.section_('Site')
config.Site.blacklist = ['T3_TW_NTU_HEP']
config.Site.storageSite = 'T3_TW_NTU_HEP'
