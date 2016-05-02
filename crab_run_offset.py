from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'offset_treemaker'
config.General.workArea = 'crab_projects/Data80x2'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_offset.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ["Fall15_25nsV2_MC.db"]
config.JobType.inputFiles = ["Fall15_25nsV2_DATA.db", "pileup_JSON_1911.txt"]
config.JobType.outputFiles = ["Offset_Data.root"]

config.section_("Data")
config.Data.inputDataset = '/ZeroBias/CMSSW_8_0_0-80X_dataRun2_relval_v0_RelVal_zb2015D-v1/RECO'
 # '/ZeroBias/Run2015D-16Dec2015-v1/AOD'
 # '/SingleNeutrino/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'
config.Data.splitting = 'LumiBased' # FileBased
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/charring/MikkoTrees/Data80x'
config.Data.publication = False
#config.Data.ignoreLocality = True
#config.Data.publishDataName = 'offset_analysis'

config.section_("Site")
#config.Site.blacklist = ['T1_US_FNAL']
#config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = "T3_US_FNALLPC"

# source /cvmfs/cms.cern.ch/crab3/crab.csh
