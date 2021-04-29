NEWCONDITIONS = False
OUTPUTSITE = 'T2_US_UCSD'
DATASET = '/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/Run3Winter20DRPremixMiniAOD-packHBTDC_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW'

# template for crab submission.
# submit_jobs.py will insert definitions above
conditionType = "def" # default
# if new L1TriggerObjects conditons have been specified with either
# a new tag or file
if NEWCONDITIONS:
    conditionType = "new_cond"

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = "hcal_" + "11_X_MC" + "_" + conditionType+"_LLP_350_10m"
config.General.transferLogs = True
config.General.transferOutputs = True

# Name of the CMSSW configuration file                                                                                                       
config.JobType.psetName = 'ntuple_maker_' + conditionType + '.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
#config.JobType.maxMemoryMB = 2500                                                                                                           
config.JobType.outputFiles = ['L1Ntuple.root']

config.Data.inputDataset = DATASET
#config.Data.userInputFiles = ['/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/10D453C2-6161-9C46-81FC-86BE05BFE75A.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/2EE8E62A-EA26-9D4E-8225-09E121869187.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/3B24D3E3-EF10-AA4F-B83E-43687B149592.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/686168B8-276C-C540-B71F-D484E12E5D9B.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/7F07E812-2A4C-1F44-B391-7293548BA25B.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/97CB9490-E89A-FB4A-A4BC-FD8E3D7E7595.root/', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/9A5BA94B-2821-9544-95E4-DBA1D82A42B8.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/9B6B9056-2B2B-2049-90B9-F4A5945E1B0E.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/AE1E95BB-61AE-504C-925B-54DB70519F28.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/C699C1FC-D398-A040-8195-3267099C7739.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/DDB82C98-2E22-D743-8C8F-B0BF6DED54A0.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/E9ECD2BD-FE90-3942-B0FA-A170D92E7AD7.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/EC854807-A260-4940-8860-DA7BE9113291.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/F4C824A6-CC41-B24A-AE34-3752F607F909.root']
config.Data.ignoreLocality = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 3
#config.Data.totalUnits = 15
config.Data.useParent = False
#config.Data.outputPrimaryDataset = 'QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL'

# This string is used to construct the output dataset name                                                                                   
config.Data.outputDatasetTag = 'Hcal' + "11_X_MC" + '_' + conditionType+"_LLP_350_10m"

# These values only make sense for processing data                                                                                           
#    Select input data based on a lumi mask                                                                                                  
#config.Data.lumiMask = LUMIMASK

#    Select input data based on run-ranges                                                                                                   
#config.Data.runRange = str(RUN)

# Where the output files will be transmitted to                                                                                              
config.Site.storageSite = OUTPUTSITE
#config.Site.whitelist = ["T3_UK_London_QMUL", "T2_FR_CCIN2P3", "T3_BG_UNI_SOFIA", "T2_CN_Beijing", "T1_RU_JINR", "T2_CH_CSCS", "T2_DE_RWTH", "T3_IT_Bologna", "T1_IT_CNAF", "T3_FR_IPNL"]
#config.Site.whitelist = [OUTPUTSITE]
