# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=ntuple_maker_def.py -n 10 --no_output --era=Run3 --data --conditions=106X_upgrade2021_realistic_v4 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=/store/user/lowang/LLP_htobbbb/step1/MH-1000_MFF-450_CTau-1000mm_step1.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RAW2DIGI',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200710_204538/0000/RelValNuGun_PU_step1_2.root'),
                            fileNames = cms.untracked.vstring('/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/10D453C2-6161-9C46-81FC-86BE05BFE75A.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/2EE8E62A-EA26-9D4E-8225-09E121869187.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/3B24D3E3-EF10-AA4F-B83E-43687B149592.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/686168B8-276C-C540-B71F-D484E12E5D9B.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/7F07E812-2A4C-1F44-B391-7293548BA25B.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/97CB9490-E89A-FB4A-A4BC-FD8E3D7E7595.root/', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/9A5BA94B-2821-9544-95E4-DBA1D82A42B8.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/9B6B9056-2B2B-2049-90B9-F4A5945E1B0E.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/AE1E95BB-61AE-504C-925B-54DB70519F28.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/C699C1FC-D398-A040-8195-3267099C7739.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/DDB82C98-2E22-D743-8C8F-B0BF6DED54A0.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/E9ECD2BD-FE90-3942-B0FA-A170D92E7AD7.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/EC854807-A260-4940-8860-DA7BE9113291.root', '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/F4C824A6-CC41-B24A-AE34-3752F607F909.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:200'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAWsimHcalTP 

#call to customisation function L1TReEmulFromRAWsimHcalTP imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAWsimHcalTP(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMU
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleGEN 

#call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process)
process = L1NtupleGEN(process)

# End of customisation functions

# Customisation from command line

process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
