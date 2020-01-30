# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=ntuple_maker_def.py -n 10 --no_output --era=Run3 --data --conditions=106X_upgrade2021_realistic_v4 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=file:/eos/cms/store/user/lowang/mh1000_pl10000_step1.root --no_exec
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
#    fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/mh1000_pl1000_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/mh350_ml160_pl500_step1.root'),
fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/mh250_ml120_pl500_step1.root'),
#                            fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/RelValNuGun_step1.root'),
#    fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_1.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_2.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_3.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_4.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_5.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_6.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_7.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_8.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_9.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_10.root'), # QCD files from Long
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2021_realistic_v4', '')

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

#call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process)

# End of customisation functions

# Customisation from command line

process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
