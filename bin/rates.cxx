// Script for calculating rate histograms
// Originally from Aaron Bundock
// Edited by Gillian Kopp to add multiplicity studies for LLP L1 trigger using depth and timing (2020)
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"


/* TODO: put errors in rates...
creates the the rates and distributions for l1 trigger objects
How to use:
1. input the number of bunches in the run (~line 35)
2. change the variables "newConditionsNtuples" and "oldConditionsNtuples" to ntuple paths
3. If good run JSON is not applied during ntuple production, modify isGoodLumiSection()

Optionally, if you want to rescale to a given instantaneous luminosity:
1. input the instantaneous luminosity of the run (~line 32) [only if we scale to 2016 nominal]
2. select whether you rescale to L=1.5e34 (~line606??...) generally have it setup to rescale
nb: for 2&3 I have provided the info in runInfoForRates.txt
*/

// configurable parameters
double numBunch = 1537; //the number of bunches colliding for the run of interest
double runLum = 0.02; // 0.44: 275783  0.58:  276363 //luminosity of the run of interest (*10^34)
double expectedLum = 1.15; //expected luminosity of 2016 runs (*10^34)

void rates(bool newConditions, const std::string& inputFileDirectory);

int main(int argc, char *argv[])
{
  bool newConditions = true;
  std::string ntuplePath("");

  if (argc != 3) {
    std::cout << "Usage: rates.exe [new/def] [path to ntuples]\n"
	      << "[new/def] indicates new or default (existing) conditions" << std::endl;
    exit(1);
  }
  else {
    std::string par1(argv[1]);
    std::transform(par1.begin(), par1.end(), par1.begin(), ::tolower);
    if(par1.compare("new") == 0) newConditions = true;
    else if(par1.compare("def") == 0) newConditions = false;
    else {
      std::cout << "First parameter must be \"new\" or \"def\"" << std::endl;
      exit(1);
    }
    ntuplePath = argv[2];
  }

  rates(newConditions, ntuplePath);

  return 0;
}

// only need to edit this section if good run JSON
// is not used during ntuple production
bool isGoodLumiSection(int lumiBlock)
{
  if (lumiBlock >= 1
      || lumiBlock <= 10000) {
    return true;
  }

  return false;
}

void rates(bool newConditions, const std::string& inputFileDirectory){
  
  bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
  bool emuOn = true;  //are we using data from emulator?

  if (hwOn==false && emuOn==false){
    std::cout << "exiting as neither hardware or emulator selected" << std::endl;
    return;
  }

  std::string inputFile(inputFileDirectory);
  inputFile += "/L1Ntuple_*.root";
  std::string outputDirectory = "emu";  //***runNumber, triggerType, version, hw/emu/both***MAKE SURE IT EXISTS
  std::string outputFilename = "rates_def.root";
  if(newConditions) outputFilename = "rates_new_cond.root";
  TFile* kk = TFile::Open( outputFilename.c_str() , "recreate");
  // if (kk!=0){
  //   cout << "TERMINATE: not going to overwrite file " << outputFilename << endl;
  //   return;
  // }


  // make trees
  std::cout << "Loading up the TChain..." << std::endl;
  TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  if (emuOn){
    treeL1emu->Add(inputFile.c_str());
  }
  TChain * treeL1hw = new TChain("l1UpgradeTree/L1UpgradeTree");
  if (hwOn){
    treeL1hw->Add(inputFile.c_str());
  }
  TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1Towemu->Add(inputFile.c_str());
  }
  TChain * treeL1CaloTPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1CaloTPemu->Add(inputFile.c_str());
  }
  TChain * eventTree = new TChain("l1EventTree/L1EventTree");
  eventTree->Add(inputFile.c_str());

  // In case you want to include PU info
  // TChain * vtxTree = new TChain("l1RecoTree/RecoTree");
  // if(binByPileUp){
  //   vtxTree->Add(inputFile.c_str());
  // }


  TChain * treeL1TPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1TPemu->Add(inputFile.c_str());
  }

  TChain * treeL1TPhw = new TChain("l1CaloTowerTree/L1CaloTowerTree");
  if (hwOn){
    treeL1TPhw->Add(inputFile.c_str());
  }

  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1hw_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1hw->SetBranchAddress("L1Upgrade", &l1hw_);
  L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
  L1Analysis::L1AnalysisCaloTPDataFormat     *l1CaloTPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1CaloTPemu->SetBranchAddress("CaloTP", &l1CaloTPemu_);
  L1Analysis::L1AnalysisEventDataFormat    *event_ = new L1Analysis::L1AnalysisEventDataFormat();
  eventTree->SetBranchAddress("Event", &event_);
  // L1Analysis::L1AnalysisRecoVertexDataFormat    *vtx_ = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  // vtxTree->SetBranchAddress("Vertex", &vtx_);

  L1Analysis::L1AnalysisCaloTPDataFormat    *l1TPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemu->SetBranchAddress("CaloTP", &l1TPemu_);
  L1Analysis::L1AnalysisCaloTPDataFormat    *l1TPhw_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPhw->SetBranchAddress("CaloTP", &l1TPhw_);


  // get number of entries
  Long64_t nentries;
  if (emuOn) nentries = treeL1emu->GetEntries();
  else nentries = treeL1hw->GetEntries();
  int goodLumiEventCount = 0;

  std::string outputTxtFilename = "output_rates/" + outputDirectory + "/extraInfo.txt";
  std::ofstream myfile; // save info about the run, including rates for a given lumi section, and number of events we used.
  myfile.open(outputTxtFilename.c_str());
  eventTree->GetEntry(0);
  myfile << "run number = " << event_->run << std::endl;

  // set parameters for histograms
  // jet bins
  int nJetBins = 400;
  float jetLo = 0.;
  float jetHi = 400.;
  float jetBinWidth = (jetHi-jetLo)/nJetBins;

  // EG bins
  int nEgBins = 300;
  float egLo = 0.;
  float egHi = 300.;
  float egBinWidth = (egHi-egLo)/nEgBins;

  // tau bins
  int nTauBins = 300;
  float tauLo = 0.;
  float tauHi = 300.;
  float tauBinWidth = (tauHi-tauLo)/nTauBins;

  // htSum bins
  int nHtSumBins = 600;
  float htSumLo = 0.;
  float htSumHi = 600.;
  float htSumBinWidth = (htSumHi-htSumLo)/nHtSumBins;

  // mhtSum bins
  int nMhtSumBins = 300;
  float mhtSumLo = 0.;
  float mhtSumHi = 300.;
  float mhtSumBinWidth = (mhtSumHi-mhtSumLo)/nMhtSumBins;

  // etSum bins
  int nEtSumBins = 600;
  float etSumLo = 0.;
  float etSumHi = 600.;
  float etSumBinWidth = (etSumHi-etSumLo)/nEtSumBins;

  // metSum bins
  int nMetSumBins = 300;
  float metSumLo = 0.;
  float metSumHi = 300.;
  float metSumBinWidth = (metSumHi-metSumLo)/nMetSumBins;

  // metHFSum bins
  int nMetHFSumBins = 300;
  float metHFSumLo = 0.;
  float metHFSumHi = 300.;
  float metHFSumBinWidth = (metHFSumHi-metHFSumLo)/nMetHFSumBins;

  // tp bins
  int nTpBins = 100;
  float tpLo = 0.;
  float tpHi = 100.;

  std::string axR = ";Threshold E_{T} (GeV);rate (Hz)";
  std::string axD = ";E_{T} (GeV);events/bin";
  std::string mult = ";Hit Multiplicity;Number of Entries";

  //make histos
  TH1F* singleJetRates_emu = new TH1F("singleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_emu = new TH1F("doubleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_emu = new TH1F("tripleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_emu = new TH1F("quadJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_emu = new TH1F("singleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_emu = new TH1F("doubleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_emu = new TH1F("singleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_emu = new TH1F("doubleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_emu = new TH1F("singleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_emu = new TH1F("doubleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_emu = new TH1F("singleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_emu = new TH1F("doubleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_emu = new TH1F("htSumRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_emu = new TH1F("mhtSumRates_emu",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_emu = new TH1F("etSumRates_emu",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_emu = new TH1F("metSumRates_emu",axR.c_str(), nMetSumBins, metSumLo, metSumHi); 
  TH1F* metHFSumRates_emu = new TH1F("metHFSumRates_emu",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  
  TH1F* singleJetRates_hw = new TH1F("singleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_hw = new TH1F("doubleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_hw = new TH1F("tripleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_hw = new TH1F("quadJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_hw = new TH1F("singleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_hw = new TH1F("doubleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_hw = new TH1F("singleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_hw = new TH1F("doubleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_hw = new TH1F("singleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_hw = new TH1F("doubleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_hw = new TH1F("singleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_hw = new TH1F("doubleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_hw = new TH1F("htSumRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_hw = new TH1F("mhtSumRates_hw",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_hw = new TH1F("etSumRates_hw",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_hw = new TH1F("metSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  TH1F* metHFSumRates_hw = new TH1F("metHFSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 

  TH1F* hcalTP_emu = new TH1F("hcalTP_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_emu = new TH1F("ecalTP_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1F* hcalTP_hw = new TH1F("hcalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_hw = new TH1F("ecalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1D * hJetEt = new TH1D("jetET",";ET;",100,0,1000);

  // 3 GeV energy cuts, scanning time cuts
  // inclusive
  TH1F * dt3GeV1nsMult_emu = new TH1F("dt3GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsMult_emu = new TH1F("dt3GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsMult_emu = new TH1F("dt3GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsMult_emu = new TH1F("dt3GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsMult_emu = new TH1F("dt3GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  // HE
  TH1F * dt3GeV1nsHEMult_emu = new TH1F("dt3GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHEMult_emu = new TH1F("dt3GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHEMult_emu = new TH1F("dt3GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHEMult_emu = new TH1F("dt3GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHEMult_emu = new TH1F("dt3GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  // HB
  TH1F * dt3GeV1nsHBMult_emu = new TH1F("dt3GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHBMult_emu = new TH1F("dt3GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBMult_emu = new TH1F("dt3GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHBMult_emu = new TH1F("dt3GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHBMult_emu = new TH1F("dt3GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt3GeV1nsJetMult_emu = new TH1F("dt3GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsJetMult_emu = new TH1F("dt3GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsJetMult_emu = new TH1F("dt3GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsJetMult_emu = new TH1F("dt3GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsJetMult_emu = new TH1F("dt3GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // HE
  TH1F * dt3GeV1nsHEJetMult_emu = new TH1F("dt3GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHEJetMult_emu = new TH1F("dt3GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHEJetMult_emu = new TH1F("dt3GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHEJetMult_emu = new TH1F("dt3GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHEJetMult_emu = new TH1F("dt3GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // HB
  TH1F * dt3GeV1nsHBJetMult_emu = new TH1F("dt3GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHBJetMult_emu = new TH1F("dt3GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBJetMult_emu = new TH1F("dt3GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHBJetMult_emu = new TH1F("dt3GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHBJetMult_emu = new TH1F("dt3GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // 2 GeV energy cuts, scanning time cut
  // inclusive
  TH1F * dt2GeV1nsMult_emu = new TH1F("dt2GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsMult_emu = new TH1F("dt2GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsMult_emu = new TH1F("dt2GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsMult_emu = new TH1F("dt2GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsMult_emu = new TH1F("dt2GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  // HE
  TH1F * dt2GeV1nsHEMult_emu = new TH1F("dt2GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHEMult_emu = new TH1F("dt2GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHEMult_emu = new TH1F("dt2GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHEMult_emu = new TH1F("dt2GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHEMult_emu = new TH1F("dt2GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  // HB
  TH1F * dt2GeV1nsHBMult_emu = new TH1F("dt2GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHBMult_emu = new TH1F("dt2GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHBMult_emu = new TH1F("dt2GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHBMult_emu = new TH1F("dt2GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHBMult_emu = new TH1F("dt2GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt2GeV1nsJetMult_emu = new TH1F("dt2GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsJetMult_emu = new TH1F("dt2GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsJetMult_emu = new TH1F("dt2GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsJetMult_emu = new TH1F("dt2GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsJetMult_emu = new TH1F("dt2GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // HE
  TH1F * dt2GeV1nsHEJetMult_emu = new TH1F("dt2GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHEJetMult_emu = new TH1F("dt2GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHEJetMult_emu = new TH1F("dt2GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHEJetMult_emu = new TH1F("dt2GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHEJetMult_emu = new TH1F("dt2GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // HB
  TH1F * dt2GeV1nsHBJetMult_emu = new TH1F("dt2GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHBJetMult_emu = new TH1F("dt2GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHBJetMult_emu = new TH1F("dt2GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHBJetMult_emu = new TH1F("dt2GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHBJetMult_emu = new TH1F("dt2GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // 1 GeV energy cuts, scanning time cuts
  // inclusive
  TH1F * dt1GeV1nsMult_emu = new TH1F("dt1GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsMult_emu = new TH1F("dt1GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsMult_emu = new TH1F("dt1GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsMult_emu = new TH1F("dt1GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsMult_emu = new TH1F("dt1GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  // HE
  TH1F * dt1GeV1nsHEMult_emu = new TH1F("dt1GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHEMult_emu = new TH1F("dt1GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHEMult_emu = new TH1F("dt1GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHEMult_emu = new TH1F("dt1GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHEMult_emu = new TH1F("dt1GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  // HB
  TH1F * dt1GeV1nsHBMult_emu = new TH1F("dt1GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHBMult_emu = new TH1F("dt1GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHBMult_emu = new TH1F("dt1GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHBMult_emu = new TH1F("dt1GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHBMult_emu = new TH1F("dt1GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt1GeV1nsJetMult_emu = new TH1F("dt1GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsJetMult_emu = new TH1F("dt1GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsJetMult_emu = new TH1F("dt1GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsJetMult_emu = new TH1F("dt1GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsJetMult_emu = new TH1F("dt1GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  // HE
  TH1F * dt1GeV1nsHEJetMult_emu = new TH1F("dt1GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHEJetMult_emu = new TH1F("dt1GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHEJetMult_emu = new TH1F("dt1GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHEJetMult_emu = new TH1F("dt1GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHEJetMult_emu = new TH1F("dt1GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  // HB
  TH1F * dt1GeV1nsHBJetMult_emu = new TH1F("dt1GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHBJetMult_emu = new TH1F("dt1GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHBJetMult_emu = new TH1F("dt1GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHBJetMult_emu = new TH1F("dt1GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHBJetMult_emu = new TH1F("dt1GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);

  // making TH2F for the energy depth plots
  TH2F * Energy_Depth = new TH2F("Energy_Depth", "TP Energy Fraction vs. Depth;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_Depth = new TH2F("Timing_Depth", "TP Timing Value vs. Depth;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE = new TH2F("Energy_DepthHE", "TP Energy Fraction vs. Depth in HE;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHE = new TH2F("Timing_DepthHE", "TP Timing Value vs. Depth in HE;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB = new TH2F("Energy_DepthHB", "TP Energy Fraction vs. Depth in HB;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHB = new TH2F("Timing_DepthHB", "TP Timing Value vs. Depth in HB;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // making TH2F for the energy depth plots for high energy TPs
  TH2F * Energy_Depth_HighE = new TH2F("Energy_Depth_HighE", "TP Energy Fraction vs. Depth for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_Depth_HighE = new TH2F("Timing_Depth_HighE", "TP Timing Value vs. Depth for TP E_{T} > 10 GeV;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE_HighE = new TH2F("Energy_DepthHE_HighE", "TP Energy Fraction vs. Depth in HE for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_DepthHE_HighE = new TH2F("Timing_DepthHE_HighE", "TP Timing Value vs. Depth in HE for TP E_{T} > 10 GeV;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB_HighE = new TH2F("Energy_DepthHB_HighE", "TP Energy Fraction vs. Depth in HB for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_DepthHB_HighE = new TH2F("Timing_DepthHB_HighE", "TP Timing Value vs. Depth in HB for TP E_{T} > 10 GeV;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // TH2F for energy depth, where matched to jets
  TH2F * Energy_Depth_Jets = new TH2F("Energy_Depth_Jets", "TP Energy Fraction vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_Depth_Jets = new TH2F("Timing_Depth_Jets", "TP Timing Value vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE_Jets = new TH2F("Energy_DepthHE_Jets", "TP Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHE_Jets = new TH2F("Timing_DepthHE_Jets", "TP Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB_Jets = new TH2F("Energy_DepthHB_Jets", "TP Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHB_Jets = new TH2F("Timing_DepthHB_Jets", "TP Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // TH2F for energy depth, where matched to jets for high energy TPs
  TH2F * Energy_Depth_Jets_HighE = new TH2F("Energy_Depth_Jets_HighE", "TP Energy Fraction vs. Depth for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_Depth_Jets_HighE = new TH2F("Timing_Depth_Jets_HighE", "TP Timing Value vs. Depth for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE_Jets_HighE = new TH2F("Energy_DepthHE_Jets_HighE", "TP Energy Fraction vs. Depth in HE for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_DepthHE_Jets_HighE = new TH2F("Timing_DepthHE_Jets_HighE", "TP Timing Value vs. Depth in HE for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB_Jets_HighE = new TH2F("Energy_DepthHB_Jets_HighE", "TP Energy Fraction vs. Depth in HB for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  //  TH2F * Timing_DepthHB_Jets_HighE = new TH2F("Timing_DepthHB_Jets_HighE", "TP Timing Value vs. Depth in HB for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // making TH1D for the ProfileX() from the TH2F of energy_depth or timing_depth plots
  TH1D * Energy_Depth_avg = new TH1D("Energy_Depth_avg", "TP Avg Energy Fraction vs. Depth;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_Depth_avg = new TH1D("Timing_Depth_avg", "TP Avg Timing Value vs. Depth;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHE_avg = new TH1D("Energy_DepthHE_avg", "TP Avg Energy Fraction vs. Depth in HE;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHE_avg = new TH1D("Timing_DepthHE_avg", "TP Avg Timing Value vs. Depth in HE;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg = new TH1D("Energy_DepthHB_avg", "TP Avg Energy Fraction vs. Depth in HB;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHB_avg = new TH1D("Timing_DepthHB_avg", "TP Avg Timing Value vs. Depth in HB;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  // TH1D for the ProfileX(), where matched to jets
  TH1D * Energy_Depth_avg_Jets = new TH1D("Energy_Depth_avg_Jets", "TP Avg Energy Fraction vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_Depth_avg_Jets = new TH1D("Timing_Depth_avg_Jets", "TP Avg Timing Value vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHE_avg_Jets = new TH1D("Energy_DepthHE_avg_Jets", "TP Avg Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHE_avg_Jets = new TH1D("Timing_DepthHE_avg_Jets", "TP Avg Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg_Jets = new TH1D("Energy_DepthHB_avg_Jets", "TP Avg Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHB_avg_Jets = new TH1D("Timing_DepthHB_avg_Jets", "TP Avg Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  // Ratio of energy in HCAL depth layers
  TH1F * Ratio_Depth = new TH1F("Ratio_Depth", "Ratio of First 2 HCAL Layers to E_{T};Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE = new TH1F("Ratio_DepthHE", "Ratio of First 2 HCAL Layers to E_{T} in HE;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB = new TH1F("Ratio_DepthHB", "Ratio of First 2 HCAL Layers to E_{T} in HB;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_Depth_Jets = new TH1F("Ratio_Depth_Jets", "Ratio of First 2 HCAL Layers to E_{T}, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE_Jets = new TH1F("Ratio_DepthHE_Jets", "Ratio of First 2 HCAL Layers to E_{T} in HE, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB_Jets = new TH1F("Ratio_DepthHB_Jets", "Ratio of First 2 HCAL Layers to E_{T} in HB, matched w/Jets;Ratio;Number of Events", 50,0,1);

  // timing values for center of barrel
  TH1F * centralTiming = new TH1F("centralTiming", "Time of arrival - TOF (central barrel iEta);Time (ns);Number of Events",50,-10,40);
  // HCAL / ECAL+HCAL energy to check cut for rates plots
  std::vector<TString> ratioStrings = {"HOvE","HOvE3","HOvE9","H3OvE3","H9OvE9"};
  TH1F * HoverEtotal = new TH1F("HoverEtotal", "HCAL energy / HCAL+ECAL energy;H/E;Number of Events",50,0,1);

  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%10000)==0) std::cout << "Done " << jentry  << " events of " << nentries << std::endl;

    //lumi break clause
    eventTree->GetEntry(jentry);
    //skip the corresponding event
    if (!isGoodLumiSection(event_->lumi)) continue;
    goodLumiEventCount++;

    //do routine for L1 emulator quantites
    if (emuOn){

      treeL1TPemu->GetEntry(jentry);
      treeL1Towemu->GetEntry(jentry);
      treeL1CaloTPemu->GetEntry(jentry);
      double tpEt(0.);
      
      for(int i=0; i < l1TPemu_->nHCALTP; i++){
	tpEt = l1TPemu_->hcalTPet[i];
	hcalTP_emu->Fill(tpEt);
      }
      for(int i=0; i < l1TPemu_->nECALTP; i++){
	tpEt = l1TPemu_->ecalTPet[i];
	ecalTP_emu->Fill(tpEt);
      }

      treeL1emu->GetEntry(jentry);
      // get jetEt*, egEt*, tauEt, htSum, mhtSum, etSum, metSum
      // ALL EMU OBJECTS HAVE BX=0...
      double jetEt_1 = 0;
      double jetEt_2 = 0;
      double jetEt_3 = 0;
      double jetEt_4 = 0;
      if (l1emu_->nJets>0) jetEt_1 = l1emu_->jetEt[0];
      if (l1emu_->nJets>1) jetEt_2 = l1emu_->jetEt[1];
      if (l1emu_->nJets>2) jetEt_3 = l1emu_->jetEt[2];
      if (l1emu_->nJets>3) jetEt_4 = l1emu_->jetEt[3];       
      double egEt_1 = 0;
      double egEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1emu_->nEGs; c++){
        if (l1emu_->egEt[c] > egEt_1){
          egEt_2 = egEt_1;
          egEt_1 = l1emu_->egEt[c];
        }
        else if (l1emu_->egEt[c] <= egEt_1 && l1emu_->egEt[c] > egEt_2){
          egEt_2 = l1emu_->egEt[c];
        }
      }

      double tauEt_1 = 0;
      double tauEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1emu_->nTaus; c++){
        if (l1emu_->tauEt[c] > tauEt_1){
          tauEt_2 = tauEt_1;
          tauEt_1 = l1emu_->tauEt[c];
        }
        else if (l1emu_->tauEt[c] <= tauEt_1 && l1emu_->tauEt[c] > tauEt_2){
          tauEt_2 = l1emu_->tauEt[c];
        }
      }

      double egISOEt_1 = 0;
      double egISOEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1emu_->nEGs; c++){
        if (l1emu_->egEt[c] > egISOEt_1 && l1emu_->egIso[c]==1){
          egISOEt_2 = egISOEt_1;
          egISOEt_1 = l1emu_->egEt[c];
        }
        else if (l1emu_->egEt[c] <= egISOEt_1 && l1emu_->egEt[c] > egISOEt_2 && l1emu_->egIso[c]==1){
          egISOEt_2 = l1emu_->egEt[c];
        }
      }

      double tauISOEt_1 = 0;
      double tauISOEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1emu_->nTaus; c++){
        if (l1emu_->tauEt[c] > tauISOEt_1 && l1emu_->tauIso[c]>0){
          tauISOEt_2 = tauISOEt_1;
          tauISOEt_1 = l1emu_->tauEt[c];
        }
        else if (l1emu_->tauEt[c] <= tauISOEt_1 && l1emu_->tauEt[c] > tauISOEt_2 && l1emu_->tauIso[c]>0){
          tauISOEt_2 = l1emu_->tauEt[c];
        }
      }

      double htSum(0.0);
      double mhtSum(0.0);
      double etSum(0.0);
      double metSum(0.0);
      double metHFSum(0.0);
      for (unsigned int c=0; c<l1emu_->nSums; c++){
          if( l1emu_->sumBx[c] != 0 ) continue;
          if( l1emu_->sumType[c] == L1Analysis::kTotalEt ) etSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kTotalHt ) htSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kMissingEt ) metSum = l1emu_->sumEt[c];
	  if( l1emu_->sumType[c] == L1Analysis::kMissingEtHF ) metHFSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kMissingHt ) mhtSum = l1emu_->sumEt[c];
      }

      int seedTowerIEta(-1);
      int seedTowerIPhi(-1);
      int nDepth(-1);
      double nCaloTPemu(0), tpEtaemu(0), tpPhiemu(0), tpEtemu(0);
      uint nJetemu(0);
      nCaloTPemu = l1CaloTPemu_->nHCALTP;
      nJetemu = l1emu_->nJets;
      //      std::cout << nJetemu << std::endl;
      // hcalTPdepth and hcalTPtiming will store the timing and depth variables from the 7 HCAL layers
      double hcalTPdepth[7] = {0};
      double hcalTPtiming[7] = {0};
      std::map<const TString, std::vector<double> > TimingVariablesAllJets;
      std::map<const TString, std::vector<double> > DepthVariablesAllJets;
      // multiplicity for all HCAL TPs in an entry
      double mult3GeV1ns(0), mult3GeV2ns(0), mult3GeV3ns(0), mult3GeV4ns(0), mult3GeV5ns(0);
      double mult3GeV1nsHE(0), mult3GeV2nsHE(0), mult3GeV3nsHE(0), mult3GeV4nsHE(0), mult3GeV5nsHE(0);
      double mult3GeV1nsHB(0), mult3GeV2nsHB(0), mult3GeV3nsHB(0), mult3GeV4nsHB(0), mult3GeV5nsHB(0);
      double mult2GeV1ns(0), mult2GeV2ns(0), mult2GeV3ns(0), mult2GeV4ns(0), mult2GeV5ns(0);
      double mult2GeV1nsHE(0), mult2GeV2nsHE(0), mult2GeV3nsHE(0), mult2GeV4nsHE(0), mult2GeV5nsHE(0);
      double mult2GeV1nsHB(0), mult2GeV2nsHB(0), mult2GeV3nsHB(0), mult2GeV4nsHB(0), mult2GeV5nsHB(0);
      double mult1GeV1ns(0), mult1GeV2ns(0), mult1GeV3ns(0), mult1GeV4ns(0), mult1GeV5ns(0);
      double mult1GeV1nsHE(0), mult1GeV2nsHE(0), mult1GeV3nsHE(0), mult1GeV4nsHE(0), mult1GeV5nsHE(0);
      double mult1GeV1nsHB(0), mult1GeV2nsHB(0), mult1GeV3nsHB(0), mult1GeV4nsHB(0), mult1GeV5nsHB(0);
      // and multiplicity for when HCAL TP is matched with Jets
      double mult3GeV1ns_Jets(0), mult3GeV2ns_Jets(0), mult3GeV3ns_Jets(0), mult3GeV4ns_Jets(0), mult3GeV5ns_Jets(0);
      double mult3GeV1nsHE_Jets(0), mult3GeV2nsHE_Jets(0), mult3GeV3nsHE_Jets(0), mult3GeV4nsHE_Jets(0), mult3GeV5nsHE_Jets(0);
      double mult3GeV1nsHB_Jets(0), mult3GeV2nsHB_Jets(0), mult3GeV3nsHB_Jets(0), mult3GeV4nsHB_Jets(0), mult3GeV5nsHB_Jets(0);
      double mult2GeV1ns_Jets(0), mult2GeV2ns_Jets(0), mult2GeV3ns_Jets(0), mult2GeV4ns_Jets(0), mult2GeV5ns_Jets(0);
      double mult2GeV1nsHE_Jets(0), mult2GeV2nsHE_Jets(0), mult2GeV3nsHE_Jets(0), mult2GeV4nsHE_Jets(0), mult2GeV5nsHE_Jets(0);
      double mult2GeV1nsHB_Jets(0), mult2GeV2nsHB_Jets(0), mult2GeV3nsHB_Jets(0), mult2GeV4nsHB_Jets(0), mult2GeV5nsHB_Jets(0);
      double mult1GeV1ns_Jets(0), mult1GeV2ns_Jets(0), mult1GeV3ns_Jets(0), mult1GeV4ns_Jets(0), mult1GeV5ns_Jets(0);
      double mult1GeV1nsHE_Jets(0), mult1GeV2nsHE_Jets(0), mult1GeV3nsHE_Jets(0), mult1GeV4nsHE_Jets(0), mult1GeV5nsHE_Jets(0);
      double mult1GeV1nsHB_Jets(0), mult1GeV2nsHB_Jets(0), mult1GeV3nsHB_Jets(0), mult1GeV4nsHB_Jets(0), mult1GeV5nsHB_Jets(0);
      
      // for H/E calculation 
      std::map<const TString, std::vector<double> > hadVariablesAllJets;
      std::map<const TString, std::vector<double> > emVariablesAllJets;
      int maxTowerEndcap = 28;
      int maxTowerBarrel = 16;
      int minTowerForHOvE = maxTowerBarrel+1;
      int maxTowerForHOvE = maxTowerEndcap;
      double towEtemu(0), towHademu(0), towEmemu(0), towEtaemu(0), towPhiemu(0), nTowemu(0);
      nTowemu = l1Towemu_->nTower;

      // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
      for(uint jetIt=0; jetIt < nJetemu && jetIt < 4; jetIt++){
	hJetEt->Fill(l1emu_->jetEt[jetIt]); // these are already in order of highest E_T
	seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
	seedTowerIEta = l1emu_->jetTowerIEta[jetIt];

	// calculate the HCAL / ECAL+HCAL energy to determine the cut for the jet rates plots. From Matthew's code
	double seedTowerHad(0), seedTowerEm(0), seedTower3x3Em(0), seedTower3x3Had(0), seedTower9x9Em(0), seedTower9x9Had(0);
	for (int towIt = 0; towIt < nTowemu; towIt++){
	  towEtemu  = l1Towemu_->iet[towIt];
	  towHademu = l1Towemu_->ihad[towIt];
	  towEmemu  = l1Towemu_->iem[towIt];
	  towEtaemu = l1Towemu_->ieta[towIt];
	  towPhiemu = l1Towemu_->iphi[towIt];
	  if (abs(towEtaemu) >= minTowerForHOvE && abs(towEtaemu) <= maxTowerForHOvE){
	    if (towEtaemu == seedTowerIEta && towPhiemu == seedTowerIPhi){
	      seedTowerHad = towHademu;
	      seedTowerEm = towEmemu;
	    }
	    for (int iSeedTowerIEta = -4; iSeedTowerIEta <= 4; ++iSeedTowerIEta){
	      for (int iSeedTowerIPhi = -4; iSeedTowerIPhi <= 4; ++iSeedTowerIPhi){
		int wrappedIPhi = seedTowerIPhi+iSeedTowerIPhi;
		if (wrappedIPhi > 72) wrappedIPhi -= 72;
		if (wrappedIPhi < 0) wrappedIPhi += 72;
		if (towEtaemu == seedTowerIEta+iSeedTowerIEta && towPhiemu == wrappedIPhi){
		  seedTower9x9Em += towEmemu;
		  seedTower9x9Had += towHademu;
		  if (abs(iSeedTowerIPhi) <= 1 && abs(iSeedTowerIEta) <= 1){
		    seedTower3x3Em += towEmemu;
		    seedTower3x3Had += towHademu;
		  }
		}
	      }
	    }
	  } // closing min max tower statement
	} // closing seed tower loop
	hadVariablesAllJets["HOvE"].push_back(seedTowerHad);
	hadVariablesAllJets["HOvE3"].push_back(seedTowerHad);
	hadVariablesAllJets["HOvE9"].push_back(seedTowerHad);
	hadVariablesAllJets["H3OvE3"].push_back(seedTower3x3Had);
	hadVariablesAllJets["H9OvE9"].push_back(seedTower9x9Had);

	emVariablesAllJets["HOvE"].push_back(seedTowerEm);
	emVariablesAllJets["HOvE3"].push_back(seedTower3x3Em);
	emVariablesAllJets["HOvE9"].push_back(seedTower9x9Em);
	emVariablesAllJets["H3OvE3"].push_back(seedTower3x3Em);
	emVariablesAllJets["H9OvE9"].push_back(seedTower9x9Em);


	// Still in the jets loop, now match HCAL TPs to L1 Jets for multiplicity studies

	if (jetIt != 0 ) continue; // only do matching to the highest energy jet
	if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to HCAL TP
	// loop over HCAL TPs to find ones that match with L1 Jet
	double maxE = 0;
	int maxE_HcalTPIt = 0;
	int maxE_iEta = 50;
	int maxE_iPhi = -1;
	// loop over HCAL TPs, and this is only for the first four L1 Jets (since these are the highest energy)
	for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	  tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions
	  tpPhiemu = l1CaloTPemu_->hcalTPiphi[HcalTPIt];
	  tpEtemu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots
	  nDepth = l1CaloTPemu_->hcalTPnDepths[HcalTPIt];

	  if (nDepth == 0) continue; // skipping events where depth = 0, since here timing = -1 and energy = 0 (invalid event)     
	 
	  // convert ieta, iphi to physical eta, phi for delta R matching to L1 Jet. From https://github.com/gk199/cms-hcal-debug/blob/PulseShape/plugins/HcalCompareUpgradeChains.cc#L915-L956
	  double Jet_eta;
	  double TP_eta;
	  if (abs(seedTowerIEta) >= 24) {
	    Jet_eta = .1695*seedTowerIEta - 1.9931*(seedTowerIEta/(abs(seedTowerIEta)));
	  }
	  else {
	    Jet_eta = .0875*seedTowerIEta - 0.0489*(seedTowerIEta/(abs(seedTowerIEta)));
	  }
          if (abs(tpEtaemu) >= 24) {
            TP_eta = .1695*tpEtaemu - 1.9931*(tpEtaemu/(abs(tpEtaemu)));
          }
          else {
            TP_eta = .0875*tpEtaemu - 0.0489*(tpEtaemu/(abs(tpEtaemu)));
          }
	  double Jet_phi;
	  double TP_phi;
	  Jet_phi = double(seedTowerIPhi)*(2.*TMath::Pi()/72);
	  if (seedTowerIPhi > 36) Jet_phi -= 2.*TMath::Pi();
          TP_phi = double(tpPhiemu)*(2.*TMath::Pi()/72);
          if (tpPhiemu > 36) TP_phi -= 2.*TMath::Pi();

	  // Delta R matching in a cone, based on converted ieta, iphi values from above to physical eta, phi values
	  double DeltaEta = Jet_eta - TP_eta;
	  double DeltaPhi = Jet_phi - TP_phi;
	  if (sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi) > 0.5 ) continue;
	  
	  // save the iterator where the TP that is matched to L1 Jet has the max energy of all matched TPs. Also save the ieta, iphi, max energy, and iterator position of this TP
	  if ( l1CaloTPemu_->hcalTPet[HcalTPIt] > maxE ) {
	    maxE = l1CaloTPemu_->hcalTPet[HcalTPIt];
	    maxE_HcalTPIt = HcalTPIt;
	    maxE_iEta = l1CaloTPemu_->hcalTPieta[HcalTPIt];
	    maxE_iPhi = l1CaloTPemu_->hcalTPiphi[HcalTPIt];
	  }
	  //	} // closing the TP loop -- move if using more than 1 HCAL TP
	  // if using more than one HCAL TP, change maxE_HcalTPIt to HcalTPIt
	  //	  if (maxE == 0) continue;

	  // Energy deposited in each depth layer for every HCAL TP (4 in HB, 7 in HE)  
	  hcalTPdepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	  hcalTPdepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	  hcalTPdepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	  hcalTPdepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	  hcalTPdepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	  hcalTPdepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	  hcalTPdepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	  // timing info for each layer, in 25 ns with resolution 0.5 ns 
	  hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];
	  hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];
	  hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];
	  hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];
	  hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];
	  hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];
	  hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];

	  /*
	  // print outs for confirming and troubleshooting energy depth
	  std::cout << "Info for when HCAL TPs are matched to L1 Jets with Delta R < 0.5 " << std::endl;
	  std::cout << "Jet iterator = " << jetIt << " and jet energy = " << l1emu_->jetEt[jetIt] << std::endl;
	  std::cout << "iPhi values = " << seedTowerIPhi << " and " << maxE_iPhi << std::endl;
	  std::cout << "iEta values = " << seedTowerIEta << " and " << maxE_iEta << std::endl;
	  std::cout << "HCAL TP iterator = " << maxE_HcalTPIt << " and TP energy = " << maxE << std::endl;
	  std::cout << hcalTPdepth[0] << ", " <<hcalTPdepth[1] << ", " << hcalTPdepth[2] << ", " << hcalTPdepth[3] << ", " << hcalTPdepth[4] << ", " << hcalTPdepth[5] << ", " << hcalTPdepth[6] << std::endl;
	  std::cout << " " << std::endl;
	  */
	
	  // filling energy and time plots for each of 7 HCAL depths.
	  for (int i = 0; i < 7; i++){
	    Energy_Depth_Jets->Fill(i+1,hcalTPdepth[i]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer
	    Timing_Depth_Jets->Fill(i+1,hcalTPtiming[i]); // raw timing value in each layer
	    if (tpEtemu > 10 ) {
	      Energy_Depth_Jets_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu); // restricting to high energy HCAL TPs
	    }
	    if (abs(tpEtaemu) < 16) {
	      Energy_DepthHB_Jets->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	      Timing_DepthHB_Jets->Fill(i+1,hcalTPtiming[i]);
	      if (tpEtemu > 10 ) {
		Energy_DepthHB_Jets_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	      }
	    }
	    if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	      Energy_DepthHE_Jets->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	      Timing_DepthHE_Jets->Fill(i+1,hcalTPtiming[i]);
	      if (tpEtemu > 10 ) {
		Energy_DepthHE_Jets_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	      }
	    }
	  }

	  // Ratio of energy in first two HCAL layers to all HCAL layers. Only consider for high energy TPs > 10 GeV	
	  if ( tpEtemu > 10 ) {
	    Ratio_Depth_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    if (abs(tpEtaemu) < 16) {
	      Ratio_DepthHB_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    }
	    if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	      Ratio_DepthHE_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    }
	  }
	  
	  // loop over HCAL depths for every HCAL TP
	  for (int depthIt = 0; depthIt < 7; depthIt++){
	    // count multiplicity of layers given a timing and energy threshold   
	    // 3 GeV energy cut
	    if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1){
	      mult3GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult3GeV2ns_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult3GeV3ns_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult3GeV4ns_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult3GeV5ns_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    // 3 GeV HB HE regions
	    if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	      mult3GeV1nsHB_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult3GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult3GeV3nsHB_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult3GeV4nsHB_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult3GeV5nsHB_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	      mult3GeV1nsHE_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult3GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult3GeV3nsHE_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult3GeV4nsHE_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult3GeV5nsHE_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    // 2 GeV energy cut
	    if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1){
	      mult2GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult2GeV2ns_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult2GeV3ns_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult2GeV4ns_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult2GeV5ns_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    // 2 GeV HB HE regions                                
	    if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	      mult2GeV1nsHB_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult2GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult2GeV3nsHB_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult2GeV4nsHB_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult2GeV5nsHB_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	      mult2GeV1nsHE_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult2GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult2GeV3nsHE_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult2GeV4nsHE_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult2GeV5nsHE_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    // 1 GeV energy cut
	    if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1){
	      mult1GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult1GeV2ns_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult1GeV3ns_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult1GeV4ns_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult1GeV5ns_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    // 1 GeV HB HE regions                                                      
	    if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	      mult1GeV1nsHB_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult1GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult1GeV3nsHB_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult1GeV4nsHB_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult1GeV5nsHB_Jets += 1;
		    }
		  }
		}
	      }
	    }
	    if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	      mult1GeV1nsHE_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2){
		mult1GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3){
		  mult1GeV3nsHE_Jets += 1;
		  if (hcalTPtiming[depthIt] > 4){
		    mult1GeV4nsHE_Jets += 1;
		    if (hcalTPtiming[depthIt] > 5){
		      mult1GeV5nsHE_Jets += 1;
		    }
		  }
		}
	      }
	    }
	  }// closing HCAL depths loop
	} // closing HCAL TP loop (all TPs)
      } // closing L1 Jets loop
      // after HCAL depth loop and L1 Jet loop fill histograms with multiplicity variables. Multiplicity counter reset on each loop iteration. These are for where HCAL TP is matched to the L1 Jet

      // 3 GeV histograms
      dt3GeV1nsJetMult_emu->Fill(mult3GeV1ns_Jets);
      dt3GeV2nsJetMult_emu->Fill(mult3GeV2ns_Jets);
      dt3GeV3nsJetMult_emu->Fill(mult3GeV3ns_Jets);
      dt3GeV4nsJetMult_emu->Fill(mult3GeV4ns_Jets);
      dt3GeV5nsJetMult_emu->Fill(mult3GeV5ns_Jets);
      dt3GeV1nsHEJetMult_emu->Fill(mult3GeV1nsHE_Jets);
      dt3GeV2nsHEJetMult_emu->Fill(mult3GeV2nsHE_Jets);
      dt3GeV3nsHEJetMult_emu->Fill(mult3GeV3nsHE_Jets);
      dt3GeV4nsHEJetMult_emu->Fill(mult3GeV4nsHE_Jets);
      dt3GeV5nsHEJetMult_emu->Fill(mult3GeV5nsHE_Jets);
      dt3GeV1nsHBJetMult_emu->Fill(mult3GeV1nsHB_Jets);
      dt3GeV2nsHBJetMult_emu->Fill(mult3GeV2nsHB_Jets);
      dt3GeV3nsHBJetMult_emu->Fill(mult3GeV3nsHB_Jets);
      dt3GeV4nsHBJetMult_emu->Fill(mult3GeV4nsHB_Jets);
      dt3GeV5nsHBJetMult_emu->Fill(mult3GeV5nsHB_Jets);
      // 2 GeV histograms
      dt2GeV1nsJetMult_emu->Fill(mult2GeV1ns_Jets);
      dt2GeV2nsJetMult_emu->Fill(mult2GeV2ns_Jets);
      dt2GeV3nsJetMult_emu->Fill(mult2GeV3ns_Jets);
      dt2GeV4nsJetMult_emu->Fill(mult2GeV4ns_Jets);
      dt2GeV5nsJetMult_emu->Fill(mult2GeV5ns_Jets);
      dt2GeV1nsHEJetMult_emu->Fill(mult2GeV1nsHE_Jets);
      dt2GeV2nsHEJetMult_emu->Fill(mult2GeV2nsHE_Jets);
      dt2GeV3nsHEJetMult_emu->Fill(mult2GeV3nsHE_Jets);
      dt2GeV4nsHEJetMult_emu->Fill(mult2GeV4nsHE_Jets);
      dt2GeV5nsHEJetMult_emu->Fill(mult2GeV5nsHE_Jets);
      dt2GeV1nsHBJetMult_emu->Fill(mult2GeV1nsHB_Jets);
      dt2GeV2nsHBJetMult_emu->Fill(mult2GeV2nsHB_Jets);
      dt2GeV3nsHBJetMult_emu->Fill(mult2GeV3nsHB_Jets);
      dt2GeV4nsHBJetMult_emu->Fill(mult2GeV4nsHB_Jets);
      dt2GeV5nsHBJetMult_emu->Fill(mult2GeV5nsHB_Jets);
      // 1 GeV histograms
      dt1GeV1nsJetMult_emu->Fill(mult1GeV1ns_Jets);
      dt1GeV2nsJetMult_emu->Fill(mult1GeV2ns_Jets);
      dt1GeV3nsJetMult_emu->Fill(mult1GeV3ns_Jets);
      dt1GeV4nsJetMult_emu->Fill(mult1GeV4ns_Jets);
      dt1GeV5nsJetMult_emu->Fill(mult1GeV5ns_Jets);
      dt1GeV1nsHEJetMult_emu->Fill(mult1GeV1nsHE_Jets);
      dt1GeV2nsHEJetMult_emu->Fill(mult1GeV2nsHE_Jets);
      dt1GeV3nsHEJetMult_emu->Fill(mult1GeV3nsHE_Jets);
      dt1GeV4nsHEJetMult_emu->Fill(mult1GeV4nsHE_Jets);
      dt1GeV5nsHEJetMult_emu->Fill(mult1GeV5nsHE_Jets);
      dt1GeV1nsHBJetMult_emu->Fill(mult1GeV1nsHB_Jets);
      dt1GeV2nsHBJetMult_emu->Fill(mult1GeV2nsHB_Jets);
      dt1GeV3nsHBJetMult_emu->Fill(mult1GeV3nsHB_Jets);
      dt1GeV4nsHBJetMult_emu->Fill(mult1GeV4nsHB_Jets);
      dt1GeV5nsHBJetMult_emu->Fill(mult1GeV5nsHB_Jets);

      // H/E for the first jet information
      HoverEtotal->Fill((hadVariablesAllJets["H3OvE3"][0])/(hadVariablesAllJets["H3OvE3"][0]+emVariablesAllJets["H3OvE3"][0]));

      // HCAL TP information when TPs are not matched to L1 Jets
      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	tpPhiemu = l1CaloTPemu_->hcalTPiphi[HcalTPIt];
	tpEtemu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	nDepth = l1CaloTPemu_->hcalTPnDepths[HcalTPIt];
	
	if (nDepth == 0) continue; // skipping events where depth = 0, since here timing = -1 and energy = 0 (invalid event)

	// Energy deposited in each depth layer for every HCAL TP (4 in HB, 7 in HE)  
	hcalTPdepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	hcalTPdepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	hcalTPdepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	hcalTPdepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	hcalTPdepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	hcalTPdepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	hcalTPdepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	// timing info for each layer, in 25 ns with resolution 0.5 ns 
	hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];
	hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];
	hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];
	hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];
	hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];
	hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];
	hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];

	for (int i = 0; i < 4; i++ ) {
	  if (tpEtaemu == 1 && hcalTPtiming[i] > -0.5 ) {
	    centralTiming->Fill(hcalTPtiming[i]-5.5);
	  }
	}

        // filling energy and time plots for each of 7 HCAL depths  
	for (int i = 0; i < 7; i++){
	  Energy_Depth->Fill(i+1,hcalTPdepth[i]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer      
	  Timing_Depth->Fill(i+1,hcalTPtiming[i]); // raw timing value in each layer  
          if (tpEtemu > 10 ) {
            Energy_Depth_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu); // energy depth for high energy HCAL TPs
          }
	  if (abs(tpEtaemu) < 16) {
	    Energy_DepthHB->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    Timing_DepthHB->Fill(i+1,hcalTPtiming[i]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHB_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    }
	  }
	  if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	    Energy_DepthHE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    Timing_DepthHE->Fill(i+1,hcalTPtiming[i]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHE_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    }
	  }
	}

	// ratio of energy in first two HCAL layers to total energy in HCAL, only for high energy TPs
	if ( tpEtemu > 10 ) {
	  Ratio_Depth->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  if (abs(tpEtaemu) < 16) {
	    Ratio_DepthHB->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  }
	  if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	    Ratio_DepthHE->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  }
	}

	// loop over HCAL depths for every HCAL TP
        for (int depthIt = 0; depthIt < nDepth-1; depthIt++){
          // count multiplicity of layers given a timing and energy threshold   
          // 3 GeV energy cut
          if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1){
            mult3GeV1ns += 1;
            if (hcalTPtiming[depthIt] > 2){
	      mult3GeV2ns += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult3GeV3ns += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult3GeV4ns += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult3GeV5ns += 1;
		  }
		}
	      }
	    }
	  }
	  // 3 GeV HB HE regions
	  if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	    mult3GeV1nsHB += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult3GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult3GeV3nsHB += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult3GeV4nsHB += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult3GeV5nsHB += 1;
		  }
		}
	      }
	    }
	  }
	  if (hcalTPdepth[depthIt] > 3 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	    mult3GeV1nsHE += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult3GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult3GeV3nsHE += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult3GeV4nsHE += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult3GeV5nsHE += 1;
		  }
		}
	      }
	    }
	  }
	  // 2 GeV energy cut
	  if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1){
	    mult2GeV1ns += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult2GeV2ns += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult2GeV3ns += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult2GeV4ns += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult2GeV5ns += 1;
		  }
		}
	      }
	    }
	  }
	  // 2 GeV HB HE regions                                
	  if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	    mult2GeV1nsHB += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult2GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult2GeV3nsHB += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult2GeV4nsHB += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult2GeV5nsHB += 1;
		  }
		}
	      }
	    }
	  }
	  if (hcalTPdepth[depthIt] > 2 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	    mult2GeV1nsHE += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult2GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult2GeV3nsHE += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult2GeV4nsHE += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult2GeV5nsHE += 1;
		  }
		}
	      }
	    }
	  }
	  // 1 GeV energy cut
	  if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1){
	    mult1GeV1ns += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult1GeV2ns += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult1GeV3ns += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult1GeV4ns += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult1GeV5ns += 1;
		  }
		}
	      }
	    }
	  }
	  // 1 GeV HB HE regions                                                      
	  if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) < 16){
	    mult1GeV1nsHB += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult1GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult1GeV3nsHB += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult1GeV4nsHB += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult1GeV5nsHB += 1;
		  }
		}
	      }
	    }
	  }
	  if (hcalTPdepth[depthIt] > 1 && hcalTPtiming[depthIt] > 1 && abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29){
	    mult1GeV1nsHE += 1;
	    if (hcalTPtiming[depthIt] > 2){
	      mult1GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3){
		mult1GeV3nsHE += 1;
		if (hcalTPtiming[depthIt] > 4){
		  mult1GeV4nsHE += 1;
		  if (hcalTPtiming[depthIt] > 5){
		    mult1GeV5nsHE += 1;
		  }
		}
	      }
	    }
	  }
	} // closing HCAL depths loop
      } // closing HCAL TP loop

      // after HCAL depth and HCAL TP loops fill the histograms with multiplicity variables. The multiplicity counter is reset on each loop iteration 
      // 3 GeV histograms
      dt3GeV1nsMult_emu->Fill(mult3GeV1ns);
      dt3GeV2nsMult_emu->Fill(mult3GeV2ns);
      dt3GeV3nsMult_emu->Fill(mult3GeV3ns);
      dt3GeV4nsMult_emu->Fill(mult3GeV4ns);
      dt3GeV5nsMult_emu->Fill(mult3GeV5ns);
      dt3GeV1nsHEMult_emu->Fill(mult3GeV1nsHE);
      dt3GeV2nsHEMult_emu->Fill(mult3GeV2nsHE);
      dt3GeV3nsHEMult_emu->Fill(mult3GeV3nsHE);
      dt3GeV4nsHEMult_emu->Fill(mult3GeV4nsHE);
      dt3GeV5nsHEMult_emu->Fill(mult3GeV5nsHE);
      dt3GeV1nsHBMult_emu->Fill(mult3GeV1nsHB);
      dt3GeV2nsHBMult_emu->Fill(mult3GeV2nsHB);
      dt3GeV3nsHBMult_emu->Fill(mult3GeV3nsHB);
      dt3GeV4nsHBMult_emu->Fill(mult3GeV4nsHB);
      dt3GeV5nsHBMult_emu->Fill(mult3GeV5nsHB);
      // 2 GeV histograms
      dt2GeV1nsMult_emu->Fill(mult2GeV1ns);
      dt2GeV2nsMult_emu->Fill(mult2GeV2ns);
      dt2GeV3nsMult_emu->Fill(mult2GeV3ns);
      dt2GeV4nsMult_emu->Fill(mult2GeV4ns);
      dt2GeV5nsMult_emu->Fill(mult2GeV5ns);
      dt2GeV1nsHEMult_emu->Fill(mult2GeV1nsHE);
      dt2GeV2nsHEMult_emu->Fill(mult2GeV2nsHE);
      dt2GeV3nsHEMult_emu->Fill(mult2GeV3nsHE);
      dt2GeV4nsHEMult_emu->Fill(mult2GeV4nsHE);
      dt2GeV5nsHEMult_emu->Fill(mult2GeV5nsHE);
      dt2GeV1nsHBMult_emu->Fill(mult2GeV1nsHB);
      dt2GeV2nsHBMult_emu->Fill(mult2GeV2nsHB);
      dt2GeV3nsHBMult_emu->Fill(mult2GeV3nsHB);
      dt2GeV4nsHBMult_emu->Fill(mult2GeV4nsHB);
      dt2GeV5nsHBMult_emu->Fill(mult2GeV5nsHB);
      // 1 GeV histograms
      dt1GeV1nsMult_emu->Fill(mult1GeV1ns);
      dt1GeV2nsMult_emu->Fill(mult1GeV2ns);
      dt1GeV3nsMult_emu->Fill(mult1GeV3ns);
      dt1GeV4nsMult_emu->Fill(mult1GeV4ns);
      dt1GeV5nsMult_emu->Fill(mult1GeV5ns);
      dt1GeV1nsHEMult_emu->Fill(mult1GeV1nsHE);
      dt1GeV2nsHEMult_emu->Fill(mult1GeV2nsHE);
      dt1GeV3nsHEMult_emu->Fill(mult1GeV3nsHE);
      dt1GeV4nsHEMult_emu->Fill(mult1GeV4nsHE);
      dt1GeV5nsHEMult_emu->Fill(mult1GeV5nsHE);
      dt1GeV1nsHBMult_emu->Fill(mult1GeV1nsHB);
      dt1GeV2nsHBMult_emu->Fill(mult1GeV2nsHB);
      dt1GeV3nsHBMult_emu->Fill(mult1GeV3nsHB);
      dt1GeV4nsHBMult_emu->Fill(mult1GeV4nsHB);
      dt1GeV5nsHBMult_emu->Fill(mult1GeV5nsHB);

      // for each bin fill according to whether our object has a larger corresponding energy
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) && mult3GeV2ns_Jets > 3 ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      } 

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) && mult3GeV2ns_Jets > 3 ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) && mult3GeV2ns_Jets > 3 ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) &&  mult3GeV2ns_Jets > 3 ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
             
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_2) >= tauLo + (bin*tauBinWidth) ) doubleISOTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_emu->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_emu->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      }
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_emu->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      }
    }// closes if 'emuOn' is true


    //do routine for L1 hardware quantities
    if (hwOn){

      treeL1TPhw->GetEntry(jentry);
      double tpEt(0.);
      
      for(int i=0; i < l1TPhw_->nHCALTP; i++){
	tpEt = l1TPhw_->hcalTPet[i];
	hcalTP_hw->Fill(tpEt);
      }
      for(int i=0; i < l1TPhw_->nECALTP; i++){
	tpEt = l1TPhw_->ecalTPet[i];
	ecalTP_hw->Fill(tpEt);
      }

      treeL1hw->GetEntry(jentry);
      // get jetEt*, egEt*, tauEt, htSum, mhtSum, etSum, metSum
      // ***INCLUDES NON_ZERO bx*** can't just read values off
      double jetEt_1 = 0;
      double jetEt_2 = 0;
      double jetEt_3 = 0;
      double jetEt_4 = 0;
      for (UInt_t c=0; c<l1hw_->nJets; c++){
        if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] > jetEt_1){
          jetEt_4 = jetEt_3;
          jetEt_3 = jetEt_2;
          jetEt_2 = jetEt_1;
          jetEt_1 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_1 && l1hw_->jetEt[c] > jetEt_2){
          jetEt_4 = jetEt_3;
          jetEt_3 = jetEt_2;      
          jetEt_2 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_2 && l1hw_->jetEt[c] > jetEt_3){
          jetEt_4 = jetEt_3;     
          jetEt_3 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_3 && l1hw_->jetEt[c] > jetEt_4){   
          jetEt_4 = l1hw_->jetEt[c];
        }
      }

      double egEt_1 = 0;
      double egEt_2 = 0;
      for (UInt_t c=0; c<l1hw_->nEGs; c++){
        if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] > egEt_1){
          egEt_2 = egEt_1;
          egEt_1 = l1hw_->egEt[c];
        }
        else if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] <= egEt_1 && l1hw_->egEt[c] > egEt_2){
          egEt_2 = l1hw_->egEt[c];
        }
      }

      double tauEt_1 = 0;
      double tauEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1hw_->nTaus; c++){
        if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] > tauEt_1){
          tauEt_1 = l1hw_->tauEt[c];
        }
        else if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] <= tauEt_1 && l1hw_->tauEt[c] > tauEt_2){
          tauEt_2 = l1hw_->tauEt[c];
        }
      }

      double egISOEt_1 = 0;
      double egISOEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1hw_->nEGs; c++){
        if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] > egISOEt_1 && l1hw_->egIso[c]==1){
          egISOEt_2 = egISOEt_1;
          egISOEt_1 = l1hw_->egEt[c];
        }
        else if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] <= egISOEt_1 && l1hw_->egEt[c] > egISOEt_2 && l1hw_->egIso[c]==1){
          egISOEt_2 = l1hw_->egEt[c];
        }
      }

      double tauISOEt_1 = 0;
      double tauISOEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1hw_->nTaus; c++){
        if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] > tauISOEt_1 && l1hw_->tauIso[c]>0){
          tauISOEt_2 = tauISOEt_1;
          tauISOEt_1 = l1hw_->tauEt[c];
        }
        else if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] <= tauISOEt_1 && l1hw_->tauEt[c] > tauISOEt_2 && l1hw_->tauIso[c]>0){
          tauISOEt_2 = l1hw_->tauEt[c];
        }
      }

      double htSum = 0;
      double mhtSum = 0;
      double etSum = 0;
      double metSum = 0;
      double metHFSum = 0;
      // HW includes -2,-1,0,1,2 bx info (hence the different numbers, could cause a seg fault if this changes)
      for (unsigned int c=0; c<l1hw_->nSums; c++){
          if( l1hw_->sumBx[c] != 0 ) continue;
          if( l1hw_->sumType[c] == L1Analysis::kTotalEt ) etSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kTotalHt ) htSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kMissingEt ) metSum = l1hw_->sumEt[c];
	  if( l1hw_->sumType[c] == L1Analysis::kMissingEtHF ) metHFSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kMissingHt ) mhtSum = l1hw_->sumEt[c];
      }

      // for each bin fill according to whether our object has a larger corresponding energy
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      } 

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
             
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_2) >= tauLo + (bin*tauBinWidth) ) doubleISOTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_hw->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_hw->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_hw->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_hw->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      } 
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_hw->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      } 
    }// closes if 'hwOn' is true
  }// closes loop through events

  //  TFile g( outputFilename.c_str() , "new");
  kk->cd();
  // normalisation factor for rate histograms (11kHz is the orbit frequency)
  double norm = 11246*(numBunch/goodLumiEventCount); // no lumi rescale
  //  double norm = 11246*(numBunch/goodLumiEventCount)*(expectedLum/runLum); //scale to nominal lumi

  if (emuOn){
    singleJetRates_emu->Scale(norm);
    doubleJetRates_emu->Scale(norm);
    tripleJetRates_emu->Scale(norm);
    quadJetRates_emu->Scale(norm);
    singleEgRates_emu->Scale(norm);
    doubleEgRates_emu->Scale(norm);
    singleTauRates_emu->Scale(norm);
    doubleTauRates_emu->Scale(norm);
    singleISOEgRates_emu->Scale(norm);
    doubleISOEgRates_emu->Scale(norm);
    singleISOTauRates_emu->Scale(norm);
    doubleISOTauRates_emu->Scale(norm);
    htSumRates_emu->Scale(norm);
    mhtSumRates_emu->Scale(norm);
    etSumRates_emu->Scale(norm);
    metSumRates_emu->Scale(norm);
    metHFSumRates_emu->Scale(norm);
    /*
    // 3 GeV
    dt3GeV1nsMult_emu->Scale(norm);
    dt3GeV2nsMult_emu->Scale(norm);
    dt3GeV3nsMult_emu->Scale(norm);
    dt3GeV4nsMult_emu->Scale(norm);
    dt3GeV5nsMult_emu->Scale(norm);
    dt3GeV1nsHEMult_emu->Scale(norm);
    dt3GeV2nsHEMult_emu->Scale(norm);
    dt3GeV3nsHEMult_emu->Scale(norm);
    dt3GeV4nsHEMult_emu->Scale(norm);
    dt3GeV5nsHEMult_emu->Scale(norm);
    dt3GeV1nsHBMult_emu->Scale(norm);
    dt3GeV2nsHBMult_emu->Scale(norm);
    dt3GeV3nsHBMult_emu->Scale(norm);
    dt3GeV4nsHBMult_emu->Scale(norm);
    dt3GeV5nsHBMult_emu->Scale(norm);
    // 2 GeV
    dt2GeV1nsMult_emu->Scale(norm);
    dt2GeV2nsMult_emu->Scale(norm);
    dt2GeV3nsMult_emu->Scale(norm);
    dt2GeV4nsMult_emu->Scale(norm);
    dt2GeV5nsMult_emu->Scale(norm);
    dt2GeV1nsHEMult_emu->Scale(norm);
    dt2GeV2nsHEMult_emu->Scale(norm);
    dt2GeV3nsHEMult_emu->Scale(norm);
    dt2GeV4nsHEMult_emu->Scale(norm);
    dt2GeV5nsHEMult_emu->Scale(norm);
    dt2GeV1nsHBMult_emu->Scale(norm);
    dt2GeV2nsHBMult_emu->Scale(norm);
    dt2GeV3nsHBMult_emu->Scale(norm);
    dt2GeV4nsHBMult_emu->Scale(norm);
    dt2GeV5nsHBMult_emu->Scale(norm);
    // 1 GeV
    dt1GeV1nsMult_emu->Scale(norm);
    dt1GeV2nsMult_emu->Scale(norm);
    dt1GeV3nsMult_emu->Scale(norm);
    dt1GeV4nsMult_emu->Scale(norm);
    dt1GeV5nsMult_emu->Scale(norm);
    dt1GeV1nsHEMult_emu->Scale(norm);
    dt1GeV2nsHEMult_emu->Scale(norm);
    dt1GeV3nsHEMult_emu->Scale(norm);
    dt1GeV4nsHEMult_emu->Scale(norm);
    dt1GeV5nsHEMult_emu->Scale(norm);
    dt1GeV1nsHBMult_emu->Scale(norm);
    dt1GeV2nsHBMult_emu->Scale(norm);
    dt1GeV3nsHBMult_emu->Scale(norm);
    dt1GeV4nsHBMult_emu->Scale(norm);
    dt1GeV5nsHBMult_emu->Scale(norm);
    // energy and timing depth
    Energy_Depth->Scale(norm);
    Timing_Depth->Scale(norm);
    Energy_DepthHE->Scale(norm);
    Timing_DepthHE->Scale(norm);
    Energy_DepthHB->Scale(norm);
    Timing_DepthHB->Scale(norm);
    Ratio_Depth->Scale(norm);
    Ratio_DepthHE->Scale(norm);
    Ratio_DepthHB->Scale(norm);
    Ratio_Depth_Jets->Scale(norm);
    Ratio_DepthHE_Jets->Scale(norm);
    Ratio_DepthHB_Jets->Scale(norm);
    */

    // TH1D as a ProfileX of the depth TH2F for timing and energy
    Energy_Depth_avg = Energy_Depth->ProfileX();
    Energy_DepthHE_avg = Energy_DepthHE->ProfileX();
    Energy_DepthHB_avg = Energy_DepthHB->ProfileX();
    Timing_Depth_avg = Timing_Depth->ProfileX();
    Timing_DepthHE_avg = Timing_DepthHE->ProfileX();
    Timing_DepthHB_avg = Timing_DepthHB->ProfileX();
    // for the ones matched with L1 jets
    Energy_Depth_avg_Jets = Energy_Depth_Jets->ProfileX();
    Energy_DepthHE_avg_Jets = Energy_DepthHE_Jets->ProfileX();
    Energy_DepthHB_avg_Jets = Energy_DepthHB_Jets->ProfileX();
    Timing_Depth_avg_Jets = Timing_Depth_Jets->ProfileX();
    Timing_DepthHE_avg_Jets = Timing_DepthHE_Jets->ProfileX();
    Timing_DepthHB_avg_Jets = Timing_DepthHB_Jets->ProfileX();

    //set the errors for the rates
    //want error -> error * sqrt(norm) ?

    hcalTP_emu->Write();
    ecalTP_emu->Write();
    singleJetRates_emu->Write();
    doubleJetRates_emu->Write();
    tripleJetRates_emu->Write();
    quadJetRates_emu->Write();
    singleEgRates_emu->Write();
    doubleEgRates_emu->Write();
    singleTauRates_emu->Write();
    doubleTauRates_emu->Write();
    singleISOEgRates_emu->Write();
    doubleISOEgRates_emu->Write();
    singleISOTauRates_emu->Write();
    doubleISOTauRates_emu->Write();
    htSumRates_emu->Write();
    mhtSumRates_emu->Write();
    etSumRates_emu->Write();
    metSumRates_emu->Write();
    metHFSumRates_emu->Write();
    // 3 GeV
    dt3GeV1nsMult_emu->Write();
    dt3GeV2nsMult_emu->Write();
    dt3GeV3nsMult_emu->Write();
    dt3GeV4nsMult_emu->Write();
    dt3GeV5nsMult_emu->Write();
    dt3GeV1nsHEMult_emu->Write();
    dt3GeV2nsHEMult_emu->Write();
    dt3GeV3nsHEMult_emu->Write();
    dt3GeV4nsHEMult_emu->Write();
    dt3GeV5nsHEMult_emu->Write();
    dt3GeV1nsHBMult_emu->Write();
    dt3GeV2nsHBMult_emu->Write();
    dt3GeV3nsHBMult_emu->Write();
    dt3GeV4nsHBMult_emu->Write();
    dt3GeV5nsHBMult_emu->Write();
    // 2 GeV
    dt2GeV1nsMult_emu->Write();
    dt2GeV2nsMult_emu->Write();
    dt2GeV3nsMult_emu->Write();
    dt2GeV4nsMult_emu->Write();
    dt2GeV5nsMult_emu->Write();
    dt2GeV1nsHEMult_emu->Write();
    dt2GeV2nsHEMult_emu->Write();
    dt2GeV3nsHEMult_emu->Write();
    dt2GeV4nsHEMult_emu->Write();
    dt2GeV5nsHEMult_emu->Write();
    dt2GeV1nsHBMult_emu->Write();
    dt2GeV2nsHBMult_emu->Write();
    dt2GeV3nsHBMult_emu->Write();
    dt2GeV4nsHBMult_emu->Write();
    dt2GeV5nsHBMult_emu->Write();
    // 1 GeV
    dt1GeV1nsMult_emu->Write();
    dt1GeV2nsMult_emu->Write();
    dt1GeV3nsMult_emu->Write();
    dt1GeV4nsMult_emu->Write();
    dt1GeV5nsMult_emu->Write();
    dt1GeV1nsHEMult_emu->Write();
    dt1GeV2nsHEMult_emu->Write();
    dt1GeV3nsHEMult_emu->Write();
    dt1GeV4nsHEMult_emu->Write();
    dt1GeV5nsHEMult_emu->Write();
    dt1GeV1nsHBMult_emu->Write();
    dt1GeV2nsHBMult_emu->Write();
    dt1GeV3nsHBMult_emu->Write();
    dt1GeV4nsHBMult_emu->Write();
    dt1GeV5nsHBMult_emu->Write();

    // and for the multiplicity counter for HCAL TPs matched with L1 Jets
    // 3 GeV
    dt3GeV1nsJetMult_emu->Write();
    dt3GeV2nsJetMult_emu->Write();
    dt3GeV3nsJetMult_emu->Write();
    dt3GeV4nsJetMult_emu->Write();
    dt3GeV5nsJetMult_emu->Write();
    dt3GeV1nsHEJetMult_emu->Write();
    dt3GeV2nsHEJetMult_emu->Write();
    dt3GeV3nsHEJetMult_emu->Write();
    dt3GeV4nsHEJetMult_emu->Write();
    dt3GeV5nsHEJetMult_emu->Write();
    dt3GeV1nsHBJetMult_emu->Write();
    dt3GeV2nsHBJetMult_emu->Write();
    dt3GeV3nsHBJetMult_emu->Write();
    dt3GeV4nsHBJetMult_emu->Write();
    dt3GeV5nsHBJetMult_emu->Write();
    // 2 GeV
    dt2GeV1nsJetMult_emu->Write();
    dt2GeV2nsJetMult_emu->Write();
    dt2GeV3nsJetMult_emu->Write();
    dt2GeV4nsJetMult_emu->Write();
    dt2GeV5nsJetMult_emu->Write();
    dt2GeV1nsHEJetMult_emu->Write();
    dt2GeV2nsHEJetMult_emu->Write();
    dt2GeV3nsHEJetMult_emu->Write();
    dt2GeV4nsHEJetMult_emu->Write();
    dt2GeV5nsHEJetMult_emu->Write();
    dt2GeV1nsHBJetMult_emu->Write();
    dt2GeV2nsHBJetMult_emu->Write();
    dt2GeV3nsHBJetMult_emu->Write();
    dt2GeV4nsHBJetMult_emu->Write();
    dt2GeV5nsHBJetMult_emu->Write();
    // 1 GeV      
    dt1GeV1nsJetMult_emu->Write();
    dt1GeV2nsJetMult_emu->Write();
    dt1GeV3nsJetMult_emu->Write();
    dt1GeV4nsJetMult_emu->Write();
    dt1GeV5nsJetMult_emu->Write();
    dt1GeV1nsHEJetMult_emu->Write();
    dt1GeV2nsHEJetMult_emu->Write();
    dt1GeV3nsHEJetMult_emu->Write();
    dt1GeV4nsHEJetMult_emu->Write();
    dt1GeV5nsHEJetMult_emu->Write();
    dt1GeV1nsHBJetMult_emu->Write();
    dt1GeV2nsHBJetMult_emu->Write();
    dt1GeV3nsHBJetMult_emu->Write();
    dt1GeV4nsHBJetMult_emu->Write();
    dt1GeV5nsHBJetMult_emu->Write();

    Energy_Depth->Write();
    Timing_Depth->Write();
    Energy_DepthHB->Write();
    Timing_DepthHB->Write();
    Energy_DepthHE->Write();
    Timing_DepthHE->Write();

    Energy_Depth_HighE->Write();
    Energy_DepthHB_HighE->Write();
    Energy_DepthHE_HighE->Write();

    Energy_Depth_Jets->Write();
    Timing_Depth_Jets->Write();
    Energy_DepthHB_Jets->Write();
    Timing_DepthHB_Jets->Write();
    Energy_DepthHE_Jets->Write();
    Timing_DepthHE_Jets->Write();

    Energy_Depth_Jets_HighE->Write();
    Energy_DepthHB_Jets_HighE->Write();
    Energy_DepthHE_Jets_HighE->Write();

    Energy_Depth_avg->Write();
    Energy_DepthHE_avg->Write();
    Energy_DepthHB_avg->Write();
    Timing_Depth_avg->Write();
    Timing_DepthHE_avg->Write();
    Timing_DepthHB_avg->Write();

    Energy_Depth_avg_Jets->Write();
    Energy_DepthHE_avg_Jets->Write();
    Energy_DepthHB_avg_Jets->Write();
    Timing_Depth_avg_Jets->Write();
    Timing_DepthHE_avg_Jets->Write();
    Timing_DepthHB_avg_Jets->Write();

    Ratio_Depth->Write();
    Ratio_DepthHE->Write();
    Ratio_DepthHB->Write();
    Ratio_Depth_Jets->Write();
    Ratio_DepthHE_Jets->Write();
    Ratio_DepthHB_Jets->Write();

    centralTiming->Write();
    HoverEtotal->Write();
  }

  if (hwOn){

    singleJetRates_hw->Scale(norm);
    doubleJetRates_hw->Scale(norm);
    tripleJetRates_hw->Scale(norm);
    quadJetRates_hw->Scale(norm);
    singleEgRates_hw->Scale(norm);
    doubleEgRates_hw->Scale(norm);
    singleTauRates_hw->Scale(norm);
    doubleTauRates_hw->Scale(norm);
    singleISOEgRates_hw->Scale(norm);
    doubleISOEgRates_hw->Scale(norm);
    singleISOTauRates_hw->Scale(norm);
    doubleISOTauRates_hw->Scale(norm);
    htSumRates_hw->Scale(norm);
    mhtSumRates_hw->Scale(norm);
    etSumRates_hw->Scale(norm);
    metSumRates_hw->Scale(norm);
    metHFSumRates_hw->Scale(norm);

    hcalTP_hw->Write();
    ecalTP_hw->Write();
    singleJetRates_hw->Write();
    doubleJetRates_hw->Write();
    tripleJetRates_hw->Write();
    quadJetRates_hw->Write();
    singleEgRates_hw->Write();
    doubleEgRates_hw->Write();
    singleTauRates_hw->Write();
    doubleTauRates_hw->Write();
    singleISOEgRates_hw->Write();
    doubleISOEgRates_hw->Write();
    singleISOTauRates_hw->Write();
    doubleISOTauRates_hw->Write();
    htSumRates_hw->Write();
    mhtSumRates_hw->Write();
    etSumRates_hw->Write();
    metSumRates_hw->Write();
    metHFSumRates_hw->Write();
  }
  myfile << "using the following ntuple: " << inputFile << std::endl;
  myfile << "number of colliding bunches = " << numBunch << std::endl;
  myfile << "run luminosity = " << runLum << std::endl;
  myfile << "expected luminosity = " << expectedLum << std::endl;
  myfile << "norm factor used = " << norm << std::endl;
  myfile << "number of good events = " << goodLumiEventCount << std::endl;
  myfile.close(); 
}//closes the function 'rates'
