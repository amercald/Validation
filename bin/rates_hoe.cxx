// Script for calculating rate histograms
// Originally from Aaron Bundock
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"

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

void rates(std::string sampleType, const std::string& inputFileDirectory);

int main(int argc, char *argv[])
{
  std::string ntuplePath("");
  std::string sampleType("");
  if (argc != 3) {
    std::cout << "Usage: rates_hoe.exe [sample type] [path to ntuples]\n" << std::endl;
    exit(1);
  }
  else {
    sampleType = argv[1];
    ntuplePath = argv[2];
  }

  rates(sampleType, ntuplePath);

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

void rates(std::string sampleType, const std::string& inputFileDirectory){
  
  bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
  bool emuOn = true;  //are we using data from emulator?

  if (hwOn==false && emuOn==false){
    std::cout << "exiting as neither hardware or emulator selected" << std::endl;
    return;
  }

  std::string inputFile(inputFileDirectory);
  //  inputFile += "/L1Ntuple_*.root";
  std::string outputDirectory = "emu";  //***runNumber, triggerType, version, hw/emu/both***MAKE SURE IT EXISTS
  std::string outputFilename = "rates_hoe_"+sampleType+".root";
  std::cout << "Will write to " << outputFilename << std::endl;
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
  L1Analysis::L1AnalysisEventDataFormat    *event_ = new L1Analysis::L1AnalysisEventDataFormat();
  eventTree->SetBranchAddress("Event", &event_);
  L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
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
  std::cout << "this nTuple has " << nentries << std::endl;
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
  float htSumHi = 1600.;
  float htSumBinWidth = (htSumHi-htSumLo)/nHtSumBins;

  // mhtSum bins
  int nMhtSumBins = 300;
  float mhtSumLo = 0.;
  float mhtSumHi = 300.;
  float mhtSumBinWidth = (mhtSumHi-mhtSumLo)/nMhtSumBins;

  // etSum bins
  int nEtSumBins = 600;
  float etSumLo = 0.;
  float etSumHi = 1600.;
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

  TH1D * hJetEt = new TH1D("jetET",";ET;",50,0,700);
  //  std::vector<TString> ratioStrings = {"HOvE","HOvE3","HOvE9","H3OvE3","H9OvE9"};
  TH1D * hJetEta = new TH1D("jetEta","jetEta",25,-5,5);
  TH1D * hNJets = new TH1D("njets","njets",20,0,20);
  TH1F* HovEtotal_1x1_emu = new TH1F("HovEtotal_1x1_emu", "HCAL energy / ECAL+HCAL energy for Jets (1x1);H/E;# Entries", 50,0,1);
  TH1F* HovEtotal_3x3_emu = new TH1F("HovEtotal_3x3_emu", "HCAL energy / ECAL+HCAL energy for Jets (3x3);H/E;# Entries", 50,0,1);
  TH1F* HovEtotal_1x1_emu_AllJets = new TH1F("HovEtotal_1x1_emu_AllJets", "HCAL energy / ECAL+HCAL energy for All Jets (1x1);H/E;# Entries", 50,0,1);
  TH1F* HovEtotal_3x3_emu_AllJets = new TH1F("HovEtotal_3x3_emu_AllJets", "HCAL energy / ECAL+HCAL energy for All Jets (3x3);H/E;# Entries", 50,0,1);
  TH1F* HEnergytotal_1x1_emu_AllJets = new TH1F("HEnergytotal_1x1_emu_AllJet", "HCAL Energy for All Jets;# Entries", 150,0,150);
  TH1F* EEnergytotal_1x1_emu_AllJets = new TH1F("EEnergytotal_1x1_emu_AllJet", "ECAL Energy for All Jets;# Entries", 150,0,150);
  TH1F* EEnergytotal_3x3_emu_AllJets = new TH1F("EEnergytotal_3x3_emu_AllJets", "Ecal Energy for All Jets (3x3);H/E;# Entries", 150,0,150);
  TH1F* HEnergytotal_3x3_emu_AllJets = new TH1F("HEnergytotal_3x3_emu_AllJets", "Hcal Energy for All Jets (3x3);H/E;# Entries", 150,0,150);
  TH2F* HEEnergytotal_1x1_emu_AllJets = new TH2F("HEEnergytotal_1x1_emu_AllJet", "HCAL vs ECAL Energy for All Jets", 100, 0, 400, 100, 0, 400);
  TH2F* HEEnergytotal_3x3_emu_AllJets = new TH2F("HEEnergytotal_3x3_emu_AllJet", "HCAL vs ECAL Energy for All Jets", 100, 0, 400, 100, 0, 400);

  TH1D * hJetEtaLeading1 = new TH1D("jetEtaLeading1", "#eta for leading jet 1", 50, -5, 5);
  TH1D * hJetEtaLeading2 = new TH1D("jetEtaLeading2", "#eta for leading jet 2", 50, -5, 5);
  TH1D * hJetEtaLeading3 = new TH1D("jetEtaLeading3", "#eta for leading jet 3", 50, -5, 5);
  TH1D * hJetEtaLeading4 = new TH1D("jetEtaLeading4", "#eta for leading jet 4", 50, -5, 5);

  TH1D * hJetETLeading1 = new TH1D("jetETLeading1", "E_{T} for leading jet 1", 14, 0, 700);
  TH1D * hJetETLeading2 = new TH1D("jetETLeading2", "E_{T} for leading jet 2", 14, 0, 700);
  TH1D * hJetETLeading3 = new TH1D("jetETLeading3", "E_{T} for leading jet 3", 14, 0, 700);
  TH1D * hJetETLeading4 = new TH1D("jetETLeading4", "E_{T} for leading jet 4", 14, 0, 700);

  TH1D * hJetET_cutHoE_1x1_Leading1 = new TH1D("jetET_cutHoE_1x1_Leading1", "E_{T} for leading jet 1", 14, 0, 700);
  TH1D * hJetET_cutHoE_1x1_Leading2 = new TH1D("jetET_cutHoE_1x1_Leading2", "E_{T} for leading jet 2", 14, 0, 700);
  TH1D * hJetET_cutHoE_1x1_Leading3 = new TH1D("jetET_cutHoE_1x1_Leading3", "E_{T} for leading jet 3", 14, 0, 700);
  TH1D * hJetET_cutHoE_1x1_Leading4 = new TH1D("jetET_cutHoE_1x1_Leading4", "E_{T} for leading jet 4", 14, 0, 700);
  TH1D * hJetET_cutHoE_3x3_Leading1 = new TH1D("jetET_cutHoE_3x3_Leading1", "E_{T} for leading jet 1", 14, 0, 700);
  TH1D * hJetET_cutHoE_3x3_Leading2 = new TH1D("jetET_cutHoE_3x3_Leading2", "E_{T} for leading jet 2", 14, 0, 700);
  TH1D * hJetET_cutHoE_3x3_Leading3 = new TH1D("jetET_cutHoE_3x3_Leading3", "E_{T} for leading jet 3", 14, 0, 700);
  TH1D * hJetET_cutHoE_3x3_Leading4 = new TH1D("jetET_cutHoE_3x3_Leading4", "E_{T} for leading jet 4", 14, 0, 700);


  TH1F* HovEtotal_1x1_emu_Leading1 = new TH1F("HovEtotal_1x1_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (1x1);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_1x1_emu_Leading2 = new TH1F("HovEtotal_1x1_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (1x1);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_1x1_emu_Leading3 = new TH1F("HovEtotal_1x1_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (1x1);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_1x1_emu_Leading4 = new TH1F("HovEtotal_1x1_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (1x1);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_3x3_emu_Leading1 = new TH1F("HovEtotal_3x3_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (3x3);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_3x3_emu_Leading2 = new TH1F("HovEtotal_3x3_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (3x3);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_3x3_emu_Leading3 = new TH1F("HovEtotal_3x3_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (3x3);H/E;# Entries", 20,0,1);
  TH1F* HovEtotal_3x3_emu_Leading4 = new TH1F("HovEtotal_3x3_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (3x3);H/E;# Entries", 20,0,1);
  //  TH1F* HovEtotal_5x5_emu_Leading1 = new TH1F("HovEtotal_5x5_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (5x5);H/E;# Entries", 20,0,1);
  //TH1F* HovEtotal_5x5_emu_Leading2 = new TH1F("HovEtotal_5x5_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (5x5);H/E;# Entries", 20,0,1);
  //TH1F* HovEtotal_5x5_emu_Leading3 = new TH1F("HovEtotal_5x5_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (5x5);H/E;# Entries", 20,0,1);
  //TH1F* HovEtotal_5x5_emu_Leading4 = new TH1F("HovEtotal_5x5_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (5x5);H/E;# Entries", 20,0,1);

  //Binning in ET
  TH2F* HovEtotal_1x1_ET_emu_Leading1 = new TH2F("HovEtotal_1x1_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (1x1);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_1x1_ET_emu_Leading2 = new TH2F("HovEtotal_1x1_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (1x1);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_1x1_ET_emu_Leading3 = new TH2F("HovEtotal_1x1_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (1x1);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_1x1_ET_emu_Leading4 = new TH2F("HovEtotal_1x1_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (1x1);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_3x3_ET_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (3x3);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_3x3_ET_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (3x3);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_3x3_ET_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (3x3);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  TH2F* HovEtotal_3x3_ET_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (3x3);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  //  TH2F* HovEtotal_5x5_ETemu_Leading1 = new TH2F("HovEtotal_5x5_ETemu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (5x5);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  //TH2F* HovEtotal_5x5_ETemu_Leading2 = new TH2F("HovEtotal_5x5_ETemu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (5x5);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  //TH2F* HovEtotal_5x5_ETemu_Leading3 = new TH2F("HovEtotal_5x5_ETemu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (5x5);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  //TH2F* HovEtotal_5x5_ETemu_Leading4 = new TH2F("HovEtotal_5x5_ETemu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (5x5);H/E;# Entries", 10, 0, 500, 20, 0, 1);
  

  TH2F* HEEnergytotal_1x1_emu_Leading1 = new TH2F("HEEnergytotal_1x1_emu_Leading1", "HCAL vs ECAL Energy for Leading 1", 100, 0, 300, 100, 0, 300);
  TH2F* HEEnergytotal_1x1_emu_Leading2 = new TH2F("HEEnergytotal_1x1_emu_Leading2", "HCAL vs ECAL Energy for Leading 2", 100, 0, 300, 100, 0, 300);
  TH2F* HEEnergytotal_1x1_emu_Leading3 = new TH2F("HEEnergytotal_1x1_emu_Leading3", "HCAL vs ECAL Energy for Leading 3", 100, 0, 300, 100, 0, 300);
  TH2F* HEEnergytotal_1x1_emu_Leading4 = new TH2F("HEEnergytotal_1x1_emu_Leading4", "HCAL vs ECAL Energy for Leading 4", 100, 0, 300, 100, 0, 300);
  TH2F* HEEnergytotal_3x3_emu_Leading1 = new TH2F("HEEnergytotal_3x3_emu_Leading1", "HCAL vs ECAL Energy for Leading 1", 100, 0, 400, 100, 0, 400);
  TH2F* HEEnergytotal_3x3_emu_Leading2 = new TH2F("HEEnergytotal_3x3_emu_Leading2", "HCAL vs ECAL Energy for Leading 2", 100, 0, 400, 100, 0, 400);
  TH2F* HEEnergytotal_3x3_emu_Leading3 = new TH2F("HEEnergytotal_3x3_emu_Leading3", "HCAL vs ECAL Energy for Leading 3", 100, 0, 400, 100, 0, 400);
  TH2F* HEEnergytotal_3x3_emu_Leading4 = new TH2F("HEEnergytotal_3x3_emu_Leading4", "HCAL vs ECAL Energy for Leading 4", 100, 0, 400, 100, 0, 400);
  //TH2F* HEEnergytotal_5x5_emu_Leading1 = new TH2F("HEEnergytotal_5x5_emu_Leading1", "HCAL vs ECAL Energy for Leading 1", 100, 0, 400, 100, 0, 400);
  //TH2F* HEEnergytotal_5x5_emu_Leading2 = new TH2F("HEEnergytotal_5x5_emu_Leading2", "HCAL vs ECAL Energy for Leading 2", 100, 0, 400, 100, 0, 400);
  //TH2F* HEEnergytotal_5x5_emu_Leading3 = new TH2F("HEEnergytotal_5x5_emu_Leading3", "HCAL vs ECAL Energy for Leading 3", 100, 0, 400, 100, 0, 400);
  //TH2F* HEEnergytotal_5x5_emu_Leading4 = new TH2F("HEEnergytotal_5x5_emu_Leading4", "HCAL vs ECAL Energy for Leading 4", 100, 0, 400, 100, 0, 400);

  /*  TProfile* HovE_ET_profile_1x1_Leading1 = new TProfile("HovE_ET_profile_1x1_Leading1", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_1x1_Leading2 = new TProfile("HovE_ET_profile_1x1_Leading2", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_1x1_Leading3 = new TProfile("HovE_ET_profile_1x1_Leading3", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_1x1_Leading4 = new TProfile("HovE_ET_profile_1x1_Leading4", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_3x3_Leading1 = new TProfile("HovE_ET_profile_3x3_Leading1", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_3x3_Leading2 = new TProfile("HovE_ET_profile_3x3_Leading2", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_3x3_Leading3 = new TProfile("HovE_ET_profile_3x3_Leading3", "", 10, 0, 500, 10, 0, 1);
  TProfile* HovE_ET_profile_3x3_Leading4 = new TProfile("HovE_ET_profile_3x3_Leading4", "", 10, 0, 500, 10, 0, 1);
  */
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

      // for each bin fill according to whether our object has a larger corresponding energy
      int seedTowerIEta(-1);
      int seedTowerIPhi(-1);
      int maxTowerEndcap = 28;
      //      int maxTowerBarrel = 16;
      int minTowerForHOvE = -999; //maxTowerBarrel+1;
      int maxTowerForHOvE = maxTowerEndcap;
      
      uint nJetemu(0);
      //towEtamu not used for now
      //      double towEtemu(0), towHademu(0), towEmemu(0), towEtaemu(0), towPhiemu(0), nTowemu(0);
      double towHademu(0), towEmemu(0), towEtaemu(0), towPhiemu(0), nTowemu(0);
      nTowemu = l1Towemu_->nTower;
      nJetemu = l1emu_->nJets;
      hNJets->Fill(nJetemu);
      //std::cout << "nTower emu = " << nTowemu << " and nJet emu = " << nJetemu << std::endl;
      std::map<const TString, std::vector<double> > hadVariablesAllJets;
      std::map<const TString, std::vector<double> > emVariablesAllJets;
      //      if (nJetemu ==0) continue;
      for(uint jetIt=0; jetIt<nJetemu; jetIt++){
	hJetEt->Fill(l1emu_->jetEt[jetIt]);
	hJetEta ->Fill(l1emu_->jetEta[jetIt]);
	seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
	seedTowerIEta = l1emu_->jetTowerIEta[jetIt];
	double seedTowerHad(0), seedTowerEm(0), seedTower3x3Em(0), seedTower3x3Had(0), seedTower5x5Em(0), seedTower5x5Had(0), seedTower9x9Em(0), seedTower9x9Had(0);
	for (int towIt = 0; towIt < nTowemu; towIt++){
	  //towEtemu  = l1Towemu_->iet[towIt];
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
		int wrappedIPhi = (seedTowerIPhi+iSeedTowerIPhi); //% 72) + 1;
		if (wrappedIPhi > 72) wrappedIPhi -= 72;
		if (wrappedIPhi < 0) wrappedIPhi += 72;
		if (towEtaemu == seedTowerIEta+iSeedTowerIEta && towPhiemu == wrappedIPhi){
		  seedTower9x9Em += towEmemu;
		  seedTower9x9Had += towHademu;
		  if (abs(iSeedTowerIPhi) <= 1 && abs(iSeedTowerIEta) <= 1){
		    seedTower3x3Em += towEmemu;
		    seedTower3x3Had += towHademu;
		  }
		  if (abs(iSeedTowerIPhi) <= 2 && abs(iSeedTowerIEta) <= 2){
		    seedTower5x5Em += towEmemu;
		    seedTower5x5Had += towHademu;
		  }
		}
	      }
	    }
	  } // closing min max tower statement
	} // closing seed tower loop
	if ( (seedTowerHad / seedTower5x5Had) <= 0.2) continue; //requirement for throwing out junk jets
	double HoE_value_1x1 = 0, HoE_value_3x3 = 0;
	if (seedTowerHad != 0 || seedTowerEm != 0) HoE_value_1x1 = seedTowerHad / (seedTowerHad + seedTowerEm);
	if (seedTower3x3Had != 0 || seedTower3x3Em != 0) HoE_value_3x3 = seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em);
	
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
	
	if (seedTowerHad > 0 && seedTowerEm > 0)
	  {
	    HovEtotal_1x1_emu_AllJets->Fill( seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_emu_AllJets->Fill( seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );
	    HEnergytotal_1x1_emu_AllJets->Fill( seedTowerHad );
	    EEnergytotal_1x1_emu_AllJets->Fill( seedTowerEm );
	    EEnergytotal_3x3_emu_AllJets->Fill( seedTower3x3Em );
	    HEnergytotal_3x3_emu_AllJets->Fill( seedTower3x3Had );
	    HEEnergytotal_1x1_emu_AllJets->Fill(seedTowerEm, seedTowerHad);
	    HEEnergytotal_3x3_emu_AllJets->Fill(seedTower3x3Em, seedTower3x3Had);
	  }
	if (jetIt == 0) 
	  {
	    hJetEtaLeading1->Fill(l1emu_->jetEta[jetIt]);
	    hJetETLeading1->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_1x1 > 0.9) hJetET_cutHoE_1x1_Leading1->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_3x3 > 0.9) hJetET_cutHoE_3x3_Leading1->Fill(l1emu_->jetEt[jetIt]);
	    HEEnergytotal_1x1_emu_Leading1->Fill(seedTowerEm, seedTowerHad);
	    HEEnergytotal_3x3_emu_Leading1->Fill(seedTower3x3Em, seedTower3x3Had);	   
	    HovEtotal_1x1_emu_Leading1->Fill(seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_emu_Leading1->Fill(seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );
	    
	    HovEtotal_1x1_ET_emu_Leading1->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_ET_emu_Leading1->Fill(l1emu_->jetEt[jetIt], seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	    // HovE_ET_profile_1x1_Leading1->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) ); 
	    //HovE_ET_profile_3x3_Leading1->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );

	  }
	else if (jetIt == 1) 
	  {
	    hJetEtaLeading2->Fill(l1emu_->jetEta[jetIt]);
	    hJetETLeading2->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_1x1 > 0.9) hJetET_cutHoE_1x1_Leading2->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_3x3 > 0.9) hJetET_cutHoE_3x3_Leading2->Fill(l1emu_->jetEt[jetIt]);
	    HEEnergytotal_1x1_emu_Leading2->Fill(seedTowerEm, seedTowerHad);
	    HEEnergytotal_3x3_emu_Leading2->Fill(seedTower3x3Em, seedTower3x3Had);	    
	    HovEtotal_1x1_emu_Leading2->Fill(seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_emu_Leading2->Fill(seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	    HovEtotal_1x1_ET_emu_Leading2->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_ET_emu_Leading2->Fill(l1emu_->jetEt[jetIt], seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	    //HovE_ET_profile_1x1_Leading2->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) ); 
	    //HovE_ET_profile_3x3_Leading2->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	  }
	else if (jetIt == 2)
	  { 
	    hJetEtaLeading3->Fill(l1emu_->jetEta[jetIt]);
	    hJetETLeading3->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_1x1 > 0.9) hJetET_cutHoE_1x1_Leading3->Fill(l1emu_->jetEt[jetIt]);
	    if (HoE_value_3x3 > 0.9) hJetET_cutHoE_3x3_Leading3->Fill(l1emu_->jetEt[jetIt]);
	    HEEnergytotal_1x1_emu_Leading3->Fill(seedTowerEm, seedTowerHad);
	    HEEnergytotal_3x3_emu_Leading3->Fill(seedTower3x3Em, seedTower3x3Had);
	    HovEtotal_1x1_emu_Leading3->Fill(seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_emu_Leading3->Fill(seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	    HovEtotal_1x1_ET_emu_Leading3->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	    HovEtotal_3x3_ET_emu_Leading3->Fill(l1emu_->jetEt[jetIt], seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	    //HovE_ET_profile_1x1_Leading3->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	    //HovE_ET_profile_3x3_Leading3->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	  }
	else if (jetIt == 3) 
	  {
	  hJetEtaLeading4->Fill(l1emu_->jetEta[jetIt]);
	  hJetETLeading4->Fill(l1emu_->jetEt[jetIt]);
	  if (HoE_value_1x1 > 0.9) hJetET_cutHoE_1x1_Leading4->Fill(l1emu_->jetEt[jetIt]);
	  if (HoE_value_3x3 > 0.9) hJetET_cutHoE_3x3_Leading4->Fill(l1emu_->jetEt[jetIt]);
	  HEEnergytotal_1x1_emu_Leading4->Fill(seedTowerEm, seedTowerHad);
	  HEEnergytotal_3x3_emu_Leading4->Fill(seedTower3x3Em, seedTower3x3Had);
	  HovEtotal_1x1_emu_Leading4->Fill(seedTowerHad / (seedTowerHad + seedTowerEm) );
	  HovEtotal_3x3_emu_Leading4->Fill(seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	  HovEtotal_1x1_ET_emu_Leading4->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) );
	  HovEtotal_3x3_ET_emu_Leading4->Fill(l1emu_->jetEt[jetIt], seedTower3x3Had / (seedTower3x3Had + seedTower3x3Em) );

	  //HovE_ET_profile_1x1_Leading4->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) ); 
	  //HovE_ET_profile_3x3_Leading4->Fill(l1emu_->jetEt[jetIt], seedTowerHad / (seedTowerHad + seedTowerEm) ); 
	  }
      } // Closing the jet loop
      HovEtotal_1x1_emu->Fill((hadVariablesAllJets["HOvE"][0])/(hadVariablesAllJets["HOvE"][0]+emVariablesAllJets["HOvE"][0]));
      HovEtotal_3x3_emu->Fill((hadVariablesAllJets["H3OvE3"][0])/(hadVariablesAllJets["H3OvE3"][0]+emVariablesAllJets["H3OvE3"][0]));
      std::vector<bool> pass_HoE(4, false);
      for (unsigned int ijet = 0; ijet < pass_HoE.size(); ijet++)
	{
	  if ((hadVariablesAllJets["H3OvE3"].at(ijet))/(hadVariablesAllJets["H3OvE3"].at(ijet)+emVariablesAllJets["H3OvE3"].at(ijet)) > 0.85) pass_HoE.at(ijet) = true;
	}
      for(int bin=0; bin<nJetBins; bin++){
        if( (pass_HoE[0]) && ( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      } 
      for(int bin=0; bin<nJetBins; bin++){
        if( (pass_HoE[1]) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
      for(int bin=0; bin<nJetBins; bin++){
        if( (pass_HoE[2]) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
      for(int bin=0; bin<nJetBins; bin++){
        if( (pass_HoE[3]) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
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
        if( (pass_HoE[0] || pass_HoE[1] || pass_HoE[2] || pass_HoE[3]) && (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_emu->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (pass_HoE[0] || pass_HoE[1] || pass_HoE[2] || pass_HoE[3]) && (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
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

    hJetEta->Write();
    hJetEt->Write();
    hNJets->Write();
    HovEtotal_1x1_emu->Write();
    HovEtotal_3x3_emu->Write();
    HovEtotal_1x1_emu_AllJets->Write();
    HovEtotal_3x3_emu_AllJets->Write();
    HEnergytotal_1x1_emu_AllJets->Write();
    HEnergytotal_3x3_emu_AllJets->Write();
    EEnergytotal_1x1_emu_AllJets->Write();
    EEnergytotal_3x3_emu_AllJets->Write();
    HEEnergytotal_1x1_emu_AllJets->Write();
    HEEnergytotal_3x3_emu_AllJets->Write();

    //Leading 4 jet histograms
    HovEtotal_1x1_emu_Leading1->Write();
    HovEtotal_3x3_emu_Leading1->Write();
    HovEtotal_1x1_emu_Leading2->Write();
    HovEtotal_3x3_emu_Leading2->Write();
    HovEtotal_1x1_emu_Leading3->Write();
    HovEtotal_3x3_emu_Leading3->Write();
    HovEtotal_1x1_emu_Leading4->Write();
    HovEtotal_3x3_emu_Leading4->Write();

    HovEtotal_1x1_ET_emu_Leading1->Write();
    HovEtotal_3x3_ET_emu_Leading1->Write();
    HovEtotal_1x1_ET_emu_Leading2->Write();
    HovEtotal_3x3_ET_emu_Leading2->Write();
    HovEtotal_1x1_ET_emu_Leading3->Write();
    HovEtotal_3x3_ET_emu_Leading3->Write();
    HovEtotal_1x1_ET_emu_Leading4->Write();
    HovEtotal_3x3_ET_emu_Leading4->Write();


    hJetEtaLeading1->Write();
    hJetEtaLeading2->Write();
    hJetEtaLeading3->Write();
    hJetEtaLeading4->Write();

    hJetETLeading1->Write();
    hJetETLeading2->Write();
    hJetETLeading3->Write();
    hJetETLeading4->Write();

    hJetET_cutHoE_1x1_Leading1->Write();
    hJetET_cutHoE_1x1_Leading2->Write();
    hJetET_cutHoE_1x1_Leading3->Write();
    hJetET_cutHoE_1x1_Leading4->Write();

    hJetET_cutHoE_3x3_Leading1->Write();
    hJetET_cutHoE_3x3_Leading2->Write();
    hJetET_cutHoE_3x3_Leading3->Write();
    hJetET_cutHoE_3x3_Leading4->Write();


    HEEnergytotal_1x1_emu_Leading1->Write();
    HEEnergytotal_3x3_emu_Leading1->Write();
    HEEnergytotal_1x1_emu_Leading2->Write();
    HEEnergytotal_3x3_emu_Leading2->Write();
    HEEnergytotal_1x1_emu_Leading3->Write();
    HEEnergytotal_3x3_emu_Leading3->Write();
    HEEnergytotal_1x1_emu_Leading4->Write();
    HEEnergytotal_3x3_emu_Leading4->Write();
    /*
    HovE_ET_profile_1x1_Leading1->Write();
    HovE_ET_profile_1x1_Leading2->Write();
    HovE_ET_profile_1x1_Leading3->Write();
    HovE_ET_profile_1x1_Leading4->Write();
    HovE_ET_profile_3x3_Leading1->Write();
    HovE_ET_profile_3x3_Leading2->Write();
    HovE_ET_profile_3x3_Leading3->Write();
    HovE_ET_profile_3x3_Leading4->Write();
    */
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
