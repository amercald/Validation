//Script for calculating %rate histograms
// Originally from Aaron Bundock
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"


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
        //    exit(1);
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
    if (lumiBlock >= 1 || lumiBlock <= 10000) {
        return true;
    }

    return false;
}

// this is a correction function for the eta phi value of the intersection of a particle not resulting from the IP       
// used since the b quark is matched to the TP based on eta phi, but the b quark results from a LLP decay so has a different vertex 
// from HcalCompareUpgradeChains intersect function https://github.com/gk199/cms-hcal-debug/blob/PulseShape/plugins/HcalCompareUpgradeChains.cc#L961
std::vector<double> intersect(double vx, double vy,double vz, double px, double py, double pz) {
    double lightSpeed = 29979245800;
    double radius = 179; // 130 for calorimeters (ECAL + HCAL)
    double length = 388; // 300 for calorimeters (ECAL + HCAL)
    double energy = sqrt(px*px + py*py + pz*pz);
    // First work out intersection with cylinder (barrel)        
    double a = (px*px + py*py)*lightSpeed*lightSpeed/(energy*energy);
    double b = 2*(vx*px + vy*py)*lightSpeed/energy;
    double c = (vx*vx + vy*vy) - radius*radius;
    double sqrt_disc = sqrt(b*b - 4*a*c);
    double tCircle1 = (-b + sqrt_disc)/(2*a);
    double tCircle2 = (-b - sqrt_disc)/(2*a);
    // If intersection in the past it doesn't count         
    if (tCircle1 < 0) tCircle1 = 1E9;
    if (tCircle2 < 0) tCircle2 = 1E9;
    // If the intsersection occurs outside the barrel length it doesn't count                       
    double zPosCircle1 = tCircle1*(pz/energy)*lightSpeed + vz;
    double zPosCircle2 = tCircle2*(pz/energy)*lightSpeed + vz;
    if (zPosCircle1 > length) tCircle1 = 1E9;
    if (zPosCircle2 > length) tCircle2 = 1E9;
    // Now work out if it intersects the endcap                      
    double tPlane1 = (length-vz)*energy/(pz*lightSpeed);
    double tPlane2 = (-length-vz)*energy/(pz*lightSpeed);
    // If intersection in the past it doesn't count                     
    if (tPlane1 < 0) tPlane1 = 1E9;
    if (tPlane2 < 0) tPlane2 = 1E9;
    double xPosPlane1 = tPlane1*(px/energy)*lightSpeed + vx;
    double yPosPlane1 = tPlane1*(py/energy)*lightSpeed + vy;
    double xPosPlane2 = tPlane2*(px/energy)*lightSpeed + vx;
    double yPosPlane2 = tPlane2*(py/energy)*lightSpeed + vy;
    // If the intsersection occurs outside the endcap radius it doesn't count     
    if (sqrt(xPosPlane1*xPosPlane1 + yPosPlane1*yPosPlane1) > radius) tPlane1 = 1E9;
    if (sqrt(xPosPlane2*xPosPlane2+yPosPlane2*yPosPlane2) > radius) tPlane2 = 1E9;
    // Find the first intersection                          
    double tInter = std::min({tCircle1,tCircle2,tPlane1,tPlane2});
    // Return 1000,1000 if not intersection with barrel or endcap             
    std::vector<double> etaphi;
    if (tInter > 1E6)
    {
        etaphi.push_back(1000);
        etaphi.push_back(1000);
        return etaphi;
    }
    // Find position of intersection                          
    double xPos = tInter*(px/energy)*lightSpeed + vx;
    double yPos = tInter*(py/energy)*lightSpeed + vy;
    double zPos = tInter*(pz/energy)*lightSpeed + vz;
    // Find eta/phi of intersection                          
    double phi = atan2(yPos,xPos); // return the arc tan in radians                                                                                                                               
    double theta = acos(zPos/sqrt(xPos*xPos + yPos*yPos + zPos*zPos));
    double eta = -log(tan(theta/2.));
    etaphi.push_back(eta);
    etaphi.push_back(phi);
    return etaphi;
}
  

double deltaPhi(double phi1, double phi2) { 
    double result = phi1 - phi2;
    if(fabs(result) > 9999) return result;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
}

double DeltaR(double phi1, double phi2, double eta1, double eta2)
{
    double phidiff = deltaPhi(phi1, phi2);
    double etadiff = eta1 - eta2;
    return sqrt(phidiff*phidiff + etadiff*etadiff);
}

double etaVal(int ieta) { // calculate eta given ieta
    double etavl;
    if (ieta <= -24){
        etavl = .1695*ieta + 1.9931;
    }
    else if (ieta <= -1){
        etavl = .0875*ieta + .0489;
    }
    else if (ieta < 24){
        etavl = .0875*ieta - .0489;
    }
    else {
        etavl = .1695*ieta - 1.9931;
    }
    return etavl;
}
double phiVal(int iphi) { // calculate phi given iphi
    double phiBins=72.;
    double phivl;
    phivl=double(iphi)*(2.*TMath::Pi()/phiBins);
    if (iphi > 36) phivl -= 2.*TMath::Pi();
    return phivl;
}

bool sort_jets (int i,int j) { return (i<j); }

void rates(std::string sampleType, const std::string& inputFileDirectory){
  
    bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
    bool emuOn = true;  //are we using data from emulator?

    if (hwOn==false && emuOn==false){
        std::cout << "exiting as neither hardware or emulator selected" << std::endl;
        return;
    }

    std::string inputFile(inputFileDirectory);
    std::string outputDirectory = "emu";  //***runNumber, triggerType, version, hw/emu/both***MAKE SURE IT EXISTS
    std::string outputFilename = "rates_"+sampleType+".root";
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
    TChain * treeL1CaloTPhw = new TChain("l1CaloTowerTree/L1CaloTowerTree");
    if (hwOn){
        treeL1CaloTPhw->Add(inputFile.c_str());
    }
    TChain * eventTree = new TChain("l1EventTree/L1EventTree");
    eventTree->Add(inputFile.c_str());
    TChain * genTree = new TChain("l1GeneratorTree/L1GenTree");
    genTree->Add(inputFile.c_str());

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
    L1Analysis::L1AnalysisGeneratorDataFormat    *generator_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
    genTree->SetBranchAddress("Generator", &generator_);
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

    TH1F * hJetEt = new TH1F("jetET",";ET;",50,0,700);
    TH1F * hJetEta = new TH1F("jetEta","jetEta",25,-5,5);
    TH1F * hNJets = new TH1F("njets","njets",20,0,20);

    TH1F * hJetEtaLeading1 = new TH1F("jetEtaLeading1", "#eta for leading jet 1", 50, -5, 5);
    TH1F * hJetEtaLeading2 = new TH1F("jetEtaLeading2", "#eta for leading jet 2", 50, -5, 5);
    TH1F * hJetEtaLeading3 = new TH1F("jetEtaLeading3", "#eta for leading jet 3", 50, -5, 5);
    TH1F * hJetEtaLeading4 = new TH1F("jetEtaLeading4", "#eta for leading jet 4", 50, -5, 5);

    TH1F * hJetETLeading1 = new TH1F("jetETLeading1", "E_{T} for leading jet 1", 14, 0, 700);
    TH1F * hJetETLeading2 = new TH1F("jetETLeading2", "E_{T} for leading jet 2", 14, 0, 700);
    TH1F * hJetETLeading3 = new TH1F("jetETLeading3", "E_{T} for leading jet 3", 14, 0, 700);
    TH1F * hJetETLeading4 = new TH1F("jetETLeading4", "E_{T} for leading jet 4", 14, 0, 700);

    TH1F* betagammaLLP = new TH1F("betagammaLLP", ";#beta#gamma; # Entries", 100, 0, 10);
    TH1F* h_nGenParticles = new TH1F("nGenParticles", ";# Generator Particles; # Entries", 11, -0.5, 10.5);
    //  TH2F* flightlength_eta_Barrel = new TH2F("flightlength_eta_Barrel", ";#eta; # Flight length", 100, -3, 3, 100, 0, 1500);
    //TH2F* flightlength_eta_Barrel = new TH2F("flightlength_eta_Barrel", ";#eta; # Flight length", 100, -3, 3, 100, 0, 1500);

    //HoE study plots

    TH1F * hJetET_cutHoE_1x1_Leading1 = new TH1F("jetET_cutHoE_1x1_Leading1", "E_{T} for leading jet 1", 14, 0, 700);
    TH1F * hJetET_cutHoE_1x1_Leading2 = new TH1F("jetET_cutHoE_1x1_Leading2", "E_{T} for leading jet 2", 14, 0, 700);
    TH1F * hJetET_cutHoE_1x1_Leading3 = new TH1F("jetET_cutHoE_1x1_Leading3", "E_{T} for leading jet 3", 14, 0, 700);
    TH1F * hJetET_cutHoE_1x1_Leading4 = new TH1F("jetET_cutHoE_1x1_Leading4", "E_{T} for leading jet 4", 14, 0, 700);
    TH1F * hJetET_cutHoE_3x3_Leading1 = new TH1F("jetET_cutHoE_3x3_Leading1", "E_{T} for leading jet 1", 14, 0, 700);
    TH1F * hJetET_cutHoE_3x3_Leading2 = new TH1F("jetET_cutHoE_3x3_Leading2", "E_{T} for leading jet 2", 14, 0, 700);
    TH1F * hJetET_cutHoE_3x3_Leading3 = new TH1F("jetET_cutHoE_3x3_Leading3", "E_{T} for leading jet 3", 14, 0, 700);
    TH1F * hJetET_cutHoE_3x3_Leading4 = new TH1F("jetET_cutHoE_3x3_Leading4", "E_{T} for leading jet 4", 14, 0, 700);

    TH1D * hJet1x1ov5x5 = new TH1D("hJet1x1ov5x5 ", "1#times1/5#times5", 50, 0, 1);

    TH1F* HovEtotal_1x1_emu = new TH1F("HovEtotal_1x1_emu", "HCAL energy / ECAL+HCAL energy for Jets (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu = new TH1F("HovEtotal_3x3_emu", "HCAL energy / ECAL+HCAL energy for Jets (3x3);H/E;# Entries", 24,0,1.2);

    TH1F* HovEtotal_1x1_emu_Barrel = new TH1F("HovEtotal_1x1_emu_Barrel", "HCAL energy / ECAL+HCAL energy for Jets (1x1) in Barrel;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Barrel = new TH1F("HovEtotal_3x3_emu_Barrel", "HCAL energy / ECAL+HCAL energy for Jets (3x3) in Barrel;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_Endcap = new TH1F("HovEtotal_1x1_emu_Endcap", "HCAL energy / ECAL+HCAL energy for Jets (1x1) in Endcap;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Endcap = new TH1F("HovEtotal_3x3_emu_Endcap", "HCAL energy / ECAL+HCAL energy for Jets (3x3) in Endcap;H/E;# Entries", 24,0,1.2);


    TH1F* HovEtotalLog_1x1_emu = new TH1F("HovEtotalLog_1x1_emu", "log(H/E) for Jets (1x1);log_{10}(H/E);# Entries", 40,-10,10);
    TH1F* HovEtotalLog_3x3_emu = new TH1F("HovEtotalLog_3x3_emu", "log(H/E) for Jets (1x1);log_{10}(H/E) (3x3);H/E;# Entries", 40,-10,10);


    TH1F* HovEtotal_1x1_emu_Leading1 = new TH1F("HovEtotal_1x1_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading2 = new TH1F("HovEtotal_1x1_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading3 = new TH1F("HovEtotal_1x1_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading4 = new TH1F("HovEtotal_1x1_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading1 = new TH1F("HovEtotal_3x3_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (3x3);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading2 = new TH1F("HovEtotal_3x3_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (3x3);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading3 = new TH1F("HovEtotal_3x3_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (3x3);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading4 = new TH1F("HovEtotal_3x3_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (3x3);H/E;# Entries", 24,0,1.2);
    //  TH1F* HovEtotal_5x5_emu_Leading1 = new TH1F("HovEtotal_5x5_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (5x5);H/E;# Entries", 24,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading2 = new TH1F("HovEtotal_5x5_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (5x5);H/E;# Entries", 24,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading3 = new TH1F("HovEtotal_5x5_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (5x5);H/E;# Entries", 24,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading4 = new TH1F("HovEtotal_5x5_emu_

    TH1F* HovEtotal_1x1_emu_GenMatchedJets = new TH1F("HovEtotal_1x1_emu_GenMatchedJets", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets = new TH1F("HovEtotal_3x3_emu_GenMatchedJets", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1);H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_GenMatchedJets_Barrel = new TH1F("HovEtotal_1x1_emu_GenMatchedJets_Barrel", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Barrel;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets_Barrel = new TH1F("HovEtotal_3x3_emu_GenMatchedJets_Barrel", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Barrel;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_1x1_emu_GenMatchedJets_Endcap = new TH1F("HovEtotal_1x1_emu_GenMatchedJets_Endcap", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Endcap;H/E;# Entries", 24,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets_Endcap = new TH1F("HovEtotal_3x3_emu_GenMatchedJets_Endcap", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Endcap;H/E;# Entries", 24,0,1.2);



    //Binning in ET
    TH2F* HovEtotal_1x1_ET_emu_Leading1 = new TH2F("HovEtotal_1x1_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (1x1);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading2 = new TH2F("HovEtotal_1x1_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (1x1);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading3 = new TH2F("HovEtotal_1x1_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (1x1);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading4 = new TH2F("HovEtotal_1x1_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (1x1);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (3x3);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (3x3);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (3x3);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (3x3);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    //  TH2F* HovEtotal_5x5_ETemu_Leading1 = new TH2F("HovEtotal_5x5_ETemu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (5x5);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading2 = new TH2F("HovEtotal_5x5_ETemu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (5x5);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading3 = new TH2F("HovEtotal_5x5_ETemu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (5x5);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading4 = new TH2F("HovEtotal_5x5_ETemu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (5x5);H/E;# Entries", 10, 0, 500, 24, 0, 1.2);
    
    //Binning in depth energy
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 24, 0, 1.2);

    TH2F* HovEtotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu = new TH2F("HovEtotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (3-4/1-2) (3x3);E_{D3-D4/D1-D2};H/(H+E)", 100, 0, 10, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu = new TH2F("HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (4/1-3) (3x3);E_{D4/D1-D3};H/(H+E)", 100, 0, 10, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu = new TH2F("HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (4-7/1-3) (3x3);E_{D4-7/D1-D3};H/(H+E)", 100, 0, 10, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_DepthRatio_Gen_Barrel_3_4_ov_1_2_emu = new TH2F("HovEtotal_3x3_DepthRatio_Gen_Barrel_3_4_ov_1_2_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (3-4/1-2) (3x3);E_{D3-D4/D1-D2};H/(H+E)", 100, 0, 10, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_DepthRatio_Gen_Barrel_4_ov_1_2_3_emu = new TH2F("HovEtotal_3x3_DepthRatio_Gen_Barrel_4_ov_1_2_3_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (4/1-3) (3x3);E_{D4/D1-D3};H/(H+E)", 100, 0, 10, 24, 0, 1.2);
    TH2F* HovEtotal_3x3_DepthRatio_Gen_Endcap_4_7_ov_1_3_emu = new TH2F("HovEtotal_3x3_DepthRatio_Gen_Endcap_4_7_ov_1_3_emu", "HCAL energy / ECAL+HCAL energy vs Depth Ratio (4-7/1-3) (3x3);E_{D4-7/D1-D3};H/(H+E)", 100, 0, 10, 24, 0, 1.2);

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


    TH1F* L1JetTPDR = new TH1F("L1JetTPDR", ";#DeltaR; # Entries", 100, 0, 5);
    TH1F* L1JetNTP = new TH1F("L1JetNTP", ";# TPs; # Entries", 20, -0.5, 19.5);

    //Energy depth study plots

    TH2F* energyDepth_Barrel = new TH2F("energyDepth_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Endcap = new TH2F("energyDepth_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_Jet1_Barrel = new TH2F("energyDepth_Jet1_Barrel", "Depth profile, leading jet 1, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet2_Barrel = new TH2F("energyDepth_Jet2_Barrel", "Depth profile, leading jet 2, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet3_Barrel = new TH2F("energyDepth_Jet3_Barrel", "Depth profile, leading jet 3, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet4_Barrel = new TH2F("energyDepth_Jet4_Barrel", "Depth profile, leading jet 4, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet1_Endcap = new TH2F("energyDepth_Jet1_Endcap", "Depth profile, leading jet 1, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet2_Endcap = new TH2F("energyDepth_Jet2_Endcap", "Depth profile, leading jet 2, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet3_Endcap = new TH2F("energyDepth_Jet3_Endcap", "Depth profile, leading jet 3, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Jet4_Endcap = new TH2F("energyDepth_Jet4_Endcap", "Depth profile, leading jet 4, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_Barrel = new TH2F("energyDepth_genMatchInclusive_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_Endcap = new TH2F("energyDepth_genMatchInclusive_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchTP_Barrel = new TH2F("energyDepth_genMatchTP_Barrel", "Depth profile, matched TPs, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchTP_Endcap = new TH2F("energyDepth_genMatchTP_Endcap", "Depth profile, matched TPs, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH1F* DeltaRLLP = new TH1F("DeltaRLLP", ";#Delta R; # Entries", 100, 0, 2);
    
    TH1F* energyDepth_Ratio_Barrel_D4 = new TH1F("energyDepth_Ratio_Barrel_D4", "Energy Depth Ratio In Barrel Layers D4" ";E_{D4}/E_{TP}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Barrel_D3_D4 = new TH1F("energyDepth_Ratio_Barrel_D3_D4", "Energy Depth Ratio In Barrel Layers D3-D4" ";E_{D3-D4}/E_{TP}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Endcap_D6_D7 = new TH1F("energyDepth_Ratio_Endcap_D6_D7", "Energy Depth Ratio In Endcap Layers D6-D7" ";E_{D6-D7/E_{TP}}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Endcap_D4_D7 = new TH1F("energyDepth_Ratio_Endcap_D4_D7", "Energy Depth Ratio In Endcap Layers D4-D7" ";E_{D4-D7/E_{TP}}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Gen_Barrel_D4 = new TH1F("energyDepth_Ratio_Gen_Barrel_D4", "Energy Depth Ratio In Barrel Layers D4 (Matched)" ";E_{D4}/E_{TP}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Gen_Barrel_D3_D4 = new TH1F("energyDepth_Ratio_Gen_Barrel_D3_D4", "Energy Depth Ratio In Barrel Layers D3-D4 (Matched)" ";E_{D3-D4}/E_{TP}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Gen_Endcap_D6_D7 = new TH1F("energyDepth_Ratio_Gen_Endcap_D6_D7", "Energy Depth Ratio In Endcap Layers D6-D7 (Matched)" ";E_{D6-D7/E_{TP}}; # Entries", 50, 0, 1);
    TH1F* energyDepth_Ratio_Gen_Endcap_D4_D7 = new TH1F("energyDepth_Ratio_Gen_Endcap_D4_D7", "Energy Depth Ratio In Endcap Layers D4-D7 (Matched)" ";E_{D4-D7/E_{TP}}; # Entries", 50, 0, 1);

      TH1F* energyDepth_NTPs_Barrel_D3_D4_ge06 = new TH1F("energyDepth_NTPs_Barrel_D3_D4_ge06", "NTPs with D3-D4 ratio ge 0.6" ";NTPs; # Entries", 10, -0.5, 9.5);
      TH1F* energyDepth_NTPs_Endcap_D4_D7_ge045 = new TH1F("energyDepth_NTPs_Endcap_D4_D7_ge045", "NTPs with D4-D7 ratio ge 0.6" ";NTPs; # Entries", 10, -0.5, 9.5);
      TH1F* energyDepth_NTPs_Gen_Barrel_D3_D4_ge06 = new TH1F("energyDepth_NTPs_Gen_Barrel_D3_D4_ge06", "NTPs with D3-D4 ratio ge 0.6 (matched)" ";NTPs; # Entries", 10, -0.5, 9.5);
      TH1F* energyDepth_NTPs_Gen_Endcap_D4_D7_ge045 = new TH1F("energyDepth_NTPs_Gen_Endcap_D4_D7_ge045", "NTPs with D4-D7 ratio ge 0.6 (matched)" ";NTPs; # Entries", 10, -0.5, 9.5);


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
            treeL1CaloTPemu->GetEntry(jentry);
            treeL1Towemu->GetEntry(jentry);
            treeL1emu->GetEntry(jentry);
            genTree->GetEntry(jentry);
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

            //Begin depth profiles part
      
            int nHCALTP = l1CaloTPemu_->nHCALTP;
            int nJetemu = l1emu_->nJets;
      
            int nGenPart = generator_->nPart;
            std::vector<bool> goodGenParticles(nGenPart, false);
            std::vector<double> genParticlesEta, genParticlesPhi;
            std::vector<bool> matchedJet(nJetemu, false);
            int nGoodGen = 0;
            int nOkayGen = 0;
            int nHad = 0;
            int nMatched = 0;

            double LLPcounter = 0;
            for(int genpart = 0; genpart < nGenPart; genpart++)
            {
                double Vz = generator_->partVz[genpart];
                double Vx = generator_->partVx[genpart];
                double Vy = generator_->partVy[genpart];
                double Pz = generator_->partPz[genpart];
                double Px = generator_->partPx[genpart];
                double Py = generator_->partPy[genpart];
                double energy = generator_->partE[genpart];
                double radius = sqrt(Vx*Vx + Vy*Vy);
                int pdgId = generator_->partId[genpart];
                bool inEndcap = abs(Vz) > 388 && abs(Vz) < 568 && radius < 568;
                bool inBarrel = abs(Vz) < 388 && radius > 179 && radius < 295;
                bool inHCAL = inEndcap || inBarrel;
                bool isQuarkorGluon = abs(pdgId) <=5  || abs(pdgId) == 21;

                bool isLLP = pdgId == 6000113;
                //if (generator_->partHardProcess[genpart] != 0) std::cout << pdgId << " " << generator_->partParent[genpart] << std::endl;

                std::vector<double> intersection = intersect(Vx, Vy, Vz, Px, Py, Pz);
                genParticlesEta.push_back(intersection.at(0));
                genParticlesPhi.push_back(intersection.at(1));
                if (isQuarkorGluon && inHCAL && generator_->partHardProcess[genpart] != 0)// && abs(intersection.at(0)) < 2.4)
                {	    
                    goodGenParticles.at(genpart) = true;
                    nGoodGen += 1;
                }
                if (inHCAL) nOkayGen += 1;
                if (isQuarkorGluon) nHad += 1;

                if (isLLP)
                {
                    double total_momentum = sqrt(Pz*Pz + Px*Px + Py*Py);
                    double betagamma = total_momentum / sqrt(energy*energy - total_momentum*total_momentum);
                    betagammaLLP->Fill(betagamma);
                    LLPcounter += 1;
                }
            }
            //std::cout << LLPcounter << std::endl;
            h_nGenParticles->Fill(nGoodGen);

            std::vector<std::vector<double>> hcalTPdepth;
            for (int TPIt = 0; TPIt < nHCALTP; TPIt++)
            {
                std::vector<double> depth_vector(7, 0);
                depth_vector.at(0) = l1CaloTPemu_->hcalTPDepth1[TPIt];
                depth_vector.at(1) = l1CaloTPemu_->hcalTPDepth2[TPIt];
                depth_vector.at(2) = l1CaloTPemu_->hcalTPDepth3[TPIt];
                depth_vector.at(3) = l1CaloTPemu_->hcalTPDepth4[TPIt];
                depth_vector.at(4) = l1CaloTPemu_->hcalTPDepth5[TPIt];
                depth_vector.at(5) = l1CaloTPemu_->hcalTPDepth6[TPIt];
                depth_vector.at(6) = l1CaloTPemu_->hcalTPDepth7[TPIt];

                hcalTPdepth.push_back(depth_vector);
            }

            std::vector<bool> matchedDirectTP(nHCALTP, false);
            //match gen particles directly to TPs
            for(int genpart = 0; genpart < nGenPart; genpart++)
            {
                for(int tpIt = 0; tpIt < nHCALTP; tpIt++)
                {
                    double TPeta = etaVal(l1CaloTPemu_->hcalTPieta[tpIt]);
                    double TPphi = phiVal(l1CaloTPemu_->hcalTPiphi[tpIt]);
                    double deltaRIt = DeltaR(genParticlesPhi.at(genpart), TPphi, genParticlesEta.at(genpart), TPeta);
		
                    if (goodGenParticles.at(genpart) && deltaRIt < 0.5)
                    {
                        matchedDirectTP.at(tpIt) = true;
                    }	       	      
                }
            }

            //match gen particles to jets
      
            for(int jetIt = 0; jetIt < nJetemu; jetIt++)
            {
                double jetEta = l1emu_->jetEta[jetIt];
                double jetPhi = l1emu_->jetPhi[jetIt];

                for(int genpart = 0; genpart < nGenPart; genpart++)
                {
                    //if (!goodGenParticles.at(genpart)) continue;
                    double genEta = generator_->partEta[genpart];
                    double genPhi = generator_->partPhi[genpart];
                    //double genEta = genParticlesEta.at(genpart);
                    //double genPhi = genParticlesPhi.at(genpart);
                    double genpartDeltaR = DeltaR(jetPhi, genPhi, jetEta, genEta);
                    if (nGoodGen > 0 && genpartDeltaR < 0.2 && generator_->partId[genpart] == 6000113) 
                    {
                        matchedJet.at(jetIt) = true;
                    }
                }
                if (matchedJet.at(jetIt)) nMatched += 1;
            }

            int nDepth = 0;
            double tpiEtaemu = 0, tpEtemu = 0, scaledEDepth = 0;
            std::vector<double> depthTPIt;
            std::vector<bool> usedTP(nHCALTP, false);

            //Fill histograms

            for(int TPIt = 0; TPIt < nHCALTP; TPIt++)
            {
                nDepth = l1CaloTPemu_->hcalTPnDepths[TPIt];
                depthTPIt = hcalTPdepth[TPIt];
                tpEtemu = l1CaloTPemu_->hcalTPet[TPIt];
                tpiEtaemu = l1CaloTPemu_->hcalTPieta[TPIt];
                if (tpEtemu < 5 || nDepth == 0) continue;
                for(int depthIt = 0; depthIt < nDepth; depthIt++)
                {
                    scaledEDepth = depthTPIt[depthIt]/tpEtemu;
                    if (abs(tpiEtaemu) < 16)
                    {
                        if (matchedDirectTP.at(TPIt)) energyDepth_genMatchTP_Barrel->Fill(depthIt+1, scaledEDepth);
                    } 
                    else if (abs(tpiEtaemu) > 16 && abs(tpiEtaemu) < 29)
                    {
                        if (matchedDirectTP.at(TPIt)) energyDepth_genMatchTP_Endcap->Fill(depthIt+1, scaledEDepth);
                    }	      
                }
            }	  
    
            //H/E study
            int seedTowerIEta(-1);
            int seedTowerIPhi(-1);
            int maxTowerEndcap = 28;
            //      int maxTowerBarrel = 16;
            int minTowerForHOvE = -999; //maxTowerBarrel+1;
            int maxTowerForHOvE = maxTowerEndcap;
      
            //towEtamu not used for now
            //      double towEtemu(0), towHademu(0), towEmemu(0), towEtaemu(0), towPhiemu(0), nTowemu(0);
            double towHademu(0), towEmemu(0), towEtaemu(0), towPhiemu(0), nTowemu(0);
            nTowemu = l1Towemu_->nTower;
            hNJets->Fill(nJetemu);
            //std::cout << "nTower emu = " << nTowemu << " and nJet emu = " << nJetemu << std::endl;
            std::vector<TH2F*> HovE_3x3_ET_AllDepthBarrel_hists{HovEtotal_3x3_ET_Depth1_Barrel_emu, HovEtotal_3x3_ET_Depth2_Barrel_emu, HovEtotal_3x3_ET_Depth3_Barrel_emu, HovEtotal_3x3_ET_Depth4_Barrel_emu, HovEtotal_3x3_ET_Depth5_Barrel_emu, HovEtotal_3x3_ET_Depth6_Barrel_emu, HovEtotal_3x3_ET_Depth7_Barrel_emu};
            std::vector<TH2F*> HovE_3x3_ET_AllDepth_GenBarrel_hists{HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu, HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu};
            std::vector<TH2F*> HovE_3x3_ET_AllDepthEndcap_hists{HovEtotal_3x3_ET_Depth1_Endcap_emu, HovEtotal_3x3_ET_Depth2_Endcap_emu, HovEtotal_3x3_ET_Depth3_Endcap_emu, HovEtotal_3x3_ET_Depth4_Endcap_emu, HovEtotal_3x3_ET_Depth5_Endcap_emu, HovEtotal_3x3_ET_Depth6_Endcap_emu, HovEtotal_3x3_ET_Depth7_Endcap_emu};
            std::vector<TH2F*> HovE_3x3_ET_AllDepth_GenEndcap_hists{HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu, HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu};

            std::map<const std::string, std::vector<double> > jetVariablesAllJets;
            std::map<const std::string, std::vector<double> > hadVariablesAllJets;
            std::map<const std::string, std::vector<double> > emVariablesAllJets;
            std::map<const std::string, std::vector<double> > jetVariablesGenMatchedJets;
            std::map<const std::string, std::vector<double> > hadVariablesGenMatchedJets;
            std::map<const std::string, std::vector<double> > emVariablesGenMatchedJets;

            std::vector< std::vector<double>> jetEnergyPerDepth;
            std::vector< std::vector<double>> jetEnergyPerDepth_Gen;
            std::vector<int> inHBorHE; // 0 means neither, 1 means in HB, 2 means in HE
            std::vector<int> inHBorHE_Gen;

            int nPassedJets(0);
            int nPassedGenMatchedJets(0);
            //      if (nJetemu ==0) continue;
            for(int jetIt=0; jetIt<nJetemu; jetIt++){
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
                                int wrappedIPhi = (seedTowerIPhi+iSeedTowerIPhi);
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
	
                if (jetIt < 4) hJet1x1ov5x5->Fill(seedTowerHad / seedTower5x5Had);

                if ( (seedTowerHad / seedTower5x5Had) > 0.1)// && seedTower3x3Had > 30)  //requirement for throwing out junk jets
                {
                    double jetEta = l1emu_->jetEta[jetIt];
                    double jetPhi = l1emu_->jetPhi[jetIt];
                    double minDR = 100;
                    int TPcounter = 0;
                    std::map<std::string, int> DepthTPcounter;
                    DepthTPcounter["BarrelD3D4_ge06"] = 0;
                    DepthTPcounter["EndcapD4D7_ge045"] = 0;
                    DepthTPcounter["BarrelD3D4_ge06_Gen"] = 0;
                    DepthTPcounter["EndcapD4D7_ge045_Gen"] = 0;
                
                    for(int TPIt = 0; TPIt < nHCALTP; TPIt++)
                    {
                        double TPeta = etaVal(l1CaloTPemu_->hcalTPieta[TPIt]);
                        double TPphi = phiVal(l1CaloTPemu_->hcalTPiphi[TPIt]);
                        double deltaRIt = DeltaR(jetPhi, TPphi, jetEta, TPeta);
                        if (/*l1CaloTPemu_->hcalTPet[TPIt] > 10 &&*/ deltaRIt < minDR) minDR = deltaRIt;
                        if (deltaRIt < 0.3 && !usedTP.at(TPIt))
                        {
                            usedTP.at(TPIt) = true;
//                            nDepth = l1CaloTPemu_->hcalTPnDepths[TPIt];
                            depthTPIt = hcalTPdepth[TPIt];
                            tpEtemu = l1CaloTPemu_->hcalTPet[TPIt];
                            tpiEtaemu = l1CaloTPemu_->hcalTPieta[TPIt];
                            if (tpEtemu < 10) continue;
                            TPcounter++;
                            if( abs(tpiEtaemu) < 16)
                            {
                                if (nPassedJets < 4)
                                {
                                    if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.6) DepthTPcounter["BarrelD3D4_ge06"] += 1;
                                    energyDepth_Ratio_Barrel_D4->Fill(depthTPIt.at(3) / tpEtemu);
                                    energyDepth_Ratio_Barrel_D3_D4->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                }
                                if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4)
                                {
                                    if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.6) DepthTPcounter["BarrelD3D4_ge06_Gen"] += 1;
                                    energyDepth_Ratio_Gen_Barrel_D4->Fill(depthTPIt.at(3) / tpEtemu);
                                    energyDepth_Ratio_Gen_Barrel_D3_D4->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);

                                }
                            }
                            else if ( abs(tpiEtaemu) > 16 && abs(tpiEtaemu) < 29)
                            {
                                if (nPassedJets < 4)
                                {
                                    if ((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.45) DepthTPcounter["EndcapD4D7_ge045"] += 1;
                                    energyDepth_Ratio_Endcap_D6_D7->Fill((depthTPIt.at(5) + depthTPIt.at(6))  / tpEtemu);
                                    energyDepth_Ratio_Endcap_D4_D7->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                }
                                if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4)
                                {
                                    if ((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.45) DepthTPcounter["EndcapD4D7_ge045"] += 1;
                                    energyDepth_Ratio_Gen_Endcap_D6_D7->Fill((depthTPIt.at(5) + depthTPIt.at(6))  / tpEtemu);
                                    energyDepth_Ratio_Gen_Endcap_D4_D7->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                }
                                
                            }   
                            for(int depthIt = 0; depthIt < 7; depthIt++)
                            {
                                scaledEDepth = depthTPIt[depthIt]/tpEtemu;
                                if (abs(tpiEtaemu) < 16)
                                {
                                    if (nPassedJets < 4)
                                    {
                                        energyDepth_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (jetIt == 0) energyDepth_Jet1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 1) energyDepth_Jet2_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 2) energyDepth_Jet3_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 3) energyDepth_Jet4_Barrel->Fill(depthIt+1, scaledEDepth);

                                        HovE_3x3_ET_AllDepthBarrel_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                        
                                    }
                                    if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4)
                                    {
                                        energyDepth_genMatchInclusive_Barrel->Fill(depthIt+1, scaledEDepth);
                                        HovE_3x3_ET_AllDepth_GenBarrel_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                } 
                                else if (abs(tpiEtaemu) > 16 && abs(tpiEtaemu) < 29)
                                {
                                    if (nPassedJets < 4)
                                    {
                                        energyDepth_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (jetIt == 0) energyDepth_Jet1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 1) energyDepth_Jet2_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 2) energyDepth_Jet3_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 3) energyDepth_Jet4_Endcap->Fill(depthIt+1, scaledEDepth);
                                        
                                        HovE_3x3_ET_AllDepthEndcap_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                    if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4) 
                                    {
                                        energyDepth_genMatchInclusive_Endcap->Fill(depthIt+1, scaledEDepth);
                                        HovE_3x3_ET_AllDepth_GenEndcap_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                }
                            }
                        }
                    }	  
                    
                    if (nPassedJets < 4) 
                    {
                        L1JetTPDR->Fill(minDR);
                        L1JetNTP->Fill(TPcounter);
                        if (TPcounter > 0) energyDepth_NTPs_Barrel_D3_D4_ge06->Fill(DepthTPcounter["BarrelD3D4_ge06"]);
                        if (TPcounter > 0) energyDepth_NTPs_Endcap_D4_D7_ge045->Fill(DepthTPcounter["EndcapD4D7_ge045"]);
                    }
                    if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4)
                    {
                        if (TPcounter > 0) energyDepth_NTPs_Gen_Barrel_D3_D4_ge06->Fill(DepthTPcounter["BarrelD3D4_ge06_Gen"]);
                        if (TPcounter > 0) energyDepth_NTPs_Gen_Endcap_D4_D7_ge045->Fill(DepthTPcounter["EndcapD4D7_ge045_Gen"]);
                    }
	
                    jetVariablesAllJets["ET"].push_back(l1emu_->jetEt[jetIt]);
                    jetVariablesAllJets["eta"].push_back(l1emu_->jetEta[jetIt]);
                    jetVariablesAllJets["ieta"].push_back(seedTowerIEta);
	    
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

                    nPassedJets++;
	    
                    if (matchedJet.at(jetIt) && nGoodGen > 0)
                    {
                        jetVariablesGenMatchedJets["ET"].push_back(l1emu_->jetEt[jetIt]);
                        jetVariablesGenMatchedJets["eta"].push_back(l1emu_->jetEta[jetIt]);
                        jetVariablesGenMatchedJets["ieta"].push_back(seedTowerIEta);
		
                        hadVariablesGenMatchedJets["HOvE"].push_back(seedTowerHad);
                        hadVariablesGenMatchedJets["HOvE3"].push_back(seedTowerHad);
                        hadVariablesGenMatchedJets["HOvE9"].push_back(seedTowerHad);
                        hadVariablesGenMatchedJets["H3OvE3"].push_back(seedTower3x3Had);
                        hadVariablesGenMatchedJets["H9OvE9"].push_back(seedTower9x9Had);
		
                        emVariablesGenMatchedJets["HOvE"].push_back(seedTowerEm);
                        emVariablesGenMatchedJets["HOvE3"].push_back(seedTower3x3Em);
                        emVariablesGenMatchedJets["HOvE9"].push_back(seedTower9x9Em);
                        emVariablesGenMatchedJets["H3OvE3"].push_back(seedTower3x3Em);
                        emVariablesGenMatchedJets["H9OvE9"].push_back(seedTower9x9Em);
		
                        nPassedGenMatchedJets++;
                    }
                }
            }

            std::vector<TH1F*> jetET_hists{hJetETLeading1, hJetETLeading2, hJetETLeading3, hJetETLeading4}; 
            std::vector<TH1F*> jetEta_hists{hJetEtaLeading1, hJetEtaLeading2, hJetEtaLeading3, hJetEtaLeading4}; 
            std::vector<TH2F*> HEEnergy_1x1_hists{HEEnergytotal_1x1_emu_Leading1, HEEnergytotal_1x1_emu_Leading2, HEEnergytotal_1x1_emu_Leading3, HEEnergytotal_1x1_emu_Leading4};      
            std::vector<TH2F*> HEEnergy_3x3_hists{HEEnergytotal_3x3_emu_Leading1, HEEnergytotal_3x3_emu_Leading2, HEEnergytotal_3x3_emu_Leading3, HEEnergytotal_3x3_emu_Leading4};      
            std::vector<TH1F*> HovE_1x1_hists{HovEtotal_1x1_emu_Leading1, HovEtotal_1x1_emu_Leading2, HovEtotal_1x1_emu_Leading3, HovEtotal_1x1_emu_Leading4};      
            std::vector<TH1F*> HovE_3x3_hists{HovEtotal_3x3_emu_Leading1, HovEtotal_3x3_emu_Leading2, HovEtotal_3x3_emu_Leading3, HovEtotal_3x3_emu_Leading4}; 
            std::vector<TH2F*> HovE_1x1_ET_hists{HovEtotal_1x1_ET_emu_Leading1, HovEtotal_1x1_ET_emu_Leading2, HovEtotal_1x1_ET_emu_Leading3, HovEtotal_1x1_ET_emu_Leading4};      
            std::vector<TH2F*> HovE_3x3_ET_hists{HovEtotal_3x3_ET_emu_Leading1, HovEtotal_3x3_ET_emu_Leading2, HovEtotal_3x3_ET_emu_Leading3, HovEtotal_3x3_ET_emu_Leading4};
/*            std::vector<TH2F*> HovE_3x3_ET_Depth1_hists{HovEtotal_3x3_ET_Depth1_emu_Leading1, HovEtotal_3x3_ET_Depth1_emu_Leading2, HovEtotal_3x3_ET_Depth1_emu_Leading3, HovEtotal_3x3_ET_Depth1_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth2_hists{HovEtotal_3x3_ET_Depth2_emu_Leading1, HovEtotal_3x3_ET_Depth2_emu_Leading2, HovEtotal_3x3_ET_Depth2_emu_Leading3, HovEtotal_3x3_ET_Depth2_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth3_hists{HovEtotal_3x3_ET_Depth3_emu_Leading1, HovEtotal_3x3_ET_Depth3_emu_Leading2, HovEtotal_3x3_ET_Depth3_emu_Leading3, HovEtotal_3x3_ET_Depth3_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth4_hists{HovEtotal_3x3_ET_Depth4_emu_Leading1, HovEtotal_3x3_ET_Depth4_emu_Leading2, HovEtotal_3x3_ET_Depth4_emu_Leading3, HovEtotal_3x3_ET_Depth4_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth5_hists{HovEtotal_3x3_ET_Depth5_emu_Leading1, HovEtotal_3x3_ET_Depth5_emu_Leading2, HovEtotal_3x3_ET_Depth5_emu_Leading3, HovEtotal_3x3_ET_Depth5_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth6_hists{HovEtotal_3x3_ET_Depth6_emu_Leading1, HovEtotal_3x3_ET_Depth6_emu_Leading2, HovEtotal_3x3_ET_Depth6_emu_Leading3, HovEtotal_3x3_ET_Depth6_emu_Leading4};
            std::vector<TH2F*> HovE_3x3_ET_Depth7_hists{HovEtotal_3x3_ET_Depth7_emu_Leading1, HovEtotal_3x3_ET_Depth7_emu_Leading2, HovEtotal_3x3_ET_Depth7_emu_Leading3, HovEtotal_3x3_ET_Depth7_emu_Leading4};
            std::vector< std::vector<TH2F*> > HovE_Depth_Hists{HovE_3x3_ET_Depth1_hists, HovE_3x3_ET_Depth2_hists, HovE_3x3_ET_Depth3_hists, HovE_3x3_ET_Depth4_hists, HovE_3x3_ET_Depth5_hists, HovE_3x3_ET_Depth6_hists, HovE_3x3_ET_Depth7_hists}; */
            
            for (int pjet = 0; pjet < nPassedJets && pjet < 4; pjet++)
            {
                jetEta_hists.at(pjet)->Fill(jetVariablesAllJets["eta"].at(pjet));
                jetET_hists.at(pjet)->Fill(jetVariablesAllJets["ET"].at(pjet));
                HEEnergy_1x1_hists.at(pjet)->Fill( emVariablesAllJets["HOvE"].at(pjet),  hadVariablesAllJets["HOvE"].at(pjet));
                HEEnergy_3x3_hists.at(pjet)->Fill( emVariablesAllJets["H3OvE3"].at(pjet),  hadVariablesAllJets["H3OvE3"].at(pjet));   
                HovE_1x1_hists.at(pjet)->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)) );
                HovE_3x3_hists.at(pjet)->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)) );
	        
                HovE_1x1_ET_hists.at(pjet)->Fill(jetVariablesAllJets["ET"].at(pjet), hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)) );
                HovE_3x3_ET_hists.at(pjet)->Fill(jetVariablesAllJets["ET"].at(pjet), hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)) );

                HovEtotal_1x1_emu->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                HovEtotal_3x3_emu->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                double logHovE = emVariablesAllJets["HOvE"].at(pjet) > 0 ? log10(hadVariablesAllJets["HOvE"].at(pjet)/emVariablesAllJets["HOvE"].at(pjet)) : 10;
                double logHovE_3x3 = emVariablesAllJets["H3OvE3"].at(pjet) > 0 ? log10(hadVariablesAllJets["H3OvE3"].at(pjet)/emVariablesAllJets["H3OvE3"].at(pjet)) : 10;	 
                HovEtotalLog_1x1_emu->Fill(logHovE);
                HovEtotalLog_3x3_emu->Fill(logHovE_3x3);
                if (abs(jetVariablesAllJets["ieta"].at(pjet)) < 16) 
                {
                    HovEtotal_1x1_emu_Barrel->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                    HovEtotal_3x3_emu_Barrel->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                }
                else if (abs(jetVariablesAllJets["ieta"].at(pjet)) > 16 && abs(jetVariablesAllJets["ieta"].at(pjet)) < 29) 
                {
                    HovEtotal_1x1_emu_Endcap->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                    HovEtotal_3x3_emu_Endcap->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                }
            }

            for (int pgjet = 0; pgjet < nPassedGenMatchedJets && pgjet < 4; pgjet++)
            {
                //std::cout << hadVariablesGenMatchedJets["HOvE"].at(pgjet) / (hadVariablesGenMatchedJets["HOvE"].at(pgjet) + emVariablesGenMatchedJets["HOvE"].at(pgjet)) << std::endl;
                HovEtotal_1x1_emu_GenMatchedJets->Fill(hadVariablesGenMatchedJets["HOvE"].at(pgjet) / (hadVariablesGenMatchedJets["HOvE"].at(pgjet) + emVariablesGenMatchedJets["HOvE"].at(pgjet)));
                HovEtotal_3x3_emu_GenMatchedJets->Fill(hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) / (hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) + emVariablesGenMatchedJets["H3OvE3"].at(pgjet)));
                if (abs(jetVariablesGenMatchedJets["ieta"].at(pgjet)) < 16) 
                {
                    HovEtotal_1x1_emu_GenMatchedJets_Barrel->Fill(hadVariablesGenMatchedJets["HOvE"].at(pgjet) / (hadVariablesGenMatchedJets["HOvE"].at(pgjet) + emVariablesGenMatchedJets["HOvE"].at(pgjet)));
                    HovEtotal_3x3_emu_GenMatchedJets_Barrel->Fill(hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) / (hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) + emVariablesGenMatchedJets["H3OvE3"].at(pgjet)));
                }
                else if (abs(jetVariablesGenMatchedJets["ieta"].at(pgjet)) > 16 && abs(jetVariablesGenMatchedJets["ieta"].at(pgjet)) < 29) 
                {
                    HovEtotal_1x1_emu_GenMatchedJets_Endcap->Fill(hadVariablesGenMatchedJets["HOvE"].at(pgjet) / (hadVariablesGenMatchedJets["HOvE"].at(pgjet) + emVariablesGenMatchedJets["HOvE"].at(pgjet)));
                    HovEtotal_3x3_emu_GenMatchedJets_Endcap->Fill(hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) / (hadVariablesGenMatchedJets["H3OvE3"].at(pgjet) + emVariablesGenMatchedJets["H3OvE3"].at(pgjet)));
                }

            }

            //      std::cout << nPassedJets << " " << nPassedGenMatchedJets << std::endl;
            std::vector<bool> pass_HoE(std::max(nPassedJets, 4) , false);
            for (unsigned int ijet = 0; ijet < hadVariablesAllJets["H3OvE3"].size(); ijet++)
            {
                if ((hadVariablesAllJets["H3OvE3"].at(ijet))/(hadVariablesAllJets["H3OvE3"].at(ijet)+emVariablesAllJets["H3OvE3"].at(ijet)) > -2) pass_HoE.at(ijet) = true;
            }

            // for each bin fill according to whether our object has a larger corresponding energy
            for(int bin=0; bin<nJetBins; bin++){
                if( pass_HoE.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
            } 

            for(int bin=0; bin<nJetBins; bin++){
                if( (pass_HoE.at(0) || pass_HoE.at(1)) || (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
            }  

            for(int bin=0; bin<nJetBins; bin++){
                if( (pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2)) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
            }  

            for(int bin=0; bin<nJetBins; bin++){
                if( (pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2) || pass_HoE.at(3)) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
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
                if( (pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2) || pass_HoE.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
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

        hJetEtaLeading1->Write();
        hJetEtaLeading2->Write();
        hJetEtaLeading3->Write();
        hJetEtaLeading4->Write();

        hJetETLeading1->Write();
        hJetETLeading2->Write();
        hJetETLeading3->Write();
        hJetETLeading4->Write();

        betagammaLLP->Write();
        DeltaRLLP->Write();
        h_nGenParticles->Write();

        hJetET_cutHoE_1x1_Leading1->Write();
        hJetET_cutHoE_1x1_Leading2->Write();
        hJetET_cutHoE_1x1_Leading3->Write();
        hJetET_cutHoE_1x1_Leading4->Write();

        hJetET_cutHoE_3x3_Leading1->Write();
        hJetET_cutHoE_3x3_Leading2->Write();
        hJetET_cutHoE_3x3_Leading3->Write();
        hJetET_cutHoE_3x3_Leading4->Write();

        hJet1x1ov5x5->Write();

        HovEtotal_1x1_emu->Write();
        HovEtotal_3x3_emu->Write();

        HovEtotalLog_1x1_emu->Write();
        HovEtotalLog_3x3_emu->Write();

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
        HovEtotal_1x1_emu_Barrel->Write();
        HovEtotal_3x3_emu_Barrel->Write();
        HovEtotal_1x1_emu_Endcap->Write();
        HovEtotal_3x3_emu_Endcap->Write();


        HovEtotal_1x1_emu_GenMatchedJets->Write();
        HovEtotal_3x3_emu_GenMatchedJets->Write();
        HovEtotal_1x1_emu_GenMatchedJets_Barrel->Write();
        HovEtotal_3x3_emu_GenMatchedJets_Barrel->Write();
        HovEtotal_1x1_emu_GenMatchedJets_Endcap->Write();
        HovEtotal_3x3_emu_GenMatchedJets_Endcap->Write();

        HEEnergytotal_1x1_emu_Leading1->Write();
        HEEnergytotal_3x3_emu_Leading1->Write();
        HEEnergytotal_1x1_emu_Leading2->Write();
        HEEnergytotal_3x3_emu_Leading2->Write();
        HEEnergytotal_1x1_emu_Leading3->Write();
        HEEnergytotal_3x3_emu_Leading3->Write();
        HEEnergytotal_1x1_emu_Leading4->Write();
        HEEnergytotal_3x3_emu_Leading4->Write();

        //energy depth plots
    
        energyDepth_Barrel->Write();
        energyDepth_Endcap->Write();

        energyDepth_Jet1_Barrel->Write();
        energyDepth_Jet2_Barrel->Write();
        energyDepth_Jet3_Barrel->Write();
        energyDepth_Jet4_Barrel->Write();
        energyDepth_Jet1_Endcap->Write();
        energyDepth_Jet2_Endcap->Write();
        energyDepth_Jet3_Endcap->Write();
        energyDepth_Jet4_Endcap->Write();

        energyDepth_genMatchInclusive_Barrel->Write();
        energyDepth_genMatchInclusive_Endcap->Write();

        energyDepth_genMatchTP_Barrel->Write();
        energyDepth_genMatchTP_Endcap->Write();

        //HovE and Depth
        HovEtotal_3x3_ET_Depth1_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth1_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth1_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth1_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth2_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth2_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth2_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth2_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth3_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth3_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth3_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth3_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth4_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth4_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth4_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth4_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth5_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth5_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth5_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth5_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth6_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth6_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth6_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth6_emu_Leading4->Write();
        HovEtotal_3x3_ET_Depth7_emu_Leading1->Write();
        HovEtotal_3x3_ET_Depth7_emu_Leading2->Write();
        HovEtotal_3x3_ET_Depth7_emu_Leading3->Write();
        HovEtotal_3x3_ET_Depth7_emu_Leading4->Write();

        HovEtotal_3x3_ET_Depth1_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth2_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth3_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth4_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth5_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth6_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth7_Barrel_emu->Write();

        HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu->Write();
        HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu->Write();

        HovEtotal_3x3_ET_Depth1_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth2_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth3_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth4_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth5_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth6_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth7_Endcap_emu->Write();

        HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu->Write();
        HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu->Write();

        HovEtotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu->Write();
        HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu->Write();
        HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu->Write();
        HovEtotal_3x3_DepthRatio_Gen_Barrel_3_4_ov_1_2_emu->Write();
        HovEtotal_3x3_DepthRatio_Gen_Barrel_4_ov_1_2_3_emu->Write();
        HovEtotal_3x3_DepthRatio_Gen_Endcap_4_7_ov_1_3_emu->Write();

        energyDepth_Ratio_Barrel_D4->Write();
        energyDepth_Ratio_Barrel_D3_D4->Write();
        energyDepth_Ratio_Endcap_D6_D7->Write();
        energyDepth_Ratio_Endcap_D4_D7->Write();
        energyDepth_Ratio_Gen_Barrel_D4->Write();
        energyDepth_Ratio_Gen_Barrel_D3_D4->Write();
        energyDepth_Ratio_Gen_Endcap_D6_D7->Write();
        energyDepth_Ratio_Gen_Endcap_D4_D7->Write();

        energyDepth_NTPs_Barrel_D3_D4_ge06->Write();
        energyDepth_NTPs_Endcap_D4_D7_ge045->Write();
        energyDepth_NTPs_Gen_Barrel_D3_D4_ge06->Write();
        energyDepth_NTPs_Gen_Endcap_D4_D7_ge045->Write();


        L1JetTPDR->Write();
        L1JetNTP->Write();

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
