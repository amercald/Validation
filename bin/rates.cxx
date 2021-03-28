//Script for calculating %rate histograms
// Originally from Aaron Bundock
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
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

double PhiFromVertex(double vx, double vy) {
    double phi = atan2(vy, vx);
    return phi;
}  

double EtaFromVertex(double vx, double vy, double vz) {
    double theta = acos( vz / sqrt(vx*vx + vy*vy + vz*vz) );
    double eta = -log(tan(theta/2.));
    return eta;
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

bool compareTP(std::pair<int, double> &a, std::pair<int, double> &b)
{
    return a.second > b.second;
}

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
    int nTpBins = 256;
    float tpLo = 0.;
    float tpHi = 128.;

    std::string axR = ";Threshold E_{T} (GeV);rate (Hz)";
    std::string axD = ";E_{T} (GeV);events/bin";

    //make histos
    TH1F* singleJetRates_emu = new TH1F("singleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_emu = new TH1F("singleJetRates_HoE_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_TP05_emu = new TH1F("singleJetRates_HoE_TP05_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_TP5_emu = new TH1F("singleJetRates_HoE_TP5_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_Ratio02_emu = new TH1F("singleJetRates_HoE_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_Ratio06_emu = new TH1F("singleJetRates_HoE_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_TP2_Ratio02_emu = new TH1F("singleJetRates_HoE_TP2_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_TP2_Ratio06_emu = new TH1F("singleJetRates_HoE_TP2_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleJetRates_HoE_TP1_Ratio06_emu = new TH1F("singleJetRates_HoE_TP1_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_emu = new TH1F("doubleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_emu = new TH1F("doubleJetRates_HoE_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_TP05_emu = new TH1F("doubleJetRates_HoE_TP05_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_TP5_emu = new TH1F("doubleJetRates_HoE_TP5_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_Ratio02_emu = new TH1F("doubleJetRates_HoE_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_Ratio06_emu = new TH1F("doubleJetRates_HoE_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_TP2_Ratio02_emu = new TH1F("doubleJetRates_HoE_TP2_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_TP2_Ratio06_emu = new TH1F("doubleJetRates_HoE_TP2_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* doubleJetRates_HoE_TP1_Ratio06_emu = new TH1F("doubleJetRates_HoE_TP1_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_emu = new TH1F("tripleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_emu = new TH1F("tripleJetRates_HoE_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_TP05_emu = new TH1F("tripleJetRates_HoE_TP05_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_TP5_emu = new TH1F("tripleJetRates_HoE_TP5_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_Ratio02_emu = new TH1F("tripleJetRates_HoE_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_Ratio06_emu = new TH1F("tripleJetRates_HoE_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_TP2_Ratio02_emu = new TH1F("tripleJetRates_HoE_TP2_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_TP2_Ratio06_emu = new TH1F("tripleJetRates_HoE_TP2_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* tripleJetRates_HoE_TP1_Ratio06_emu = new TH1F("tripleJetRates_HoE_TP1_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_emu = new TH1F("quadJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_emu = new TH1F("quadJetRates_HoE_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_TP05_emu = new TH1F("quadJetRates_HoE_TP05_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_TP5_emu = new TH1F("quadJetRates_HoE_TP5_emu", axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_Ratio02_emu = new TH1F("quadJetRates_HoE_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_Ratio06_emu = new TH1F("quadJetRates_HoE_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_TP2_Ratio02_emu = new TH1F("quadJetRates_HoE_TP2_Ratio02_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_TP2_Ratio06_emu = new TH1F("quadJetRates_HoE_TP2_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* quadJetRates_HoE_TP1_Ratio06_emu = new TH1F("quadJetRates_HoE_TP1_Ratio06_emu",axR.c_str(), nJetBins, jetLo, jetHi);
    TH1F* singleEgRates_emu = new TH1F("singleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
    TH1F* doubleEgRates_emu = new TH1F("doubleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
    TH1F* singleTauRates_emu = new TH1F("singleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
    TH1F* doubleTauRates_emu = new TH1F("doubleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
    TH1F* singleISOEgRates_emu = new TH1F("singleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
    TH1F* doubleISOEgRates_emu = new TH1F("doubleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
    TH1F* singleISOTauRates_emu = new TH1F("singleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
    TH1F* doubleISOTauRates_emu = new TH1F("doubleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
    TH1F* htSumRates_emu = new TH1F("htSumRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_emu = new TH1F("htSumRates_HoE_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP05_emu = new TH1F("htSumRates_HoE_TP05_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP1_emu = new TH1F("htSumRates_HoE_TP1_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_emu = new TH1F("htSumRates_HoE_TP2_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP3_emu = new TH1F("htSumRates_HoE_TP3_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP4_emu = new TH1F("htSumRates_HoE_TP4_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP5_emu = new TH1F("htSumRates_HoE_TP5_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP05_ORHT360_emu = new TH1F("htSumRates_HoE_TP05_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP1_ORHT360_emu = new TH1F("htSumRates_HoE_TP1_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_ORHT360_emu = new TH1F("htSumRates_HoE_TP2_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP3_ORHT360_emu = new TH1F("htSumRates_HoE_TP3_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP4_ORHT360_emu = new TH1F("htSumRates_HoE_TP4_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP5_ORHT360_emu = new TH1F("htSumRates_HoE_TP5_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio02_emu = new TH1F("htSumRates_HoE_Ratio02_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio04_emu = new TH1F("htSumRates_HoE_Ratio04_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio06_emu = new TH1F("htSumRates_HoE_Ratio06_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio08_emu = new TH1F("htSumRates_HoE_Ratio08_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio02_ORHT360_emu = new TH1F("htSumRates_HoE_Ratio02_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio04_ORHT360_emu = new TH1F("htSumRates_HoE_Ratio04_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio06_ORHT360_emu = new TH1F("htSumRates_HoE_Ratio06_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_Ratio08_ORHT360_emu = new TH1F("htSumRates_HoE_Ratio08_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio02_emu = new TH1F("htSumRates_HoE_TP2_Ratio02_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio04_emu = new TH1F("htSumRates_HoE_TP2_Ratio04_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio06_emu = new TH1F("htSumRates_HoE_TP2_Ratio06_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio08_emu = new TH1F("htSumRates_HoE_TP2_Ratio08_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio02_ORHT360_emu = new TH1F("htSumRates_HoE_TP2_Ratio02_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio04_ORHT360_emu = new TH1F("htSumRates_HoE_TP2_Ratio04_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio06_ORHT360_emu = new TH1F("htSumRates_HoE_TP2_Ratio06_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    TH1F* htSumRates_HoE_TP2_Ratio08_ORHT360_emu = new TH1F("htSumRates_HoE_TP2_Ratio08_ORHT360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);

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

    TH1F* hcalTP_Barrel_emu = new TH1F("hcalTP_Barrel_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
    TH1F* hcalTP_Endcap_emu = new TH1F("hcalTP_Endcap_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

    TH1F* hcalTP_nearL1Jet_emu = new TH1F("hcalTP_nearL1Jet_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
    TH1F* hcalTP_nearL1Jet_Gen_emu = new TH1F("hcalTP_nearL1Jet_Gen_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

    TH1F* hcalTP_nearL1Jet_Endcap_emu = new TH1F("hcalTP_nearL1Jet_Endcap_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
    TH1F* hcalTP_nearL1Jet_Gen_Endcap_emu = new TH1F("hcalTP_nearL1Jet_Gen_Endcap_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

    TH1F* hcalTP_nearL1Jet_Barrel_emu = new TH1F("hcalTP_nearL1Jet_Barrel_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
    TH1F* hcalTP_nearL1Jet_Gen_Barrel_emu = new TH1F("hcalTP_nearL1Jet_Gen_Barrel_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

    TH1F* hHTSum_emu = new TH1F("hHTSum_emu", ";H_{T}; # Entries", 240, 0, 1200);
    TH1F* hHTSum_Gen_emu = new TH1F("hHTSum_Gen_emu", ";H_{T}; # Entries", 240, 0, 1200);

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

    TH1F* HovEtotal_1x1_emu = new TH1F("HovEtotal_1x1_emu", "HCAL energy / ECAL+HCAL energy for Jets (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu = new TH1F("HovEtotal_3x3_emu", "HCAL energy / ECAL+HCAL energy for Jets (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_HETge50_emu = new TH1F("HovEtotal_3x3_HETge50_emu", "HCAL energy / ECAL+HCAL energy for Jets (3x3);H/E;# Entries", 48,0,1.2);

    TH1F* HovEtotal_1x1_emu_Barrel = new TH1F("HovEtotal_1x1_emu_Barrel", "HCAL energy / ECAL+HCAL energy for Jets (1x1) in Barrel;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Barrel = new TH1F("HovEtotal_3x3_emu_Barrel", "HCAL energy / ECAL+HCAL energy for Jets (3x3) in Barrel;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_Endcap = new TH1F("HovEtotal_1x1_emu_Endcap", "HCAL energy / ECAL+HCAL energy for Jets (1x1) in Endcap;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Endcap = new TH1F("HovEtotal_3x3_emu_Endcap", "HCAL energy / ECAL+HCAL energy for Jets (3x3) in Endcap;H/E;# Entries", 48,0,1.2);


    TH1F* HovEtotalLog_1x1_emu = new TH1F("HovEtotalLog_1x1_emu", "log(H/E) for Jets (1x1);log_{10}(H/E);# Entries", 40,-10,10);
    TH1F* HovEtotalLog_3x3_emu = new TH1F("HovEtotalLog_3x3_emu", "log(H/E) for Jets (1x1);log_{10}(H/E) (3x3);H/E;# Entries", 40,-10,10);


    TH1F* HovEtotal_1x1_emu_Leading1 = new TH1F("HovEtotal_1x1_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading2 = new TH1F("HovEtotal_1x1_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading3 = new TH1F("HovEtotal_1x1_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_Leading4 = new TH1F("HovEtotal_1x1_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading1 = new TH1F("HovEtotal_3x3_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading2 = new TH1F("HovEtotal_3x3_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading3 = new TH1F("HovEtotal_3x3_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_Leading4 = new TH1F("HovEtotal_3x3_emu_Leading4", "HCAL energy / ECAL+HCAL energy for Leading Jet 4 (3x3);H/E;# Entries", 48,0,1.2);
    //  TH1F* HovEtotal_5x5_emu_Leading1 = new TH1F("HovEtotal_5x5_emu_Leading1", "HCAL energy / ECAL+HCAL energy for Leading Jet 1 (5x5);H/E;# Entries", 48,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading2 = new TH1F("HovEtotal_5x5_emu_Leading2", "HCAL energy / ECAL+HCAL energy for Leading Jet 2 (5x5);H/E;# Entries", 48,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading3 = new TH1F("HovEtotal_5x5_emu_Leading3", "HCAL energy / ECAL+HCAL energy for Leading Jet 3 (5x5);H/E;# Entries", 48,0,1.2);
    //TH1F* HovEtotal_5x5_emu_Leading4 = new TH1F("HovEtotal_5x5_emu_

    TH1F* HovEtotal_1x1_emu_GenMatchedJets = new TH1F("HovEtotal_1x1_emu_GenMatchedJets", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets = new TH1F("HovEtotal_3x3_emu_GenMatchedJets", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_HETge50_GenMatchedJets = new TH1F("HovEtotal_3x3_emu_HETge50_GenMatchedJets", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (3x3);H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_GenMatchedJets_Barrel = new TH1F("HovEtotal_1x1_emu_GenMatchedJets_Barrel", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Barrel;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets_Barrel = new TH1F("HovEtotal_3x3_emu_GenMatchedJets_Barrel", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (3x3) in Barrel;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_1x1_emu_GenMatchedJets_Endcap = new TH1F("HovEtotal_1x1_emu_GenMatchedJets_Endcap", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (1x1) in Endcap;H/E;# Entries", 48,0,1.2);
    TH1F* HovEtotal_3x3_emu_GenMatchedJets_Endcap = new TH1F("HovEtotal_3x3_emu_GenMatchedJets_Endcap", "HCAL energy / ECAL+HCAL energy for Gen Matched Jets 1 (3x3) in Endcap;H/E;# Entries", 48,0,1.2);



    //Binning in ET
    TH2F* HovEtotal_1x1_ET_emu_Leading1 = new TH2F("HovEtotal_1x1_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (1x1);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading2 = new TH2F("HovEtotal_1x1_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (1x1);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading3 = new TH2F("HovEtotal_1x1_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (1x1);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_1x1_ET_emu_Leading4 = new TH2F("HovEtotal_1x1_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (1x1);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (3x3);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (3x3);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (3x3);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (3x3);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    //  TH2F* HovEtotal_5x5_ETemu_Leading1 = new TH2F("HovEtotal_5x5_ETemu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 1 (5x5);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading2 = new TH2F("HovEtotal_5x5_ETemu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 2 (5x5);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading3 = new TH2F("HovEtotal_5x5_ETemu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 3 (5x5);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    //TH2F* HovEtotal_5x5_ETemu_Leading4 = new TH2F("HovEtotal_5x5_ETemu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET for Leading Jet 4 (5x5);H/E;# Entries", 10, 0, 500, 48, 0, 1.2);
    
    //Binning in depth energy
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth1_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth1_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth2_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth3_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth4_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth5_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth6_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading1 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading1", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading2 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading2", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading3 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading3", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_emu_Leading4 = new TH2F("HovEtotal_3x3_ET_Depth7_emu_Leading4", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 for Leading Jet 4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);

    TH2F* HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth1 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth2 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth3 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth4 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth5 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth6 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);
    TH2F* HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu = new TH2F("HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu", "HCAL energy / ECAL+HCAL energy vs ET_Depth7 (3x3);TP E_{T};H/(H+E)", 60, 0, 60, 48, 0, 1.2);

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

    TH2F* HEEnergytotal_1x1_emu_GenMatched = new TH2F("HEEnergytotal_1x1_emu_GenMatched", "HCAL vs ECAL Energy for Matched Jet", 100, 0, 300, 100, 0, 300);
    TH2F* HEEnergytotal_3x3_emu_GenMatched = new TH2F("HEEnergytotal_3x3_emu_GenMatched", "HCAL vs ECAL Energy for Matched Jet", 100, 0, 300, 100, 0, 300);

    //Energy depth study plots

    TH2F* energyDepth_Barrel = new TH2F("energyDepth_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_Endcap = new TH2F("energyDepth_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge1_Barrel = new TH2F("energyDepth_TPge1_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge1_Endcap = new TH2F("energyDepth_TPge1_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge2_Barrel = new TH2F("energyDepth_TPge2_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge2_Endcap = new TH2F("energyDepth_TPge2_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge3_Barrel = new TH2F("energyDepth_TPge3_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge3_Endcap = new TH2F("energyDepth_TPge3_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);


    TH2F* energyDepth_TPge5_Barrel = new TH2F("energyDepth_TPge5_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge5_Endcap = new TH2F("energyDepth_TPge5_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);


    TH2F* energyDepth_HT120_Barrel = new TH2F("energyDepth_HT120_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_HT120_Endcap = new TH2F("energyDepth_HT120_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_L1_Barrel = new TH2F("energyDepth_L1_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_L1_Endcap = new TH2F("energyDepth_L1_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge5_L1_Barrel = new TH2F("energyDepth_TPge5_L1_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge5_L1_Endcap = new TH2F("energyDepth_TPge5_L1_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_HT120_L1_Barrel = new TH2F("energyDepth_HT120_L1_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_HT120_L1_Endcap = new TH2F("energyDepth_HT120_L1_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge5_LowRatio_Barrel = new TH2F("energyDepth_TPge5_LowRatio_Barrel", "Depth profile, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge5_LowRatio_Endcap = new TH2F("energyDepth_TPge5_LowRatio_Endcap", "Depth profile, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_LowRatio_Barrel = new TH2F("energyDepth_LowRatio_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_LowRatio_Endcap = new TH2F("energyDepth_LowRatio_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPge5_HoEcut_Barrel = new TH2F("energyDepth_TPge5_HoEcut_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_TPge5_HoEcut_Endcap = new TH2F("energyDepth_TPge5_HoEcut_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_HoEcut_Barrel = new TH2F("energyDepth_HoEcut_Barrel", "Depth profile, inclusive, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_HoEcut_Endcap = new TH2F("energyDepth_HoEcut_Endcap", "Depth profile, inclusive, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPE_Barrel = new TH2F("energyDepth_TPE_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_Endcap = new TH2F("energyDepth_TPE_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_TPge5_Barrel = new TH2F("energyDepth_TPE_TPge5_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_TPge5_Endcap = new TH2F("energyDepth_TPE_TPge5_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);

    TH2F* energyDepth_TPE_HoEcut_Barrel = new TH2F("energyDepth_TPE_HoEcut_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_HoEcut_Endcap = new TH2F("energyDepth_TPE_HoEcut_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);

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

    TH2F* energyDepth_genMatchInclusive_TPge5_Barrel = new TH2F("energyDepth_genMatchInclusive_TPge5_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_TPge5_Endcap = new TH2F("energyDepth_genMatchInclusive_TPge5_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_HT120_Barrel = new TH2F("energyDepth_genMatchInclusive_HT120_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_HT120_Endcap = new TH2F("energyDepth_genMatchInclusive_HT120_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel = new TH2F("energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap = new TH2F("energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_HoEcut_Barrel = new TH2F("energyDepth_genMatchInclusive_HoEcut_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_HoEcut_Endcap = new TH2F("energyDepth_genMatchInclusive_HoEcut_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_L1_Barrel = new TH2F("energyDepth_genMatchInclusive_L1_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_L1_Endcap = new TH2F("energyDepth_genMatchInclusive_L1_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_TPge5_L1_Barrel = new TH2F("energyDepth_genMatchInclusive_TPge5_L1_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_TPge5_L1_Endcap = new TH2F("energyDepth_genMatchInclusive_TPge5_L1_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_HT120_L1_Barrel = new TH2F("energyDepth_genMatchInclusive_HT120_L1_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_HT120_L1_Endcap = new TH2F("energyDepth_genMatchInclusive_HT120_L1_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_TPge5_LowRatio_Barrel = new TH2F("energyDepth_genMatchInclusive_TPge5_LowRatio_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_TPge5_LowRatio_Endcap = new TH2F("energyDepth_genMatchInclusive_TPge5_LowRatio_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_genMatchInclusive_LowRatio_Barrel = new TH2F("energyDepth_genMatchInclusive_LowRatio_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchInclusive_LowRatio_Endcap = new TH2F("energyDepth_genMatchInclusive_LowRatio_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);


    TH2F* energyDepth_genMatchTP_Barrel = new TH2F("energyDepth_genMatchTP_Barrel", "Depth profile, matched TPs, in Barrel", 8, -0.5, 7.5, 60, 0, 1.2);
    TH2F* energyDepth_genMatchTP_Endcap = new TH2F("energyDepth_genMatchTP_Endcap", "Depth profile, matched TPs, in Endcap", 8, -0.5, 7.5, 60, 0, 1.2);

    TH2F* energyDepth_TPE_genMatchInclusive_Barrel = new TH2F("energyDepth_TPE_genMatchInclusive_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_genMatchInclusive_Endcap = new TH2F("energyDepth_TPE_genMatchInclusive_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);

    TH2F* energyDepth_TPE_genMatchInclusive_TPge5_Barrel = new TH2F("energyDepth_TPE_genMatchInclusive_TPge5_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_genMatchInclusive_TPge5_Endcap = new TH2F("energyDepth_TPE_genMatchInclusive_TPge5_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);

    TH2F* energyDepth_TPE_genMatchInclusive_HoEcut_Barrel = new TH2F("energyDepth_TPE_genMatchInclusive_HoEcut_Barrel", "Depth profile, matched jets, in Barrel", 8, -0.5, 7.5, 60, 0, 10);
    TH2F* energyDepth_TPE_genMatchInclusive_HoEcut_Endcap = new TH2F("energyDepth_TPE_genMatchInclusive_HoEcut_Endcap", "Depth profile, matched jets, in Endcap", 8, -0.5, 7.5, 60, 0, 10);

    TH1F* TPEt_HB = new TH1F("TPEt_HB", "E_{T} of TPs with Low Depth Energy Ratio in HB" ";TP E_{T} (GeV); Entries", 20, 0, 10);
    TH1F* TPEt_HE = new TH1F("TPEt_HE", "E_{T} of TPs with Low Depth Energy Ratio in HE" ";TP E_{T} (GeV); Entries", 20, 0, 10);

    TH1F* TPEt_Matched_HB = new TH1F("TPEt_Matched_HB", "E_{T} of TPs in HB" ";TP E_{T} (GeV); Entries", 20, 0, 10);
    TH1F* TPEt_Matched_HE = new TH1F("TPEt_Matched_HE", "E_{T} of TPs in HE" ";TP E_{T} (GeV); Entries", 20, 0, 10);

    TH1F* TPEt_LowRatio_HB = new TH1F("TPEt_LowRatio_HB", "E_{T} of TPs with Low Depth Energy Ratio in HB" ";TP E_{T} (GeV); Entries", 20, 0, 10);
    TH1F* TPEt_LowRatio_HE = new TH1F("TPEt_LowRatio_HE", "E_{T} of TPs with Low Depth Energy Ratio in HE" ";TP E_{T} (GeV); Entries", 20, 0, 10);

    TH1F* TPEt_LowRatio_Matched_HB = new TH1F("TPEt_LowRatio_Matched_HB", "E_{T} of TPs with Low Depth Energy Ratio in HB" ";TP E_{T} (GeV); Entries", 20, 0, 10);
    TH1F* TPEt_LowRatio_Matched_HE = new TH1F("TPEt_LowRatio_Matched_HE", "E_{T} of TPs with Low Depth Energy Ratio in HE" ";TP E_{T} (GeV); Entries", 20, 0, 10);

    TH1F* DeltaRLLP = new TH1F("DeltaRLLP", ";#Delta R; # Entries", 100, 0, 2);
    
    TH2F* energyDepth_NTPs_HBD4_HED47_Max = new TH2F("energyDepth_NTPs_HBD4_HED47_Max", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);
    TH2F* energyDepth_NTPs_HBD4_HED47_Max_Gen = new TH2F("energyDepth_NTPs_HBD4_HED47_Max_Gen", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);
    TH2F* energyDepth_NTPs_HBD4_HED47_Max_HT120 = new TH2F("energyDepth_NTPs_HBD4_HED47_Max_HT120", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);
    TH2F* energyDepth_NTPs_HBD4_HED47_Max_Gen_HT120 = new TH2F("energyDepth_NTPs_HBD4_HED47_Max_Gen_HT120", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);
    TH2F* energyDepth_NTPs_HBD4_HED47_Max_HT360 = new TH2F("energyDepth_NTPs_HBD4_HED47_Max_HT360", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);
    TH2F* energyDepth_NTPs_HBD4_HED47_Max_Gen_HT360 = new TH2F("energyDepth_NTPs_HBD4_HED47_Max_Gen_HT360", "Max NTPs above threshold" ";Threshold (GeV); # TPs", 50, 0, 10, 10, -0.5, 9.5);

    TH1F* energyDepth_HT120_DepthEnergy1_HB = new TH1F("energyDepth_HT120_DepthEnergy1_HB", "Energy in Depth Layer 1 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy1_HE = new TH1F("energyDepth_HT120_DepthEnergy1_HE", "Energy in Depth Layer 1 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy2_HB = new TH1F("energyDepth_HT120_DepthEnergy2_HB", "Energy in Depth Layer 2 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy2_HE = new TH1F("energyDepth_HT120_DepthEnergy2_HE", "Energy in Depth Layer 2 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy3_HB = new TH1F("energyDepth_HT120_DepthEnergy3_HB", "Energy in Depth Layer 3 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy3_HE = new TH1F("energyDepth_HT120_DepthEnergy3_HE", "Energy in Depth Layer 3 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy4_HB = new TH1F("energyDepth_HT120_DepthEnergy4_HB", "Energy in Depth Layer 4 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy4_HE = new TH1F("energyDepth_HT120_DepthEnergy4_HE", "Energy in Depth Layer 4 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy5_HE = new TH1F("energyDepth_HT120_DepthEnergy5_HE", "Energy in Depth Layer 5 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy6_HE = new TH1F("energyDepth_HT120_DepthEnergy6_HE", "Energy in Depth Layer 6 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_DepthEnergy7_HE = new TH1F("energyDepth_HT120_DepthEnergy7_HE", "Energy in Depth Layer 7 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy1_HB = new TH1F("energyDepth_HT120_Gen_DepthEnergy1_HB", "Energy in Depth Layer 1 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy1_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy1_HE", "Energy in Depth Layer 1 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy2_HB = new TH1F("energyDepth_HT120_Gen_DepthEnergy2_HB", "Energy in Depth Layer 2 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy2_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy2_HE", "Energy in Depth Layer 2 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy3_HB = new TH1F("energyDepth_HT120_Gen_DepthEnergy3_HB", "Energy in Depth Layer 3 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy3_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy3_HE", "Energy in Depth Layer 3 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy4_HB = new TH1F("energyDepth_HT120_Gen_DepthEnergy4_HB", "Energy in Depth Layer 4 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy4_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy4_HE", "Energy in Depth Layer 4 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy5_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy5_HE", "Energy in Depth Layer 5 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy6_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy6_HE", "Energy in Depth Layer 6 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_HT120_Gen_DepthEnergy7_HE = new TH1F("energyDepth_HT120_Gen_DepthEnergy7_HE", "Energy in Depth Layer 7 in HE" ";Energy (GeV); Entries", 100, 0, 10);

    TH1F* energyDepth_DepthEnergy1_HB = new TH1F("energyDepth_DepthEnergy1_HB", "Energy in Depth Layer 1 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy1_HE = new TH1F("energyDepth_DepthEnergy1_HE", "Energy in Depth Layer 1 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy2_HB = new TH1F("energyDepth_DepthEnergy2_HB", "Energy in Depth Layer 2 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy2_HE = new TH1F("energyDepth_DepthEnergy2_HE", "Energy in Depth Layer 2 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy3_HB = new TH1F("energyDepth_DepthEnergy3_HB", "Energy in Depth Layer 3 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy3_HE = new TH1F("energyDepth_DepthEnergy3_HE", "Energy in Depth Layer 3 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy4_HB = new TH1F("energyDepth_DepthEnergy4_HB", "Energy in Depth Layer 4 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy4_HE = new TH1F("energyDepth_DepthEnergy4_HE", "Energy in Depth Layer 4 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy5_HE = new TH1F("energyDepth_DepthEnergy5_HE", "Energy in Depth Layer 5 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy6_HE = new TH1F("energyDepth_DepthEnergy6_HE", "Energy in Depth Layer 6 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_DepthEnergy7_HE = new TH1F("energyDepth_DepthEnergy7_HE", "Energy in Depth Layer 7 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy1_HB = new TH1F("energyDepth_Gen_DepthEnergy1_HB", "Energy in Depth Layer 1 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy1_HE = new TH1F("energyDepth_Gen_DepthEnergy1_HE", "Energy in Depth Layer 1 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy2_HB = new TH1F("energyDepth_Gen_DepthEnergy2_HB", "Energy in Depth Layer 2 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy2_HE = new TH1F("energyDepth_Gen_DepthEnergy2_HE", "Energy in Depth Layer 2 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy3_HB = new TH1F("energyDepth_Gen_DepthEnergy3_HB", "Energy in Depth Layer 3 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy3_HE = new TH1F("energyDepth_Gen_DepthEnergy3_HE", "Energy in Depth Layer 3 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy4_HB = new TH1F("energyDepth_Gen_DepthEnergy4_HB", "Energy in Depth Layer 4 in HB" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy4_HE = new TH1F("energyDepth_Gen_DepthEnergy4_HE", "Energy in Depth Layer 4 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy5_HE = new TH1F("energyDepth_Gen_DepthEnergy5_HE", "Energy in Depth Layer 5 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy6_HE = new TH1F("energyDepth_Gen_DepthEnergy6_HE", "Energy in Depth Layer 6 in HE" ";Energy (GeV); Entries", 100, 0, 10);
    TH1F* energyDepth_Gen_DepthEnergy7_HE = new TH1F("energyDepth_Gen_DepthEnergy7_HE", "Energy in Depth Layer 7 in HE" ";Energy (GeV); Entries", 100, 0, 10);

    TH1F* energyDepth_Ratio_Gen_Depth4_HB = new TH1F("energyDepth_Ratio_Gen_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_Gen_Depth4_7_HE = new TH1F("energyDepth_Ratio_Gen_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_Depth4_HB = new TH1F("energyDepth_Ratio_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_Depth4_7_HE = new TH1F("energyDepth_Ratio_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_Gen_Depth3_4_HB = new TH1F("energyDepth_Ratio_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_Gen_Depth2_7_HE = new TH1F("energyDepth_Ratio_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_Depth3_4_HB = new TH1F("energyDepth_Ratio_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_Depth2_7_HE = new TH1F("energyDepth_Ratio_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
 
    TH1F* energyDepth_Ratio_HT120_Gen_Depth4_HB = new TH1F("energyDepth_Ratio_HT120_Gen_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_Gen_Depth4_7_HE = new TH1F("energyDepth_Ratio_HT120_Gen_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_HT120_Depth4_HB = new TH1F("energyDepth_Ratio_HT120_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_Depth4_7_HE = new TH1F("energyDepth_Ratio_HT120_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_HT120_Gen_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_Gen_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_HT120_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_HT120_TP1_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_TP1_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP1_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_TP1_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP2_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_TP2_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP2_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_TP2_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP3_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_TP3_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP3_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_TP3_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP4_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_TP4_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP4_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_TP4_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP5_Depth3_4_HB = new TH1F("energyDepth_Ratio_HT120_TP5_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_HT120_TP5_Depth2_7_HE = new TH1F("energyDepth_Ratio_HT120_TP5_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);


    TH1F* energyDepth_Ratio_TPge10_Gen_Depth4_HB = new TH1F("energyDepth_Ratio_TPge10_Gen_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_TPge10_Gen_Depth4_7_HE = new TH1F("energyDepth_Ratio_TPge10_Gen_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_TPge10_Depth4_HB = new TH1F("energyDepth_Ratio_TPge10_DepthEnergy4_HB", "Ratio of Energy in HB Depth Layers 3 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_TPge10_Depth4_7_HE = new TH1F("energyDepth_Ratio_TPge10_DepthEnergy4_7_HE", "Ratio of Energy in HE Depth Layers 4 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_TPge10_Gen_Depth3_4_HB = new TH1F("energyDepth_Ratio_TPge10_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_TPge10_Gen_Depth2_7_HE = new TH1F("energyDepth_Ratio_TPge10_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);

    TH1F* energyDepth_Ratio_TPge10_Depth3_4_HB = new TH1F("energyDepth_Ratio_TPge10_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH1F* energyDepth_Ratio_TPge10_Depth2_7_HE = new TH1F("energyDepth_Ratio_TPge10_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2);
    TH2F* energyDepth_Ratio_IEta_Gen_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_Gen_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);

    TH2F* energyDepth_Ratio_IEta_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);


    TH2F* energyDepth_Ratio_IEta_HT120_Gen_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_HT120_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_HT120_Gen_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_HT120_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);

    TH2F* energyDepth_Ratio_IEta_HT120_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_HT120_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_HT120_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_HT120_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);


    TH2F* energyDepth_Ratio_IEta_TPge10_Gen_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_TPge10_Gen_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_TPge10_Gen_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_TPge10_Gen_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);

    TH2F* energyDepth_Ratio_IEta_TPge10_Depth3_4_HB = new TH2F("energyDepth_Ratio_IEta_TPge10_Depth3_4_HB", "Ratio of Energy in HB Depth Layers 2 - 4 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 16, 0.5, 16.5);
    TH2F* energyDepth_Ratio_IEta_TPge10_Depth2_7_HE = new TH2F("energyDepth_Ratio_IEta_TPge10_Depth2_7_HE", "Ratio of Energy in HE Depth Layers 2 - 7 to all Layers" ";Depth Energy Ratio; Entries", 120, 0, 1.2, 13, 16.5, 29.5);


    //Efficiencies for JetIDs
    TH1F * effJetID_HoE_DepthFB_TPE = new TH1F("effJetID_HoE_DepthFB_TPE", "", 14, 0, 14);
    TH1F * effJetID_HoE_DepthFB_TPE_Gen = new TH1F("effJetID_HoE_DepthFB_TPE_Gen", "", 14, 0, 14);

    TH1F * effJetID_HoE_DepthFB_HTscan = new TH1F("effJetID_HoE_DepthFB_HTscan", "", 32, 0, 32);
    TH1F * effJetID_HoE_DepthFB_Gen_HTscan = new TH1F("effJetID_HoE_DepthFB_Gen_HTscan", "", 32, 0, 32);

    TH1F * effJetID_HoE_DepthFB_Ratio = new TH1F("effJetID_HoE_DepthFB_Ratio", "", 19, 0, 19);
    TH1F * effJetID_HoE_DepthFB_Ratio_Gen = new TH1F("effJetID_HoE_DepthFB_Ratio_Gen", "", 19, 0, 19);

    TH1F * effJetID_HoE_DepthFB_Jetscan = new TH1F("effJetID_HoE_DepthFB_Jetscan", "", 32, 0, 32);
    TH1F * effJetID_HoE_DepthFB_Gen_Jetscan = new TH1F("effJetID_HoE_DepthFB_Gen_Jetscan", "", 32, 0, 32);

    // Gen Matching verification plots
    TH1F * hJetGenPartDR_LLPdaught = new TH1F("hJetGenPartDR_LLPdaugh",";#DeltaR;",100,0,5);
    TH1F * hJetGenPartDR_LLP = new TH1F("hJetGenPartDR_LLP",";#DeltaR;",100,0,5);
    TH1F * hJetGenPartDR_LLP_inHCAL = new TH1F("hJetGenPartDR_LLP_inHCAL",";#DeltaR;",100,0,5);
    TH1F * hJetGenPartDRfromVertex_LLP_inHCAL = new TH1F("hJetGenPartDRfromVertex_LLP_inHCAL",";#DeltaR;",100,0,5);

    TH1F * hLLP_vertex_DR = new TH1F("hLLP_vertex_DR",";#DeltaR;",100,0,5);

    TH1F * hNMatchedLLP_inHCAL_DR02 = new TH1F("hNMatchedLLP_inHCAL_DR02", "# LLP in HCAL Matched (DR < 0.2) to Jet;# LLP;",5, -0.5, 4.5);
    TH1F * hNMatchedLLP_inHCAL_DR05 = new TH1F("hNMatchedLLP_inHCAL_DR05", "# LLP in HCAL Matched (DR < 0.5) to Jet;# LLP;",5, -0.5, 4.5);
    TH1F * hNMatchedLLP_DR02 = new TH1F("hNMatchedLLP_DR02", "# LLPs Matched (DR < 0.2) to Jet;# LLP;",5, -0.5, 4.5);
    TH1F * hNMatchedLLP_DR05 = new TH1F("hNMatchedLLP_DR05", "# LLPs Matched (DR < 0.5) to Jet;# LLP;",5, -0.5, 4.5);
    TH1F * hNLLPdaughts = new TH1F("hnLLPdaughts",";# LLP daughters;",5, -0.5, 4.5);
    TH1F * hNLLPdaughts_pteta = new TH1F("hnLLPdaughts_pteta",";# LLP daughters;",5, -0.5, 4.5);
    TH1F * hNLLPdaughts_inHCAL = new TH1F("hnLLPdaughts_inHCAL",";# LLP daughters in HCAL;",5, -0.5, 4.5);
    TH1F * hNMatched_LLPdaught_DR02 = new TH1F("hNMatched_LLPdaught_DR02","# LLP Daughters in HCAL Matched (DR < 0.2) to Jet;# LLP Daughters;",5, -0.5, 4.5);
    TH1F * hNMatched_LLPdaught_DR05 = new TH1F("hNMatched_LLPdaught_DR05","# LLP Daughters (DR < 0.5) to Jet;# LLP Daughters;",5, -0.5, 4.5);
    TH1F * hfracMatched_LLPdaught_DR02 = new TH1F("hfracMatched_LLPdaught_DR02","Fraction of LLP Daughters in HCAL Matched (DR < 0.2) to Jet;Fraction Matched;",5,0,1.25);
    TH1F * hfracMatched_LLPdaught_DR05 = new TH1F("hfracMatched_LLPdaught_DR05","Fraction of LLP Daughters in HCAL Matched (DR < 0.5) to Jet;Fraction Matched;",5,0,1.25);

    TH1F* fracLLP_Separated = new TH1F("fracLLP_Separated", ";HCAL Region; # Entries", 3, 0, 3);
    TH1F* betagammaLLP = new TH1F("betagammaLLP", ";#beta#gamma; # Entries", 100, 0, 10);
    TH1F* velocityLLP = new TH1F("velocityLLP", ";v (m/s); # Entries", 100, 0, 10000000);
    TH1F* betaLLP = new TH1F("betaLLP", ";#beta; # Entries", 100, 0, 1);
    TH1F* vertex_HB_LLPD = new TH1F("vertex_HB_LLPD", ";Radius (cm); # Entries", 100, 0, 500);
    TH1F* vertex_HE_LLPD = new TH1F("vertex_HE_LLPD", ";V_{z} (cm); # Entries", 100, 0, 800);
    //text file for writing some events

    std::ofstream TPtextfile;
    TPtextfile.open("jet_tp_info_"+sampleType+".txt");
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
            double tpEt(0.), tpiEta(0.);
      
            for(int i=0; i < l1TPemu_->nHCALTP; i++){
                tpEt = l1TPemu_->hcalTPet[i];
                tpiEta = l1TPemu_->hcalTPieta[i];
                hcalTP_emu->Fill(tpEt);
                if (abs(tpiEta) <= 16) hcalTP_Barrel_emu->Fill(tpEt);
                if (abs(tpiEta) >= 16 && abs(tpiEta) <= 29) hcalTP_Endcap_emu->Fill(tpEt);
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

//            if (htSum < 120) continue;
            //Define good gen particles i.e. LLP daughters decaying in HCAL
            int nHCALTP = l1CaloTPemu_->nHCALTP;
            int nJetemu = l1emu_->nJets;
      
            int nGenPart = generator_->nPart;
            std::vector<bool> goodGenParticles(nGenPart, false);
            std::vector<bool> goodGenParticles_HB(nGenPart, false);
            std::vector<bool> goodGenParticles_HE(nGenPart, false);
            std::vector<double> genParticlesEta, genParticlesPhi;

            int nGoodGen = 0;
            int nOkayGen = 0;
            int nHad = 0;
            int nMatched = 0;
            int nLLPdaught = 0;
            int nLLPdaught_pteta = 0;

            double LLPcounter = 0;
            std::pair<int, std::string> first_vertex_LLP (-1, "None");
            std::pair<int, std::string> second_vertex_LLP (-1, "None");

            std::vector<int> LLPD_Idx;

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
                bool inEndcap = abs(Vz) > 388 && abs(Vz) < 568 && radius < 295;
                bool inBarrel = abs(Vz) < 388 && radius > 179 && radius < 295;
                bool inHCAL = inEndcap || inBarrel;
                bool isQuarkorGluon = abs(pdgId) <=5  || abs(pdgId) == 21;
                bool isLLP = pdgId == 6000113;
                //if (generator_->partHardProcess[genpart] != 0) std::cout << pdgId << " " << generator_->partParent[genpart] << std::endl;
                std::vector<double> intersection = intersect(Vx, Vy, Vz, Px, Py, Pz);
                genParticlesEta.push_back(intersection.at(0));
                genParticlesPhi.push_back(intersection.at(1));

//                std::cout << Vx << " " << Vy << " " << Vz << " " << PhiFromVertex(Vx, Vy) << " " << EtaFromVertex(Vx, Vy, Vz) << std::endl; 
                if (isQuarkorGluon && inHCAL && generator_->partHardProcess[genpart] != 0 && generator_->partParent[genpart] == 6000113)// && abs(intersection.at(0)) < 2.4)
                {	    
                    goodGenParticles.at(genpart) = true;
                    nGoodGen += 1;
                    if (inBarrel) goodGenParticles_HB.at(genpart) = true;
                    if (inEndcap) goodGenParticles_HE.at(genpart) = true;
                    LLPD_Idx.push_back(genpart);
                }
                if (inHCAL) nOkayGen += 1;
                if (isQuarkorGluon) nHad += 1;
                if (isQuarkorGluon && generator_->partHardProcess[genpart] != 0 && generator_->partParent[genpart] == 6000113) 
                {
                    nLLPdaught += 1;
                    if (generator_->partVz[genpart-1] != Vz)
                    { 
                        if (abs(generator_->partEta[genpart]) < 1.3)
                        {
                            vertex_HB_LLPD->Fill(radius);
                            fracLLP_Separated->Fill(0.5);
                        }
                        else if (abs(generator_->partEta[genpart]) > 1.3 && abs(generator_->partEta[genpart]) < 3)
                        {
                            vertex_HE_LLPD->Fill(Vz);
                            fracLLP_Separated->Fill(1.5);
                        }
                        else
                        {
                            fracLLP_Separated->Fill(2.5);
                        }
                    }
                }
                if (isQuarkorGluon && generator_->partHardProcess[genpart] != 0 && generator_->partParent[genpart] == 6000113 && abs(intersection.at(0)) < 3 && generator_->partPt[genpart] > 20) nLLPdaught_pteta += 1;
                if (isLLP)
                {
                    double total_momentum = sqrt(Pz*Pz + Px*Px + Py*Py);
                    double lightspeed = 29979245800; // cm/s
                    double beta = sqrt(Pz*Pz + Px*Px + Py*Py)/energy;
                    double gamma = 1./sqrt(1. - beta*beta);
                    double betagamma = beta*gamma;
                    double velocity = beta*lightspeed;
                    betagammaLLP->Fill(betagamma);
                    velocityLLP->Fill(velocity/100);
                    betaLLP->Fill(beta);
                    LLPcounter += 1;                    
                }
            }
//            std::cout << first_vertex_LLP.first << " " << first_vertex_LLP.second << " " << second_vertex_LLP.first << " " << second_vertex_LLP.second << std::endl;
            //std::cout << LLPcounter << std::endl;
            h_nGenParticles->Fill(nGoodGen);            
            //make a more easily iterable vector of depths
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
            std::vector<bool> matchedJet(nJetemu, false);
            std::vector<bool> matchedJet_HB(nJetemu, false);
            std::vector<bool> matchedJet_HE(nJetemu, false);

            std::vector<bool> matchedJet_noVertex(nJetemu, false);
            
            for(int jetIt = 0; jetIt < nJetemu; jetIt++)
            {
                double jetEta = l1emu_->jetEta[jetIt];
                double jetPhi = l1emu_->jetPhi[jetIt];
                bool nearLLP = false;
                bool nearVertex = false;

                for(int genpart = 0; genpart < nGenPart; genpart++)
                {
                    double genEta = generator_->partEta[genpart];
                    double genPhi = generator_->partPhi[genpart];
                    double Vz = generator_->partVz[genpart];
                    double Vx = generator_->partVx[genpart];
                    double Vy = generator_->partVy[genpart];
                    double vertexPhi = PhiFromVertex(Vx, Vy);
                    double vertexEta = EtaFromVertex(Vx, Vy, Vz);
                    double genpartDeltaR = DeltaR(jetPhi, genPhi, jetEta, genEta);
                    double vertexDeltaR  = DeltaR(jetPhi, vertexPhi, jetEta, vertexEta);
                    if (nGoodGen > 0 && genpartDeltaR < 0.2 && generator_->partId[genpart] == 6000113) 
                    {
                        nearLLP = true;
                    }
                    if (goodGenParticles.at(genpart) && vertexDeltaR < 0.2)
                    {
                        nearVertex = true;
                        if (goodGenParticles_HB.at(genpart)) matchedJet_HB.at(jetIt) = true;
                        if (goodGenParticles_HE.at(genpart)) matchedJet_HE.at(jetIt) = true;
                    } 
                }
                if (nearLLP && nearVertex) matchedJet.at(jetIt) = true;
                if (matchedJet.at(jetIt)) nMatched += 1;
                if (nearLLP) matchedJet_noVertex.at(jetIt) = true;
            }
            //FOR GEN MATCHING STUDY

            int nMatchedGoodGen_DR02 = 0, nMatchedGoodGen_DR05 = 0;
            int nMatchedLLP_DR02 = 0, nMatchedLLP_DR05 = 0, nMatchedLLP_inHCAL_DR02 = 0,  nMatchedLLP_inHCAL_DR05 = 0;
            for(int genpart = 0; genpart < nGenPart; genpart++)
            {                
                if (abs(genParticlesEta.at(genpart)) > 3 || generator_->partPt[genpart] < 20) continue;
                double genEta = generator_->partEta[genpart];
                double genPhi = generator_->partPhi[genpart];
                double Vz = generator_->partVz[genpart];
                double Vx = generator_->partVx[genpart];
                double Vy = generator_->partVy[genpart];
                double vertexPhi = PhiFromVertex(Vx, Vy);
                double vertexEta = EtaFromVertex(Vx, Vy, Vz);
                double minDR_LLP = 100, minDR_LLP_inHCAL = 100, minDRfromVertex_LLP_inHCAL = 100, minDR_LLPdaught = 100;
                for(int jetIt = 0; jetIt < nJetemu; jetIt++)
                {
                    double jetEta = l1emu_->jetEta[jetIt];
                    double jetPhi = l1emu_->jetPhi[jetIt];
                    if (abs(jetEta) > 3 || l1emu_->jetEt[jetIt] < 20) continue;
                    double genpartDeltaR = DeltaR(jetPhi, genPhi, jetEta, genEta);
                    double genpartDeltaR_corrected = DeltaR(jetPhi, genParticlesPhi.at(genpart), jetEta, genParticlesEta.at(genpart)); 
                    double vertexDeltaR  = DeltaR(jetPhi, vertexPhi, jetEta, vertexEta);

                    if (generator_->partId[genpart] == 6000113)
                    {
                        if (genpartDeltaR < minDR_LLP) minDR_LLP = genpartDeltaR;
                        if (nGoodGen > 0 && genpartDeltaR < minDR_LLP_inHCAL)
                        { 
                            minDR_LLP_inHCAL = genpartDeltaR;
                        }
                    }
                    if (goodGenParticles.at(genpart))
                    {
                        if (matchedJet_noVertex.at(jetIt) && vertexDeltaR < minDRfromVertex_LLP_inHCAL) minDRfromVertex_LLP_inHCAL = vertexDeltaR;
                        if (genpartDeltaR_corrected < minDR_LLPdaught) minDR_LLPdaught = genpartDeltaR_corrected;
                    }
                }

                if (minDR_LLPdaught < 0.2) nMatchedGoodGen_DR02 += 1;
                if (minDR_LLPdaught < 0.5) nMatchedGoodGen_DR05 += 1;
                if (minDR_LLP < 0.2) nMatchedLLP_DR02 += 1;
                if (minDR_LLP < 0.5) nMatchedLLP_DR05 += 1;
                if (minDR_LLP_inHCAL < 0.2) nMatchedLLP_inHCAL_DR02 += 1;
                if (minDR_LLP_inHCAL < 0.5) nMatchedLLP_inHCAL_DR05 += 1;
                

                hJetGenPartDR_LLPdaught->Fill(minDR_LLPdaught);
                hJetGenPartDR_LLP->Fill(minDR_LLP);
                hJetGenPartDR_LLP_inHCAL->Fill(minDR_LLP_inHCAL);
                hJetGenPartDRfromVertex_LLP_inHCAL->Fill(minDRfromVertex_LLP_inHCAL);

                if (generator_->partId[genpart] == 6000113)
                {
                    double min_DR_LLP_vertex = 100;
                    for (uint LLPD = 0; LLPD < LLPD_Idx.size(); LLPD++)
                    { 
                        int indx = LLPD_Idx.at(LLPD);
                        double DR_LLP_vertex = DeltaR(genEta, EtaFromVertex(generator_->partVx[indx], generator_->partVy[indx], generator_->partVz[indx]), genPhi,  PhiFromVertex(generator_->partVx[indx], generator_->partVy[indx]));
                        if (DR_LLP_vertex < min_DR_LLP_vertex) min_DR_LLP_vertex = DR_LLP_vertex;
                    }
                    hLLP_vertex_DR->Fill(min_DR_LLP_vertex);
                }
            }


            if (nGoodGen > 0 ) hNMatchedLLP_inHCAL_DR02->Fill(nMatchedLLP_DR02);
            if (nGoodGen > 0 ) hNMatchedLLP_inHCAL_DR05->Fill(nMatchedLLP_DR05);
            hNMatchedLLP_DR02->Fill(nMatchedLLP_DR02);
            hNMatchedLLP_DR05->Fill(nMatchedLLP_DR05);
            if (nGoodGen > 0 ) hNMatched_LLPdaught_DR02->Fill(nMatchedGoodGen_DR02);
            if (nGoodGen > 0 ) hNMatched_LLPdaught_DR05->Fill(nMatchedGoodGen_DR05);
            if (nGoodGen > 0) hfracMatched_LLPdaught_DR02->Fill(nMatchedGoodGen_DR02 / nGoodGen);
            if (nGoodGen > 0) hfracMatched_LLPdaught_DR05->Fill(nMatchedGoodGen_DR05 / nGoodGen);
            hNLLPdaughts_inHCAL->Fill(nGoodGen);
            hNLLPdaughts->Fill(nLLPdaught);
            hNLLPdaughts_pteta->Fill(nLLPdaught_pteta);

        
            int nDepth = 0;
            double tpiEtaemu = 0, tpEtemu = 0, scaledEDepth = 0;
            std::vector<double> depthTPIt;
            std::vector<bool> usedTP(nHCALTP, false);

            //Fill energy depth histograms for matching gen particles to TPs

            for(int TPIt = 0; TPIt < nHCALTP; TPIt++)
            {
                nDepth = l1CaloTPemu_->hcalTPnDepths[TPIt];
                depthTPIt = hcalTPdepth[TPIt];
                tpEtemu = l1CaloTPemu_->hcalTPet[TPIt];
                tpiEtaemu = l1CaloTPemu_->hcalTPieta[TPIt];
                if (nDepth == 0) continue;
                for(int depthIt = 0; depthIt < nDepth; depthIt++)
                {
                    scaledEDepth = depthTPIt[depthIt]/tpEtemu;
                    if (abs(tpiEtaemu) <= 16)
                    {
                        if (matchedDirectTP.at(TPIt)) energyDepth_genMatchTP_Barrel->Fill(depthIt+1, scaledEDepth);
                    } 
                    else if (abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29)
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
            std::vector<bool> GoodJets_pt20(nJetemu, false);
            std::map< std::string, std::vector<bool> > LLPJetTags;
            LLPJetTags["HB02_HE02"] = GoodJets_pt20;
            LLPJetTags["HB04_HE04"] = GoodJets_pt20;
            LLPJetTags["HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["HB08_HE08"] = GoodJets_pt20;
            LLPJetTags["TP2_HB02_HE02"] = GoodJets_pt20;
            LLPJetTags["TP2_HB04_HE04"] = GoodJets_pt20;
            LLPJetTags["TP2_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP2_HB08_HE08"] = GoodJets_pt20;
            LLPJetTags["TP05_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP1_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP2_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP3_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP4_HB06_HE06"] = GoodJets_pt20;
            LLPJetTags["TP5_HB06_HE06"] = GoodJets_pt20;
            std::vector< std::vector<double>> jetEnergyPerDepth;
            std::vector< std::vector<double>> jetEnergyPerDepth_Gen;

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
                if ( (seedTowerHad / seedTower5x5Had) > 0.1 && abs(l1emu_->jetEta[jetIt]) < 3 && l1emu_->jetEt[jetIt] > 20 ) GoodJets_pt20.at(jetIt) = true;                
                if ( (seedTowerHad / seedTower5x5Had) > 0.1 && abs(l1emu_->jetEta[jetIt]) < 3 && seedTower3x3Had/ (seedTower3x3Had + seedTower3x3Em) > 0.9 /*&& l1emu_->jetEt[jetIt] > 20*/ )// && seedTower3x3Had > 30)  //requirement for throwing out junk jets
                {

                    double jetEta = l1emu_->jetEta[jetIt];
                    double jetPhi = l1emu_->jetPhi[jetIt];
                    for(int TPIt = 0; TPIt < nHCALTP; TPIt++)
                    {
                        double TPeta = etaVal(l1CaloTPemu_->hcalTPieta[TPIt]);
                        double TPphi = phiVal(l1CaloTPemu_->hcalTPiphi[TPIt]);
                        double deltaRIt = DeltaR(jetPhi, TPphi, jetEta, TPeta);
                        if (deltaRIt < 0.5 && !usedTP.at(TPIt) && l1CaloTPemu_->hcalTPet[TPIt] > 0.5 )
                        {
//                            usedTP.at(TPIt) = true;
//                            nDepth = l1CaloTPemu_->hcalTPnDepths[TPIt];
                            depthTPIt = hcalTPdepth[TPIt];
                            tpEtemu = l1CaloTPemu_->hcalTPet[TPIt];
                            tpiEtaemu = l1CaloTPemu_->hcalTPieta[TPIt];

                            if (nPassedJets < 4) hcalTP_nearL1Jet_emu->Fill(tpEtemu);
                            if (nPassedJets < 4 && abs(tpiEtaemu) <= 16) hcalTP_nearL1Jet_Barrel_emu->Fill(tpEtemu);
                            if (nPassedJets < 4 && abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29) hcalTP_nearL1Jet_Endcap_emu->Fill(tpEtemu);
                            if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4 && abs(tpiEtaemu) <= 16) hcalTP_nearL1Jet_Gen_Barrel_emu->Fill(tpEtemu);
                            if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4 && abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29) hcalTP_nearL1Jet_Gen_Endcap_emu->Fill(tpEtemu);

                            if( abs(tpiEtaemu) <= 16)
                            {
                                if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.2) LLPJetTags["HB02_HE02"].at(jetIt) = true;
                                if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.4) LLPJetTags["HB04_HE04"].at(jetIt) = true;
                                if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.6) LLPJetTags["HB06_HE06"].at(jetIt) = true;
                                if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.8) LLPJetTags["HB08_HE08"].at(jetIt) = true;
                                double HB_TPE_thresh = 2;
                                if (tpEtemu > HB_TPE_thresh && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.2) LLPJetTags["TP2_HB02_HE02"].at(jetIt) = true;
                                if (tpEtemu > HB_TPE_thresh && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.4) LLPJetTags["TP2_HB04_HE04"].at(jetIt) = true;
                                if (tpEtemu > HB_TPE_thresh && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.6) LLPJetTags["TP2_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > HB_TPE_thresh && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > 0.8) LLPJetTags["TP2_HB08_HE08"].at(jetIt) = true;
                                double HB_depthratio_thresh = 0.6;
                                if ((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP05_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 1 && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP1_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 2 && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP2_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 3 && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP3_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 4 && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP4_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 5 && (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu > HB_depthratio_thresh) LLPJetTags["TP5_HB06_HE06"].at(jetIt) = true;

                                if (nPassedJets < 4)
                                {
                                    energyDepth_DepthEnergy1_HB->Fill(depthTPIt.at(0));
                                    energyDepth_DepthEnergy2_HB->Fill(depthTPIt.at(1));
                                    energyDepth_DepthEnergy3_HB->Fill(depthTPIt.at(2));
                                    energyDepth_DepthEnergy4_HB->Fill(depthTPIt.at(3));
                                    energyDepth_Ratio_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);

                                    energyDepth_Ratio_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                    energyDepth_Ratio_IEta_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));

                                    if (htSum > 120)
                                    {
                                        energyDepth_HT120_DepthEnergy1_HB->Fill(depthTPIt.at(0));
                                        energyDepth_HT120_DepthEnergy2_HB->Fill(depthTPIt.at(1));
                                        energyDepth_HT120_DepthEnergy3_HB->Fill(depthTPIt.at(2));
                                        energyDepth_HT120_DepthEnergy4_HB->Fill(depthTPIt.at(3));

                                        energyDepth_Ratio_HT120_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);
                                        energyDepth_Ratio_HT120_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        if (tpEtemu > 1) energyDepth_Ratio_HT120_TP1_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        if (tpEtemu > 2) energyDepth_Ratio_HT120_TP2_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        if (tpEtemu > 3) energyDepth_Ratio_HT120_TP3_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        if (tpEtemu > 4) energyDepth_Ratio_HT120_TP4_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        if (tpEtemu > 5) energyDepth_Ratio_HT120_TP5_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);

                                        energyDepth_Ratio_IEta_HT120_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    if (tpEtemu >= 10) 
                                    {
                                        energyDepth_Ratio_TPge10_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);
                                        energyDepth_Ratio_TPge10_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);

                                        energyDepth_Ratio_IEta_TPge10_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    
                                }
                                if (nPassedGenMatchedJets < 4 && matchedJet.at(jetIt))
                                {
                                    energyDepth_Gen_DepthEnergy1_HB->Fill(depthTPIt.at(0));
                                    energyDepth_Gen_DepthEnergy2_HB->Fill(depthTPIt.at(1));
                                    energyDepth_Gen_DepthEnergy3_HB->Fill(depthTPIt.at(2));
                                    energyDepth_Gen_DepthEnergy4_HB->Fill(depthTPIt.at(3));
                                    energyDepth_Ratio_Gen_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);
                                    energyDepth_Ratio_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                    energyDepth_Ratio_IEta_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));
                                    if (htSum > 120)
                                    {
                                        energyDepth_HT120_Gen_DepthEnergy1_HB->Fill(depthTPIt.at(0));
                                        energyDepth_HT120_Gen_DepthEnergy2_HB->Fill(depthTPIt.at(1));
                                        energyDepth_HT120_Gen_DepthEnergy3_HB->Fill(depthTPIt.at(2));
                                        energyDepth_HT120_Gen_DepthEnergy4_HB->Fill(depthTPIt.at(3));

                                        energyDepth_Ratio_HT120_Gen_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);
                                        energyDepth_Ratio_HT120_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);

                                        energyDepth_Ratio_IEta_HT120_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    if (tpEtemu >= 10) 
                                    {
                                        energyDepth_Ratio_TPge10_Gen_Depth4_HB->Fill((depthTPIt.at(3)) / tpEtemu);
                                        energyDepth_Ratio_TPge10_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu);
                                        energyDepth_Ratio_IEta_TPge10_Gen_Depth3_4_HB->Fill((depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    
                                }

                            }
                            else if ( abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29)
                            {
                                if ( (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.2) LLPJetTags["HB02_HE02"].at(jetIt) = true;
                                if ( (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.4) LLPJetTags["HB04_HE04"].at(jetIt) = true;
                                if ( (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.6) LLPJetTags["HB06_HE06"].at(jetIt) = true;
                                if ( (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.8) LLPJetTags["HB08_HE08"].at(jetIt) = true;
                                double HE_TPE_thresh = 2;
                                if (tpEtemu > HE_TPE_thresh && (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.2) LLPJetTags["TP2_HB02_HE02"].at(jetIt) = true;
                                if (tpEtemu > HE_TPE_thresh && (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.4) LLPJetTags["TP2_HB04_HE04"].at(jetIt) = true;
                                if (tpEtemu > HE_TPE_thresh && (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.6) LLPJetTags["TP2_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > HE_TPE_thresh && (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > 0.8) LLPJetTags["TP2_HB08_HE08"].at(jetIt) = true;
                                double HE_depthratio_thresh = 0.6;
                                if ((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP05_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 1 && (depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP1_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 2 && (depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP2_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 3 && (depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP3_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 4 && (depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP4_HB06_HE06"].at(jetIt) = true;
                                if (tpEtemu > 5 && (depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu > HE_depthratio_thresh) LLPJetTags["TP5_HB06_HE06"].at(jetIt) = true;
                
                                if (nPassedJets < 4)
                                {
                                    energyDepth_DepthEnergy1_HE->Fill(depthTPIt.at(0));
                                    energyDepth_DepthEnergy2_HE->Fill(depthTPIt.at(1));
                                    energyDepth_DepthEnergy3_HE->Fill(depthTPIt.at(2));
                                    energyDepth_DepthEnergy4_HE->Fill(depthTPIt.at(3));
                                    energyDepth_DepthEnergy5_HE->Fill(depthTPIt.at(4));
                                    energyDepth_DepthEnergy6_HE->Fill(depthTPIt.at(5));
                                    energyDepth_DepthEnergy7_HE->Fill(depthTPIt.at(6));

                                    energyDepth_Ratio_Depth4_7_HE->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                    energyDepth_Ratio_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);

                                    energyDepth_Ratio_IEta_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    if (htSum > 120)
                                    {
                                        energyDepth_HT120_DepthEnergy1_HE->Fill(depthTPIt.at(0));
                                        energyDepth_HT120_DepthEnergy2_HE->Fill(depthTPIt.at(1));
                                        energyDepth_HT120_DepthEnergy3_HE->Fill(depthTPIt.at(2));
                                        energyDepth_HT120_DepthEnergy4_HE->Fill(depthTPIt.at(3));
                                        energyDepth_HT120_DepthEnergy5_HE->Fill(depthTPIt.at(4));
                                        energyDepth_HT120_DepthEnergy6_HE->Fill(depthTPIt.at(5));
                                        energyDepth_HT120_DepthEnergy7_HE->Fill(depthTPIt.at(6));

                                        energyDepth_Ratio_HT120_Depth4_7_HE->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_HT120_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        if (tpEtemu > 1) energyDepth_Ratio_HT120_TP1_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        if (tpEtemu > 2) energyDepth_Ratio_HT120_TP2_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        if (tpEtemu > 3) energyDepth_Ratio_HT120_TP3_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        if (tpEtemu > 4) energyDepth_Ratio_HT120_TP4_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        if (tpEtemu > 5) energyDepth_Ratio_HT120_TP5_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_IEta_HT120_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    }   
                                    if (tpEtemu >= 10) 
                                    {
                                        energyDepth_Ratio_TPge10_Depth4_7_HE->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_TPge10_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_IEta_TPge10_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    
                                }
                                if (nPassedGenMatchedJets < 4 && matchedJet.at(jetIt))
                                {
                                    energyDepth_Gen_DepthEnergy1_HE->Fill(depthTPIt.at(0));
                                    energyDepth_Gen_DepthEnergy2_HE->Fill(depthTPIt.at(1));
                                    energyDepth_Gen_DepthEnergy3_HE->Fill(depthTPIt.at(2));
                                    energyDepth_Gen_DepthEnergy4_HE->Fill(depthTPIt.at(3));
                                    energyDepth_Gen_DepthEnergy5_HE->Fill(depthTPIt.at(4));
                                    energyDepth_Gen_DepthEnergy6_HE->Fill(depthTPIt.at(5));
                                    energyDepth_Gen_DepthEnergy7_HE->Fill(depthTPIt.at(6));
                                
                                    energyDepth_Ratio_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                    energyDepth_Ratio_IEta_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    if (htSum > 120)
                                    {
                                        energyDepth_HT120_Gen_DepthEnergy1_HE->Fill(depthTPIt.at(0));
                                        energyDepth_HT120_Gen_DepthEnergy2_HE->Fill(depthTPIt.at(1));
                                        energyDepth_HT120_Gen_DepthEnergy3_HE->Fill(depthTPIt.at(2));
                                        energyDepth_HT120_Gen_DepthEnergy4_HE->Fill(depthTPIt.at(3));
                                        energyDepth_HT120_Gen_DepthEnergy5_HE->Fill(depthTPIt.at(4));
                                        energyDepth_HT120_Gen_DepthEnergy6_HE->Fill(depthTPIt.at(5));
                                        energyDepth_HT120_Gen_DepthEnergy7_HE->Fill(depthTPIt.at(6));
                                        
                                        energyDepth_Ratio_HT120_Gen_Depth4_7_HE->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_HT120_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_IEta_HT120_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    if (tpEtemu >= 10) 
                                    {
                                        energyDepth_Ratio_TPge10_Gen_Depth4_7_HE->Fill((depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_TPge10_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu);
                                        energyDepth_Ratio_IEta_TPge10_Gen_Depth2_7_HE->Fill((depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu, abs(tpiEtaemu));
                                    }
                                    
                                }
                            }   
                            for(int depthIt = 0; depthIt < 7; depthIt++)
                            {                                
                                scaledEDepth = depthTPIt[depthIt]/tpEtemu;
                                bool is_low_ratio_HB = (depthTPIt.at(2) + depthTPIt.at(3)) / tpEtemu < 0.00001;
                                bool is_low_ratio_HE = (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu < 0.00001;
                                if (abs(tpiEtaemu) <= 16)
                                {
                                    if (nPassedJets < 4)
                                    {
                                        energyDepth_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0) energyDepth_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0 && tpEtemu > 1) energyDepth_TPge5_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0 && htSum > 120) energyDepth_HT120_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (is_low_ratio_HB) energyDepth_LowRatio_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0) TPEt_HB->Fill(tpEtemu);
                                        if (is_low_ratio_HB) TPEt_LowRatio_HB->Fill(tpEtemu);
                                        if (tpEtemu > 5 && is_low_ratio_HB) energyDepth_TPge5_LowRatio_Barrel->Fill(depthIt+1, scaledEDepth);
                                        energyDepth_TPE_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975)  energyDepth_TPE_HoEcut_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPE_TPge5_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPge5_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (htSum > 120) energyDepth_HT120_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (tpEtemu > 5 && seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_TPge5_HoEcut_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_HoEcut_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (jetIt == 0) energyDepth_Jet1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 1) energyDepth_Jet2_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 2) energyDepth_Jet3_Barrel->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 3) energyDepth_Jet4_Barrel->Fill(depthIt+1, scaledEDepth);

                                        HovE_3x3_ET_AllDepthBarrel_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                        
                                    }
                                    if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4)
                                    {
                                        energyDepth_genMatchInclusive_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0) energyDepth_genMatchInclusive_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0 && tpEtemu > 1) energyDepth_genMatchInclusive_TPge5_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0 && htSum > 120) energyDepth_genMatchInclusive_HT120_L1_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (is_low_ratio_HB) energyDepth_genMatchInclusive_LowRatio_Barrel->Fill(depthIt+1, scaledEDepth);
                                        TPEt_Matched_HB->Fill(tpEtemu);
                                        if (is_low_ratio_HB) TPEt_LowRatio_Matched_HB->Fill(tpEtemu);
                                        TPEt_Matched_HB->Fill(tpEtemu);
                                        if (tpEtemu > 5 && is_low_ratio_HB) energyDepth_genMatchInclusive_TPge5_LowRatio_Barrel->Fill(depthIt+1, scaledEDepth);
                                        energyDepth_TPE_genMatchInclusive_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975)  energyDepth_TPE_genMatchInclusive_HoEcut_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPE_genMatchInclusive_TPge5_Barrel->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_genMatchInclusive_TPge5_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (htSum > 120) energyDepth_genMatchInclusive_HT120_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (tpEtemu > 5 && seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel->Fill(depthIt+1, scaledEDepth);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_genMatchInclusive_HoEcut_Barrel->Fill(depthIt+1, scaledEDepth);
                                        HovE_3x3_ET_AllDepth_GenBarrel_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                } 
                                else if (abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29)
                                {
                                    if (nPassedJets < 4)
                                    {
                                        energyDepth_Endcap->Fill(depthIt+1, scaledEDepth);
                                        energyDepth_TPE_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (nPassedJets == 0) energyDepth_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0 && tpEtemu > 1) energyDepth_TPge5_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedJets == 0 && htSum > 120) energyDepth_HT120_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (is_low_ratio_HE) energyDepth_LowRatio_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (tpEtemu > 5 && is_low_ratio_HE) energyDepth_TPge5_LowRatio_Endcap->Fill(depthIt+1, scaledEDepth);
                                        TPEt_HE->Fill(tpEtemu);
                                        if (is_low_ratio_HE) TPEt_LowRatio_HE->Fill(tpEtemu);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975)  energyDepth_TPE_HoEcut_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPE_TPge5_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPge5_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (htSum > 120) energyDepth_HT120_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (tpEtemu > 5 && seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_TPge5_HoEcut_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_HoEcut_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (jetIt == 0) energyDepth_Jet1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 1) energyDepth_Jet2_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 2) energyDepth_Jet3_Endcap->Fill(depthIt+1, scaledEDepth);
                                        else if (jetIt == 3) energyDepth_Jet4_Endcap->Fill(depthIt+1, scaledEDepth);
                                        
                                        HovE_3x3_ET_AllDepthEndcap_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                    if (matchedJet.at(jetIt) && nPassedGenMatchedJets < 4) 
                                    {
                                        energyDepth_genMatchInclusive_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (tpEtemu > 5) energyDepth_genMatchInclusive_TPge5_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (htSum > 120) energyDepth_genMatchInclusive_HT120_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0) energyDepth_genMatchInclusive_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0 && tpEtemu > 1) energyDepth_genMatchInclusive_TPge5_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (nPassedGenMatchedJets == 0 && htSum > 120) energyDepth_genMatchInclusive_HT120_L1_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (is_low_ratio_HE) energyDepth_genMatchInclusive_LowRatio_Endcap->Fill(depthIt+1, scaledEDepth);
                                        TPEt_Matched_HE->Fill(tpEtemu);
                                        if (is_low_ratio_HE) TPEt_LowRatio_Matched_HE->Fill(tpEtemu);
                                        if (tpEtemu > 5 && is_low_ratio_HE) energyDepth_genMatchInclusive_TPge5_LowRatio_Endcap->Fill(depthIt+1, scaledEDepth);
                                        energyDepth_TPE_genMatchInclusive_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975)  energyDepth_TPE_genMatchInclusive_HoEcut_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);
                                        if (tpEtemu > 5) energyDepth_TPE_genMatchInclusive_TPge5_Endcap->Fill(depthIt+1, depthTPIt[depthIt]);

                                        if (tpEtemu > 5 && seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap->Fill(depthIt+1, scaledEDepth);
                                        if (seedTower3x3Had/(seedTower3x3Had / seedTower3x3Em) > 0.975) energyDepth_genMatchInclusive_HoEcut_Endcap->Fill(depthIt+1, scaledEDepth);
                                        HovE_3x3_ET_AllDepth_GenEndcap_hists.at(depthIt)->Fill(depthTPIt[depthIt], seedTower3x3Had/(seedTower3x3Had + seedTower3x3Em));
                                    }
                                }
                            }//close depth loop
                        }//close TP near L1 jet requirement
                    }//close TP loop
                    nPassedJets++;
                    if (matchedJet.at(jetIt)) nPassedGenMatchedJets++;
                }// end jet requirement 
	
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

            }

            std::vector<TH1F*> jetET_hists{hJetETLeading1, hJetETLeading2, hJetETLeading3, hJetETLeading4}; 
            std::vector<TH1F*> jetEta_hists{hJetEtaLeading1, hJetEtaLeading2, hJetEtaLeading3, hJetEtaLeading4}; 
            std::vector<TH2F*> HEEnergy_1x1_hists{HEEnergytotal_1x1_emu_Leading1, HEEnergytotal_1x1_emu_Leading2, HEEnergytotal_1x1_emu_Leading3, HEEnergytotal_1x1_emu_Leading4};      
            std::vector<TH2F*> HEEnergy_3x3_hists{HEEnergytotal_3x3_emu_Leading1, HEEnergytotal_3x3_emu_Leading2, HEEnergytotal_3x3_emu_Leading3, HEEnergytotal_3x3_emu_Leading4};      
            std::vector<TH1F*> HovE_1x1_hists{HovEtotal_1x1_emu_Leading1, HovEtotal_1x1_emu_Leading2, HovEtotal_1x1_emu_Leading3, HovEtotal_1x1_emu_Leading4};      
            std::vector<TH1F*> HovE_3x3_hists{HovEtotal_3x3_emu_Leading1, HovEtotal_3x3_emu_Leading2, HovEtotal_3x3_emu_Leading3, HovEtotal_3x3_emu_Leading4}; 
            std::vector<TH2F*> HovE_1x1_ET_hists{HovEtotal_1x1_ET_emu_Leading1, HovEtotal_1x1_ET_emu_Leading2, HovEtotal_1x1_ET_emu_Leading3, HovEtotal_1x1_ET_emu_Leading4};      
            std::vector<TH2F*> HovE_3x3_ET_hists{HovEtotal_3x3_ET_emu_Leading1, HovEtotal_3x3_ET_emu_Leading2, HovEtotal_3x3_ET_emu_Leading3, HovEtotal_3x3_ET_emu_Leading4};
            for (int pjet = 0; pjet < nJetemu && pjet < 4; pjet++)
            {

                if (!GoodJets_pt20.at(pjet)) continue;
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
                if (hadVariablesAllJets["H3OvE3"].at(pjet) > 50) HovEtotal_3x3_HETge50_emu->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                double logHovE = emVariablesAllJets["HOvE"].at(pjet) > 0 ? log10(hadVariablesAllJets["HOvE"].at(pjet)/emVariablesAllJets["HOvE"].at(pjet)) : 10;
                double logHovE_3x3 = emVariablesAllJets["H3OvE3"].at(pjet) > 0 ? log10(hadVariablesAllJets["H3OvE3"].at(pjet)/emVariablesAllJets["H3OvE3"].at(pjet)) : 10;	 
                HovEtotalLog_1x1_emu->Fill(logHovE);
                HovEtotalLog_3x3_emu->Fill(logHovE_3x3);
                if (abs(jetVariablesAllJets["ieta"].at(pjet)) <= 16) 
                {
                    HovEtotal_1x1_emu_Barrel->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                    HovEtotal_3x3_emu_Barrel->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                }
                else if (abs(jetVariablesAllJets["ieta"].at(pjet)) > 16 && abs(jetVariablesAllJets["ieta"].at(pjet)) <= 29) 
                {
                    HovEtotal_1x1_emu_Endcap->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                    HovEtotal_3x3_emu_Endcap->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                }
                if (matchedJet.at(pjet)) 
                {
                    HovEtotal_1x1_emu_GenMatchedJets->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                    HovEtotal_3x3_emu_GenMatchedJets->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                    HEEnergytotal_1x1_emu_GenMatched->Fill( emVariablesAllJets["HOvE"].at(pjet),  hadVariablesAllJets["HOvE"].at(pjet));
                    HEEnergytotal_3x3_emu_GenMatched->Fill( emVariablesAllJets["H3OvE3"].at(pjet),  hadVariablesAllJets["H3OvE3"].at(pjet));
                    if (hadVariablesAllJets["H3OvE3"].at(pjet) > 50) HovEtotal_3x3_emu_HETge50_GenMatchedJets->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                    if (abs(jetVariablesAllJets["ieta"].at(pjet)) <= 16) 
                    {
                        HovEtotal_1x1_emu_GenMatchedJets_Barrel->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                        HovEtotal_3x3_emu_GenMatchedJets_Barrel->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                    }
                    else if (abs(jetVariablesAllJets["ieta"].at(pjet)) > 16 && abs(jetVariablesAllJets["ieta"].at(pjet)) <= 29) 
                    {
                        HovEtotal_1x1_emu_GenMatchedJets_Endcap->Fill(hadVariablesAllJets["HOvE"].at(pjet) / (hadVariablesAllJets["HOvE"].at(pjet) + emVariablesAllJets["HOvE"].at(pjet)));
                        HovEtotal_3x3_emu_GenMatchedJets_Endcap->Fill(hadVariablesAllJets["H3OvE3"].at(pjet) / (hadVariablesAllJets["H3OvE3"].at(pjet) + emVariablesAllJets["H3OvE3"].at(pjet)));
                    }
                }
            }
            
            std::vector<bool> pass_HoE(4, false);
            std::vector<bool> pass_HoE_DFB_02(4, false), pass_HoE_DFB_04(4, false), pass_HoE_DFB_06(4, false), pass_HoE_DFB_08(4, false);
            std::vector<bool> pass_HoE_TP2_DFB_02(4, false), pass_HoE_TP2_DFB_04(4, false), pass_HoE_TP2_DFB_06(4, false), pass_HoE_TP2_DFB_08(4, false);
            std::vector<bool> pass_HoE_DFBTP05(4, false), pass_HoE_DFBTP1(4, false), pass_HoE_DFBTP2(4, false), pass_HoE_DFBTP3(4, false), pass_HoE_DFBTP4(4, false), pass_HoE_DFBTP5(4, false);
            for (int ijet = 0; ijet < nJetemu && ijet < 4; ijet++)
            {
              if ((hadVariablesAllJets["H3OvE3"].at(ijet))/(hadVariablesAllJets["H3OvE3"].at(ijet)+emVariablesAllJets["H3OvE3"].at(ijet)) > -1) pass_HoE.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP05_HB06_HE06"].at(ijet)) pass_HoE_DFBTP05.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP1_HB06_HE06"].at(ijet)) pass_HoE_DFBTP1.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP2_HB06_HE06"].at(ijet)) pass_HoE_DFBTP2.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP3_HB06_HE06"].at(ijet)) pass_HoE_DFBTP3.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP4_HB06_HE06"].at(ijet)) pass_HoE_DFBTP4.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP5_HB06_HE06"].at(ijet)) pass_HoE_DFBTP5.at(ijet) = true;

              if (pass_HoE.at(ijet) && LLPJetTags["HB02_HE02"].at(ijet)) pass_HoE_DFB_02.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["HB04_HE04"].at(ijet)) pass_HoE_DFB_04.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["HB06_HE06"].at(ijet)) pass_HoE_DFB_06.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["HB08_HE08"].at(ijet)) pass_HoE_DFB_08.at(ijet) = true;

              if (pass_HoE.at(ijet) && LLPJetTags["TP2_HB02_HE02"].at(ijet)) pass_HoE_TP2_DFB_02.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP2_HB04_HE04"].at(ijet)) pass_HoE_TP2_DFB_04.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP2_HB06_HE06"].at(ijet)) pass_HoE_TP2_DFB_06.at(ijet) = true;
              if (pass_HoE.at(ijet) && LLPJetTags["TP2_HB08_HE08"].at(ijet)) pass_HoE_TP2_DFB_08.at(ijet) = true;
            }

            hHTSum_emu->Fill(htSum);
            effJetID_HoE_DepthFB_TPE->Fill(0); 
            if (htSum > 360) effJetID_HoE_DepthFB_TPE->Fill(1); 
            if (htSum > 120) effJetID_HoE_DepthFB_TPE->Fill(2); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(3); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(4); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(5); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(6); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(7); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(8);
            if ((htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(9); 
            if ((htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(10); 
            if ((htSum > 120 && (pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(11); 
            if ((htSum > 120 && (pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(12); 
            if ((htSum > 120 && (pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(13); 
            if ((htSum > 120 && (pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) ) ) effJetID_HoE_DepthFB_TPE->Fill(14); 

            effJetID_HoE_DepthFB_HTscan->Fill(0);
            if (htSum > 120) effJetID_HoE_DepthFB_HTscan->Fill(1);
            if (htSum > 180) effJetID_HoE_DepthFB_HTscan->Fill(2);
            if (htSum > 240) effJetID_HoE_DepthFB_HTscan->Fill(3);
            if (htSum > 300) effJetID_HoE_DepthFB_HTscan->Fill(4);            
            if (htSum > 360) effJetID_HoE_DepthFB_HTscan->Fill(5);
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)))) effJetID_HoE_DepthFB_HTscan->Fill(6);
            if (htSum > 360 || (htSum > 180 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(7);
            if (htSum > 360 || (htSum > 240 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(8);
            if (htSum > 360 || (htSum > 300 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(9);
            if (htSum > 360 || (htSum > 360 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(10);
            if ((htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)))) effJetID_HoE_DepthFB_HTscan->Fill(11);
            if ((htSum > 180 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(12);
            if ((htSum > 240 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(13);
            if ((htSum > 300 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(14);
            if ((htSum > 360 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(15);

            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)))) effJetID_HoE_DepthFB_HTscan->Fill(16);
            if (htSum > 360 || (htSum > 180 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(17);
            if (htSum > 360 || (htSum > 240 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(18);
            if (htSum > 360 || (htSum > 300 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(19);
            if (htSum > 360 || (htSum > 360 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(20);
            if ((htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)))) effJetID_HoE_DepthFB_HTscan->Fill(21);
            if ((htSum > 180 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(22);
            if ((htSum > 240 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(23);
            if ((htSum > 300 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(24);
            if ((htSum > 360 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_HTscan->Fill(25);

//            if (htSum > 120 && (pass_HoE_DFBTP05.at(0)4|| pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) effJetID_HoE_DepthFB->Fill(4); 
            effJetID_HoE_DepthFB_Ratio->Fill(0);
            if (htSum > 360) effJetID_HoE_DepthFB_Ratio->Fill(1); 
            if (htSum > 120) effJetID_HoE_DepthFB_Ratio->Fill(2); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(3); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(4); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(5); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(6); 
            if ((htSum > 120 && (pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(7); 
            if ((htSum > 120 && (pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(8); 
            if ((htSum > 120 && (pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(9); 
            if ((htSum > 120 && (pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(10); 

            if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(11); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(12); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(13); 
            if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(14); 
            if ((htSum > 120 && (pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(15); 
            if ((htSum > 120 && (pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(16); 
            if ((htSum > 120 && (pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(17); 
            if ((htSum > 120 && (pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio->Fill(18); 

            effJetID_HoE_DepthFB_Jetscan->Fill(0);
            if (jetEt_3 > 20) effJetID_HoE_DepthFB_Jetscan->Fill(1);
            if (jetEt_3 > 30) effJetID_HoE_DepthFB_Jetscan->Fill(2);
            if (jetEt_3 > 40) effJetID_HoE_DepthFB_Jetscan->Fill(3);
            if (jetEt_3 > 50) effJetID_HoE_DepthFB_Jetscan->Fill(4);
            if (jetEt_3 > 60) effJetID_HoE_DepthFB_Jetscan->Fill(5);
            if (jetEt_3 > 20 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(6);
            if (jetEt_3 > 30 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(7);
            if (jetEt_3 > 40 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(8);
            if (jetEt_3 > 50 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(9);
            if (jetEt_3 > 60 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(10);
            if (jetEt_3 > 20 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(11);
            if (jetEt_3 > 30 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(12);
            if (jetEt_3 > 40 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(13);
            if (jetEt_3 > 50 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(14);
            if (jetEt_3 > 60 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Jetscan->Fill(15);
            if (jetEt_4 > 20) effJetID_HoE_DepthFB_Jetscan->Fill(16);
            if (jetEt_4 > 30) effJetID_HoE_DepthFB_Jetscan->Fill(17);
            if (jetEt_4 > 40) effJetID_HoE_DepthFB_Jetscan->Fill(18);
            if (jetEt_4 > 50) effJetID_HoE_DepthFB_Jetscan->Fill(19);
            if (jetEt_4 > 60) effJetID_HoE_DepthFB_Jetscan->Fill(20);
            if (jetEt_4 > 20 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(21);
            if (jetEt_4 > 30 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(22);
            if (jetEt_4 > 40 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(23);
            if (jetEt_4 > 50 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(24);
            if (jetEt_4 > 60 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(25);
            if (jetEt_4 > 20 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(26);
            if (jetEt_4 > 30 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(27);
            if (jetEt_4 > 40 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(28);
            if (jetEt_4 > 50 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(29);
            if (jetEt_4 > 60 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Jetscan->Fill(30);

            
            if (nGoodGen > 0) 
            {
                hHTSum_Gen_emu->Fill(htSum);
                effJetID_HoE_DepthFB_TPE_Gen->Fill(0);
                if (htSum > 360) effJetID_HoE_DepthFB_TPE_Gen->Fill(1); 
                if (htSum > 120) effJetID_HoE_DepthFB_TPE_Gen->Fill(2); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(3); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(4); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(5); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(6); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(7); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(8); 
                if ((htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(9); 
                if ((htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(10); 
                if ((htSum > 120 && (pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(11); 
                if ((htSum > 120 && (pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(12); 
                if ((htSum > 120 && (pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(13); 
                if ((htSum > 120 && (pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) ) ) effJetID_HoE_DepthFB_TPE_Gen->Fill(14); 


                effJetID_HoE_DepthFB_Gen_HTscan->Fill(0);
                if (htSum > 120) effJetID_HoE_DepthFB_Gen_HTscan->Fill(1);
                if (htSum > 180) effJetID_HoE_DepthFB_Gen_HTscan->Fill(2);
                if (htSum > 240) effJetID_HoE_DepthFB_Gen_HTscan->Fill(3);
                if (htSum > 300) effJetID_HoE_DepthFB_Gen_HTscan->Fill(4);            
                if (htSum > 360) effJetID_HoE_DepthFB_Gen_HTscan->Fill(5);
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)))) effJetID_HoE_DepthFB_Gen_HTscan->Fill(6);
                if (htSum > 360 || (htSum > 180 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(7);
                if (htSum > 360 || (htSum > 240 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(8);
                if (htSum > 360 || (htSum > 300 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(9);
                if (htSum > 360 || (htSum > 360 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(10);
                if ((htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)))) effJetID_HoE_DepthFB_Gen_HTscan->Fill(11);
                if ((htSum > 180 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(12);
                if ((htSum > 240 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(13);
                if ((htSum > 300 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(14);
                if ((htSum > 360 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(15);
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)))) effJetID_HoE_DepthFB_Gen_HTscan->Fill(16);
                if (htSum > 360 || (htSum > 180 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(17);
                if (htSum > 360 || (htSum > 240 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(18);
                if (htSum > 360 || (htSum > 300 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(19);
                if (htSum > 360 || (htSum > 360 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(20);
                if ((htSum > 120 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)))) effJetID_HoE_DepthFB_Gen_HTscan->Fill(21);
                if ((htSum > 180 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(22);
                if ((htSum > 240 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(23);
                if ((htSum > 300 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(24);
                if ((htSum > 360 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) )) effJetID_HoE_DepthFB_Gen_HTscan->Fill(25);


            //if ((htSum > 120 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) ) ) effJetID_HoE_DepthFB_Gen->Fill(4); 
                effJetID_HoE_DepthFB_Ratio_Gen->Fill(0);
                if (htSum > 360) effJetID_HoE_DepthFB_Ratio_Gen->Fill(1); 
                if (htSum > 120) effJetID_HoE_DepthFB_Ratio_Gen->Fill(2); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(3); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(4); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(5); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(6); 
                if ((htSum > 120 && (pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(7); 
                if ((htSum > 120 && (pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(8); 
                if ((htSum > 120 && (pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(9); 
                if ((htSum > 120 && (pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(10); 

                if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(11); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(12); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(13); 
                if (htSum > 360 || (htSum > 120 && (pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(14); 
                if ((htSum > 120 && (pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(15); 
                if ((htSum > 120 && (pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(16); 
                if ((htSum > 120 && (pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(17); 
                if ((htSum > 120 && (pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) ) ) effJetID_HoE_DepthFB_Ratio_Gen->Fill(18); 

                effJetID_HoE_DepthFB_Gen_Jetscan->Fill(0);
                if (jetEt_3 > 20) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(1);
                if (jetEt_3 > 30) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(2);
                if (jetEt_3 > 40) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(3);
                if (jetEt_3 > 50) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(4);
                if (jetEt_3 > 60) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(5);
                if (jetEt_3 > 20 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(6);
                if (jetEt_3 > 30 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(7);
                if (jetEt_3 > 40 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(8);
                if (jetEt_3 > 50 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(9);
                if (jetEt_3 > 60 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(10);
                if (jetEt_3 > 20 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(11);
                if (jetEt_3 > 30 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(12);
                if (jetEt_3 > 40 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(13);
                if (jetEt_3 > 50 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(14);
                if (jetEt_3 > 60 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(15);
                if (jetEt_4 > 20) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(16);
                if (jetEt_4 > 30) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(17);
                if (jetEt_4 > 40) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(18);
                if (jetEt_4 > 50) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(19);
                if (jetEt_4 > 60) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(20);
                if (jetEt_4 > 20 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(21);
                if (jetEt_4 > 30 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(22);
                if (jetEt_4 > 40 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(23);
                if (jetEt_4 > 50 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(24);
                if (jetEt_4 > 60 && (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(25);
                if (jetEt_4 > 20 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(26);
                if (jetEt_4 > 30 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(27);
                if (jetEt_4 > 40 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(28);
                if (jetEt_4 > 50 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(29);
                if (jetEt_4 > 60 && (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3))) effJetID_HoE_DepthFB_Gen_Jetscan->Fill(30);


            }

            // for each bin fill according to whether our object has a larger corresponding energy
            for(int bin=0; bin<nJetBins; bin++){
                if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( pass_HoE.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( pass_HoE_DFBTP05.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_TP05_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( pass_HoE_DFBTP5.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_TP5_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if(  pass_HoE_DFB_02.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if(  pass_HoE_DFB_06.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if(  pass_HoE_TP2_DFB_02.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_TP2_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if(  pass_HoE_TP2_DFB_06.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_TP2_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if(  pass_HoE_DFBTP1.at(0) && (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_HoE_TP1_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV


            } 

            for(int bin=0; bin<nJetBins; bin++){
                if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE.at(0) || pass_HoE.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_TP05_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_TP5_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_TP2_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_TP2_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( (pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1)) && (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_HoE_TP1_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV


            }  

            for(int bin=0; bin<nJetBins; bin++){
                if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_TP05_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_TP5_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_02.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_TP2_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_02.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_TP2_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) ) && (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_HoE_TP1_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV

            }  

            for(int bin=0; bin<nJetBins; bin++){
                if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2) || pass_HoE.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_TP05_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_TP5_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_TP2_Ratio02_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_TP2_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
                if( ( pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3) ) && (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_HoE_TP1_Ratio06_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV

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
                if( ((pass_HoE.at(0) || pass_HoE.at(1) || pass_HoE.at(2) || pass_HoE.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( htSum > 360 || ((pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP05_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP1_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( htSum > 360 || ((pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( htSum > 360 || ((pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP3_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( htSum > 360 || ((pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP4_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                     
                if( htSum > 360 || ((pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP5_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                

                if( ((pass_HoE_DFBTP05.at(0) || pass_HoE_DFBTP05.at(1) || pass_HoE_DFBTP05.at(2) || pass_HoE_DFBTP05.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP05_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_DFBTP1.at(0) || pass_HoE_DFBTP1.at(1) || pass_HoE_DFBTP1.at(2) || pass_HoE_DFBTP1.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP1_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( ((pass_HoE_DFBTP2.at(0) || pass_HoE_DFBTP2.at(1) || pass_HoE_DFBTP2.at(2) || pass_HoE_DFBTP2.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( ((pass_HoE_DFBTP3.at(0) || pass_HoE_DFBTP3.at(1) || pass_HoE_DFBTP3.at(2) || pass_HoE_DFBTP3.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP3_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                
                if( ((pass_HoE_DFBTP4.at(0) || pass_HoE_DFBTP4.at(1) || pass_HoE_DFBTP4.at(2) || pass_HoE_DFBTP4.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP4_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                     
                if( ((pass_HoE_DFBTP5.at(0) || pass_HoE_DFBTP5.at(1) || pass_HoE_DFBTP5.at(2) || pass_HoE_DFBTP5.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP5_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV                

                if( htSum > 360 || ((pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio02_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio04_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio06_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio08_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           

                if( ((pass_HoE_DFB_02.at(0) || pass_HoE_DFB_02.at(1) || pass_HoE_DFB_02.at(2) || pass_HoE_DFB_02.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio02_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_DFB_04.at(0) || pass_HoE_DFB_04.at(1) || pass_HoE_DFB_04.at(2) || pass_HoE_DFB_04.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio04_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_DFB_06.at(0) || pass_HoE_DFB_06.at(1) || pass_HoE_DFB_06.at(2) || pass_HoE_DFB_06.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio06_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_DFB_08.at(0) || pass_HoE_DFB_08.at(1) || pass_HoE_DFB_08.at(2) || pass_HoE_DFB_08.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_Ratio08_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           

                if( htSum > 360 || ((pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio02_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio04_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio06_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( htSum > 360 || ((pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio08_ORHT360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           

                if( ((pass_HoE_TP2_DFB_02.at(0) || pass_HoE_TP2_DFB_02.at(1) || pass_HoE_TP2_DFB_02.at(2) || pass_HoE_TP2_DFB_02.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio02_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_TP2_DFB_04.at(0) || pass_HoE_TP2_DFB_04.at(1) || pass_HoE_TP2_DFB_04.at(2) || pass_HoE_TP2_DFB_04.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio04_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_TP2_DFB_06.at(0) || pass_HoE_TP2_DFB_06.at(1) || pass_HoE_TP2_DFB_06.at(2) || pass_HoE_TP2_DFB_06.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio06_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           
                if( ((pass_HoE_TP2_DFB_08.at(0) || pass_HoE_TP2_DFB_08.at(1) || pass_HoE_TP2_DFB_08.at(2) || pass_HoE_TP2_DFB_08.at(3)) && (htSum) >= htSumLo+(bin*htSumBinWidth) )) htSumRates_HoE_TP2_Ratio08_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV           

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


//            if (jentry > 1000) break; //JUST FOR TEXT DUMP STUDY REMOVE AFTER!!
/*            
            //continuing gen study; printing jet info
            int jetcounter = 0;
            int textspace = 10;
            for (int jetIt = 0; jetIt < nJetemu; jetIt++)
            {
//                if (jentry > 500) break;
                double jetEta = l1emu_->jetEta[jetIt];
                double jetPhi = l1emu_->jetPhi[jetIt];
                double jetEt = l1emu_->jetEt[jetIt];
                if (abs(jetEta) < 3 && jetEt > 20 && jetcounter < 4)
                {
                    TPtextfile << std::left << "Jet ET: " << std::setw(textspace) << jetEt << "Eta: " << std::setw(textspace) << jetEta << "Phi: " << std::setw(textspace) << jetPhi << " H: " << std::setw(textspace) << hadVariablesAllJets["H3OvE3"].at(jetIt) << " E: " << std::setw(textspace) << emVariablesAllJets["H3OvE3"].at(jetIt) << std::endl;
                        
                    if (matchedJet.at(jetIt))
                    {
                        for(int genpart = 0; genpart < nGenPart; genpart++)
                        {
                            double genEta = generator_->partEta[genpart];
                            double genPhi = generator_->partPhi[genpart];
                            double genPt = generator_->partPt[genpart];
                            double Vz = generator_->partVz[genpart];
                            double Vx = generator_->partVx[genpart];
                            double Vy = generator_->partVy[genpart];
                            double genEta_corrected = genParticlesEta.at(genpart);
                            double genPhi_corrected = genParticlesPhi.at(genpart);
                            double DeltaRtoJet = DeltaR(jetPhi, genPhi, jetEta, genEta);
                            double DeltaRtoJet_corrected = DeltaR(jetPhi, genPhi_corrected, jetEta, genEta_corrected);
                            double vertexPhi = PhiFromVertex(Vx, Vy);
                            double vertexEta = EtaFromVertex(Vx, Vy, Vz);
                            double vertexR = sqrt(Vx*Vx + Vy*Vy);
                            double genp_vertex_DR = DeltaR(vertexPhi, genPhi, vertexEta, genEta);
                            double genp_vertex_DR_corrected = DeltaR(vertexPhi, genPhi_corrected, vertexEta, genEta_corrected);

                            if (generator_->partId[genpart] == 6000113)
                            {
                                TPtextfile << std::left  << std::setprecision(3) <<  "LLP PT: " << std::setw(textspace) << genPt  << "Eta: " << std::setw(textspace) <<  genEta  << "Phi " << std::setw(textspace) <<  genPhi  << "DR from jet: " << std::setw(textspace) <<  DeltaRtoJet << std::endl;  
                            }
                            if (generator_->partParent[genpart] == 6000113)
                            {
                                TPtextfile << std::left << std::setprecision(3)  << "LLPD PT: " << std::setw(textspace) << genPt << "Vertex Phi: " << std::setw(textspace) << vertexPhi  << "Vertex Eta: " << std::setw(textspace) << vertexEta << "Vertex R: " << std::setw(textspace) << vertexR  << "Vertex Z: " << std::setw(textspace) << Vz << "Eta: " << std::setw(textspace) <<  genEta  << "Phi " << std::setw(textspace) <<  genPhi  << "DR to jet: " << std::setw(textspace) <<  DeltaRtoJet << "DR to vertex: " << std::setw(textspace) << genp_vertex_DR <<  "Eta*: " << std::setw(textspace) <<  genEta_corrected  << "Phi*: " << std::setw(textspace) <<  genPhi_corrected  << "DR* to jet: " << std::setw(textspace) <<  DeltaRtoJet_corrected << "DR* to vertex: " << std::setw(textspace) << genp_vertex_DR_corrected << std::endl;
                            }
                        }
                    }
                    std::vector< std::pair<int, double>> TP_energy_sorted;
                    for(int TPIt = 0; TPIt < nHCALTP; TPIt++)
                    {
                        double TPeta = etaVal(l1CaloTPemu_->hcalTPieta[TPIt]);
                        double TPphi = phiVal(l1CaloTPemu_->hcalTPiphi[TPIt]);
                        double deltaRIt = DeltaR(jetPhi, TPphi, jetEta, TPeta);
                        if (deltaRIt < 0.5 ) TP_energy_sorted.push_back(std::make_pair(TPIt, l1CaloTPemu_->hcalTPet[TPIt]));
                    }
                    std::sort(TP_energy_sorted.begin(), TP_energy_sorted.end(), compareTP);
                    for (uint TPIt = 0; TPIt < TP_energy_sorted.size(); TPIt++)
                    {

                        double TPeta = etaVal(l1CaloTPemu_->hcalTPieta[TP_energy_sorted.at(TPIt).first]);
                        double TPphi = phiVal(l1CaloTPemu_->hcalTPiphi[TP_energy_sorted.at(TPIt).first]);
                        double tpEtemu = l1CaloTPemu_->hcalTPet[TP_energy_sorted.at(TPIt).first];
                        double tpiEtaemu = l1CaloTPemu_->hcalTPieta[TP_energy_sorted.at(TPIt).first];
                        double tpiPhiemu = l1CaloTPemu_->hcalTPCaliphi[TP_energy_sorted.at(TPIt).first];
                        double deltaRIt = DeltaR(jetPhi, TPphi, jetEta, TPeta);
                        int closest_jet = -1;
                        double min_DR = 100;
                        for (int jetIt2 = 0; jetIt2 < nJetemu; jetIt2++)
                        {
                            double deltaRIt2 = DeltaR(l1emu_->jetPhi[jetIt2], TPphi, l1emu_->jetEta[jetIt2], TPeta);
                            if (deltaRIt2 < min_DR ) 
                            {
                                closest_jet = jetIt2;
                                min_DR = deltaRIt2;
                            }
                        }
                        TPtextfile << std::left << std::setprecision(3) <<  "TP ET: " << std::setw(textspace) << tpEtemu <<  "TP IEta: " << std::setw(textspace) << tpiEtaemu << "TP IPhi: " << std::setw(textspace) << tpiPhiemu <<  "TP Eta: " << std::setw(textspace) << TPeta << "TP Phi: " << std::setw(textspace) << TPphi << " Delta R: " << std::setw(textspace) << deltaRIt << " Closest jet: " << std::setw(textspace) << closest_jet + 1;
                        for (int depthIt = 0; depthIt < 7; depthIt++)
                        {
                            TPtextfile << std::left << "D" << depthIt + 1 << ": " << std::setw(textspace) << hcalTPdepth.at(TP_energy_sorted.at(TPIt).first).at(depthIt);
                        }
                        std::vector<double> depthTPIt = hcalTPdepth.at(TP_energy_sorted.at(TPIt).first);
                        double ratio = -1;
                        if (abs(tpiEtaemu) <= 16) ratio = (depthTPIt.at(2) + depthTPIt.at(3))/tpEtemu;
                        else if (abs(tpiEtaemu) > 16 && abs(tpiEtaemu) <= 29) ratio = (depthTPIt.at(1) + depthTPIt.at(2) + depthTPIt.at(3) + depthTPIt.at(4) + depthTPIt.at(5) + depthTPIt.at(6)) / tpEtemu;
                        TPtextfile << "Ratio: " << std::setw(textspace) << ratio << std::endl;
                    }
                    jetcounter++;
                    TPtextfile << "-----------------------------------------------" << std::endl;
                }
                
            }
            TPtextfile << "******************************************************" << std::endl;
*/            


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
        singleJetRates_HoE_emu->Scale(norm);
        singleJetRates_HoE_TP05_emu->Scale(norm);
        singleJetRates_HoE_TP5_emu->Scale(norm);
        singleJetRates_HoE_Ratio02_emu->Scale(norm);
        singleJetRates_HoE_Ratio06_emu->Scale(norm);
        singleJetRates_HoE_TP2_Ratio02_emu->Scale(norm);
        singleJetRates_HoE_TP2_Ratio06_emu->Scale(norm);
        singleJetRates_HoE_TP1_Ratio06_emu->Scale(norm);
        doubleJetRates_emu->Scale(norm);
        doubleJetRates_HoE_emu->Scale(norm);
        doubleJetRates_HoE_TP05_emu->Scale(norm);
        doubleJetRates_HoE_TP5_emu->Scale(norm);
        doubleJetRates_HoE_Ratio02_emu->Scale(norm);
        doubleJetRates_HoE_Ratio06_emu->Scale(norm);
        doubleJetRates_HoE_TP2_Ratio02_emu->Scale(norm);
        doubleJetRates_HoE_TP2_Ratio06_emu->Scale(norm);
        doubleJetRates_HoE_TP1_Ratio06_emu->Scale(norm);
        tripleJetRates_emu->Scale(norm);
        tripleJetRates_HoE_emu->Scale(norm);
        tripleJetRates_HoE_TP05_emu->Scale(norm);
        tripleJetRates_HoE_TP5_emu->Scale(norm);
        tripleJetRates_HoE_Ratio02_emu->Scale(norm);
        tripleJetRates_HoE_Ratio06_emu->Scale(norm);
        tripleJetRates_HoE_TP2_Ratio02_emu->Scale(norm);
        tripleJetRates_HoE_TP2_Ratio06_emu->Scale(norm);
        tripleJetRates_HoE_TP1_Ratio06_emu->Scale(norm);
        quadJetRates_emu->Scale(norm);
        quadJetRates_HoE_emu->Scale(norm);
        quadJetRates_HoE_TP05_emu->Scale(norm);
        quadJetRates_HoE_TP5_emu->Scale(norm);
        quadJetRates_HoE_Ratio02_emu->Scale(norm);
        quadJetRates_HoE_Ratio06_emu->Scale(norm);
        quadJetRates_HoE_TP2_Ratio02_emu->Scale(norm);
        quadJetRates_HoE_TP2_Ratio06_emu->Scale(norm);
        quadJetRates_HoE_TP1_Ratio06_emu->Scale(norm);
        singleEgRates_emu->Scale(norm);
        doubleEgRates_emu->Scale(norm);
        singleTauRates_emu->Scale(norm);
        doubleTauRates_emu->Scale(norm);
        singleISOEgRates_emu->Scale(norm);
        doubleISOEgRates_emu->Scale(norm);
        singleISOTauRates_emu->Scale(norm);
        doubleISOTauRates_emu->Scale(norm);
        htSumRates_emu->Scale(norm);
        htSumRates_HoE_emu->Scale(norm);
        htSumRates_HoE_TP05_emu->Scale(norm);
        htSumRates_HoE_TP1_emu->Scale(norm);
        htSumRates_HoE_TP2_emu->Scale(norm);
        htSumRates_HoE_TP3_emu->Scale(norm);
        htSumRates_HoE_TP4_emu->Scale(norm);
        htSumRates_HoE_TP5_emu->Scale(norm);
        htSumRates_HoE_TP05_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP1_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP2_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP3_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP4_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP5_ORHT360_emu->Scale(norm);
        htSumRates_HoE_Ratio02_emu->Scale(norm);
        htSumRates_HoE_Ratio04_emu->Scale(norm);
        htSumRates_HoE_Ratio06_emu->Scale(norm);
        htSumRates_HoE_Ratio08_emu->Scale(norm);
        htSumRates_HoE_Ratio02_ORHT360_emu->Scale(norm);
        htSumRates_HoE_Ratio04_ORHT360_emu->Scale(norm);
        htSumRates_HoE_Ratio06_ORHT360_emu->Scale(norm);
        htSumRates_HoE_Ratio08_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio02_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio04_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio06_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio08_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio02_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio04_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio06_ORHT360_emu->Scale(norm);
        htSumRates_HoE_TP2_Ratio08_ORHT360_emu->Scale(norm);

        mhtSumRates_emu->Scale(norm);
        etSumRates_emu->Scale(norm);
        metSumRates_emu->Scale(norm);
        metHFSumRates_emu->Scale(norm);

        //set the errors for the rates
        //want error -> error * sqrt(norm) ?

        hcalTP_emu->Write();
        hcalTP_Barrel_emu->Write();
        hcalTP_Endcap_emu->Write();
        ecalTP_emu->Write();
        singleJetRates_emu->Write();
        singleJetRates_HoE_emu->Write();
        singleJetRates_HoE_TP05_emu->Write();
        singleJetRates_HoE_TP5_emu->Write();
        singleJetRates_HoE_Ratio02_emu->Write();
        singleJetRates_HoE_Ratio06_emu->Write();
        singleJetRates_HoE_TP2_Ratio02_emu->Write();
        singleJetRates_HoE_TP2_Ratio06_emu->Write();
        singleJetRates_HoE_TP1_Ratio06_emu->Write();
        doubleJetRates_emu->Write();
        doubleJetRates_HoE_emu->Write();
        doubleJetRates_HoE_TP05_emu->Write();
        doubleJetRates_HoE_TP5_emu->Write();
        doubleJetRates_HoE_Ratio02_emu->Write();
        doubleJetRates_HoE_Ratio06_emu->Write();
        doubleJetRates_HoE_TP2_Ratio02_emu->Write();
        doubleJetRates_HoE_TP2_Ratio06_emu->Write();
        doubleJetRates_HoE_TP1_Ratio06_emu->Write();
        tripleJetRates_emu->Write();
        tripleJetRates_HoE_emu->Write();
        tripleJetRates_HoE_TP05_emu->Write();
        tripleJetRates_HoE_TP5_emu->Write();
        tripleJetRates_HoE_Ratio02_emu->Write();
        tripleJetRates_HoE_Ratio06_emu->Write();
        tripleJetRates_HoE_TP2_Ratio02_emu->Write();
        tripleJetRates_HoE_TP2_Ratio06_emu->Write();
        tripleJetRates_HoE_TP1_Ratio06_emu->Write();
        quadJetRates_emu->Write();
        quadJetRates_HoE_emu->Write();
        quadJetRates_HoE_TP05_emu->Write();
        quadJetRates_HoE_TP5_emu->Write();
        quadJetRates_HoE_Ratio02_emu->Write();
        quadJetRates_HoE_Ratio06_emu->Write();
        quadJetRates_HoE_TP2_Ratio02_emu->Write();
        quadJetRates_HoE_TP2_Ratio06_emu->Write();
        quadJetRates_HoE_TP1_Ratio06_emu->Write();
        singleEgRates_emu->Write();
        doubleEgRates_emu->Write();
        singleTauRates_emu->Write();
        doubleTauRates_emu->Write();
        singleISOEgRates_emu->Write();
        doubleISOEgRates_emu->Write();
        singleISOTauRates_emu->Write();
        doubleISOTauRates_emu->Write();
        htSumRates_emu->Write();
        htSumRates_HoE_emu->Write();
        htSumRates_HoE_TP05_emu->Write();
        htSumRates_HoE_TP1_emu->Write();
        htSumRates_HoE_TP2_emu->Write();
        htSumRates_HoE_TP3_emu->Write();
        htSumRates_HoE_TP4_emu->Write();
        htSumRates_HoE_TP5_emu->Write();
        htSumRates_HoE_TP05_ORHT360_emu->Write();
        htSumRates_HoE_TP1_ORHT360_emu->Write();
        htSumRates_HoE_TP2_ORHT360_emu->Write();
        htSumRates_HoE_TP3_ORHT360_emu->Write();
        htSumRates_HoE_TP4_ORHT360_emu->Write();
        htSumRates_HoE_TP5_ORHT360_emu->Write();
        htSumRates_HoE_Ratio02_emu->Write();
        htSumRates_HoE_Ratio04_emu->Write();
        htSumRates_HoE_Ratio06_emu->Write();
        htSumRates_HoE_Ratio08_emu->Write();
        htSumRates_HoE_Ratio02_ORHT360_emu->Write();
        htSumRates_HoE_Ratio04_ORHT360_emu->Write();
        htSumRates_HoE_Ratio06_ORHT360_emu->Write();
        htSumRates_HoE_Ratio08_ORHT360_emu->Write();
        htSumRates_HoE_TP2_Ratio02_emu->Write();
        htSumRates_HoE_TP2_Ratio04_emu->Write();
        htSumRates_HoE_TP2_Ratio06_emu->Write();
        htSumRates_HoE_TP2_Ratio08_emu->Write();
        htSumRates_HoE_TP2_Ratio02_ORHT360_emu->Write();
        htSumRates_HoE_TP2_Ratio04_ORHT360_emu->Write();
        htSumRates_HoE_TP2_Ratio06_ORHT360_emu->Write();
        htSumRates_HoE_TP2_Ratio08_ORHT360_emu->Write();

        mhtSumRates_emu->Write();
        etSumRates_emu->Write();
        metSumRates_emu->Write();
        metHFSumRates_emu->Write();

        hJetEta->Write();
        hJetEt->Write();
        hNJets->Write();

        hHTSum_emu->Write();
        hHTSum_Gen_emu->Write();

        hJetEtaLeading1->Write();
        hJetEtaLeading2->Write();
        hJetEtaLeading3->Write();
        hJetEtaLeading4->Write();

        hJetETLeading1->Write();
        hJetETLeading2->Write();
        hJetETLeading3->Write();
        hJetETLeading4->Write();


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
        HovEtotal_3x3_HETge50_emu->Write();

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
        HovEtotal_3x3_emu_HETge50_GenMatchedJets->Write();
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

        HEEnergytotal_1x1_emu_GenMatched->Write();
        HEEnergytotal_3x3_emu_GenMatched->Write();

        //energy depth plots
    
        energyDepth_Barrel->Write();
        energyDepth_Endcap->Write();

        energyDepth_TPge5_Barrel->Write();
        energyDepth_TPge5_Endcap->Write();

        energyDepth_HT120_Barrel->Write();
        energyDepth_HT120_Endcap->Write();

        energyDepth_L1_Barrel->Write();
        energyDepth_L1_Endcap->Write();

        energyDepth_TPge5_L1_Barrel->Write();
        energyDepth_TPge5_L1_Endcap->Write();

        energyDepth_HT120_L1_Barrel->Write();
        energyDepth_HT120_L1_Endcap->Write();

        energyDepth_LowRatio_Barrel->Write();
        energyDepth_LowRatio_Endcap->Write();

        energyDepth_TPge5_LowRatio_Barrel->Write();
        energyDepth_TPge5_LowRatio_Endcap->Write();

        energyDepth_TPE_TPge5_Barrel->Write();
        energyDepth_TPE_TPge5_Endcap->Write();

        energyDepth_TPge5_HoEcut_Barrel->Write();
        energyDepth_TPge5_HoEcut_Endcap->Write();

        energyDepth_HoEcut_Barrel->Write();
        energyDepth_HoEcut_Endcap->Write();

        energyDepth_TPE_Barrel->Write();
        energyDepth_TPE_Endcap->Write();

        energyDepth_TPE_HoEcut_Barrel->Write();
        energyDepth_TPE_HoEcut_Endcap->Write();

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
        energyDepth_genMatchInclusive_TPge5_Barrel->Write();
        energyDepth_genMatchInclusive_TPge5_Endcap->Write();
        energyDepth_genMatchInclusive_HT120_Barrel->Write();
        energyDepth_genMatchInclusive_HT120_Endcap->Write();
        energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel->Write();
        energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap->Write();
        energyDepth_genMatchInclusive_HoEcut_Barrel->Write();
        energyDepth_genMatchInclusive_HoEcut_Endcap->Write();

        energyDepth_genMatchInclusive_L1_Barrel->Write();
        energyDepth_genMatchInclusive_L1_Endcap->Write();
        energyDepth_genMatchInclusive_TPge5_L1_Barrel->Write();
        energyDepth_genMatchInclusive_TPge5_L1_Endcap->Write();
        energyDepth_genMatchInclusive_HT120_L1_Barrel->Write();
        energyDepth_genMatchInclusive_HT120_L1_Endcap->Write();


        energyDepth_genMatchInclusive_LowRatio_Barrel->Write();
        energyDepth_genMatchInclusive_LowRatio_Endcap->Write();
        energyDepth_genMatchInclusive_TPge5_LowRatio_Barrel->Write();
        energyDepth_genMatchInclusive_TPge5_LowRatio_Endcap->Write();

        energyDepth_TPE_genMatchInclusive_Barrel->Write();
        energyDepth_TPE_genMatchInclusive_Endcap->Write();

        energyDepth_TPE_genMatchInclusive_TPge5_Barrel->Write();
        energyDepth_TPE_genMatchInclusive_TPge5_Endcap->Write();

        energyDepth_TPE_genMatchInclusive_HoEcut_Barrel->Write();
        energyDepth_TPE_genMatchInclusive_HoEcut_Endcap->Write();

        TPEt_HB->Write();
        TPEt_HE->Write();

        TPEt_Matched_HB->Write();
        TPEt_Matched_HE->Write();

        TPEt_LowRatio_HB->Write();
        TPEt_LowRatio_HE->Write();

        TPEt_LowRatio_Matched_HB->Write();
        TPEt_LowRatio_Matched_HE->Write();

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

        energyDepth_NTPs_HBD4_HED47_Max->Write();
        energyDepth_NTPs_HBD4_HED47_Max_Gen->Write();
        energyDepth_NTPs_HBD4_HED47_Max_HT120->Write();
        energyDepth_NTPs_HBD4_HED47_Max_Gen_HT120->Write();                
        energyDepth_NTPs_HBD4_HED47_Max_HT360->Write();
        energyDepth_NTPs_HBD4_HED47_Max_Gen_HT360->Write();                

        energyDepth_HT120_DepthEnergy1_HB->Write();
        energyDepth_HT120_DepthEnergy1_HE->Write();
        energyDepth_HT120_DepthEnergy2_HB->Write();
        energyDepth_HT120_DepthEnergy2_HE->Write();
        energyDepth_HT120_DepthEnergy3_HB->Write();
        energyDepth_HT120_DepthEnergy3_HE->Write();
        energyDepth_HT120_DepthEnergy4_HB->Write();
        energyDepth_HT120_DepthEnergy4_HE->Write();
        energyDepth_HT120_DepthEnergy5_HE->Write();
        energyDepth_HT120_DepthEnergy6_HE->Write();
        energyDepth_HT120_DepthEnergy7_HE->Write();

        energyDepth_HT120_Gen_DepthEnergy1_HB->Write();
        energyDepth_HT120_Gen_DepthEnergy1_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy2_HB->Write();
        energyDepth_HT120_Gen_DepthEnergy2_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy3_HB->Write();
        energyDepth_HT120_Gen_DepthEnergy3_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy4_HB->Write();
        energyDepth_HT120_Gen_DepthEnergy4_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy5_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy6_HE->Write();
        energyDepth_HT120_Gen_DepthEnergy7_HE->Write();

        energyDepth_DepthEnergy1_HB->Write();
        energyDepth_DepthEnergy1_HE->Write();
        energyDepth_DepthEnergy2_HB->Write();
        energyDepth_DepthEnergy2_HE->Write();
        energyDepth_DepthEnergy3_HB->Write();
        energyDepth_DepthEnergy3_HE->Write();
        energyDepth_DepthEnergy4_HB->Write();
        energyDepth_DepthEnergy4_HE->Write();
        energyDepth_DepthEnergy5_HE->Write();
        energyDepth_DepthEnergy6_HE->Write();
        energyDepth_DepthEnergy7_HE->Write();

        energyDepth_Gen_DepthEnergy1_HB->Write();
        energyDepth_Gen_DepthEnergy1_HE->Write();
        energyDepth_Gen_DepthEnergy2_HB->Write();
        energyDepth_Gen_DepthEnergy2_HE->Write();
        energyDepth_Gen_DepthEnergy3_HB->Write();
        energyDepth_Gen_DepthEnergy3_HE->Write();
        energyDepth_Gen_DepthEnergy4_HB->Write();
        energyDepth_Gen_DepthEnergy4_HE->Write();
        energyDepth_Gen_DepthEnergy5_HE->Write();
        energyDepth_Gen_DepthEnergy6_HE->Write();
        energyDepth_Gen_DepthEnergy7_HE->Write();

        energyDepth_Ratio_Gen_Depth4_HB->Write();
        energyDepth_Ratio_Gen_Depth4_7_HE->Write();

        energyDepth_Ratio_Depth4_HB->Write();
        energyDepth_Ratio_Depth4_7_HE->Write();

        energyDepth_Ratio_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_Depth3_4_HB->Write();
        energyDepth_Ratio_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_Gen_Depth4_HB->Write();
        energyDepth_Ratio_HT120_Gen_Depth4_7_HE->Write();
        
        energyDepth_Ratio_HT120_Depth4_HB->Write();
        energyDepth_Ratio_HT120_Depth4_7_HE->Write();

        energyDepth_Ratio_HT120_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_HT120_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_TP1_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_TP1_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_TP2_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_TP2_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_TP3_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_TP3_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_TP4_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_TP4_Depth2_7_HE->Write();

        energyDepth_Ratio_HT120_TP5_Depth3_4_HB->Write();
        energyDepth_Ratio_HT120_TP5_Depth2_7_HE->Write();

        energyDepth_Ratio_TPge10_Gen_Depth4_HB->Write();
        energyDepth_Ratio_TPge10_Gen_Depth4_7_HE->Write();
        
        energyDepth_Ratio_TPge10_Depth4_HB->Write();
        energyDepth_Ratio_TPge10_Depth4_7_HE->Write();

        energyDepth_Ratio_TPge10_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_TPge10_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_TPge10_Depth3_4_HB->Write();
        energyDepth_Ratio_TPge10_Depth2_7_HE->Write();

        energyDepth_Ratio_IEta_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_IEta_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_Depth2_7_HE->Write();

        energyDepth_Ratio_IEta_HT120_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_HT120_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_IEta_HT120_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_HT120_Depth2_7_HE->Write();

        energyDepth_Ratio_IEta_TPge10_Gen_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_TPge10_Gen_Depth2_7_HE->Write();
        
        energyDepth_Ratio_IEta_TPge10_Depth3_4_HB->Write();
        energyDepth_Ratio_IEta_TPge10_Depth2_7_HE->Write();

        
        hcalTP_nearL1Jet_emu->Write();
        hcalTP_nearL1Jet_Gen_emu->Write();
        hcalTP_nearL1Jet_Barrel_emu->Write();
        hcalTP_nearL1Jet_Gen_Barrel_emu->Write();
        hcalTP_nearL1Jet_Endcap_emu->Write();
        hcalTP_nearL1Jet_Gen_Endcap_emu->Write();

        effJetID_HoE_DepthFB_TPE->Scale(1/effJetID_HoE_DepthFB_TPE->GetBinContent(1));
        effJetID_HoE_DepthFB_TPE->Write();

        effJetID_HoE_DepthFB_TPE_Gen->Scale(1/effJetID_HoE_DepthFB_TPE_Gen->GetBinContent(1));
        effJetID_HoE_DepthFB_TPE_Gen->Write();

        effJetID_HoE_DepthFB_HTscan->Scale(1/effJetID_HoE_DepthFB_HTscan->GetBinContent(1));
        effJetID_HoE_DepthFB_HTscan->Write();

        effJetID_HoE_DepthFB_Gen_HTscan->Scale(1/effJetID_HoE_DepthFB_Gen_HTscan->GetBinContent(1));
        effJetID_HoE_DepthFB_Gen_HTscan->Write();

        effJetID_HoE_DepthFB_Jetscan->Scale(1/effJetID_HoE_DepthFB_Jetscan->GetBinContent(1));
        effJetID_HoE_DepthFB_Jetscan->Write();

        effJetID_HoE_DepthFB_Gen_Jetscan->Scale(1/effJetID_HoE_DepthFB_Gen_Jetscan->GetBinContent(1));
        effJetID_HoE_DepthFB_Gen_Jetscan->Write();

        effJetID_HoE_DepthFB_Ratio->Scale(1/effJetID_HoE_DepthFB_Ratio->GetBinContent(1));
        effJetID_HoE_DepthFB_Ratio->Write();

        effJetID_HoE_DepthFB_Ratio_Gen->Scale(1/effJetID_HoE_DepthFB_Ratio_Gen->GetBinContent(1));
        effJetID_HoE_DepthFB_Ratio_Gen->Write();

        
        //FOR GEN MATCHING STUDY
        hJetGenPartDR_LLPdaught->Write();
        hJetGenPartDR_LLP->Write();
        hJetGenPartDR_LLP_inHCAL->Write();
        hJetGenPartDRfromVertex_LLP_inHCAL->Write();

        hLLP_vertex_DR->Write();

        hNMatchedLLP_inHCAL_DR02->Write();
        hNMatchedLLP_inHCAL_DR05->Write();
        hNMatchedLLP_DR02->Write();
        hNMatchedLLP_DR05->Write();
        hNMatched_LLPdaught_DR02->Write();
        hNMatched_LLPdaught_DR05->Write();
        
        hfracMatched_LLPdaught_DR02->Write();
        hfracMatched_LLPdaught_DR05->Write();
        hNLLPdaughts_inHCAL->Write();
        hNLLPdaughts->Write();
        hNLLPdaughts_pteta->Write();
        betagammaLLP->Write();
        betaLLP->Write();
        velocityLLP->Write();

//        double total_LLP = vertex_HB_LLPD->Integral() + vertex_HE_LLPD->Integral();
//        vertex_HB_LLPD->Scale(1/total_LLP);
//        vertex_HE_LLPD->Scale(1/total_LLP);
        fracLLP_Separated->GetXaxis()->SetBinLabel(1, "HB");
        fracLLP_Separated->GetXaxis()->SetBinLabel(2, "HE");
        fracLLP_Separated->GetXaxis()->SetBinLabel(3, "Neither");
        
        fracLLP_Separated->Write();
        vertex_HB_LLPD->Write();
        vertex_HE_LLPD->Write();

        TPtextfile.close();
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
