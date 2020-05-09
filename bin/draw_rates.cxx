// From HCAL Trigger Validation framework, edited by Gillian Kopp for LLP trigger multiplicity studies (2020)

 #include "PhysicsTools/Utilities/macros/setTDRStyle.C"
 #include "TCanvas.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TProfile.h"
 #include "TFile.h"
 #include "TLegend.h"
 #include "TROOT.h"
 #include "TGraph.h"
 #include "TMarker.h"
 #include <map>
 #include <string>
 #include <vector>
 #include <iostream>

 int main()
 {
   // include comparisons between HW and data TPs
   bool includeHW = false;
   int rebinFactor = 1;

   setTDRStyle();
   gROOT->ForceStyle();

   // default, then new conditions
   // files for L1 rates
   //   std::vector<std::string> filenames = {"rates_def_nugunTDC.root", "rates_new_cond_nugunTDC_1L1Jet.root"};
   std::vector<std::string> filenames = {"rates_def_nugunTDC.root", "rates_new_cond_nugunTDC_4L1Jets.root"};
   //  std::vector<std::string> filenames = {"rates_new_cond_QCD.root", "rates_new_cond_pl1000.root"};
   std::vector<std::string> rateTypes = {"singleJet", "doubleJet", "tripleJet", "quadJet",
					 "singleJetGlobal", "doubleJetGlobal", "tripleJetGlobal", "quadJetGlobal",
					 "singleEg", "singleISOEg", "doubleEg", "doubleISOEg",
					 "singleTau", "singleISOTau", "doubleTau", "doubleISOTau",
					 "htSum", "etSum", "metSum", "metHFSum",
					 "htSumGlobal", "etSumGlobal"};
   // files for multiplicity overlay plots
   //   std::vector<std::string> mult_filenames = {"rates_new_cond_pl10000_1L1Jet.root", "rates_new_cond_pl1000_1L1Jet.root", "rates_new_cond_pl500_1L1Jet.root", "rates_new_cond_QCD_1L1Jet.root"};
   std::vector<std::string> mult_filenames = {"rates_new_cond_pl500_noPU_4L1Jets.root", "rates_new_cond_pl1000_4L1Jets.root", "rates_new_cond_pl500_4L1Jets.root", "rates_new_cond_QCD_4L1Jets.root"};
   // these are multiplicity files (used for timescan and overlay plots) as well as ratio type plts
   std::vector<std::string> multTypes = {//"dt3GeV1ns","dt3GeV2ns","dt3GeV3ns","dt3GeV4ns","dt3GeV5ns",
					 //					 "dt3GeV1nsHE","dt3GeV2nsHE","dt3GeV3nsHE","dt3GeV4nsHE","dt3GeV5nsHE",
					 "dt3GeV1nsHB","dt3GeV2nsHB","dt3GeV3nsHB","dt3GeV4nsHB","dt3GeV5nsHB",
					 //					 "dt2GeV1ns","dt2GeV2ns","dt2GeV3ns","dt2GeV4ns","dt2GeV5ns",
					 //					 "dt2GeV1nsHE","dt2GeV2nsHE","dt2GeV3nsHE","dt2GeV4nsHE","dt2GeV5nsHE",
					 "dt2GeV1nsHB","dt2GeV2nsHB","dt2GeV3nsHB","dt2GeV4nsHB","dt2GeV5nsHB",
					 //					 "dt1GeV1ns","dt1GeV2ns","dt1GeV3ns","dt1GeV4ns","dt1GeV5ns",
					 //					 "dt1GeV1nsHE","dt1GeV2nsHE","dt1GeV3nsHE","dt1GeV4nsHE","dt1GeV5nsHE",
					 "dt1GeV1nsHB","dt1GeV2nsHB","dt1GeV3nsHB","dt1GeV4nsHB","dt1GeV5nsHB",
					 //					 "dt3GeV1nsJet","dt3GeV2nsJet","dt3GeV3nsJet","dt3GeV4nsJet","dt3GeV5nsJet",
					 //					 "dt3GeV1nsHEJet","dt3GeV2nsHEJet","dt3GeV3nsHEJet","dt3GeV4nsHEJet","dt3GeV5nsHEJet",
					 "dt3GeV1nsHBJet","dt3GeV2nsHBJet","dt3GeV3nsHBJet","dt3GeV4nsHBJet","dt3GeV5nsHBJet",
                                         "dt0GeV1nsHBJet","dt0GeV2nsHBJet","dt0GeV3nsHBJet","dt0GeV4nsHBJet","dt0GeV5nsHBJet",
					 //					 "dt3GeV3nsHBnearJet","dt3GeV3nsHBnearJet1","dt3GeV3nsHBnearJet2","dt3GeV3nsHBnearJet3","dt3GeV3nsHBnearJet4",
					 "dt3GeV3nsHBJet1","dt3GeV3nsHBJet2","dt3GeV3nsHBJet3","dt3GeV3nsHBJet4",
					 "dt3GeV3nsHBQuadJet","dt3GeV3nsHBTripleJet","dt3GeV3nsHBDoubleJet","dt3GeV3nsHBSingleJet",
                                         "dt0GeV5nsHBJet1","dt0GeV5nsHBJet2","dt0GeV5nsHBJet3","dt0GeV5nsHBJet4",
                                         "dt0GeV5nsHBQuadJet","dt0GeV5nsHBTripleJet","dt0GeV5nsHBDoubleJet","dt0GeV5nsHBSingleJet",
					 //					 "dt2GeV1nsJet","dt2GeV2nsJet","dt2GeV3nsJet","dt2GeV4nsJet","dt2GeV5nsJet",
					 //					 "dt2GeV1nsHEJet","dt2GeV2nsHEJet","dt2GeV3nsHEJet","dt2GeV4nsHEJet","dt2GeV5nsHEJet",
					 "dt2GeV1nsHBJet","dt2GeV2nsHBJet","dt2GeV3nsHBJet","dt2GeV4nsHBJet","dt2GeV5nsHBJet",
					 //					 "dt1GeV1nsJet","dt1GeV2nsJet","dt1GeV3nsJet","dt1GeV4nsJet","dt1GeV5nsJet",
					 //					 "dt1GeV1nsHEJet","dt1GeV2nsHEJet","dt1GeV3nsHEJet","dt1GeV4nsHEJet","dt1GeV5nsHEJet",
					 "dt1GeV1nsHBJet","dt1GeV2nsHBJet","dt1GeV3nsHBJet","dt1GeV4nsHBJet","dt1GeV5nsHBJet",
					 "Ratio_Depth", "Ratio_DepthHE", "Ratio_DepthHB","Ratio_Depth_Jets", "Ratio_DepthHE_Jets", "Ratio_DepthHB_Jets","centralTiming",
					 "DepthVariable",
					 "DelayedHitFraction2GeV1ns","DelayedHitFraction2GeV2ns","DelayedHitFraction2GeV3ns","DelayedHitFraction2GeV4ns","DelayedHitFraction2GeV5ns","DelayedHitFraction2GeV6ns","DelayedHitFraction2GeV7ns",
					 "DelayedHitFraction3GeV1ns","DelayedHitFraction3GeV2ns","DelayedHitFraction3GeV3ns","DelayedHitFraction3GeV4ns","DelayedHitFraction3GeV5ns","DelayedHitFraction3GeV6ns","DelayedHitFraction3GeV7ns",
					 "DelayedHitFraction4GeV1ns","DelayedHitFraction4GeV2ns","DelayedHitFraction4GeV3ns","DelayedHitFraction4GeV4ns","DelayedHitFraction4GeV5ns","DelayedHitFraction4GeV6ns","DelayedHitFraction4GeV7ns",
					 "dt1GeVcaloT1","dt1GeVcaloT2","dt1GeVcaloT3","dt1GeVcaloT4",
					 "dt2GeVcaloT1","dt2GeVcaloT2","dt2GeVcaloT3","dt2GeVcaloT4",
					 "dt3GeVcaloT1","dt3GeVcaloT2","dt3GeVcaloT3","dt3GeVcaloT4"};

   std::vector<std::string> EDepthTypes = {"Timing_Depth","Timing_DepthHE","Timing_DepthHB","Energy_Depth_HighE","Energy_DepthHE_HighE","Energy_DepthHB_HighE","Timing_Depth_Jets","Timing_DepthHE_Jets","Timing_DepthHB_Jets","Energy_Depth_Jets_HighE","Energy_DepthHE_Jets_HighE","Energy_DepthHB_Jets_HighE"};

   std::vector<std::string> RatioTypes = {"Ratio_Depth", "Ratio_DepthHE", "Ratio_DepthHB","Ratio_Depth_Jets", "Ratio_DepthHE_Jets", "Ratio_DepthHB_Jets"};

  std::map<std::string, int> histColor;
  histColor["singleJet"] = histColor["singleJetGlobal"] = histColor["singleEg"] = histColor["singleTau"] = histColor["etSum"] = histColor["etSumGlobal"] = histColor["metSum"] = histColor["dt3GeV1ns"] = histColor["dt3GeV1nsHE"] =histColor["dt3GeV1nsHB"] = histColor["dt2GeV1ns"] = histColor["dt2GeV1nsHE"] = histColor["dt2GeV1nsHB"] = histColor["dt1GeV1ns"] = histColor["dt1GeV1nsHE"] = histColor["dt1GeV1nsHB"] = histColor["dt3GeV1nsJet"] = histColor["dt3GeV1nsHEJet"] =histColor["dt3GeV1nsHBJet"] = histColor["dt2GeV1nsJet"] = histColor["dt2GeV1nsHEJet"] = histColor["dt2GeV1nsHBJet"] = histColor["dt1GeV1nsJet"] = histColor["dt1GeV1nsHEJet"] = histColor["dt1GeV1nsHBJet"] = kRed;
  histColor["doubleJet"] = histColor["doubleJetGlobal"] = histColor["singleISOEg"] = histColor["singleISOTau"] = histColor["htSum"] = histColor["htSumGlobal"] = histColor["metHFSum"] = histColor["dt3GeV2ns"] = histColor["dt3GeV2nsHE"] = histColor["dt3GeV2nsHB"] = histColor["dt2GeV2ns"] = histColor["dt2GeV2nsHE"] = histColor["dt2GeV2nsHB"] = histColor["dt1GeV2ns"] = histColor["dt1GeV2nsHE"] = histColor["dt1GeV2nsHB"] = histColor["dt3GeV2nsJet"] = histColor["dt3GeV2nsHEJet"] = histColor["dt3GeV2nsHBJet"] = histColor["dt2GeV2nsJet"] = histColor["dt2GeV2nsHEJet"] = histColor["dt2GeV2nsHBJet"] = histColor["dt1GeV2nsJet"] = histColor["dt1GeV2nsHEJet"] = histColor["dt1GeV2nsHBJet"] = kBlue;
  histColor["tripleJet"] = histColor["tripleJetGlobal"] = histColor["doubleEg"] = histColor["doubleTau"] = histColor["dt3GeV3ns"] = histColor["dt3GeV3nsHE"] = histColor["dt3GeV3nsHB"] = histColor["dt2GeV3ns"] = histColor["dt1GeV3ns"] =histColor["dt2GeV3nsHE"] = histColor["dt1GeV3nsHE"] = histColor["dt2GeV3nsHB"] = histColor["dt1GeV3nsHB"] = histColor["dt3GeV3nsJet"] = histColor["dt3GeV3nsHEJet"] = histColor["dt3GeV3nsHBJet"] = histColor["dt0GeV1nsHBJet"] =histColor["dt0GeV2nsHBJet"] =histColor["dt0GeV3nsHBJet"] =histColor["dt0GeV4nsHBJet"] =histColor["dt0GeV5nsHBJet"] = histColor["dt3GeV3nsHBJet1"] = histColor["dt3GeV3nsHBJet2"]= histColor["dt3GeV3nsHBJet3"] = histColor["dt3GeV3nsHBJet4"] = histColor["dt2GeV3nsJet"] = histColor["dt1GeV3nsJet"] =histColor["dt2GeV3nsHEJet"] = histColor["dt1GeV3nsHEJet"] = histColor["dt2GeV3nsHBJet"] = histColor["dt1GeV3nsHBJet"] =histColor["dt0GeV5nsHBJet1"] = histColor["dt0GeV5nsHBJet2"]= histColor["dt0GeV5nsHBJet3"] = histColor["dt0GeV5nsHBJet4"] = kGreen+1;
  histColor["quadJet"] = histColor["quadJetGlobal"] = histColor["doubleISOEg"] = histColor["doubleISOTau"] = histColor["dt3GeV4ns"] = histColor["dt3GeV4nsHE"] = histColor["dt3GeV4nsHB"] = histColor["dt2GeV4ns"] = histColor["dt1GeV4ns"] = histColor["dt2GeV4nsHE"] = histColor["dt1GeV4nsHE"] = histColor["dt2GeV4nsHB"] = histColor["dt1GeV4nsHB"] = histColor["dt3GeV4nsJet"] = histColor["dt3GeV4nsHEJet"] = histColor["dt3GeV4nsHBJet"] = histColor["dt2GeV4nsJet"] = histColor["dt1GeV4nsJet"] = histColor["dt2GeV4nsHEJet"] = histColor["dt1GeV4nsHEJet"] = histColor["dt2GeV4nsHBJet"] = histColor["dt1GeV4nsHBJet"] = kBlack;
  histColor["dt3GeV5ns"] = histColor["dt3GeV5nsHE"] = histColor["dt3GeV5nsHB"] = histColor["dt2GeV5ns"] = histColor["dt1GeV5ns"]  = histColor["dt2GeV5nsHE"] = histColor["dt1GeV5nsHE"] = histColor["dt2GeV5nsHB"] = histColor["dt1GeV5nsHB"] = histColor["dt3GeV5nsJet"] = histColor["dt3GeV5nsHEJet"] = histColor["dt3GeV5nsHBJet"] = histColor["dt2GeV5nsJet"] = histColor["dt1GeV5nsJet"]  = histColor["dt2GeV5nsHEJet"] = histColor["dt1GeV5nsHEJet"] = histColor["dt2GeV5nsHBJet"] = histColor["dt1GeV5nsHBJet"] = kCyan;

  std::map<std::string, TH1F*> rateHists_def;
  std::map<std::string, TH1F*> rateHists_new_cond;
  std::map<std::string, TH1F*> rateHists_hw;
  std::map<std::string, TH1F*> rateHistsRatio;
  
  std::map<std::string, TH1F*> multHists_QCD;
  std::map<std::string, TH1F*> multHists_LLP10000;
  std::map<std::string, TH1F*> multHists_LLP1000;
  std::map<std::string, TH1F*> multHists_LLP500;
  std::map<std::string, TH1F*> multHistsRatio;
  std::map<std::string, TH1F*> multHists_hw;

  std::map<std::string, TH2F*> energy_depth_QCD;
  std::map<std::string, TH2F*> energy_depth_LLP10000;
  std::map<std::string, TH2F*> energy_depth_LLP1000;
  std::map<std::string, TH2F*> energy_depth_LLP500;

  std::map<std::string, TH1D*> energy_profile_QCD_overlay;
  std::map<std::string, TH1D*> energy_profile_LLP500_overlay;
  std::map<std::string, TH1D*> energy_profile_LLP1000_overlay;
  std::map<std::string, TH1D*> energy_profile_LLP10000_overlay;

  int SJet60GeV = 0;
  int htSum350GeV = 0;
  int SJet60GeV_l = 0;
  int htSum350GeV_l = 0;

  std::vector<TFile*> files;
  for(auto file : filenames) {
    files.push_back(TFile::Open(file.c_str()));
  }
  // making rate plots for current and new conditions - filling histograms from the root files
  for(auto rateType : rateTypes) {
    std::string histName(rateType);
    std::string histNameHw(histName);
    histName += "Rates_emu";
    //    std::cout <<histName<< std::endl;
    histNameHw += "Rates_hw";
    rateHists_def[rateType] = dynamic_cast<TH1F*>(files.at(0)->Get(histName.c_str()));
    rateHists_hw[rateType] = dynamic_cast<TH1F*>(files.at(0)->Get(histNameHw.c_str()));
    rateHists_new_cond[rateType] = dynamic_cast<TH1F*>(files.at(1)->Get(histName.c_str())); 
    rateHists_def[rateType]->Rebin(rebinFactor);
    rateHists_hw[rateType]->Rebin(rebinFactor);
    rateHists_new_cond[rateType]->Rebin(rebinFactor);
    rateHists_def[rateType]->SetLineColor(histColor[rateType]);
    rateHists_hw[rateType]->SetLineColor(histColor[rateType]);
    rateHists_new_cond[rateType]->SetLineColor(histColor[rateType]);
    // for the rate and efficiency plots
    if ( rateType == "singleJetGlobal" ) {
      int xval = rateHists_new_cond[rateType]->GetXaxis()->FindBin(60); // get x value of bin of interest
      SJet60GeV = rateHists_new_cond[rateType]->GetBinContent(xval);
    }
    if ( rateType == "htSumGlobal" ) {
      int xval = rateHists_new_cond[rateType]->GetXaxis()->FindBin(350);
      htSum350GeV = rateHists_new_cond[rateType]->GetBinContent(xval);
    }
    if ( rateType == "singleJet" ) {
      int xval = rateHists_new_cond[rateType]->GetXaxis()->FindBin(60); // get x value of bin of interest   
      SJet60GeV_l = rateHists_new_cond[rateType]->GetBinContent(xval);
    }
    if ( rateType == "htSum" ) {
      int xval = rateHists_new_cond[rateType]->GetXaxis()->FindBin(350);
      htSum350GeV_l = rateHists_new_cond[rateType]->GetBinContent(xval);
    }
    TString name(rateHists_new_cond[rateType]->GetName());
    name += "_ratio";
    if(includeHW) {
      rateHistsRatio[rateType] = dynamic_cast<TH1F*>(rateHists_def[rateType]->Clone(name));
      rateHistsRatio[rateType]->Divide(rateHists_hw[rateType]);
    }
    else {
      rateHistsRatio[rateType] = dynamic_cast<TH1F*>(rateHists_new_cond[rateType]->Clone(name));
      rateHistsRatio[rateType]->Divide(rateHists_def[rateType]);
    }
    rateHistsRatio[rateType]->SetMinimum(-0.2);    // -0.5 for singleJet  // previously 0.6
    rateHistsRatio[rateType]->SetMaximum(1.4);    // 80 for singleJet // previously 1.4
    if ((rateType == "singleJet") || (rateType ==  "doubleJet") || (rateType ==  "tripleJet") || (rateType ==  "quadJet")) {
      rateHistsRatio[rateType]->SetMaximum(0.4);
      rateHistsRatio[rateType]->SetMinimum(0);
    }
    if ((rateType ==  "htSum") || (rateType ==  "etSum")) {
      rateHistsRatio[rateType]->SetMaximum(0.04);
      rateHistsRatio[rateType]->SetMinimum(0);
    }
    rateHistsRatio[rateType]->SetLineWidth(2);    
  }

  std::vector<TFile*> mult_files;
  for(auto file : mult_filenames) {
    mult_files.push_back(TFile::Open(file.c_str()));
  }
  // opening the files for multiplicity and energy ratio plots - filling histograms from root files
  for(auto multType : multTypes) {
    std::string histName(multType);
    std::string histNameHw(histName);
    if (multType.substr(0,2) == "dt" ) {
      histName += "Mult_emu"; // multiplicity plots all have this appended to name in histogram, example "dt3GeV2nsMult_emu". Ratio plots do not have this added, nor does central timing
    }
    multHists_QCD[multType]  = dynamic_cast<TH1F*>(mult_files.at(3)->Get(histName.c_str()));
    multHists_LLP10000[multType] = dynamic_cast<TH1F*>(mult_files.at(0)->Get(histName.c_str()));
    multHists_LLP1000[multType] = dynamic_cast<TH1F*>(mult_files.at(1)->Get(histName.c_str()));
    multHists_LLP500[multType] = dynamic_cast<TH1F*>(mult_files.at(2)->Get(histName.c_str()));
    
    multHists_QCD[multType]->Rebin(rebinFactor);
    multHists_LLP10000[multType]->Rebin(rebinFactor);
    multHists_LLP1000[multType]->Rebin(rebinFactor);
    multHists_LLP500[multType]->Rebin(rebinFactor);

    multHists_QCD[multType]->SetLineColor(histColor[multType]);
    multHists_LLP10000[multType]->SetLineColor(histColor[multType]);
    multHists_LLP1000[multType]->SetLineColor(histColor[multType]);
    multHists_LLP500[multType]->SetLineColor(histColor[multType]);
  }

  // get histograms for energy depth similar for multiplicty. Done in separate loop because these are TH2F instead of TH1F
  for(auto EDepthType : EDepthTypes) {
    std::string histName(EDepthType);
    energy_depth_QCD[EDepthType]  = dynamic_cast<TH2F*>(mult_files.at(3)->Get(histName.c_str()));
    energy_depth_LLP10000[EDepthType] = dynamic_cast<TH2F*>(mult_files.at(0)->Get(histName.c_str()));
    energy_depth_LLP1000[EDepthType] = dynamic_cast<TH2F*>(mult_files.at(1)->Get(histName.c_str()));
    energy_depth_LLP500[EDepthType] = dynamic_cast<TH2F*>(mult_files.at(2)->Get(histName.c_str()));
  }

  for(auto pair : rateHists_new_cond) pair.second->SetLineWidth(2);
  for(auto pair : rateHists_hw) pair.second->SetLineStyle(kDashed);
  for(auto pair : rateHists_def) pair.second->SetLineStyle(kDotted);

  for(auto pair : multHists_LLP10000) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLP1000) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLP500) pair.second->SetLineWidth(2);
  for(auto pair : multHists_QCD) pair.second->SetLineWidth(2);

  // jet, eg, tau, energy rate plot types
  std::vector<std::string> jetPlots = {"singleJet", "doubleJet", "tripleJet", "quadJet"};
  std::vector<std::string> jetPlotsGlobal = {"singleJetGlobal", "doubleJetGlobal", "tripleJetGlobal", "quadJetGlobal"};
  std::vector<std::string> egPlots = {"singleEg", "singleISOEg", "doubleEg", "doubleISOEg"};
  std::vector<std::string> tauPlots = {"singleTau", "singleISOTau", "doubleTau", "doubleISOTau"};
  std::vector<std::string> scalarSumPlots = {"etSum", "htSum"};
  std::vector<std::string> scalarSumPlotsGlobal = {"etSumGlobal", "htSumGlobal"};
  std::vector<std::string> vectorSumPlots = {"metSum", "metHFSum"};
  // multiplicity plot types 
  //  std::vector<std::string> multPlots3GeV = {"dt3GeV5ns","dt3GeV4ns","dt3GeV3ns","dt3GeV2ns","dt3GeV1ns"};
  //  std::vector<std::string> multPlots3GeVHE = {"dt3GeV5nsHE","dt3GeV4nsHE","dt3GeV3nsHE","dt3GeV2nsHE","dt3GeV1nsHE"};
  std::vector<std::string> multPlots3GeVHB = {"dt3GeV5nsHB","dt3GeV4nsHB","dt3GeV3nsHB","dt3GeV2nsHB","dt3GeV1nsHB"};
  //  std::vector<std::string> multPlots2GeV = {"dt2GeV5ns","dt2GeV4ns","dt2GeV3ns","dt2GeV2ns","dt2GeV1ns"};
  //  std::vector<std::string> multPlots2GeVHE = {"dt2GeV5nsHE","dt2GeV4nsHE","dt2GeV3nsHE","dt2GeV2nsHE","dt2GeV1nsHE"};
  std::vector<std::string> multPlots2GeVHB = {"dt2GeV5nsHB","dt2GeV4nsHB","dt2GeV3nsHB","dt2GeV2nsHB","dt2GeV1nsHB"};
  //  std::vector<std::string> multPlots1GeV = {"dt1GeV5ns","dt1GeV4ns","dt1GeV3ns","dt1GeV2ns","dt1GeV1ns"}; 
  //  std::vector<std::string> multPlots1GeVHE = {"dt1GeV5nsHE","dt1GeV4nsHE","dt1GeV3nsHE","dt1GeV2nsHE","dt1GeV1nsHE"};
  std::vector<std::string> multPlots1GeVHB = {"dt1GeV5nsHB","dt1GeV4nsHB","dt1GeV3nsHB","dt1GeV2nsHB","dt1GeV1nsHB"};
  // multiplicity plot types for matched with L1 Jets
  //  std::vector<std::string> multPlots3GeV_Jet = {"dt3GeV5nsJet","dt3GeV4nsJet","dt3GeV3nsJet","dt3GeV2nsJet","dt3GeV1nsJet"};
  //  std::vector<std::string> multPlots3GeVHE_Jet = {"dt3GeV5nsHEJet","dt3GeV4nsHEJet","dt3GeV3nsHEJet","dt3GeV2nsHEJet","dt3GeV1nsHEJet"};
  std::vector<std::string> multPlots3GeVHB_Jet = {"dt3GeV5nsHBJet","dt3GeV4nsHBJet","dt3GeV3nsHBJet","dt3GeV2nsHBJet","dt3GeV1nsHBJet"};
  std::vector<std::string> multPlots0GeVHB_Jet = {"dt0GeV5nsHBJet","dt0GeV4nsHBJet","dt0GeV3nsHBJet","dt0GeV2nsHBJet","dt0GeV1nsHBJet"};
  //  std::vector<std::string> multPlots3GeVHB_Jet_near = {"dt3GeV3nsHBnearJet","dt3GeV3nsHBnearJet1","dt3GeV3nsHBnearJet2","dt3GeV3nsHBnearJet3","dt3GeV3nsHBnearJet4"}; // from associating HCAL TPs to one of four L1 jets and then DR restrictions
  std::vector<std::string> multPlots3GeVHB_Jet_L1DRcone = {"dt3GeV3nsHBJet1","dt3GeV3nsHBJet2","dt3GeV3nsHBJet3","dt3GeV3nsHBJet4"}; // just DR restrictions around a L1 jet
  std::vector<std::string> multPlots0GeVHB_Jet_L1DRcone = {"dt0GeV5nsHBJet1","dt0GeV5nsHBJet2","dt0GeV5nsHBJet3","dt0GeV5nsHBJet4"}; // just DR restrictions around a L1 jet           
  //  std::vector<std::string> multPlots2GeV_Jet = {"dt2GeV5nsJet","dt2GeV4nsJet","dt2GeV3nsJet","dt2GeV2nsJet","dt2GeV1nsJet"};
  //  std::vector<std::string> multPlots2GeVHE_Jet = {"dt2GeV5nsHEJet","dt2GeV4nsHEJet","dt2GeV3nsHEJet","dt2GeV2nsHEJet","dt2GeV1nsHEJet"};
  std::vector<std::string> multPlots2GeVHB_Jet = {"dt2GeV5nsHBJet","dt2GeV4nsHBJet","dt2GeV3nsHBJet","dt2GeV2nsHBJet","dt2GeV1nsHBJet"};
  //  std::vector<std::string> multPlots1GeV_Jet = {"dt1GeV5nsJet","dt1GeV4nsJet","dt1GeV3nsJet","dt1GeV2nsJet","dt1GeV1nsJet"};
  //  std::vector<std::string> multPlots1GeVHE_Jet = {"dt1GeV5nsHEJet","dt1GeV4nsHEJet","dt1GeV3nsHEJet","dt1GeV2nsHEJet","dt1GeV1nsHEJet"};
  std::vector<std::string> multPlots1GeVHB_Jet = {"dt1GeV5nsHBJet","dt1GeV4nsHBJet","dt1GeV3nsHBJet","dt1GeV2nsHBJet","dt1GeV1nsHBJet"};
  
  std::vector<TCanvas*> canvases;
  std::vector<TPad*> pad1;
  std::vector<TPad*> pad2;
  std::map<std::string, std::vector<std::string> > plots;
  plots["jet"] = jetPlots;
  plots["jetGlobal"] = jetPlotsGlobal;
  plots["eg"] = egPlots;
  plots["tau"] = tauPlots;
  plots["scalarSum"] = scalarSumPlots;
  plots["scalarSumGlobal"] = scalarSumPlotsGlobal;
  plots["vectorSum"] = vectorSumPlots;

  std::map<std::string, std::vector<std::string> > mult_plots;
  //  mult_plots["3GeV_timescan"] = multPlots3GeV;
  //  mult_plots["3GeV_timescanHE"] = multPlots3GeVHE;
  mult_plots["3GeV_timescanHB"]= multPlots3GeVHB;
  //  mult_plots["2GeV_timescan"] = multPlots2GeV;
  //  mult_plots["2GeV_timescanHE"] = multPlots2GeVHE;
  mult_plots["2GeV_timescanHB"]= multPlots2GeVHB;
  //  mult_plots["1GeV_timescan"] = multPlots1GeV;  
  //  mult_plots["1GeV_timescanHE"] = multPlots1GeVHE;
  mult_plots["1GeV_timescanHB"]= multPlots1GeVHB;
  //  mult_plots["3GeV_timescan_Jet"] = multPlots3GeV_Jet;
  //  mult_plots["3GeV_timescanHE_Jet"] = multPlots3GeVHE_Jet;
  mult_plots["3GeV_timescanHB_Jet"]= multPlots3GeVHB_Jet;
  mult_plots["0GeV_timescanHB_Jet"]= multPlots0GeVHB_Jet;
  //  mult_plots["3GeV_timescanHB_Jet_near"]= multPlots3GeVHB_Jet_near;
  mult_plots["multPlots3GeVHB_Jet_L1DRcone"]= multPlots3GeVHB_Jet_L1DRcone;
  mult_plots["multPlots0GeVHB_Jet_L1DRcone"]= multPlots0GeVHB_Jet_L1DRcone;
  //  mult_plots["2GeV_timescan"] = multPlots2GeV_Jet;
  //  mult_plots["2GeV_timescanHE_Jet"] = multPlots2GeVHE_Jet;
  mult_plots["2GeV_timescanHB_Jet"]= multPlots2GeVHB_Jet;
  //  mult_plots["1GeV_timescan_Jet"] = multPlots1GeV_Jet;
  //  mult_plots["1GeV_timescanHE_Jet"] = multPlots1GeVHE_Jet;
  mult_plots["1GeV_timescanHB_Jet"]= multPlots1GeVHB_Jet;
  // looping through all plot collections (jets, eg, tau, scalar, vector)
  for(auto iplot : plots) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.3, 1, 1));
    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad2.push_back(new TPad("pad2", "pad2", 0, 0, 1, 0.3));
    pad2.back()->SetGrid();
    pad2.back()->Draw();    
    pad1.back()->cd();
    rateHists_def[iplot.second.front()]->GetYaxis()->SetTitle("rate (Hz, unnormalized)");
    rateHists_def[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.55, 0.9 - 0.1*iplot.second.size(), 0.95, 0.93);
    for(auto hist : iplot.second) {
      rateHists_def[hist]->Draw("hist same");
      rateHists_def[hist]->GetYaxis()->SetRangeUser(1000, 100000000); // setting the range of the Y axis to show low rates
      if(includeHW) rateHists_hw[hist]->Draw("hist same");
      rateHists_new_cond[hist]->Draw("hist same");
      TString name(rateHists_def[hist]->GetName());
      TString nameHw(rateHists_hw[hist]->GetName());
      leg->AddEntry(rateHists_def[hist], name + " (current)", "L");
      if(includeHW) leg->AddEntry(rateHists_hw[hist], name + " (hw)", "L");
      leg->AddEntry(rateHists_new_cond[hist], name + " (new)", "L"); 
    }
    leg->SetBorderSize(0);
    leg->Draw();
    
    pad2.back()->cd();
    rateHistsRatio[iplot.second.front()]->Draw("hist");
    if(includeHW) rateHistsRatio[iplot.second.front()]->GetYaxis()->SetTitle("Current/HW");
    else rateHistsRatio[iplot.second.front()]->GetYaxis()->SetTitle("New/Current");
    for(auto hist : iplot.second) {
      rateHistsRatio[hist]->Draw("hist same");
    }

    if(includeHW) canvases.back()->Print(Form("plots/%sRates_hw.pdf", iplot.first.c_str()));
    else canvases.back()->Print(Form("plots/%sRates_emu.pdf", iplot.first.c_str()));
  }

  // multiplicity plot loop
  // QCD
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_QCD[iplot.second.front()]->Scale(1./multHists_QCD[iplot.second.front()]->Integral());
    multHists_QCD[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) { // associative array is list of pairs, access by first entry. Second is actual name / value to access
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV. Don't do for HB
      if (hist.substr(0,3) == "dt2" && hist.substr(9,2) != "HB" && hist.substr(hist.length()-3) != "Jet") {
	multHists_QCD[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here. Don't do for barrel plots as the x axis is later restricted to smaller region
      if (hist.substr(0,3) == "dt1" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet") {
        multHists_QCD[hist]->Rebin(rebinFactor*4);
      }
      int yMax = 0;
      yMax = multHists_QCD[hist]->GetMaximum();
      multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,1.2*yMax);
      multHists_QCD[hist]->Scale(1./multHists_QCD[hist]->Integral());
      multHists_QCD[hist]->Draw("hist same");
      TString name(multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_QCD[hist], name(6,3) + " ", "L");
      multHists_QCD[hist]->SetTitle("Multiplicity for QCD, timing scan at " + name(2,4)+", TP matched w/"+ name(9,3));
      if (name(9,3) != "Jet" ) multHists_QCD[hist]->SetTitle("Multiplicity for QCD, timing scan at " + name(2,4));
      multHists_QCD[hist]->GetYaxis()->CenterTitle(true);
      if ( name(9,2) == "HE" || name(9,2) == "HB" ){
        multHists_QCD[hist]->SetTitle("Multiplicity for QCD, timing scan at " + name(2,4) + " in " + name(9,2)+", TP matched w/" + name(11,3));
	if ( name(11,3) != "Jet" ) multHists_QCD[hist]->SetTitle("Multiplicity for QCD, timing scan at " + name(2,4) + " in " + name(9,2));
	if (hist.substr(0,3) == "dt0" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,100);
        if (hist.substr(0,3) == "dt1" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,100);
        if (hist.substr(0,3) == "dt2" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,70);
        if (hist.substr(0,3) == "dt3" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,35);
      }
      //      if ( name(11,3) == "Jet" || name(9,3) == "Jet" ) { // reduce x axis range for plots where HCAL TP has been matched to L1 Jet                     
	//	multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,25);
      //      }
      multHists_QCD[hist]->GetXaxis()->SetLabelSize(0.03);
      multHists_QCD[hist]->GetYaxis()->SetLabelSize(0.03);
      multHists_QCD[hist]->GetXaxis()->SetTitleSize(0.04);
      multHists_QCD[hist]->GetYaxis()->SetTitleSize(0.04);
      multHists_QCD[hist]->GetXaxis()->SetTitleOffset(1.2);
      multHists_QCD[hist]->GetYaxis()->SetTitleOffset(1.5);
    }
    leg->SetBorderSize(0);
    leg->Draw();
    //    canvases.back()->Print(Form("plots/%s_Mult_emu_QCD.pdf", iplot.first.c_str()));
  }

  // LLP 10000
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP10000[iplot.second.front()]->Scale(1./multHists_LLP10000[iplot.second.front()]->Integral());
    //multHists_LLP10000[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV                                
      if (hist.substr(0,3) == "dt2" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet") {
        multHists_LLP10000[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here                       
      if (hist.substr(0,3) == "dt1" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet") {
        multHists_LLP10000[hist]->Rebin(rebinFactor*4);
      }
      int yMax = 0;
      yMax = multHists_LLP10000[hist]->GetMaximum();
      multHists_LLP10000[hist]->GetYaxis()->SetRangeUser(0,1.2*yMax);
      multHists_LLP10000[hist]->Scale(1./multHists_LLP10000[hist]->Integral());
      //multHists_LLP10000[hist]->Draw("hist same");
      TString name(multHists_LLP10000[hist]->GetName());
      leg->AddEntry(multHists_LLP10000[hist], name(6,3) + " ", "L");
      multHists_LLP10000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=10m, timing scan at " + name(2,4)+", TP matched w/" + name(9,3));
      if (name(9,3) != "Jet") multHists_LLP10000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=10m, timing scan at " + name(2,4));
      multHists_LLP10000[hist]->GetYaxis()->CenterTitle(true);
      if ( name(9,2) == "HE" || name(9,2) == "HB" ){
	multHists_LLP10000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=10m, timing scan at " + name(2,4) + " in " + name(9,2)+", TP matched w/" + name(11,3)); 
	if (name(11,3) != "Jet" ) multHists_LLP10000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=10m, timing scan at " + name(2,4) + " in " + name(9,2));
        if (hist.substr(0,3) == "dt0" && name(9,2) == "HB" ) multHists_LLP10000[hist]->GetXaxis()->SetRangeUser(0,100);
	if (hist.substr(0,3) == "dt1" && name(9,2) == "HB" ) multHists_LLP10000[hist]->GetXaxis()->SetRangeUser(0,100);
	if (hist.substr(0,3) == "dt2" && name(9,2) == "HB" ) multHists_LLP10000[hist]->GetXaxis()->SetRangeUser(0,70);
	if (hist.substr(0,3) == "dt3" && name(9,2) == "HB" ) multHists_LLP10000[hist]->GetXaxis()->SetRangeUser(0,35);
      }
      //      if ( name(11,3) == "Jet" || name(9,3) == "Jet" ) { // reduce x axis range for plots where HCAL TP has been matched to L1 Jet      
      //	multHists_LLP10000[hist]->GetXaxis()->SetRangeUser(0,25);
      //      }
      multHists_LLP10000[hist]->GetXaxis()->SetLabelSize(0.03);
      multHists_LLP10000[hist]->GetYaxis()->SetLabelSize(0.03);
      multHists_LLP10000[hist]->GetXaxis()->SetTitleSize(0.04);
      multHists_LLP10000[hist]->GetYaxis()->SetTitleSize(0.04);
      multHists_LLP10000[hist]->GetXaxis()->SetTitleOffset(1.2);
      multHists_LLP10000[hist]->GetYaxis()->SetTitleOffset(1.5);
    }
    leg->SetBorderSize(0);
    leg->Draw();
    //    canvases.back()->Print(Form("plots/%s_Mult_emu_LLP10000.pdf", iplot.first.c_str()));
  }

  // LLP 1000
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP1000[iplot.second.front()]->Scale(1./multHists_LLP1000[iplot.second.front()]->Integral());
    multHists_LLP1000[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV
      if ( hist.substr(0,3) == "dt2" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet") {
        multHists_LLP1000[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here
      if ( hist.substr(0,3) == "dt1" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet" ) {
        multHists_LLP1000[hist]->Rebin(rebinFactor*4);
      }
      int yMax = 0;
      yMax = multHists_LLP1000[hist]->GetMaximum();
      multHists_LLP1000[hist]->GetYaxis()->SetRangeUser(0,1.2*yMax);
      multHists_LLP1000[hist]->Scale(1./multHists_LLP1000[hist]->Integral());
      multHists_LLP1000[hist]->Draw("hist same");
      TString name(multHists_LLP1000[hist]->GetName());
      leg->AddEntry(multHists_LLP1000[hist], name(6,3) + " ", "L");
      multHists_LLP1000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=1m, timing scan at " + name(2,4)+", TP matched w/" + name(9,3));
      if (name(9,3) != "Jet" ) multHists_LLP1000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=1m, timing scan at " + name(2,4));
      multHists_LLP1000[hist]->GetYaxis()->CenterTitle(true);
      if ( name(9,2) == "HE" || name(9,2) == "HB" ){
        multHists_LLP1000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=1m, timing scan at " + name(2,4) + " in " + name(9,2)+", TP matched w/" + name(11,3));
	if (name(11,3) != "Jet" ) multHists_LLP1000[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=1m, timing scan at " + name(2,4) + " in " + name(9,2));
        if (hist.substr(0,3) == "dt0" && name(9,2) == "HB" ) multHists_LLP1000[hist]->GetXaxis()->SetRangeUser(0,100);
	if (hist.substr(0,3) == "dt1" && name(9,2) == "HB" ) multHists_LLP1000[hist]->GetXaxis()->SetRangeUser(0,100);
        if (hist.substr(0,3) == "dt2" && name(9,2) == "HB" ) multHists_LLP1000[hist]->GetXaxis()->SetRangeUser(0,70);
        if (hist.substr(0,3) == "dt3" && name(9,2) == "HB" ) multHists_LLP1000[hist]->GetXaxis()->SetRangeUser(0,35);
      }
      //      if ( name(11,3) == "Jet" || name(9,3) == "Jet" ) { // reduce x axis range for plots where HCAL TP has been matched to L1 Jet                             
      //	multHists_LLP1000[hist]->GetXaxis()->SetRangeUser(0,25);
      //      }
      multHists_LLP1000[hist]->GetXaxis()->SetLabelSize(0.03);
      multHists_LLP1000[hist]->GetYaxis()->SetLabelSize(0.03);
      multHists_LLP1000[hist]->GetXaxis()->SetTitleSize(0.04);
      multHists_LLP1000[hist]->GetYaxis()->SetTitleSize(0.04);
      multHists_LLP1000[hist]->GetXaxis()->SetTitleOffset(1.2);
      multHists_LLP1000[hist]->GetYaxis()->SetTitleOffset(1.5);
    }
    leg->SetBorderSize(0);
    leg->Draw();
    //    canvases.back()->Print(Form("plots/%s_Mult_emu_LLP1000.pdf", iplot.first.c_str()));
  }

  // LLP 500                                                                
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP500[iplot.second.front()]->Scale(1./multHists_LLP500[iplot.second.front()]->Integral());
    multHists_LLP500[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV
      if ( hist.substr(0,3) == "dt2" && hist.substr(9,2) != "HB"  && hist.substr(hist.length()-3) != "Jet" ) {
        multHists_LLP500[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here
      if ( hist.substr(0,3) =="dt1" && hist.substr(9,2) != "HB" && hist.substr(hist.length()-3) != "Jet" ) {
        multHists_LLP500[hist]->Rebin(rebinFactor*4);
      }
      int yMax = 0;
      yMax = multHists_LLP500[hist]->GetMaximum();
      multHists_LLP500[hist]->GetYaxis()->SetRangeUser(0,1.2*yMax);
      multHists_LLP500[hist]->Scale(1./multHists_LLP500[hist]->Integral());
      multHists_LLP500[hist]->Draw("hist same");
      TString name(multHists_LLP500[hist]->GetName());
      leg->AddEntry(multHists_LLP500[hist], name(6,3) + " ", "L");
      multHists_LLP500[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=0.5m, timing scan at " + name(2,4)+", TP matched w/" + name(9,3));
      if (name(9,3)!="Jet")  multHists_LLP500[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=0.5m, timing scan at " + name(2,4));
      multHists_LLP500[hist]->GetYaxis()->CenterTitle(true);
      if ( name(9,2) == "HE" || name(9,2) == "HB" ){
        multHists_LLP500[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=0.5m, timing scan at " + name(2,4) + " in region " + name(9,2)+", TP matched w/" + name(11,3));
	if (name(11,3) != "Jet" ) multHists_LLP500[hist]->SetTitle("Multiplicity for LLP c#scale[1.2]{#tau}=0.5m, timing scan at " + name(2,4) + " in region " + name(9,2));
        if (hist.substr(0,3) == "dt0" && name(9,2) == "HB" ) multHists_LLP500[hist]->GetXaxis()->SetRangeUser(0,100);
	if (hist.substr(0,3) == "dt1" && name(9,2) == "HB" ) multHists_LLP500[hist]->GetXaxis()->SetRangeUser(0,100);
        if (hist.substr(0,3) == "dt2" && name(9,2) == "HB" ) multHists_LLP500[hist]->GetXaxis()->SetRangeUser(0,70);
        if (hist.substr(0,3) == "dt3" && name(9,2) == "HB" ) multHists_LLP500[hist]->GetXaxis()->SetRangeUser(0,35);
      }
      //      if ( name(11,3) == "Jet" || name(9,3) == "Jet" ) { // reduce x axis range for plots where HCAL TP has been matched to L1 Jet              
      //	multHists_LLP500[hist]->GetXaxis()->SetRangeUser(0,25);
      //      }
      multHists_LLP500[hist]->GetXaxis()->SetLabelSize(0.03);
      multHists_LLP500[hist]->GetYaxis()->SetLabelSize(0.03);
      multHists_LLP500[hist]->GetXaxis()->SetTitleSize(0.04);
      multHists_LLP500[hist]->GetYaxis()->SetTitleSize(0.04);
      multHists_LLP500[hist]->GetXaxis()->SetTitleOffset(1.2);
      multHists_LLP500[hist]->GetYaxis()->SetTitleOffset(1.5);
    }
    leg->SetBorderSize(0);
    leg->Draw();
    //    canvases.back()->Print(Form("plots/%s_Mult_emu_LLP500.pdf", iplot.first.c_str()));
  }

  // overlay LLP and QCD at a single energy / timing cut value
  for (auto hist : multTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_QCD[hist]->SetLineColor(kBlack);
    int RatioyMax = 0; // setting max for ratio histograms
    RatioyMax = multHists_LLP1000[hist]->GetMaximum();
    multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,1.2*RatioyMax);
    TString name (multHists_QCD[hist]->GetName());
    double x1 = 0.15; // values for legend position for ratio plots (this moves it to the top left)
    double x2 = 0.5;
    if (name(0,5) != "Ratio") {
      x1 = 0.55; // for the legend position when not on a ratio plot (this moves it to the top right)
      x2 = 0.95;
      int yMax = 0;
      yMax = multHists_QCD[hist]->GetMaximum();
      multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,6*yMax);
      if (name(11,4) == "Jet1" ) multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,2*yMax); // 3.5*
      if (name(11,4) == "Jet2" ) multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,2*yMax); // 3.3*
      if (name(11,4) == "Jet3" ) multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,2*yMax); // 3*
      if (name(11,4) == "Jet4" ) multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,1.5*yMax); // 1.5*
    }
    if (name(0,11) == "DelayedHitF" ) {
      x1 = 0.25;
      x2 = 0.6;
    }
    multHists_QCD[hist]->SetFillStyle(3005); // this is the grey shading on QCD plots
    multHists_QCD[hist]->Scale(1./multHists_QCD[hist]->Integral());
    multHists_QCD[hist]->Draw("hist pfc");
    TLegend *leg = new TLegend(x1, 0.6, x2, 0.9);
    multHists_LLP500[hist]->SetLineColor(kBlue);
    multHists_LLP500[hist]->Scale(1./multHists_LLP500[hist]->Integral());
    multHists_LLP500[hist]->Draw("hist same");
    multHists_LLP1000[hist]->SetLineColor(kGreen+1);
    multHists_LLP1000[hist]->Scale(1./multHists_LLP1000[hist]->Integral());
    multHists_LLP1000[hist]->Draw("hist same");
    multHists_LLP10000[hist]->SetLineColor(kRed); // not using LLP with ctau = 10m, very far out
    multHists_LLP10000[hist]->Scale(1./multHists_LLP10000[hist]->Integral());
    //multHists_LLP10000[hist]->Draw("hist same");
    leg->AddEntry(multHists_QCD[hist],"QCD", "F");
    leg->AddEntry(multHists_LLP500[hist],"LLP, c#scale[1.2]{#tau}=0.5m", "L");
    leg->AddEntry(multHists_LLP1000[hist], "LLP, c#scale[1.2]{#tau}=1m", "L");
    //    leg->AddEntry(multHists_LLP10000[hist], "LLP, c#scale[1.2]{#tau}=0.5m, noPU", "L");
    multHists_QCD[hist]->GetYaxis()->CenterTitle(true);
    multHists_QCD[hist]->SetTitle("Multiplicity at " + name(2,4) + " and " + name(6,3)+", TP matched w/" + name (9,3)); // general name "mult at 3 GeV and 1 ns, TP matched w/jet"
    if (name(9,3) != "Jet" ) multHists_QCD[hist]->SetTitle("Multiplicity at " + name(2,4) + " and " + name(6,3)); // "mult at 3 GeV and 1 ns"
    if (name(0,6) == "DepthV" ) multHists_QCD[hist]->SetTitle("Depth Variable Multiplicity");
    if (name(0,11) == "DelayedHitF" ) {
      multHists_QCD[hist]->SetTitle("Delayed Hit Fraction " + name(18,7));
      multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,1);
    }
    if ( name(9,2) == "HE" || name(9,2) == "HB" ){ // setting range and name of histograms for HCAL barrel and endcap regions
      multHists_QCD[hist]->SetTitle("Multiplicity at " + name(2,4) + " and " + name(6,3) + " in " + name(9,2)+", TP matched w/" + name(11,3));
      if (name(11,3) != "Jet" ) multHists_QCD[hist]->SetTitle("Multiplicity at " + name(2,4) + " and " + name(6,3) + " in " + name(9,2));
      if (name(9,5) == "HBSin" || name(9,5) == "HBDou" || name(9,5) == "HBTri" ) multHists_QCD[hist]->SetTitle(name(11,6) + " Jet Hit Multiplicity at " + name(2,4) + " and " + name(6,3));
      if (name(9,5) == "HBQua" ) multHists_QCD[hist]->SetTitle(name(11,4) + " Jet Hit Multiplicity at " + name(2,4) + " and " + name(6,3));
      if (hist.substr(0,3) == "dt0" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,100);
      if (hist.substr(0,3) == "dt1" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,100);
      if (hist.substr(0,3) == "dt2" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,70);
      if (hist.substr(0,3) == "dt3" && name(9,2) == "HB" ) multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,35);
      if ( (name(11,4) == "Jet1") || (name(11,4) == "Jet2") || (name(11,4) == "Jet3") || (name(11,4) == "Jet4") ) {
        multHists_QCD[hist]->SetTitle("Multiplicity at " + name(2,4) + " and " + name(6,3) + ", TP in DR cone of L1 Jet #" + name(14,1)); // "mult at 3 GeV and 1 ns, TP in DR cone of L1 jet 1" 
	multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,25);
      }
    }
    if ( hist == "centralTiming" ) multHists_QCD[hist]->SetTitle("Time of arrival - TOF (central barrel iEta)");
    if (hist.substr(6,4) == "calo"){ // setting histogram name for calo tower ieta energy scan
      if (hist.substr(10,2) == "T1") multHists_QCD[hist]->SetTitle("Multiplicity of cells above " + name(2,4) + ", CaloTower 1 (iEta 1-4)");
      if (hist.substr(10,2) == "T2") multHists_QCD[hist]->SetTitle("Multiplicity of cells above " + name(2,4) + ", CaloTower 2 (iEta 5-8)");
      if (hist.substr(10,2) == "T3") multHists_QCD[hist]->SetTitle("Multiplicity of cells above " + name(2,4) + ", CaloTower 3 (iEta 9-12)");
      if (hist.substr(10,2) == "T4") multHists_QCD[hist]->SetTitle("Multiplicity of cells above " + name(2,4) + ", CaloTower 4 (iEta 13-16)");
      if (hist.substr(0,3) == "dt1") multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,50);
      if (hist.substr(0,3) == "dt2") multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,25);
      if (hist.substr(0,3) == "dt3") multHists_QCD[hist]->GetXaxis()->SetRangeUser(0,15);
    }
    if (name(0,5) == "Ratio" ) {
      multHists_QCD[hist]->SetTitle("Ratio of First 2 HCAL Layers to E_{T} " + name(11,7));
    }
    multHists_QCD[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_QCD[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay.pdf", hist.substr(2).c_str()));
  }

  for (auto hist : EDepthTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    //    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    TH1D *energy_profile_QCD = energy_depth_QCD[hist]->ProfileX("QCD");
    energy_profile_QCD->Draw("hist");
    double y1 = 0;
    double y2 = 0;
    if (hist.substr(0,6) == "Energy") {
      y1 = 0.65;
      y2 = 0.93;
      energy_profile_QCD->GetYaxis()->SetRangeUser(0,1);
      energy_profile_QCD->GetYaxis()->SetTitle("Average Energy Fraction");
    }
    if (hist.substr(0,6) == "Timing") {
      y1 = 0.15;
      y2 = 0.35;
      energy_profile_QCD->GetYaxis()->SetRangeUser(0,10);
      energy_profile_QCD->GetYaxis()->SetTitle("Average Timing (ns)");
    }
    energy_profile_QCD->GetYaxis()->CenterTitle(true);
    TLegend *leg = new TLegend(0.55, y1, 0.95, y2);
    energy_profile_QCD->SetLineColor(kBlack);
    energy_profile_QCD->SetLineWidth(2);
    energy_profile_QCD->SetFillStyle(3005); // this is the grey shading on QCD plots
    energy_profile_QCD->SetFillColor(kGray+2);
    energy_profile_QCD->SetDirectory(0);
    energy_profile_QCD->Draw("ehist same");
    TH1D *energy_profile_LLP500 = energy_depth_LLP500[hist]->ProfileX("LLP500");
    energy_profile_LLP500->SetLineColor(kBlue);
    energy_profile_LLP500->SetDirectory(0);
    energy_profile_LLP500->SetLineWidth(2);
    energy_profile_LLP500->Draw("ehist same");
    TH1D *energy_profile_LLP1000 = energy_depth_LLP1000[hist]->ProfileX("LLP1000");
    energy_profile_LLP1000->SetLineColor(kGreen+1);
    energy_profile_LLP1000->SetDirectory(0);
    energy_profile_LLP1000->SetLineWidth(2);
    energy_profile_LLP1000->Draw("ehist same");
    TH1D *energy_profile_LLP10000 = energy_depth_LLP10000[hist]->ProfileX("LLP10000");
    energy_profile_LLP10000->SetLineColor(kRed);
    energy_profile_LLP10000->SetDirectory(0);
    energy_profile_LLP10000->SetLineWidth(2);
    //    energy_profile_LLP10000->Draw("ehist same");
    leg->AddEntry(energy_profile_QCD,"QCD","F");
    leg->AddEntry(energy_profile_LLP500,"LLP, c#scale[1.2]{#tau}=0.5m", "L");
    leg->AddEntry(energy_profile_LLP1000, "LLP, c#scale[1.2]{#tau}=1m", "L");
    //    leg->AddEntry(energy_profile_LLP10000, "LLP, c#scale[1.2]{#tau}=0.5m, noPU", "L");
    energy_profile_QCD->GetXaxis()->SetLabelSize(0.03);
    energy_profile_QCD->GetYaxis()->SetLabelSize(0.03);
    energy_profile_QCD->GetXaxis()->SetTitleSize(0.04);
    energy_profile_QCD->GetYaxis()->SetTitleSize(0.04);
    energy_profile_QCD->GetXaxis()->SetTitleOffset(1.2);
    energy_profile_QCD->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw("same");
    canvases.back()->Print(Form("plots/%s_LLPQCD_overlay.pdf", hist.c_str()));
  }

  for (auto hist : EDepthTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    TH1D *energy_profile_QCD = energy_depth_QCD[hist]->ProfileX("QCD");
    if (hist.substr(0,6) =="Energy") {
      energy_profile_QCD->SetMaximum(1);
    }
    else energy_profile_QCD->GetYaxis()->SetRangeUser(0,12);
    energy_profile_QCD->SetLineColor(kBlack);
    energy_profile_QCD->Draw("ehist");
    energy_profile_QCD->SetTitle("TP Energy Fraction vs. Depth for QCD");
    //    canvases.back()->Print(Form("plots/%s_QCD.pdf", hist.c_str()));
  }
  for (auto hist : EDepthTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    TH1D *energy_profile_LLP500 = energy_depth_LLP500[hist]->ProfileX("LLP500");
    if (hist.substr(0,6) =="Energy") {
      energy_profile_LLP500->SetMaximum(1);
    }
    else energy_profile_LLP500->GetYaxis()->SetRangeUser(0,12);
    energy_profile_LLP500->SetLineColor(kBlue);
    energy_profile_LLP500->Draw("ehist");
    energy_profile_LLP500->SetTitle("TP Energy Fraction vs. Depth for LLP c#scale[1.2]{#tau}=0.5m");
    //    canvases.back()->Print(Form("plots/%s_LLP500.pdf", hist.c_str()));
  }
  for (auto hist : EDepthTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    TH1D *energy_profile_LLP1000 = energy_depth_LLP1000[hist]->ProfileX("LLP1000");
    if (hist.substr(0,6) =="Energy") {
      energy_profile_LLP1000->SetMaximum(1);
    }
    else energy_profile_LLP1000->GetYaxis()->SetRangeUser(0,12);
    energy_profile_LLP1000->SetLineColor(kGreen+1);
    energy_profile_LLP1000->Draw("ehist");
    energy_profile_LLP1000->SetTitle("TP Energy Fraction vs. Depth for LLP c#scale[1.2]{#tau}=1m");
    //    canvases.back()->Print(Form("plots/%s_LLP1000.pdf", hist.c_str()));
  }
  for (auto hist : EDepthTypes ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    TH1D *energy_profile_LLP10000 = energy_depth_LLP10000[hist]->ProfileX("LLP10000");
    if (hist.substr(0,6) == "Energy") {
      energy_profile_LLP10000->SetMaximum(1);
    }
    else energy_profile_LLP10000->GetYaxis()->SetRangeUser(0,12);
    energy_profile_LLP10000->SetLineColor(kRed);
    //    energy_profile_LLP10000->Draw("ehist");
    energy_profile_LLP10000->SetTitle("TP Energy Fraction vs. Depth for LLP c#scale[1.2]{#tau}=10m");
    //    canvases.back()->Print(Form("plots/%s_LLP10000.pdf", hist.c_str()));
  }

  std::cout << "Global multiplicity: single jet rate = " << SJet60GeV << " and htSum rate = " << htSum350GeV << std::endl;
  std::cout << "Jet matched multiplicity: single jet rate = " << SJet60GeV_l << " and htSum rate = " << htSum350GeV_l << std::endl;

  Double_t EffPl500_120[6], EffPl500_350[6], EffPl500Mh350_120[6], EffPl500Mh350_350[6], EffQCD_120[6], EffQCD_350[6], singleJetRate[6], quadJetRate[6], htSum120Rate[6], htSum350Rate[6];
  Double_t EffPl500Global[6], EffPl500Mh350Global[6], EffQCDGlobal[6], singleJetRateGlobal[6], quadJetRateGlobal[6], htSumRate120Global[6], htSumRate350Global[6];
  double EffPl500Mh1000_htSum120_Global, EffPl500Mh1000_htSum350_Global, EffPl500Mh1000_htSum430_Global, EffQCD_htSum120_Global, EffQCD_htSum350_Global, EffQCD_htSum430_Global, Original_htSumRate120, Original_htSumRate350, Original_htSumRate430;
  double EffPl500Mh350_htSum120_Global, EffPl500Mh350_htSum350_Global, EffPl500Mh350_htSum430_Global;

  // L1 JET MATCHED
  // neutrino gun rate at 60 GeV for single Jet, in kHz
  singleJetRate[0] = 184;
  singleJetRate[1] = 28;
  singleJetRate[2] = 7;
  singleJetRate[3] = 0;
  singleJetRate[4] = 0;
  singleJetRate[5] = 0;
  // neutrino gun rate at 60 GeV for quad Jet
  quadJetRate[0] = 14; // in kHz
  quadJetRate[1] = 7;
  quadJetRate[2] = 0;
  quadJetRate[3] = 0;
  quadJetRate[4] = 0;
  quadJetRate[5] = 0;
  // neutrino gun rate for htSum at 120 GeV
  htSum120Rate[0] = 312;
  htSum120Rate[1] = 28;
  htSum120Rate[2] = 7;
  htSum120Rate[3] = 0;
  htSum120Rate[4] = 0;
  htSum120Rate[5] = 0;
  // neutrino gun rate for htSum at 350 GeV
  htSum350Rate[0] = 78;
  htSum350Rate[1] = 21;
  htSum350Rate[2] = 0;
  htSum350Rate[3] = 0;
  htSum350Rate[4] = 0;
  htSum350Rate[5] = 0;
  // signal efficiency for pl 500, mh1000 GeV. Efficiency = % of events passing timing mult cuts AND ht > 120
  EffPl500_120[0] = 0.9075;
  EffPl500_120[1] = 0.8575;
  EffPl500_120[2] = 0.8065;
  EffPl500_120[3] = 0.7535;
  EffPl500_120[4] = 0.6855;
  EffPl500_120[5] = 0.629;
  // signal efficiency for pl 500, mh1000 GeV. Efficiency = % of events passing timing mult cuts AND ht > 350 
  EffPl500_350[0] = 0.9065;
  EffPl500_350[1] = 0.857;
  EffPl500_350[2] = 0.8065;
  EffPl500_350[3] = 0.7535;
  EffPl500_350[4] = 0.6855;
  EffPl500_350[5] = 0.629;
  // signal efficiency for pl 500, mh350 GeV. Efficiency = % of events passing timing mult cuts AND ht > 120
  EffPl500Mh350_120[0] = 0.5225;
  EffPl500Mh350_120[1] = 0.389;
  EffPl500Mh350_120[2] = 0.285;
  EffPl500Mh350_120[3] = 0.219;
  EffPl500Mh350_120[4] = 0.1665;
  EffPl500Mh350_120[5] = 0.1255;
  // signal efficiency for pl 500, mh350 GeV. Efficiency = % of events passing timing mult cuts AND ht > 350
  EffPl500Mh350_350[0] = 0.522;
  EffPl500Mh350_350[1] = 0.3695;
  EffPl500Mh350_350[2] = 0.2745;
  EffPl500Mh350_350[3] = 0.2135;
  EffPl500Mh350_350[4] = 0.162;
  EffPl500Mh350_350[5] = 0.1235;
  // background efficiency for QCD. Efficiency = % of events passing timing mult cuts AND ht > 120
  EffQCD_120[0] = 0.275; //mult3GeV3nsHB_Jets > 1
  EffQCD_120[1] = 0.202; //mult3GeV3nsHB_Jets > 2
  EffQCD_120[2] = 0.1505; //mult3GeV3nsHB_Jets > 3 
  EffQCD_120[3] = 0.1165; //mult3GeV3nsHB_Jets > 4
  EffQCD_120[4] = 0.0945; //mult3GeV3nsHB_Jets > 5
  EffQCD_120[5] = 0.0755; //mult3GeV3nsHB_Jets > 6 
  // background efficiency for QCD. Efficiency = % of events passing timing mult cuts AND ht > 350
  EffQCD_350[0] = 0.261; //mult3GeV3nsHB_Jets > 1
  EffQCD_350[1] = 0.199; //mult3GeV3nsHB_Jets > 2
  EffQCD_350[2] = 0.1505; //mult3GeV3nsHB_Jets > 3
  EffQCD_350[3] = 0.1165; //mult3GeV3nsHB_Jets > 4 
  EffQCD_350[4] = 0.0945; //mult3GeV3nsHB_Jets > 5
  EffQCD_350[5] = 0.0755; //mult3GeV3nsHB_Jets > 6  

  // GLOBAL
  // neutrino gun rate at 60 GeV for single Jet
  singleJetRateGlobal[0] = 3782; // in kHz        
  singleJetRateGlobal[1] = 1270;
  singleJetRateGlobal[2] = 269;
  singleJetRateGlobal[3] = 56;
  singleJetRateGlobal[4] = 14;
  singleJetRateGlobal[5] = 0;
  // neutrino gun rate at 60 GeV for quad Jet
  quadJetRateGlobal[0] = 283; // in kHz
  quadJetRateGlobal[1] = 141;
  quadJetRateGlobal[2] = 28;
  quadJetRateGlobal[3] = 7;
  quadJetRateGlobal[4] = 7;
  quadJetRateGlobal[5] = 0;
  // neutrino gun rate for htSum at 120 GeV
  htSumRate120Global[0] = 5607;
  htSumRate120Global[1] = 1760;
  htSumRate120Global[2] = 347;
  htSumRate120Global[3] = 85;
  htSumRate120Global[4] = 14;
  htSumRate120Global[5] = 0;
  // neutrino gun rate for htSum at 350 GeV    
  htSumRate350Global[0] = 1547;
  htSumRate350Global[1] = 574;
  htSumRate350Global[2] = 127;
  htSumRate350Global[3] = 42;
  htSumRate350Global[4] = 7;
  htSumRate350Global[5] = 0;
  // signal efficiency for pl 500, mh1000
  EffPl500Global[0] = 0.96;
  EffPl500Global[1] = 0.9205;
  EffPl500Global[2] = 0.8765;
  EffPl500Global[3] = 0.8365;
  EffPl500Global[4] = 0.781;
  EffPl500Global[5] = 0.726;
  // signal efficiency for pl500, mh350
  EffPl500Mh350Global[0] = 0.7305;
  EffPl500Mh350Global[1] = 0.583;
  EffPl500Mh350Global[2] = 0.448;
  EffPl500Mh350Global[3] = 0.3365;
  EffPl500Mh350Global[4] = 0.2545;
  EffPl500Mh350Global[5] = 0.191;
  // background efficiency for QCD
  EffQCDGlobal[0] = 0.505; //mult3GeV3nsHB > 1
  EffQCDGlobal[1] = 0.32; //mult3GeV3nsHB > 2
  EffQCDGlobal[2] = 0.2175; //mult3GeV3nsHB > 3
  EffQCDGlobal[3] = 0.159; //mult3GeV3nsHB > 4  
  EffQCDGlobal[4] = 0.1225; //mult3GeV3nsHB > 5
  EffQCDGlobal[5] = 0.1; //mult3GeV3nsHB > 6

  // comparison points for htSum rates and efficiencies
  EffPl500Mh1000_htSum120_Global = 0.9995;
  EffPl500Mh1000_htSum350_Global = 0.9975;
  EffPl500Mh1000_htSum430_Global = 0.996;
  EffPl500Mh350_htSum120_Global = 0.9955;
  EffPl500Mh350_htSum350_Global = 0.849;
  EffPl500Mh350_htSum430_Global = 0.7085;
  EffQCD_htSum120_Global = 0.965;
  EffQCD_htSum350_Global = 0.656;
  EffQCD_htSum430_Global = 0.5285;
  Original_htSumRate120 = 24117;
  Original_htSumRate350 = 4776;
  Original_htSumRate430 = 1795;

  // ************* Efficiency vs. Rate for LLP mh=1000 GeV, ct=500 mm ************
  TGraph *gr_sing_LLP = new TGraph (6, EffPl500_120, singleJetRate);
  TGraph *gr_quad_LLP = new TGraph (6, EffPl500_120, quadJetRate);
  TGraph *gr_350_LLP = new TGraph (6, EffPl500_350, htSum350Rate);
  TGraph *gr_120_LLP = new TGraph (6, EffPl500_120, htSum120Rate);
  TGraph *gr_global_sing_LLP = new TGraph (6, EffPl500Global, singleJetRateGlobal);
  TGraph *gr_global_quad_LLP = new TGraph (6, EffPl500Global, quadJetRateGlobal);
  TGraph *gr_global_350_LLP = new TGraph (6, EffPl500Global, htSumRate350Global);
  TGraph *gr_global_120_LLP = new TGraph (6, EffPl500Global, htSumRate120Global);
  // ************** Efficiency vs. Rate for LLP mh=350 GeV, ct=500 mm *************
  TGraph *gr_sing_LLPMh350 = new TGraph (6, EffPl500Mh350_120, singleJetRate);
  TGraph *gr_quad_LLPMh350 = new TGraph (6, EffPl500Mh350_120, quadJetRate);
  TGraph *gr_350_LLPMh350 = new TGraph (6, EffPl500Mh350_350, htSum350Rate);
  TGraph *gr_120_LLPMh350 = new TGraph (6, EffPl500Mh350_120, htSum120Rate);
  TGraph *gr_global_sing_LLPMh350 = new TGraph (6, EffPl500Mh350Global, singleJetRateGlobal);
  TGraph *gr_global_quad_LLPMh350 = new TGraph (6, EffPl500Mh350Global, quadJetRateGlobal);
  TGraph *gr_global_350_LLPMh350 = new TGraph (6, EffPl500Mh350Global, htSumRate350Global);
  TGraph *gr_global_120_LLPMh350 = new TGraph (6, EffPl500Mh350Global, htSumRate120Global);
  // ************** Markers for comparison with current benchmark performance ******
  TMarker *m_QCD_htSum120 = new TMarker (EffQCD_htSum120_Global, Original_htSumRate120, 21);
  TMarker *m_QCD_htSum350 = new TMarker (EffQCD_htSum350_Global, Original_htSumRate350, 21);
  TMarker *m_QCD_htSum430 = new TMarker (EffQCD_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLP_htSum120 = new TMarker (EffPl500Mh1000_htSum120_Global, Original_htSumRate120, 21);
  TMarker *m_LLP_htSum350 = new TMarker (EffPl500Mh1000_htSum350_Global, Original_htSumRate350, 21);
  TMarker *m_LLP_htSum430 = new TMarker (EffPl500Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350_htSum120 = new TMarker (EffPl500Mh350_htSum120_Global, Original_htSumRate120, 21);
  TMarker *m_LLPMh350_htSum350 = new TMarker (EffPl500Mh350_htSum350_Global, Original_htSumRate350, 21);
  TMarker *m_LLPMh350_htSum430 = new TMarker (EffPl500Mh350_htSum430_Global, Original_htSumRate430, 21);

  // Jet rate, LLP mh=1000 efficiency regional
  TCanvas *c1_sing_quad = new TCanvas("c1_sing_quad","Graph Draw Options",200,10,600,400);
  gr_sing_LLP->SetLineColor(4); // blue
  gr_sing_LLP->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP->Draw("AC*");
  gr_sing_LLP->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP->SetLineColor(2); // red
  gr_quad_LLP->Draw("C*");
  auto legend = new TLegend(0.15,0.7,0.5,0.9);
  legend->AddEntry(gr_sing_LLP,"Single Jet Rate at 60 GeV");
  legend->AddEntry(gr_quad_LLP,"Quad Jet Rate at 60 GeV");
  legend->Draw();
  gr_sing_LLP->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh1000.pdf");

  // Jet rate, LLP mh=350 efficiency regional
  TCanvas *c1_sing_quadMh350 = new TCanvas("c1_sing_quadMh350","Graph Draw Options",200,10,600,400);
  gr_sing_LLPMh350->SetLineColor(4); // blue
  gr_sing_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLPMh350->Draw("AC*");
  gr_sing_LLPMh350->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLPMh350->SetLineColor(2); // red 
  gr_quad_LLPMh350->Draw("C*");
  auto legendMh350 = new TLegend(0.15,0.7,0.5,0.9);
  legendMh350->AddEntry(gr_sing_LLPMh350,"Single Jet Rate at 60 GeV");
  legendMh350->AddEntry(gr_quad_LLPMh350,"Quad Jet Rate at 60 GeV");
  legendMh350->Draw();
  gr_sing_LLPMh350->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quadMh350->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh350.pdf");

  // htSum rate, LLP mh=1000 efficiency, regional
  TCanvas *c1_350_120 = new TCanvas("c1_350_120","Graph Draw Options",200,10,600,400);
  gr_120_LLP->SetLineColor(4); // blue 
  gr_120_LLP->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP->Draw("AL*");
  gr_120_LLP->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_350_LLP->SetLineColor(2); // red    
  //  gr_350_LLP->Draw("L*");
  m_LLP_htSum120->SetMarkerStyle(21);
  m_LLP_htSum120->SetMarkerColor(4);
  m_LLP_htSum350->SetMarkerStyle(21);
  m_LLP_htSum350->SetMarkerColor(2);
  m_LLP_htSum430->SetMarkerStyle(21);
  m_LLP_htSum430->SetMarkerColor(3);
  auto legend_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend_350_120->AddEntry(gr_120_LLP,"H_{T}>120 GeV, with timing cuts");
  //  legend_350_120->AddEntry(gr_350_LLP,"H_{T}>350 GeV, with timing cuts");
  //  legend_350_120->AddEntry(m_LLP_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend_350_120->AddEntry(m_LLP_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_350_120->Draw();
  gr_120_LLP->GetXaxis()->SetLimits(0.,1.);
  //  m_LLP_htSum350->Draw();
  m_LLP_htSum430->Draw();
  gr_120_LLP->GetHistogram()->SetMinimum(1.);
  gr_120_LLP->GetHistogram()->SetMaximum(10000.);
  c1_350_120->SetLogy(); 
  c1_350_120->SetGrid();
  c1_350_120->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh1000.pdf");

  // htSum rate, LLP mh=350 efficiency, regional
  TCanvas *c1_350_120Mh350 = new TCanvas("c1_350_120Mh350","Graph Draw Options",200,10,600,400);
  gr_120_LLPMh350->SetLineColor(4); // blue
  gr_120_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_120_LLPMh350->Draw("AL*");
  gr_120_LLPMh350->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_350_LLPMh350->SetLineColor(2); // red
  //  gr_350_LLPMh350->Draw("L*");
  m_LLPMh350_htSum120->SetMarkerStyle(21);
  m_LLPMh350_htSum120->SetMarkerColor(4);
  m_LLPMh350_htSum350->SetMarkerStyle(21);
  m_LLPMh350_htSum350->SetMarkerColor(2);
  m_LLPMh350_htSum430->SetMarkerStyle(21);
  m_LLPMh350_htSum430->SetMarkerColor(3);
  auto legend_350_120Mh350 = new TLegend(0.15,0.7,0.5,0.9);
  legend_350_120Mh350->AddEntry(gr_120_LLPMh350,"H_{T}>120 GeV, with timing cuts");
  //  legend_350_120Mh350->AddEntry(gr_350_LLPMh350,"H_{T}>350 GeV, with timing cuts");
  //  legend_350_120Mh350->AddEntry(m_LLPMh350_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend_350_120Mh350->AddEntry(m_LLPMh350_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_350_120Mh350->Draw();
  gr_120_LLPMh350->GetXaxis()->SetLimits(0.,1.);
  //  m_LLPMh350_htSum350->Draw();
  m_LLPMh350_htSum430->Draw();
  gr_120_LLPMh350->GetHistogram()->SetMinimum(1.);
  gr_120_LLPMh350->GetHistogram()->SetMaximum(10000.);
  c1_350_120Mh350->SetGrid();
  c1_350_120Mh350->SetLogy();
  c1_350_120Mh350->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh350.pdf");

  // Jet rate, LLP mh=1000 efficiency, global
  TCanvas *c1_global_sing_quad = new TCanvas("c1_global_sing_quad","Graph Draw Options",200,10,600,400);
  gr_global_sing_LLP->SetLineColor(4); // blue                                            
  gr_global_sing_LLP->GetHistogram()->SetMinimum(-5.);
  gr_global_quad_LLP->GetHistogram()->SetMinimum(-5.);
  gr_global_sing_LLP->Draw("AC*");
  gr_global_sing_LLP->SetTitle("Jet Rate vs. Signal Efficiency for Global Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_quad_LLP->SetLineColor(2); // red        
  gr_global_quad_LLP->Draw("C*");
  auto legend_global_sing_quad = new TLegend(0.15,0.7,0.5,0.9);
  legend_global_sing_quad->AddEntry(gr_global_sing_LLP,"Single Jet Rate at 60 GeV");
  legend_global_sing_quad->AddEntry(gr_global_quad_LLP,"Quad Jet Rate at 60 GeV");
  legend_global_sing_quad->Draw();
  gr_global_sing_LLP->GetXaxis()->SetLimits(0.,1.);
  c1_global_sing_quad->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh1000_Global.pdf");

  // Jet rate, LLP mh=350 efficiency, global   
  TCanvas *c1_global_sing_quadMh350 = new TCanvas("c1_global_sing_quadMh350","Graph Draw Options",200,10,600,400);
  gr_global_sing_LLPMh350->SetLineColor(4); // blue  
  gr_global_sing_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_global_quad_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_global_sing_LLPMh350->Draw("AC*");
  gr_global_sing_LLPMh350->SetTitle("Jet Rate vs. Signal Efficiency for Global Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_quad_LLPMh350->SetLineColor(2); // red 
  gr_global_quad_LLPMh350->Draw("C*");
  auto legend_global_sing_quadMh350 = new TLegend(0.15,0.7,0.5,0.9);
  legend_global_sing_quadMh350->AddEntry(gr_global_sing_LLPMh350,"Single Jet Rate at 60 GeV");
  legend_global_sing_quadMh350->AddEntry(gr_global_quad_LLPMh350,"Quad Jet Rate at 60 GeV");
  legend_global_sing_quadMh350->Draw();
  gr_global_sing_LLPMh350->GetXaxis()->SetLimits(0.,1.);
  c1_global_sing_quadMh350->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh350_Global.pdf");

  // htSum rate, LLP mh=1000 efficiency, global   
  TCanvas *c1_global_350_120 = new TCanvas("c1_global_350_120","Graph Draw Options",200,10,600,400);
  gr_global_120_LLP->SetLineColor(4); // blue
  gr_global_120_LLP->GetHistogram()->SetMinimum(-5.);
  gr_global_350_LLP->GetHistogram()->SetMinimum(-5.);
  gr_global_120_LLP->Draw("AC*");
  gr_global_120_LLP->SetTitle("htSum Rate vs. Signal Efficiency for Global Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_350_LLP->SetLineColor(2); // red              
  //  gr_global_350_LLP->Draw("C*");
  m_LLP_htSum120->SetMarkerStyle(21);
  m_LLP_htSum120->SetMarkerColor(4);
  m_LLP_htSum350->SetMarkerStyle(21);
  m_LLP_htSum350->SetMarkerColor(2);
  m_LLP_htSum430->SetMarkerStyle(21);
  m_LLP_htSum430->SetMarkerColor(3);
  m_LLP_htSum120->Draw();
  //  m_LLP_htSum350->Draw();
  m_LLP_htSum430->Draw();
  auto legend_global_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend_global_350_120->AddEntry(gr_global_120_LLP,"H_{T}>120 GeV, with timing cuts");
  //  legend_global_350_120->AddEntry(gr_global_350_LLP,"H_{T}>350 GeV, with timing cuts");
  //  legend_global_350_120->AddEntry(m_LLP_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend_global_350_120->AddEntry(m_LLP_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_global_350_120->Draw();
  gr_global_120_LLP->GetXaxis()->SetLimits(0.,1.);
  c1_global_350_120->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh1000_Global.pdf");

  // htSum rate, LLP mh=350 efficiency, global
  TCanvas *c1_global_350_120Mh350 = new TCanvas("c1_global_350_120Mh350","Graph Draw Options",200,10,600,400);
  gr_global_120_LLPMh350->SetLineColor(4); // blue
  gr_global_120_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_global_350_LLPMh350->GetHistogram()->SetMinimum(-5.);
  gr_global_120_LLPMh350->Draw("AC*");
  gr_global_120_LLPMh350->SetTitle("htSum Rate vs. Signal Efficiency for Global Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_350_LLPMh350->SetLineColor(2); // red
  //  gr_global_350_LLPMh350->Draw("C*");
  m_LLPMh350_htSum120->SetMarkerStyle(21);
  m_LLPMh350_htSum120->SetMarkerColor(4);
  m_LLPMh350_htSum350->SetMarkerStyle(21);
  m_LLPMh350_htSum350->SetMarkerColor(2);
  m_LLPMh350_htSum430->SetMarkerStyle(21);
  m_LLPMh350_htSum430->SetMarkerColor(3);
  m_LLPMh350_htSum120->Draw();
  //  m_LLPMh350_htSum350->Draw();
  m_LLPMh350_htSum430->Draw();
  auto legend_global_350_120Mh350 = new TLegend(0.15,0.7,0.5,0.9);
  legend_global_350_120Mh350->AddEntry(gr_global_120_LLPMh350,"H_{T}>120 GeV, with timing cuts");
  //  legend_global_350_120Mh350->AddEntry(gr_global_350_LLPMh350,"H_{T}>350 GeV, with timing cuts");
  //  legend_global_350_120Mh350->AddEntry(m_LLPMh350_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend_global_350_120Mh350->AddEntry(m_LLPMh350_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_global_350_120Mh350->Draw();
  gr_global_120_LLPMh350->GetXaxis()->SetLimits(0.,1.);
  c1_global_350_120Mh350->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh350_Global.pdf");

  // ************ Background (QCD) Efficiency ***********************
  TGraph *gr_sing_QCD = new TGraph (6, EffQCD_120, singleJetRate);
  TGraph *gr_quad_QCD = new TGraph (6, EffQCD_120, quadJetRate);
  TGraph *gr_350_QCD = new TGraph (6, EffQCD_350, htSum350Rate);
  TGraph *gr_120_QCD = new TGraph (6, EffQCD_120, htSum120Rate);
  TGraph *gr_global_sing_QCD = new TGraph (6, EffQCDGlobal, singleJetRateGlobal);
  TGraph *gr_global_quad_QCD = new TGraph (6, EffQCDGlobal, quadJetRateGlobal);
  TGraph *gr_global_350_QCD = new TGraph (6, EffQCDGlobal, htSumRate350Global);
  TGraph *gr_global_120_QCD = new TGraph (6, EffQCDGlobal, htSumRate120Global);

  TCanvas *c2_sing_quad = new TCanvas("c2_sing_quad","Graph Draw Options",200,10,600,400);
  gr_sing_QCD->SetLineColor(4); // blue  
  gr_sing_QCD->GetHistogram()->SetMinimum(-5.);
  gr_quad_QCD->GetHistogram()->SetMinimum(-5.);
  gr_sing_QCD->Draw("AC*");
  gr_sing_QCD->SetTitle("Jet Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_QCD->SetLineColor(2); // red  
  gr_quad_QCD->Draw("C*");
  auto legend2 = new TLegend(0.55,0.7,0.9,0.9);
  legend2->AddEntry(gr_sing_QCD,"Single Jet Rate at 60 GeV");
  legend2->AddEntry(gr_quad_QCD,"Quad Jet Rate at 60 GeV");
  legend2->Draw();
  gr_sing_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_sing_quad->SaveAs("plots/NuGun_JetRates_vs_BackgroundEff.pdf");

  TCanvas *c2_350_120 = new TCanvas("c2_350_120","Graph Draw Options",200,10,600,400);
  gr_120_QCD->SetLineColor(4); // blue  
  gr_120_QCD->GetHistogram()->SetMinimum(-5.);
  gr_350_QCD->GetHistogram()->SetMinimum(-5.);
  gr_120_QCD->Draw("AL*");
  gr_120_QCD->SetTitle("htSum Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_350_QCD->SetLineColor(2); // red
  //  gr_350_QCD->Draw("L*");
  m_QCD_htSum120->SetMarkerStyle(21);
  m_QCD_htSum120->SetMarkerColor(4);
  m_QCD_htSum350->SetMarkerStyle(21);
  m_QCD_htSum350->SetMarkerColor(2);
  m_QCD_htSum430->SetMarkerStyle(21);
  m_QCD_htSum430->SetMarkerColor(3);
  m_QCD_htSum120->Draw();
  //  m_QCD_htSum350->Draw();
  m_QCD_htSum430->Draw();
  gr_120_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_350_120->SetGrid();
  auto legend2_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend2_350_120->AddEntry(gr_120_QCD,"H_{T}>120 GeV, with timing cuts");
  //  legend2_350_120->AddEntry(gr_350_QCD,"H_{T}>350 GeV, with timing cuts");
  //  legend2_350_120->AddEntry(m_QCD_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend2_350_120->AddEntry(m_QCD_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend2_350_120->Draw();
  c2_350_120->SetLogy();
  gr_120_QCD->GetHistogram()->SetMinimum(1.);
  gr_120_QCD->GetHistogram()->SetMaximum(10000.);
  c2_350_120->SaveAs("plots/NuGun_htSumRates_vs_BackgroundEff.pdf");

  TCanvas *c2_global_sing_quad = new TCanvas("c2_global_sing_quad","Graph Draw Options",200,10,600,400);
  gr_global_sing_QCD->SetLineColor(4); // blue   
  gr_global_sing_QCD->GetHistogram()->SetMinimum(-5.);
  gr_global_quad_QCD->GetHistogram()->SetMinimum(-5.);
  gr_global_sing_QCD->Draw("AC*");
  gr_global_sing_QCD->SetTitle("Jet Rate vs. Background Efficiency for Global Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_quad_QCD->SetLineColor(2); // red    
  gr_global_quad_QCD->Draw("C*");
  auto legend2_global_sing_quad = new TLegend(0.15,0.7,0.5,0.9);
  legend2_global_sing_quad->AddEntry(gr_global_sing_QCD,"Single Jet Rate at 60 GeV");
  legend2_global_sing_quad->AddEntry(gr_global_quad_QCD,"Quad Jet Rate at 60 GeV");
  legend2_global_sing_quad->Draw();
  gr_global_sing_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_global_sing_quad->SaveAs("plots/NuGun_JetRates_vs_BackgroundEff_Global.pdf");

  TCanvas *c2_global_350_120 = new TCanvas("c2_global_120","Graph Draw Options",200,10,600,400);
  gr_global_120_QCD->SetLineColor(4); // blue    
  gr_global_120_QCD->GetHistogram()->SetMinimum(-5.);
  gr_global_350_QCD->GetHistogram()->SetMinimum(-5.);
  gr_global_120_QCD->Draw("AC*");
  gr_global_120_QCD->GetXaxis()->SetLimits(0.,1.);
  gr_global_120_QCD->SetTitle("htSum Rate vs. Background Efficiency for Global Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_global_350_QCD->SetLineColor(2); // red     
  //  gr_global_350_QCD->Draw("C*");
  m_QCD_htSum120->SetMarkerStyle(21);
  m_QCD_htSum120->SetMarkerColor(4);
  m_QCD_htSum350->SetMarkerStyle(21);
  m_QCD_htSum350->SetMarkerColor(2);
  m_QCD_htSum430->SetMarkerStyle(21);
  m_QCD_htSum430->SetMarkerColor(3);
  m_QCD_htSum120->Draw();
  //  m_QCD_htSum350->Draw();
  m_QCD_htSum430->Draw();
  auto legend2_global_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend2_global_350_120->AddEntry(gr_global_120_QCD,"H_{T}>120 GeV, with timing cuts");
  //  legend2_global_350_120->AddEntry(gr_global_350_QCD,"H_{T}>350 GeV, with timing cuts");
  //  legend2_global_350_120->AddEntry(m_LLP_htSum350,"H_{T}>350 GeV, no timing cuts");
  legend2_global_350_120->AddEntry(m_LLP_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend2_global_350_120->Draw();
  gr_global_120_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_global_350_120->SaveAs("plots/NuGun_htSumRates_vs_BackgroundEff_Global.pdf");

  return 0;
}
