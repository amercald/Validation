import ROOT
ROOT.TH1.AddDirectory(0)

path = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/rates_studies/"

nocut_dict = {"name" : "rate_none",
            "color" : ROOT.kBlack,
            "legendlabel" : "No requirement"}
cut09_dict = {"name" : "rate_09",
               "color" : ROOT.kGreen+2,
               "legendlabel" : "H(H+E) > 0.9"}
cut095_dict = {"name" : "rate_09",
               "color" : ROOT.kBlue,
               "legendlabel" : "H/(H+E) > 0.95"}
cut085_dict = {"name" : "rate_085",
               "color" : ROOT.kRed,
               "legendlabel" : "H/(H+E) > 0.85"}
cut099_dict = {"name" : "rate_099",
               "color" : ROOT.kMagenta,
               "legendlabel" : "H/(H+E) > 0.99"}
cut09_3x3_dict = {"name" : "rate_09_3x3",
               "color" : ROOT.kGreen+2,
               "legendlabel" : "3x3 H/(H+E) > 0.9"}
cut095_3x3_dict = {"name" : "rate_095_3x3",
               "color" : ROOT.kBlue,
               "legendlabel" : "3x3 H/(H+E) > 0.95"}
cut085_3x3_dict = {"name" : "rate_085_3x3",
               "color" : ROOT.kRed,
               "legendlabel" : "3x3 H/(H+E) > 0.85"}
cut099_3x3_dict = {"name" : "rate_099_3x3",
               "color" : ROOT.kMagenta,
               "legendlabel" : "H/(H+E) > 0.99"}

singlejet_dict = {"hist" : "singleJetRates_emu",
                  "current_thresh" : 180}

doublejet_dict = {"hist" : "doubleJetRates_emu",
                  "current_thresh" : 150}

def analyzerates():
    jet_list = [singlejet_dict, doublejet_dict]
    cut_list = [nocut_dict, cut085_dict, cut09_dict, cut095_dict, cut085_3x3_dict, cut09_3x3_dict, cut095_3x3_dict]
    with open('HovE_jet_thresholds.txt', 'w') as txt_file:
        for jet in jet_list:
            for cut in cut_list:
                file = ROOT.TFile.Open(path+"rates_hoe_"+cut["name"]+".root")
                txt_file.write(cut["legendlabel"]+"\n")
                hist = file.Get(jet["hist"])
                if cut == nocut_dict:
                    current_rate = hist.GetBinContent(jet["current_thresh"])
                new_thresh = hist.FindLastBinAbove(current_rate, 1)
                txt_file.write(jet["hist"]+":: Current rate: "+str(current_rate)+", New Threshold: "+str(new_thresh)+" GeV \n\n")

def main():
    analyzerates()

if __name__ == "__main__":
    main()
