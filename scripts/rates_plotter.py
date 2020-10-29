import ROOT

outpath = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.AddDirectory(0)
ROOT.TProfile.AddDirectory(0)
#ROOT.TProfile.SetDefaultSumw2()
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
               "legendlabel" : "H/(H+E) > 0.9"}
cut095_3x3_dict = {"name" : "rate_095_3x3",
               "color" : ROOT.kBlue,
               "legendlabel" : "H/(H+E) > 0.95"}
cut085_3x3_dict = {"name" : "rate_085_3x3",
               "color" : ROOT.kRed,
               "legendlabel" : "H/(H+E) > 0.85"}
cut099_3x3_dict = {"name" : "rate_099_3x3",
               "color" : ROOT.kMagenta,
               "legendlabel" : "H/(H+E) > 0.99"}

file_list = [nocut_dict, cut09_dict, cut085_dict, cut095_dict]
file_list_3x3 = [nocut_dict, cut09_3x3_dict, cut085_3x3_dict, cut095_3x3_dict]


def plotRatesHoEcut(histname, currentlist):
    
    towersize = "1x1"
    if currentlist == file_list_3x3: towersize = "3x3"
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
    pad1.SetTopMargin(magicMargins["T"])
    pad1.SetBottomMargin(magicMargins["B"])
    pad1.SetLeftMargin(magicMargins["L"])
    pad1.SetRightMargin(magicMargins["R"])    
#    pad1.SetGridy()
    pad1.SetTicks()
    pad1.SetLogy()
    pad1.Draw()

    pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
    pad2.SetTopMargin(magicMargins["T"])
    pad2.SetBottomMargin(magicMargins["B"])
    pad2.SetLeftMargin(magicMargins["L"])
    pad2.SetRightMargin(magicMargins["R"])    
#    pad2.SetGridy()
    pad2.SetTicks()
#    pad2.SetLogy()
    pad2.Draw()
    
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    hist_list = []
    ratiohist_list = []
    print("plotting "+histname)
    bfile = ROOT.TFile.Open(path+"rates_hoe_"+nocut_dict["name"]+".root")
    bhist = bfile.Get("singleJetRates_emu")

    for filename in currentlist:
#        print("using rates_depth_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_hoe_"+filename["name"]+".root")
        hist = file.Get(histname)
        hist.SetLineColor(filename["color"])
        hist.SetMarkerColor(filename["color"])
        hist.SetMarkerStyle(2)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(2)
        hist.Rebin(1)
        hist.SetStats(0)
        hist.GetXaxis().SetRangeUser(0, 200)
        if filename != nocut_dict:
            ratiohist = hist.Clone("h2")
            ratiohist.Divide(bhist)
            ratiohist.GetYaxis().SetTitle("cut/no cut")
            ratiohist_list.append(ratiohist)

        hist.GetYaxis().SetRangeUser(10**3 - 1, 10**7)
        hist.GetYaxis().SetTitle("Rate (Hz)")
        hist.GetXaxis().SetTitle("Jet Threshold")
        legend.AddEntry(hist, filename["legendlabel"], "l")
        hist_list.append(hist)

    pad1.cd()
    hcounter = 0   
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("hist")
        else:
            ihist.Draw("hist same")
        hcounter += 1
        legend.Draw("same")
    pad2.cd()
    hcounter = 0   
    for ihist in ratiohist_list:
        if hcounter == 0:
            ihist.Draw("hist")
        else:
            ihist.Draw("hist same")
        hcounter += 1
    c1.SaveAs(outpath+histname+"_"+towersize+".pdf")
    del c1

if __name__ == "__main__":
    
    plotRatesHoEcut("singleJetRates_emu", file_list)
    plotRatesHoEcut("singleJetRates_emu", file_list_3x3)

