import ROOT
import numpy as np

outpath = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/compare_plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
#ROOT.TH1.SetDefaultSumw2()
#ROOT.TProfile.SetDefaultSumw2()
path = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/rates_for_comparison/"

NuGun_11_0_X_PU = {"filename" : "NuGun_11_0_X_PU", "color" : ROOT.kRed, "legendlabel" : "RelValNuGun_11_0_X"}
NuGun_11_2_X_PU = {"filename" : "NuGun_11_2_X_PU", "color" : ROOT.kBlue, "legendlabel" : "RelValNuGun_11_2_X"}
NuGun_10_6_X_PU = {"filename" : "NuGun_10_6_X_PU", "color" : ROOT.kGreen, "legendlabel" : "RelValNuGun_10_6_X"}
NuGun_11_0_X_PT2_20 = {"filename" : "NuGun_11_0_X_PT2_20", "color" : ROOT.kBlack, "legendlabel" : "NuGun PT 2-20 11_0_X"}

def plot(filelist, histname, log, rebin, xrange, yrange, overflow):
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    logname = ""
    if log :
        ROOT.gPad.SetLogy()
        logname = "_log_"
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    hcounter = 0
    for f in filelist:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        hist = file.Get(histname)
        if hcounter == 0:
            framehist = hist.Clone()
            framehist.Reset("ICE")
            framehist.SetStats(0)
            framehist.Rebin(rebin)
            if (xrange[0] != -1 and xrange[1] != -1): framehist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            if (yrange[0] != -1 and yrange[1] != -1): framehist.GetYaxis().SetRangeUser(yrange[0], yrange[1])
            framehist.GetYaxis().SetTitle("A.U.")
            framehist.Draw()
#        if (xrange[0] != -1 and xrange[1] != -1): hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        hist.Rebin(rebin)
        if overflow and xrange[1] != -1:
            overflowbin = hist.FindBin(xrange[1])
            totaloverflow = hist.Integral(overflowbin, hist.GetNbinsX())
            normalizedoverflow = totaloverflow / hist.Integral()
            hist.SetBinContent(overflowbin - 2, totaloverflow)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.DrawNormalized("same")
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        hcounter += 1
    legend.Draw("same")
    c1.SaveAs(outpath+histname+logname+".pdf")
def main():
    #do stuff here
    filelist = [NuGun_10_6_X_PU, NuGun_11_0_X_PU, NuGun_11_2_X_PU, NuGun_11_0_X_PT2_20]
#    filelist = [NuGun_11_0_X_PU]
    plot(filelist, "hHTSum_emu", True, 4, [-1, -1], [10**-3, 1], False)
    plot(filelist, "hNJets_emu", True, 1, [-1, -1], [10**-3, 0.6], False) 
    plot(filelist, "hNJets_HT120_emu", True, 1, [-1, -1], [10**-3, 0.6], False) 
    plot(filelist, "hNJets_Good_emu", False, 1, [-1, -1], [10**-3, 0.6], False) 
    plot(filelist, "hNJets_Good_HT120_emu", False, 1, [-1, -1], [10**-3, 0.6], False) 
    plot(filelist, "hJetET_emu", True, 2, [0, 200], [10**-3, 1], False)
    plot(filelist, "hJetET_Leading1_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_Leading2_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_Leading3_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_Leading4_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT120_Leading1_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT120_Leading2_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT120_Leading3_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT120_Leading4_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT360_Leading1_emu", True, 2, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT360_Leading2_emu", True, 8, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT360_Leading3_emu", True, 8, [0, 200], [10**-3, 1], False) 
    plot(filelist, "hJetET_HT360_Leading4_emu", True, 8, [0, 200], [10**-3, 1], False) 

    plot(filelist, "hJetEta_emu", False, 1, [0, 200], [0, 0.5], False) 
    plot(filelist, "hJetEta_Leading1_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_Leading2_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_Leading3_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_Leading4_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT120_Leading1_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT120_Leading2_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT120_Leading3_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT120_Leading4_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT360_Leading1_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT360_Leading2_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT360_Leading3_emu", False, 1, [-1, -1], [0, 0.5], False) 
    plot(filelist, "hJetEta_HT360_Leading4_emu", False, 1, [-1, -1], [0, 0.5], False) 

    plot(filelist, "hcalTP_Barrel_emu", True, 1, [-1, -1], [10**-6, 1], False) 
    plot(filelist, "hcalTP_Endcap_emu", True, 1, [-1, -1], [10**-6, 1], False) 

    plot(filelist, "hcalTP_Barrel_HT120_emu", True, 1, [-1, -1], [10**-6, 1], False) 
    plot(filelist, "hcalTP_Endcap_HT120_emu", True, 1, [-1, -1], [10**-6, 1], False) 

    plot(filelist, "hcalTP_Barrel_HT360_emu", False, 1, [0, 10], [10**-6, 1], True) 
    plot(filelist, "hcalTP_Endcap_HT360_emu", False, 1, [0, 10], [10**-6, 1], True) 

    plot(filelist, "hcalTP_Barrel_emu", False, 1, [0, 10], [10**-6, 1], True) 
    plot(filelist, "hcalTP_Endcap_emu", False, 1, [0, 10], [10**-6, 1], True) 

    plot(filelist, "hcalTP_Barrel_HT120_emu", False, 1, [0, 10], [10**-6, 1], True) 
    plot(filelist, "hcalTP_Endcap_HT120_emu", False, 1, [0, 10], [10**-6, 1], True) 

    plot(filelist, "hcalTP_Barrel_HT360_emu", False, 1, [0, 10], [10**-6, 1], True) 
    plot(filelist, "hcalTP_Endcap_HT360_emu", False, 1, [0, 10], [10**-6, 1], True) 


    plot(filelist, "singleJetRates_emu", True, 1, [-1, -1], [10**-6, 1], False) 
    plot(filelist, "doubleJetRates_emu", True, 1, [-1, -1], [10**-6, 1], False) 
    plot(filelist, "tripleJetRates_emu", True, 1, [-1, -1], [10**-6, 1], False) 
    plot(filelist, "quadJetRates_emu", True, 1, [-1, -1], [10**-6, 1], False) 

if __name__ == "__main__":
    main()
