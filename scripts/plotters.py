import ROOT

class rootPlotter():

    def __init__(self):
        self.path = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/result_rates/"
        self.outpath = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/plots/"
        self.magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
        self.XCANVAS = 2400; self.YCANVAS = 2400
        ROOT.gROOT.SetBatch(True)
        ROOT.TH1.AddDirectory(0)

        self.QCD_dict = {"filename" : "QCD", "color" : ROOT.kBlack, "legendlabel" : "QCD"}
        self.NuGun_dict = {"filename" : "RelValNuGun", "color" : ROOT.kBlack, "legendlabel" : "NuGun"}
        
        self.LLP_m1000_500_dict = {"filename" : "LLP_MH1000_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
        self.LLP_m1000_1000_dict = {"filename" : "LLP_MH1000_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}
        
        self.LLP_m350_500_dict = {"filename" : "LLP_MH350_Ctau500", "color" : ROOT.kRed, "marker" : 2, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 500 mm"}
        self.LLP_m350_1000_dict = {"filename" : "LLP_MH350_Ctau1000", "color" : ROOT.kBlue, "marker" : 3, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 1000 mm"}
        
        self.LLP_m250_500_dict = {"filename" : "LLP_MH250_Ctau500", "color" : ROOT.kRed, "marker" : 4, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 500 mm"}
        self.LLP_m250_1000_dict = {"filename" : "LLP_MH250_Ctau1000", "color" : ROOT.kBlue, "marker": 5, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 1000 mm"}


    def plot1DNormalized(self, filelist, histname, bghistname, appendname = "", log = False, rebin = 1, xlabel = "", xrange = [-1, -1], yrange = [-1, -1], overflow = False):
        c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), self.XCANVAS, self.YCANVAS);
        ROOT.gPad.SetTopMargin(self.magicMargins["T"])
        ROOT.gPad.SetBottomMargin(self.magicMargins["B"])
        ROOT.gPad.SetLeftMargin(self.magicMargins["L"])
        ROOT.gPad.SetRightMargin(self.magicMargins["R"])
        #ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()
        ypos = 1 - self.magicMargins["T"] / 2
        l = ROOT.TLatex(self.magicMargins["L"], ypos, "CMS #scale[0.8]{#it{Preliminary}}" )
        l.SetTextAlign(12)
        l.SetNDC()
        l.SetTextSize(0.90 * self.magicMargins["T"])
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        logname = ""
        if log :
            ROOT.gPad.SetLogy()
            logname = "_log_"
        legend = ROOT.TLegend(0.4, 0.55, 0.8, 0.92)
        hcounter = 0
        for f in filelist:
            file = ROOT.TFile.Open(self.path+"rates_"+f["filename"]+".root")
            if "NuGun" or "QCD" in f["filename"]:
                hist = file.Get(bghistname)
            else:
                hist = file.Get(histname)
                hist.SetFillColor(16)
                hist.SetFillStyle(3004)
            if hcounter == 0:
                framehist = hist.Clone()
                framehist.Reset("ICE")
                framehist.SetStats(0)
                framehist.Rebin(rebin)
                if (xrange[0] != -1 and xrange[1] != -1): framehist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
                if (yrange[0] != -1 and yrange[1] != -1): framehist.GetYaxis().SetRangeUser(yrange[0], yrange[1])
                framehist.GetYaxis().SetTitle("A.U.")
                if xlabel: framehist.GetXaxis().SetTitle(xlabel)
                framehist.Draw()
            hist.Rebin(rebin)
            hist.Scale(1/hist.Integral())
            if overflow and xrange[1] != -1:
                overflowbin = hist.FindBin(xrange[1])
                totaloverflow = hist.Integral(overflowbin, hist.GetNbinsX())
                hist.SetBinContent(overflowbin - 1, totaloverflow)
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(1.7)
            hist.SetLineWidth(2)
            hist.SetMarkerColor(f["color"])
            hist.SetLineColor(f["color"])
            hist.Draw("h same")
            legendentr = f["legendlabel"]
            legend.AddEntry(hist, legendentr, "l")
            hcounter += 1
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.025)
        legend.SetHeader("#bf{CMS #scale[0.8]{#it{Preliminary}}}")
        legend.Draw("same")
#        l.Draw()
        c1.SaveAs(self.outpath+histname+appendname+logname+".pdf")
