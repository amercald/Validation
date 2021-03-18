import ROOT
import numpy as np

class rootPlotter():

    def __init__(self):
        self.path = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/result_rates/"
        self.outpath = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/plots/"
        self.magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
        self.XCANVAS = 2400; self.YCANVAS = 2400
        ROOT.gROOT.SetBatch(True)
        ROOT.TH1.AddDirectory(0)

        self.QCD_dict = {"filename" : "QCD", "color" : ROOT.kBlack, "legendlabel" : "QCD"}
        self.RelValNuGun_dict = {"filename" : "RelValNuGun", "color" : ROOT.kBlack, "legendlabel" : "NuGun"}
        self.NuGun_dict = {"filename" : "NuGun", "color" : ROOT.kBlack, "legendlabel" : "NuGun"}
        
        self.LLP_m1000_500_dict = {"filename" : "LLP_MH1000_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
        self.LLP_m1000_1000_dict = {"filename" : "LLP_MH1000_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}
        
        self.LLP_m350_500_dict = {"filename" : "LLP_MH350_Ctau500", "color" : ROOT.kRed, "marker" : 2, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 500 mm"}
        self.LLP_m350_1000_dict = {"filename" : "LLP_MH350_Ctau1000", "color" : ROOT.kBlue, "marker" : 3, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 1000 mm"}
        
        self.LLP_m250_500_dict = {"filename" : "LLP_MH250_Ctau500", "color" : ROOT.kRed, "marker" : 4, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 500 mm"}
        self.LLP_m250_1000_dict = {"filename" : "LLP_MH250_Ctau1000", "color" : ROOT.kBlue, "marker": 5, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 1000 mm"}

        self.LLP_m350_80_500_dict = {"filename" : "LLP_MH350_80_Ctau500", "color" : ROOT.kRed, "marker" : 2, "legendlabel" : "LLP Jets, c#tau = 0.5 m"}
        self.LLP_m350_80_1000_dict = {"filename" : "LLP_MH350_80_Ctau1000", "color" : ROOT.kMagenta+2, "marker" : 3, "legendlabel" : "LLP Jets, c#tau = 1 m"}
        
        self.LLP_m250_60_500_dict = {"filename" : "LLP_MH250_60_Ctau500", "color" : ROOT.kRed, "marker" : 4, "legendlabel" : "LLP Jets, c#tau = 0.5 m"}
        self.LLP_m250_60_1000_dict = {"filename" : "LLP_MH250_60_Ctau1000", "color" : ROOT.kMagenta+2, "marker": 5, "legendlabel" : "LLP Jets, c#tau = 1 m"}


    def plot1DNormalized(self, filelist, histname, bghistname, appendname = "", title = "", log = False, rebin = 1, xlabel = "", xrange = [-1, -1], yrange = [-1, -1], overflow = False):
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
        if log:
            ROOT.gPad.SetLogy()
            logname = "_log_"
        legend = ROOT.TLegend(0.4, 0.55, 0.8, 0.92)
        hcounter = 0
        for f in filelist:
            file = ROOT.TFile.Open(self.path+"rates_"+f["filename"]+".root")
            if "NuGun" in f["filename"] or "QCD" in f["filename"]:
                hist = file.Get(bghistname)
                hist.SetFillColor(16)
                hist.SetFillStyle(3002)
            else:
                hist = file.Get(histname)

            if "hHTSum" in histname:
                int_above_360 = hist.Integral(72, 240) / hist.Integral()
                int_above_120 = hist.Integral(24, 240) / hist.Integral()
                print("INTEGRALS: ", int_above_360, int_above_120)
            if title:
                hist.SetTitle(title)
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
            if (hist.Integral() != 0):
                hist.Scale(1/hist.Integral())
            if overflow and xrange[1] != -1:
                overflowbin = hist.FindBin(xrange[1])
                totaloverflow = hist.Integral(overflowbin, hist.GetNbinsX())
                hist.SetBinContent(overflowbin - 1, totaloverflow)
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(1.7)
            hist.SetLineWidth(3)
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

    def plotProfiles(self, filelist, histname, bghistname, appendname = "", xlabel = "", xrange = [-1, -1], yrange = [-1, -1]):
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
        hlist = []
        legend = ROOT.TLegend(0.4, 0.55, 0.8, 0.92)
        hcounter = 0
        for f in filelist:
            file = ROOT.TFile.Open(self.path+"rates_"+f["filename"]+".root")
            if "NuGun" or "QCD" in f["filename"]:
                hist = file.Get(bghistname)
                hist.SetFillColor(16)
                hist.SetFillStyle(3004)
                
            else:
                hist = file.Get(histname)
            if hcounter == 0:
                framehist = hist.Clone()
                framehist.Reset("ICE")
                framehist.SetStats(0)
                if (xrange[0] != -1 and xrange[1] != -1): framehist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
                if (yrange[0] != -1 and yrange[1] != -1): framehist.GetYaxis().SetRangeUser(yrange[0], yrange[1])
                framehist.GetYaxis().SetTitle("A.U.")
                if xlabel: framehist.GetXaxis().SetTitle(xlabel)
                framehist.Draw()
            profilehist = hist.ProfileX()
            profilehist.SetMarkerStyle(20)
            profilehist.SetMarkerSize(1.7)
            profilehist.SetLineWidth(2)
            profilehist.SetMarkerColor(f["color"])
            profilehist.SetLineColor(f["color"])
#            profilehist.Draw("h same")
            legendentr = f["legendlabel"]
            legend.AddEntry(profilehist, legendentr, "l")
            hcounter += 1
            hlist.append(profilehist)

        for h in hlist:
            h.Draw("h same")
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.025)
        legend.SetHeader("#bf{CMS #scale[0.8]{#it{Preliminary}}}")
        legend.Draw("same")
#        l.Draw()
        c1.SaveAs(self.outpath+histname+appendname+".pdf")

    def plot1BinProjections(self, filelist, histname, bghistname, appendname = "", title = "", log = False, rebin = 1, xlabel = "", xrange = [-1, -1], yrange = [-1, -1], overflow = False):
        
        if "HB" in histname:
            ietamin = 1
            ietamax = 16
        if "HE" in histname:
            ietamin = 17
            ietamax = 29
            
        for ieta in range(ietamin, ietamax+1):
            c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), self.XCANVAS, self.YCANVAS);
            ROOT.gPad.SetTopMargin(self.magicMargins["T"])
            ROOT.gPad.SetBottomMargin(self.magicMargins["B"])
            ROOT.gPad.SetLeftMargin(self.magicMargins["L"])
            ROOT.gPad.SetRightMargin(self.magicMargins["R"])
            #ROOT.gPad.SetGridy()
            ROOT.gPad.SetTicks()
            legend = ROOT.TLegend(0.4, 0.55, 0.8, 0.92)
            hcounter = 0
            logname = ""
            if log:
                ROOT.gPad.SetLogy()
                logname = "_log_"
            for f in filelist:
                file = ROOT.TFile.Open(self.path+"rates_"+f["filename"]+".root")
                if "NuGun" in f["filename"] or "QCD" in f["filename"]:
                    hist = file.Get(bghistname)
                else:
                    hist = file.Get(histname)
                    
                if "HB" in histname:
                    ietabin = ieta
                if "HE" in histname:
                    ietabin = ieta - 16
                ietahist = hist.ProjectionX(histname+"_"+str(ieta)+f["filename"], ietabin, ietabin)
                if title:
                    ietahist.SetTitle(title)
                if hcounter == 0:
                    frameietahist = ietahist.Clone()
                    frameietahist.Reset("ICE")
                    frameietahist.SetStats(0)
                    frameietahist.Rebin(rebin)
                    if (xrange[0] != -1 and xrange[1] != -1): frameietahist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
                    if (yrange[0] != -1 and yrange[1] != -1): frameietahist.GetYaxis().SetRangeUser(yrange[0], yrange[1])
                    frameietahist.GetYaxis().SetTitle("A.U.")
                    if xlabel: frameietahist.GetXaxis().SetTitle(xlabel)
                    frameietahist.GetXaxis().SetLabelSize(0.05)
                    frameietahist.GetXaxis().SetTitleOffset(1.2)
                    frameietahist.GetYaxis().SetLabelSize(0.05)
                    frameietahist.GetYaxis().SetTitleOffset(1.2)

                    frameietahist.Draw()
                ietahist.Rebin(rebin)
                if (ietahist.Integral() != 0):
                    ietahist.Scale(1/ietahist.Integral())
                if overflow and xrange[1] != -1:
                    overflowbin = ietahist.FindBin(xrange[1])
                    totaloverflow = ietahist.Integral(overflowbin, ietahist.GetNbinsX())
                    ietahist.SetBinContent(overflowbin - 1, totaloverflow)
                ietahist.SetMarkerStyle(20)
                ietahist.SetMarkerSize(1.7)
                ietahist.SetLineWidth(3)
                ietahist.SetMarkerColor(f["color"])
                ietahist.SetLineColor(f["color"])

                if "NuGun" in f["filename"] or "QCD" in f["filename"]:
                    ietahist.SetFillColor(16)
                    ietahist.SetFillStyle(3002)
                ietahist.Draw("h same")
                legendentr = f["legendlabel"]
                legend.AddEntry(ietahist, legendentr, "l")
                hcounter += 1

            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.045)
#            legend.SetHeader("#bf{CMS #scale[0.8]{#it{Preliminary}}}")
            legend.SetHeader("|i#eta| = "+str(ieta))
            legend.Draw("same")
            #        l.Draw()
            c1.SaveAs(self.outpath+histname+"_"+str(ieta)+appendname+logname+".pdf")
            del c1

    def plotRatesVsEff(self, rateff_list, background, signal_list, effhistname, title, appendname, legendentry = ""):
        
        c1 = ROOT.TCanvas("%s"%(effhistname), "%s"%(effhistname), self.XCANVAS, self.YCANVAS);
        ROOT.gPad.SetTopMargin(self.magicMargins["T"])
        ROOT.gPad.SetBottomMargin(self.magicMargins["B"])
        ROOT.gPad.SetLeftMargin(self.magicMargins["L"])
        ROOT.gPad.SetRightMargin(self.magicMargins["R"])
        #ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()
        ROOT.gPad.SetLogy()
        ROOT.gPad.SetGrid()
        multigraph = ROOT.TMultiGraph()
        multigraph.SetMinimum(1)
        multigraph.SetMaximum(10**10)
        multigraph.GetXaxis().SetRangeUser(0, 1)
        multigraph.GetXaxis().SetTitle("LLP Signal Efficiency")
        multigraph.GetYaxis().SetTitle("Rate (Hz)")
        multigraph.SetTitle(title)
        legend = ROOT.TLegend(0.15, 0.7, 0.8, 0.92)

        bgfile = ROOT.TFile.Open(self.path+"rates_"+background["filename"]+".root")
        for sg in signal_list:

            rate_list = []
            rate_360_def = 0
            rate_120_def = 0
            eff_list = []
            eff_360_def = 0
            eff_120_def = 0
            for rateff in rateff_list:
                ratehist = bgfile.Get(rateff["rateshist"])
                rate = ratehist.GetBinContent(rateff["ratebin"])

                sgfile = ROOT.TFile.Open(self.path+"rates_"+sg["filename"]+".root")
                effhist = sgfile.Get(effhistname)
                eff = effhist.GetBinContent(rateff["effbin"])
                if rateff["rateshist"] == "htSumRates_emu":
                    if rateff["ratebin"] == 361:
                        rate_360_def = rate
                        eff_360_def = eff
                    if rateff["ratebin"] == 121:
                        rate_120_def = rate
                        eff_120_def = eff

                else:
                    rate_list.append(rate)
                    eff_list.append(eff)

            print(rate_list, eff_list)

            rateff_graph_360_def = ROOT.TGraph(1, np.array([eff_360_def]), np.array([rate_360_def]))
            rateff_graph_360_def.SetMarkerStyle(22)
            rateff_graph_360_def.SetMarkerSize(7)
            rateff_graph_360_def.SetMarkerColor(sg["color"])
            multigraph.Add(rateff_graph_360_def)
            legend.AddEntry(rateff_graph_360_def, sg["legendlabel"]+"; HT > 360", "lp")

            rateff_graph_120_def = ROOT.TGraph(1, np.array([eff_120_def]), np.array([rate_120_def]))
            rateff_graph_120_def.SetMarkerStyle(21)
            rateff_graph_120_def.SetMarkerSize(7)
            rateff_graph_120_def.SetMarkerColor(sg["color"])
            multigraph.Add(rateff_graph_120_def)
            legend.AddEntry(rateff_graph_120_def, sg["legendlabel"]+"; HT > 120", "lp")

            rate_arr = np.array(rate_list)
            eff_arr = np.array(eff_list)
            rateff_graph = ROOT.TGraph(len(rate_list), eff_arr, rate_arr)
            rateff_graph.SetMarkerSize(7)
            rateff_graph.SetMarkerStyle(31)
            rateff_graph.SetMarkerColor(sg["color"])
            rateff_graph.SetLineColor(sg["color"])

            legend.AddEntry(rateff_graph, sg["legendlabel"]+";"+legendentry, "lp")
            multigraph.Add(rateff_graph)


        multigraph.Draw("ALP")
        ROOT.gPad.Modified(); ROOT.gPad.Update()
        multigraph.GetXaxis().SetLimits(0., 1.);
        legend.SetBorderSize(0)
#        legend.SetFillStyle(0)
        legend.SetTextSize(0.025)
        legend.Draw("same")
        c1.SaveAs(self.outpath+"rateff"+appendname+".pdf")


    def plotJetRates(self, histname, title):
        c1 = ROOT.TCanvas("c1", "c1", self.XCANVAS, self.YCANVAS)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
        pad1.SetTopMargin(self.magicMargins["T"])
        pad1.SetBottomMargin(self.magicMargins["B"])
        pad1.SetLeftMargin(self.magicMargins["L"])
        pad1.SetRightMargin(self.magicMargins["R"])    
        pad1.SetGrid()
        pad1.SetTicks()
        pad1.SetLogy()
        pad1.Draw()

        pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
        pad2.SetTopMargin(0.01)
        pad2.SetBottomMargin(self.magicMargins["B"])
        pad2.SetLeftMargin(self.magicMargins["L"])
        pad2.SetRightMargin(self.magicMargins["R"])    
        pad2.SetGrid()
        pad2.SetTicks()
        pad2.Draw()

        legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)

        rate_no_cut = {"histname" : "_emu",
                       "legendentry" : "No LLP ID",
                       "color" : ROOT.kBlack }

        rate_HoE_cut = {"histname" : "_HoE_emu",
                       "legendentry" : "H/(H+E) > 0.9",
                       "color" : ROOT.kGreen+2 }

        rate_05_cut = {"histname" : "_HoE_TP05_emu",
                       "legendentry" : "LLP ID E_{T,TP} > 0.5 GeV",
                       "color" : ROOT.kRed }

        rate_5_cut = {"histname" : "_HoE_TP5_emu",
                       "legendentry" : "LLP ID E_{T,TP} > 5 GeV",
                       "color" : ROOT.kBlue }
                   

        file = ROOT.TFile.Open(self.path+"rates_"+self.NuGun_dict["filename"]+".root")

        rate_list = [rate_no_cut, rate_HoE_cut, rate_05_cut] 
        ratehist_list = []
        ratiohist_list = []
        for cut in rate_list:
            print(histname+cut["histname"])
            ratehist = file.Get(histname+cut["histname"])
            ratehist.SetLineColor(cut["color"])
            ratehist.SetLineWidth(2)
            ratehist.SetStats(0)
            ratehist_list.append(ratehist)

            if cut != rate_no_cut:
                ratiohist = ratehist.Clone("h2")
                ratiohist.Divide(file.Get(histname+rate_no_cut["histname"]))
                ratiohist.GetYaxis().SetTitle("Tag / No tag")
                ratiohist_list.append(ratiohist)

            legend.AddEntry(ratehist, cut["legendentry"], "l")
        
        pad1.cd()
        hcounter = 0
        for h in ratehist_list:
            if hcounter == 0:
                h.GetYaxis().SetRangeUser(1, 10**8)
                h.SetTitle(title)
                h.Draw("hist")
            else:
                h.Draw("hist same")
            hcounter += 1
        legend.Draw("same")
        pad2.cd()
        hcounter = 0
        for h in ratiohist_list:
            if hcounter == 0:
                h.GetYaxis().SetRangeUser(0, 1.4)
                h.GetYaxis().SetTitleSize(0.075)
                h.GetYaxis().SetLabelSize(0.066)
                h.GetXaxis().SetTitleSize(0.07)
                h.GetXaxis().SetLabelSize(0.066)
                h.GetXaxis().SetTitle("")

                h.Draw("hist")
            else:
                h.Draw("hist same")
            hcounter += 1


        c1.SaveAs(self.outpath+histname+"_HovE_FB.pdf")
