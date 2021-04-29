import ROOT
import numpy as np

outpath = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
#ROOT.TH1.SetDefaultSumw2()
#ROOT.TProfile.SetDefaultSumw2()
path = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/result_rates/"
#variableplots = ["jetEta", "jetET", "jetEtaLeading1", "jetEtaLeading2", "jetEtaLeading3", "jetEtaLeading4", "njets", "jetETLeading1", "jetETLeading2", "jetETLeading3", "jetETLeading4"]
variableplots = ["hJet1x1ov5x5", "DeltaRLLP", "nGenParticles"]

QCD_dict = {"filename" : "QCD", "color" : ROOT.kBlack, "legendlabel" : "QCD"}
NuGun_dict = {"filename" : "NuGun", "color" : ROOT.kBlack, "legendlabel" : "NuGun"}

LLP_m1000_500_dict = {"filename" : "LLP_MH1000_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
LLP_m1000_1000_dict = {"filename" : "LLP_MH1000_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}

LLP_m350_500_dict = {"filename" : "LLP_MH350_Ctau500", "color" : ROOT.kRed, "marker" : 2, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 500 mm"}
LLP_m350_1000_dict = {"filename" : "LLP_MH350_Ctau1000", "color" : ROOT.kBlue, "marker" : 3, "legendlabel" : "LLP m_{H} = 350 GeV, c#tau = 1000 mm"}

LLP_m250_500_dict = {"filename" : "LLP_MH250_Ctau500", "color" : ROOT.kRed, "marker" : 4, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 500 mm"}
LLP_m250_1000_dict = {"filename" : "LLP_MH250_Ctau1000", "color" : ROOT.kBlue, "marker": 5, "legendlabel" : "LLP m_{H} = 250 GeV, c#tau = 1000 mm"}

file_list_1000 = { "list" : [LLP_m1000_1000_dict, LLP_m1000_500_dict, QCD_dict,], "name" : "MH_1000"}
file_list_350 = { "list" : [LLP_m350_1000_dict,LLP_m350_500_dict, QCD_dict], "name" : "MH_350"}
file_list_250 = { "list" : [LLP_m250_1000_dict, LLP_m250_500_dict, QCD_dict], "name" : "MH_250"}

file_list_NuGun_1000 = { "list" : [LLP_m1000_1000_dict, LLP_m1000_500_dict, NuGun_dict,], "name" : "MH_1000"}
file_list_NuGun_350 = { "list" : [LLP_m350_1000_dict,LLP_m350_500_dict, NuGun_dict], "name" : "MH_350"}
file_list_NuGun_250 = { "list" : [LLP_m250_1000_dict, LLP_m250_500_dict, NuGun_dict], "name" : "MH_250"}

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
hist_dict = {}

def plotHEEnergy(histname, filename):
    print("opening "+path+"rates_"+filename["filename"]+".root")
    file = ROOT.TFile.Open(path+"rates_"+filename["filename"]+".root")
    print("getting "+histname)
    hist = file.Get(histname)
    HCALhist = hist.ProjectionX()
    ECALhist = hist.ProjectionY()

    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    #ROOT.gPad.SetLogz()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    hist.SetStats(0)
    hist.GetXaxis().SetTitle("ECAL Energy")
    hist.GetYaxis().SetTitle("HCAL Energy")
    if filename["filename"] == "NuGun":
        hist.GetXaxis().SetRangeUser(0, 100)
        hist.GetYaxis().SetRangeUser(0, 100)
    
    hist.Draw("colz")

    c1.SaveAs(outpath+histname+"_2D_"+filename["filename"]+".pdf")
    del c1

    c2 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()

    HCALhist.SetMarkerStyle(20)
    HCALhist.SetMarkerSize(1.7)
    HCALhist.SetLineWidth(3)
    
    HCALhist.GetXaxis().SetTitle("ECAL Energy")
    HCALhist.Draw("h")
    
    c2.SaveAs(outpath+histname+"_ECAL_"+filename["filename"]+".pdf")
    del c2

    c3 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()    

    ECALhist.SetMarkerStyle(20)
    ECALhist.SetMarkerSize(1.7)
    ECALhist.SetLineWidth(3)

    
    ECALhist.GetXaxis().SetTitle("HCAL Energy")
    ECALhist.Draw("h")
    
    c3.SaveAs(outpath+histname+"_HCAL_"+filename["filename"]+".pdf")

def plotHoE(histname, filename):
    file = ROOT.TFile.Open(path+"rates_"+filename["filename"]+".root")
    hist = file.Get(histname)
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.7)
    hist.SetLineWidth(3)
    hist.GetYaxis().SetRangeUser(1, 10**3)
    
    hist.GetXaxis().SetTitle("H/(H+E)")
    #hist.SetStats(0)
    hist.Draw("h")    
    c1.SaveAs(outpath+histname+"_"+filename["filename"]+".pdf")

def plotHoEShape(histname, filename):
    file = ROOT.TFile.Open(path+"rates_"+filename+".root")
    hist = file.Get(histname)
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.7)
    hist.SetLineWidth(3)
    hist.GetYaxis().SetRangeUser(1, 10**3)
    
    hist.GetXaxis().SetTitle("H/(H+E)")
    #hist.SetStats(0)
    hist.DrawNormalized("h")    
    c1.SaveAs(outpath+histname+"_"+filename+"_Normalized.pdf")


def plotHoEGenSame(histname, bghistname, file_list_dict):

    file_list = file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.15, 0.71, 0.45, 0.92)
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        if "qcd" in f["filename"].lower() or "nugun" in f["filename"].lower():
            hist = file.Get(bghistname)
        else:
            hist = file.Get(histname)
    
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.GetXaxis().SetTitle("H/(H+E)")
        hist.SetStats(0)
        hist.GetYaxis().SetLimits(10**-5, 1000)
        hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        if hcounter == 0:
            hist.DrawNormalized("h")
        else:
            hist.DrawNormalized("h same")
        hcounter += 1
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")


def plotHoESame(histname, file_list_dict):

    file_list = file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.15, 0.71, 0.45, 0.92)
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        hist = file.Get(histname)
    
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.GetXaxis().SetTitle("H/(H+E)")
        hist.SetStats(0)

        hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        if hcounter == 0:
            hist.DrawNormalized("h")
            hist.GetYaxis().SetLimits(10**-5, 100)    
        else:
            hist.DrawNormalized("h same")
        hcounter += 1

    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")

    
def plotvariable(filename):
    file = ROOT.TFile.Open(path+"rates_"+filename+".root")
    for var in variableplots:
        print("plotting variable "+var)
        c4 = ROOT.TCanvas("%s"%(var), "%s"%(var), XCANVAS, YCANVAS);
        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()

        varhist = file.Get(var)
        xtitle = ""
        if var == "jetEta":
            xtitle = "#eta"
        elif var == "jetET":
            xtitle = "E_{T}"
        elif var == "njets":
            xtitle = "N_{#mathrm{jets}}"

        varhist.SetMarkerStyle(20)
        varhist.SetMarkerSize(1.7)
        varhist.SetLineWidth(3)

        varhist.GetXaxis().SetTitle(xtitle)
        varhist.Draw("h")
        c4.SaveAs(outpath+var+"_"+filename+".pdf")
        del c4

def plotvariablesame(histname, bghistname, file_list_dict, log, rebin, ymax, xrange):
    file_list =file_list_dict["list"]
    mh = file_list_dict["name"]
    
    c5 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    log_name = ""
    if log : 
        ROOT.gPad.SetLogy()
        log_name = "Log"
#        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
        
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        if "qcd" in f["filename"].lower() or "nugun" in f["filename"].lower():
            varhist = file.Get(bghistname)
        else:
            varhist = file.Get(histname)
        varhist.SetStats(0)
        varhist.SetMarkerStyle(20)
        varhist.SetMarkerSize(1.7)
        varhist.SetLineWidth(3)
        varhist.SetMarkerColor(f["color"])
        varhist.SetLineColor(f["color"])
        varhist.Rebin(rebin)
        if (xrange[0] != -1): varhist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if hcounter == 0:
            nbins = varhist.GetNbinsX()
            if xrange[0] == -1: xrange[0] = varhist.GetBinLowEdge(1)
            if xrange[1] == -1: xrange[1] = varhist.GetBinLowEdge(varhist.GetMaximumBin()) + varhist.GetBinWidth(varhist.GetMaximumBin())
            frame = ROOT.TH1F("frame", "", nbins, xrange[0], xrange[1])
            if ymax != -1: frame.SetMaximum(ymax)
            frame.SetStats(0)
            frame.GetXaxis().SetTitle(varhist.GetXaxis().GetTitle())
            frame.Draw()
        varhist.DrawNormalized("h same")
        legendentr = f["legendlabel"]
        legend.AddEntry(varhist, legendentr, "l")
        hcounter += 1
    legend.Draw("same")
    c5.SaveAs(outpath+histname+"_"+mh+"_"+log_name+"_All.pdf")
    del c5
    del legend

def plotETRatio(histname, histcutname):
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.15, 0.77, 0.4, 0.92)
    for f in range(len(file_list)):
        file = ROOT.TFile.Open(path+"rates_"+file_list[f]+".root")
        histnocut = file.Get(histname)    
        hist = file.Get(histcutname)
        hist.Divide(histnocut)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(colors[f])
        hist.SetLineColor(colors[f])
        hist.GetXaxis().SetTitle("E_{T} (H/(H+E) > 0.9) / E_{T}")
        hist.SetStats(0)
        
        
        legendentr = file_list[f]
        if file_list[f][:3] == "LLP":
            legendentr = "LLP c#tau = "+file_list[f][4:]
        legend.AddEntry(hist, legendentr, "l")
        if f == 0:
            hist.Draw("h")    
        else:
            hist.Draw("h same")
    legend.Draw("same")
    c1.SaveAs(outpath+histcutname+"_ratio_All.pdf")

    

def plotHoEvsET(histname, filename):
    print("opening "+path+"rates_"+filename["filename"]+".root")
    file = ROOT.TFile.Open(path+"rates_"+filename["filename"]+".root")
    hist = file.Get(histname)

    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    #ROOT.gPad.SetLogz()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    #hist.SetStats(0)
    hist.GetXaxis().SetTitle("E_{T}")
    hist.GetYaxis().SetTitle("H/(H+E)")
    
    hist.Draw("colz")    
    c1.SaveAs(outpath+histname+"_2D_"+filename["filename"]+".pdf")
    del c1
    #Draw Profile Histogram on top
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogz()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    profileHist = hist.ProfileX()
    profileHist.SetLineColor(ROOT.kMagenta)
    profileHist.SetMarkerStyle(20)
    profileHist.SetMarkerSize(1.7)
    profileHist.SetLineWidth(3)
    profileHist.SetMarkerColor(ROOT.kMagenta)
    hist.SetStats(0)
    hist.Draw("colz")
    profileHist.Draw("spread same")
    c1.SaveAs(outpath+histname+"_Profile_"+filename["filename"]+".pdf")

    

    if not histname in hist_dict.keys():
        hist_dict[histname] = {}
    
    for xbin in range(hist.GetNbinsX()):
        if not xbin in hist_dict[histname].keys():
            hist_dict[histname][xbin] = []
        binSize = 50
        histETbin = hist.ProjectionY("",xbin, xbin+1)
        histETbin.SetTitle("H/(H+E) ET from"+str(xbin*binSize)+" - "+str((xbin+1)*binSize)+" "+histname+" "+filename["filename"])
        histETbin.GetXaxis().SetTitle("H/(H+E)")
        hist_dict[histname][xbin].append(histETbin)

def plot2DwithProfile(histname, filename, log):
    print("opening "+path+"rates_"+filename["filename"]+".root")
    file = ROOT.TFile.Open(path+"rates_"+filename["filename"]+".root")
    hist = file.Get(histname)

    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    if log: ROOT.gPad.SetLogz()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    hist.DrawNormalized("colz")    
#    c1.SaveAs(outpath+histname+"_2D_"+filename["filename"]+".pdf")
#    del c1
    #Draw Profile Histogram on top
#    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
#    ROOT.gPad.SetLogz()
#    ROOT.gPad.SetTopMargin(magicMargins["T"])
#    ROOT.gPad.SetBottomMargin(magicMargins["B"])
#    ROOT.gPad.SetLeftMargin(magicMargins["L"])
#    ROOT.gPad.SetRightMargin(magicMargins["R"])
#    ROOT.gPad.SetGridy()
#    ROOT.gPad.SetTicks()
    profileHist = hist.ProfileX()
    profileHist.SetLineColor(ROOT.kMagenta)
    profileHist.SetMarkerStyle(20)
    profileHist.SetMarkerSize(1.7)
    profileHist.SetLineWidth(3)
    profileHist.SetMarkerColor(ROOT.kMagenta)
    hist.SetStats(0)
    hist.DrawNormalized("colz")
    profileHist.Draw("spread same")
    c1.SaveAs(outpath+histname+"_Profile_"+filename["filename"]+".pdf")


def plotHoEvsDepthSame(histname, file_list_dict):

    file_list = file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    #ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    hist_list = []
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.15, 0.71, 0.45, 0.92)
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        initialhist = file.Get(histname)
        hist = initialhist.ProfileX()
        if "depthratio" in histname.lower():
            hist.RebinX(5)
            hist.GetXaxis().SetRangeUser(0,4)
        else:
            hist.RebinX(4)

        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
#        hist.GetXaxis().SetTitle("TP E_{T}")
        hist.GetYaxis().SetTitle("H/(H+E)")
        hist.SetStats(0)
        hist.GetYaxis().SetRangeUser(0, 1.5)

#        hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        hist_list.append(hist)
    
    hcounter = 0
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("h")
        else:
            ihist.Draw("h same")
        hcounter += 1

    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")

def plotHoEvsDepthGenSame(histname, bghistname, file_list_dict):

    file_list = file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    #ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    hist_list = []
    
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.15, 0.71, 0.45, 0.92)
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        if "qcd" in f["filename"].lower() or "nugun" in f["filename"].lower():
            initialhist = file.Get(bghistname)
        else:
            initialhist = file.Get(histname)
            
        hist = initialhist.ProfileX()
        hist.RebinX(4)

        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
#        hist.GetXaxis().SetTitle("TP E_{T}")
        hist.GetYaxis().SetTitle("H/(H+E)")
        hist.SetStats(0)
        hist.GetYaxis().SetRangeUser(0, 1.5)

#        hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        hist_list.append(hist)
    
    hcounter = 0
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("h")
        else:
            ihist.Draw("h same")
        hcounter += 1

    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")

def plotDepthRatiosame(histname, file_list_dict, rebin):
    file_list =file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogy()
    #        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        hist = file.Get(histname)
        hist.SetStats(0)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.Rebin(rebin)
 #       hist.GetXaxis().SetRangeUser(0, 4)
        if hcounter == 0:
            hist.DrawNormalized("h")
        else:
            hist.DrawNormalized("h same")
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        hcounter += 1
        legend.Draw("same")
        c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")
    del c1
    del legend


def plotDepthRatioGensame(histname, bghistname, file_list_dict, rebin):
    file_list =file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogy()
    #        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        if "qcd" in f["filename"].lower() or "nugun" in f["filename"].lower(): 
            hist = file.Get(bghistname)
        else:
            hist = file.Get(histname)
        hist.SetStats(0)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(2)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.Rebin(rebin)
#        hist.GetXaxis().SetRangeUser(0, 4)
        if hcounter == 0:
            hist.DrawNormalized("h")
        else:
            hist.DrawNormalized("h same")
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        hcounter += 1
        legend.Draw("same")
        c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")
    del c1
    del legend

def plotDepthEnergyGensame(histname, bghistname, file_list_dict, rebin):
    file_list =file_list_dict["list"]
    mh = file_list_dict["name"]
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    #ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogy()
    #        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    
    hist_list = []
    hcounter = 0
    for f in file_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
        if "qcd" in f["filename"].lower() or "nugun" in f["filename"].lower(): 
            hist = file.Get(bghistname)
        else:
            hist = file.Get(histname)
        hist.SetStats(0)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(2)
        hist.SetMarkerColor(f["color"])
        hist.SetLineColor(f["color"])
        hist.Rebin(rebin)
#        hist.GetXaxis().SetRangeUser(0, 4)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        if hcounter == 0:
            hist.DrawNormalized("h")
        else:
            hist.DrawNormalized("h same")
        hcounter += 1
        legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mh+"_All.pdf")
    del c1
    del legend

def signaleff_Rate(effname, ratename, signal_list, background):
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])

    ROOT.gPad.SetTicks()

    multigraph = ROOT.TMultiGraph()
    multigraph.SetMinimum(1)
    multigraph.SetMaximum(10**7)
    multigraph.GetXaxis().SetTitle("Signal Efficiency")
    multigraph.GetYaxis().SetTitle("Rate (Hz)")


    eff_arr = []
    rate_arr = []

    bgfile = ROOT.TFile.Open(path+"rates_"+background["filename"]+".root")
    rateplot = bgfile.Get(ratename)
    rate_120 = rateplot.GetBinContent(120)
    rate_360 = rateplot.GetBinContent(360)

    for f in signal_list:
        file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")

        effplot = file.Get(effname)
            
        eff_360 = effplot.GetEfficiency(1)
        eff_120 = effplot.GetEfficiency(2)
        graph_120 = ROOT.TGraph(1, np.array([eff_120]), np.array([rate_120]))
        graph_120.SetMarkerSize(7)
        graph_120.SetMarkerStyle(22)
        graph_120.SetMarkerColor(f["marker"])
        graph_120.GetXaxis().SetTitle("Signal Efficiency")
        graph_120.GetYaxis().SetTitle("Rate (Hz)")
        graph_120.SetTitle(f["legendlabel"])
        graph_120.GetXaxis().SetRangeUser(0, 1)
        graph_360 = ROOT.TGraph(1, np.array([eff_360]), np.array([rate_360]))
        print(np.array([eff_360]), np.array([rate_360]))
        graph_360.SetMarkerSize(7)
        graph_360.SetMarkerStyle(21)
        graph_360.SetMarkerColor(f["marker"])
        graph_360.GetXaxis().SetTitle("Signal Efficiency")
        graph_360.GetYaxis().SetTitle("Rate (Hz)")
        graph_360.SetTitle(f["legendlabel"])
        graph_360.GetXaxis().SetRangeUser(0, 1)
        multigraph.Add(graph_120)
        multigraph.Add(graph_360)

    multigraph.GetXaxis().SetLimits(0, 1)
    multigraph.GetXaxis().SetRangeUser(0, 1)

    multigraph.Draw("AP")
    c1.BuildLegend()
    multigraph.GetXaxis().SetLimits(0, 1)
    multigraph.GetXaxis().SetRangeUser(0, 1)
    ROOT.gPad.Modified()

    c1.SaveAs(outpath+"rates_"+effname+".pdf")

        
if __name__ == "__main__":
    #file_list_list = [file_list_1000, file_list_350, file_list_250]
    #file_list_list = [file_list_250, file_list_350]
    file_list_list = [file_list_NuGun_250, file_list_NuGun_350]

    for l in file_list_list:
        for f in l["list"]:
            #plotHEEnergy("HEEnergytotal_1x1_emu_AllJet", f)
            #plotHEEnergy("HEEnergytotal_3x3_emu_AllJet", f)
            #plotvariable(f)
#            plotHoE("HovEtotal_1x1_emu_AllJets", f)
#            plotHoE("HovEtotal_3x3_emu_AllJets", f)
#            plotHoE("HovEtotal_1x1_emu", f)
#            plotHoE("HovEtotal_3x3_emu", f)
#            plotHoEShape("HovEtotal_1x1_emu", f)
#            plotHoEShape("HovEtotal_3x3_emu", f)
            #plotHEEnergy("HEEnergytotal_1x1_emu_Leading1", f)
            #plotHEEnergy("HEEnergytotal_3x3_emu_Leading1", f)        
#            plotHEEnergy("HEEnergytotal_1x1_emu_Leading2", f)
#            plotHEEnergy("HEEnergytotal_3x3_emu_Leading2", f)        
#            plotHEEnergy("HEEnergytotal_1x1_emu_Leading3", f)
#            plotHEEnergy("HEEnergytotal_3x3_emu_Leading3", f)
#            plotHEEnergy("HEEnergytotal_1x1_emu_Leading4", f)
 #           plotHEEnergy("HEEnergytotal_3x3_emu_Leading4", f)                
#            
            #plotHoE("HovEtotal_1x1_emu_Leading1", f)
            #plotHoE("HovEtotal_3x3_emu_Leading1", f)        
            #plotHoE("HovEtotal_1x1_emu_Leading2", f)
            #plotHoE("HovEtotal_3x3_emu_Leading2", f)        
            #plotHoE("HovEtotal_1x1_emu_Leading3", f)
            #plotHoE("HovEtotal_3x3_emu_Leading3", f)
            #plotHoE("HovEtotal_1x1_emu_Leading4", f)
            #plotHoE("HovEtotal_3x3_emu_Leading4", f)                
            
            plotHoEvsET("HovEtotal_1x1_ET_emu_Leading1", f)
            plotHoEvsET("HovEtotal_3x3_ET_emu_Leading1", f)        
            plotHoEvsET("HovEtotal_1x1_ET_emu_Leading2", f)
            plotHoEvsET("HovEtotal_3x3_ET_emu_Leading2", f)        
            plotHoEvsET("HovEtotal_1x1_ET_emu_Leading3", f)
            plotHoEvsET("HovEtotal_3x3_ET_emu_Leading3", f)
            plotHoEvsET("HovEtotal_1x1_ET_emu_Leading4", f)
            plotHoEvsET("HovEtotal_3x3_ET_emu_Leading4", f)

#            plot2DwithProfile("HovEtotal_3x3_ET_Depth1_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth2_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth3_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth4_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth5_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth6_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth7_Barrel_emu", f)
            
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu", f)

#            plot2DwithProfile("HovEtotal_3x3_ET_Depth1_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth2_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth3_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth4_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth5_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth6_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth7_Endcap_emu", f)
            
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu", f)
#            plot2DwithProfile("HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu", f)

#            plot2DwithProfile("energyDepth_NTPs_Barrel_D4_ge", f)
#            plot2DwithProfile("energyDepth_NTPs_Barrel_D2_le", f)
#            plot2DwithProfile("energyDepth_NTPs_Gen_Barrel_D4_ge", f)
#            plot2DwithProfile("energyDepth_NTPs_Gen_Barrel_D2_le", f)

#            plot2DwithProfile("energyDepth_NTPs_Endcap_D4_D7_OR_ge", f)
#            plot2DwithProfile("energyDepth_NTPs_Endcap_D1_le", f)
#            plot2DwithProfile("energyDepth_NTPs_Gen_Endcap_D4_D7_OR_ge", f)
#            plot2DwithProfile("energyDepth_NTPs_Gen_Endcap_D1_le", f)

            #Plot2dwithprofile("Hovetotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu", f)
            #plot2DwithProfile("HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu", f)
            #plot2DwithProfile("HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu", f)
            #plot2DwithProfile("HovEtotal_3x3_DepthRatio_Gen_Barrel_3_4_ov_1_2_emu", f)
            #plot2DwithProfile("HovEtotal_3x3_DepthRatio_Gen_Barrel_4_ov_1_2_3_emu", f)
            #plot2DwithProfile("HovEtotal_3x3_DepthRatio_Gen_Endcap_4_7_ov_1_3_emu", f)

            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max", f, True)
            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max_Gen", f, True)
            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max_HT120", f, True)
            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max_Gen_HT120", f, True)
            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max_HT360", f, True)
            plot2DwithProfile("energyDepth_NTPs_HBD4_HED47_Max_Gen_HT360", f, True)


            plot2DwithProfile("JetEt_Thresh_MaxJet", f, False)
            plot2DwithProfile("JetEt_Thresh_MaxJet_Gen", f, False)
            plot2DwithProfile("JetEt_Thresh_MaxJet_HT120", f, False)
            plot2DwithProfile("JetEt_Thresh_MaxJet_Gen_HT120", f, False)
            plot2DwithProfile("JetEt_Thresh_MaxJet_HT360", f, False)
            plot2DwithProfile("JetEt_Thresh_MaxJet_Gen_HT360", f, False)






#    draw_hist_dict()
    for l in file_list_list:
        plotvariablesame("L1JetTPDR", "L1JetTPDR", l, False, 1, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut0_5", "L1JetTPDR_cut0_5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut1", "L1JetTPDR_cut1", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut1_5", "L1JetTPDR_cut1_5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut2", "L1JetTPDR_cut2", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut2_5", "L1JetTPDR_cut2_5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut3", "L1JetTPDR_cut3", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut3_5", "L1JetTPDR_cut3_5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut4", "L1JetTPDR_cut4", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut4_5", "L1JetTPDR_cut4_5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_cut5", "L1JetTPDR_cut5", l, False, 2, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen", "L1JetTPDR", l, False, 1, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut0_5", "L1JetTPDR_cut0_5", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut1", "L1JetTPDR_cut1", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut1_5", "L1JetTPDR_cut1_5", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut2", "L1JetTPDR_cut2", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut2_5", "L1JetTPDR_cut2_5", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut3", "L1JetTPDR_cut3", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut3_5", "L1JetTPDR_cut3_5", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut4", "L1JetTPDR_cut4", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut4_5", "L1JetTPDR_cut4_5", l, False, 4, 0.8, [0, 5])
        plotvariablesame("L1JetTPDR_Gen_cut5", "L1JetTPDR_cut5", l, False, 4, 0.8, [0, 5])

        plotvariablesame("L1JetNTP", "L1JetNTP", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut0_5", "L1JetNTP_cut0_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut1", "L1JetNTP_cut1", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut1_5", "L1JetNTP_cut1_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut2", "L1JetNTP_cut2", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut2_5", "L1JetNTP_cut2_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut3", "L1JetNTP_cut3", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut3_5", "L1JetNTP_cut3_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut4", "L1JetNTP_cut4", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut4_5", "L1JetNTP_cut4_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_cut5", "L1JetNTP_cut5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut0_5", "L1JetNTP_cut0_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut1", "L1JetNTP_cut1", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut1_5", "L1JetNTP_cut1_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut2", "L1JetNTP_cut2", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut2_5", "L1JetNTP_cut2_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut3", "L1JetNTP_cut3", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut3_5", "L1JetNTP_cut3_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut4", "L1JetNTP_cut4", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut4_5", "L1JetNTP_cut4_5", l, True, 1, 1, [0, 20])
        plotvariablesame("L1JetNTP_Gen_cut5", "L1JetNTP_cut5", l, True, 1, 1, [0, 20])
        plotvariablesame("hJet1x1ov5x5", "hJet1x1ov5x5", l, True, 1, 1, [-1, -1])
        plotvariablesame("jetET", "jetET", l, True, 2, -1, [-1, -1])
        plotvariablesame("hcalTP_nearL1Jet_emu", "hcalTP_nearL1Jet_emu", l, False, 1, -1, [0, 10])
        plotvariablesame("hcalTP_nearL1Jet_Gen_emu", "hcalTP_nearL1Jet_emu", l, False, 1, -1, [0, 10])
        plotvariablesame("hcalTP_emu", "hcalTP_emu", l, False, 1, -1, [0, 10])
        plotvariablesame("jetET", "jetET", l, True, 2, -1, [-1, -1])
        plotvariablesame("hcalTP_nearL1Jet_Overflowge10_emu", "hcalTP_nearL1Jet_Overflowge10_emu", l, False, 1, 0.6, [0, 10])
        plotvariablesame("hcalTP_nearL1Jet_Gen_Overflowge10_emu", "hcalTP_nearL1Jet_Overflowge10_emu", l, False, 1, 0.6, [0, 10])
        plotvariablesame("hcalTP_Overflowge10_emu", "hcalTP_Overflowge10_emu", l, False, 1, 0.6, [0, 10])
        plotvariablesame("hcalTP_nearL1Jet_Overflowge10_emu", "hcalTP_nearL1Jet_Overflowge10_emu", l, True, 1, 1, [0, 10])
        plotvariablesame("hcalTP_nearL1Jet_Gen_Overflowge10_emu", "hcalTP_nearL1Jet_Overflowge10_emu", l, True, 1, 1, [0, 10])
        plotvariablesame("hcalTP_Overflowge10_emu", "hcalTP_Overflowge10_emu", l, True, 1, 1, [0, 10])

#        plotvariablesame("jetET_Leading4", "jetET_Leading4", l, True, 2, [-1, -1])
#        plotvariablesame("jetET_Leading4_Matched", "jetET_Leading4", l, True, 2, [-1, -1])
#        plotvariablesame("njets", "njets", l, True, 1, -1, [-1, -1])
        plotvariablesame("jetET_L4_Gen", "jetET_L4", l, False, 1, 0.3, [-1, -1])
        plotvariablesame("jetET_L4_Gen_HT120", "jetET_L4_HT120", l, False, 1, 0.3, [-1, -1])
        plotvariablesame("jetET_L4_Gen_HT360", "jetET_L4_HT360", l, False, 1, 0.3, [-1, -1])
        plotHoESame("HovEtotal_1x1_emu_Leading1", l)
        plotHoESame("HovEtotal_1x1_emu_Leading2", l)
        plotHoESame("HovEtotal_1x1_emu_Leading3", l)
        plotHoESame("HovEtotal_1x1_emu_Leading4", l)        
        plotHoESame("HovEtotalLog_1x1_emu", l)
        plotHoEGenSame("HovEtotal_1x1_emu_GenMatchedJets", "HovEtotal_1x1_emu", l)
        plotHoEGenSame("HovEtotal_1x1_emu_GenMatchedJets_Barrel", "HovEtotal_1x1_emu_Barrel", l)
        plotHoEGenSame("HovEtotal_1x1_emu_GenMatchedJets_Endcap", "HovEtotal_1x1_emu_Endcap", l)
        plotHoESame("HovEtotal_3x3_emu_Leading1", l)
        plotHoESame("HovEtotal_3x3_emu_Leading2", l)
        plotHoESame("HovEtotal_3x3_emu_Leading3", l)
        plotHoESame("HovEtotal_3x3_emu_Leading4", l)
        plotHoESame("HovEtotalLog_3x3_emu", l)
        plotHoEGenSame("HovEtotal_3x3_emu_GenMatchedJets", "HovEtotal_3x3_emu", l)
        plotHoEGenSame("HovEtotal_3x3_emu_GenMatchedJets_Barrel", "HovEtotal_3x3_emu_Barrel", l)
        plotHoEGenSame("HovEtotal_3x3_emu_GenMatchedJets_Endcap", "HovEtotal_3x3_emu_Endcap", l)

        plotHoESame("HovEtotal_1x1_emu", l)
        plotHoESame("HovEtotal_3x3_emu", l)
        plotHoESame("HovEtotal_1x1_emu_Barrel", l)
        plotHoESame("HovEtotal_3x3_emu_Barrel", l)
        plotHoESame("HovEtotal_1x1_emu_Endcap", l)
        plotHoESame("HovEtotal_3x3_emu_Endcap", l)

#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth1_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth2_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth3_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth4_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth5_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth6_Barrel_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth7_Barrel_emu", l)

#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth1_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth1_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth2_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth2_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth3_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth3_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth4_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth4_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth5_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth5_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth6_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth6_Barrel_emu",  l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth7_Gen_Barrel_emu", "HovEtotal_3x3_ET_Depth7_Barrel_emu",  l)

 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth1_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth2_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth3_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth4_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth5_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth6_Endcap_emu", l)
 #       plotHoEvsDepthSame("HovEtotal_3x3_ET_Depth7_Endcap_emu", l)

 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth1_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth1_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth2_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth2_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth3_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth3_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth4_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth4_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth5_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth5_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth6_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth6_Endcap_emu",  l)
 #       plotHoEvsDepthGenSame("HovEtotal_3x3_ET_Depth7_Gen_Endcap_emu", "HovEtotal_3x3_ET_Depth7_Endcap_emu",  l)
   #     plotDepthRatiosame("energyDepth_Ratio_Barrel_D4", l, 5)
   #     plotDepthRatiosame("energyDepth_Ratio_Barrel_D3_D4", l, 5)
   #     plotDepthRatiosame("energyDepth_Ratio_Endcap_D6_D7", l, 5)
   #     plotDepthRatiosame("energyDepth_Ratio_Endcap_D4_D7", l, 5)
   #     plotDepthRatioGensame("energyDepth_Ratio_Gen_Barrel_D4", "energyDepth_Ratio_Barrel_D4",  l, 5)
   #     plotDepthRatioGensame("energyDepth_Ratio_Gen_Barrel_D3_D4", "energyDepth_Ratio_Barrel_D3_D4", l, 5)
   #     plotDepthRatioGensame("energyDepth_Ratio_Gen_Endcap_D6_D7", "energyDepth_Ratio_Endcap_D6_D7", l, 5)
   #     plotDepthRatioGensame("energyDepth_Ratio_Gen_Endcap_D4_D7", "energyDepth_Ratio_Endcap_D4_D7", l, 5)

#        plotDepthRatiosame("energyDepth_NTPs_Barrel_D3_D4_ge06", l, 1)
#        plotDepthRatiosame("energyDepth_NTPs_Barrel_D4_ge06", l, 1)
#        plotDepthRatiosame("energyDepth_NTPs_Endcap_D4_D7_ge045", l, 1)
#        plotDepthRatioGensame("energyDepth_NTPs_Gen_Barrel_D3_D4_ge06", "energyDepth_NTPs_Barrel_D3_D4_ge06", l, 1)
#        plotDepthRatioGensame("energyDepth_NTPs_Gen_Barrel_D4_ge06", "energyDepth_NTPs_Barrel_D4_ge06", l, 1)
#E        plotDepthRatioGensame("energyDepth_NTPs_Gen_Endcap_D4_D7_ge045", "energyDepth_NTPs_Endcap_D4_D7_ge045", l, 1)

  #      plotDepthRatiosame("energyDepth_fracTPs_Barrel_D3_D4_ge06", l, 1)
  #      plotDepthRatiosame("energyDepth_fracTPs_Endcap_D4_D7_ge045", l, 1)
  #      plotDepthRatioGensame("energyDepth_fracTPs_Gen_Barrel_D3_D4_ge06", "energyDepth_fracTPs_Barrel_D3_D4_ge06", l, 1)
  #      plotDepthRatioGensame("energyDepth_fracTPs_Gen_Endcap_D4_D7_ge045", "energyDepth_fracTPs_Endcap_D4_D7_ge045", l, 1)

#        plotDepthEnergyGensame("energyDepth_Depth1_Gen_Barrel_emu", "energyDepth_Depth1_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth2_Gen_Barrel_emu", "energyDepth_Depth2_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth3_Gen_Barrel_emu", "energyDepth_Depth3_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth4_Gen_Barrel_emu", "energyDepth_Depth4_Barrel_emu",  l, 1)

#        plotDepthEnergyGensame("energyDepth_Depth1_Gen_Endcap_emu", "energyDepth_Depth1_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth2_Gen_Endcap_emu", "energyDepth_Depth2_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth3_Gen_Endcap_emu", "energyDepth_Depth3_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth4_Gen_Endcap_emu", "energyDepth_Depth4_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth5_Gen_Endcap_emu", "energyDepth_Depth5_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth6_Gen_Endcap_emu", "energyDepth_Depth6_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth7_Gen_Endcap_emu", "energyDepth_Depth7_Endcap_emu",  l, 1)

#        plotDepthEnergyGensame("energyDepth_Depth1_Barrel_emu", "energyDepth_Depth1_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth2_Barrel_emu", "energyDepth_Depth2_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth3_Barrel_emu", "energyDepth_Depth3_Barrel_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth4_Barrel_emu", "energyDepth_Depth4_Barrel_emu",  l, 1)

#        plotDepthEnergyGensame("energyDepth_Depth1_Endcap_emu", "energyDepth_Depth1_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth2_Endcap_emu", "energyDepth_Depth2_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth3_Endcap_emu", "energyDepth_Depth3_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth4_Endcap_emu", "energyDepth_Depth4_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth5_Endcap_emu", "energyDepth_Depth5_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth6_Endcap_emu", "energyDepth_Depth6_Endcap_emu",  l, 1)
#        plotDepthEnergyGensame("energyDepth_Depth7_Endcap_emu", "energyDepth_Depth7_Endcap_emu",  l, 1)


#        plotHoEvsDepthSame("HovEtotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu", l)
#        plotHoEvsDepthSame("HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu", l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_DepthRatio_Gen_Barrel_3_4_ov_1_2_emu", "HovEtotal_3x3_DepthRatio_Barrel_3_4_ov_1_2_emu", l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_DepthRatio_Gen_Barrel_4_ov_1_2_3_emu", "HovEtotal_3x3_DepthRatio_Barrel_4_ov_1_2_3_emu", l)
#        plotHoEvsDepthGenSame("HovEtotal_3x3_DepthRatio_Gen_Endcap_4_7_ov_1_3_emu", "HovEtotal_3x3_DepthRatio_Endcap_4_7_ov_1_3_emu", l)

#        plotHoESame("HovEtotal_1x1_emu_AllJets", l)
#        plotHoESame("HovEtotal_3x3_emu_AllJets", l)   

#    plotETRatio("jetETLeading1", "jetET_cutHoE_1x1_Leading1")
#    plotETRatio("jetETLeading1", "jetET_cutHoE_3x3_Leading1")
#    plotETRatio("jetETLeading2", "jetET_cutHoE_1x1_Leading2")
#    plotETRatio("jetETLeading2", "jetET_cutHoE_3x3_Leading2")
#    plotETRatio("jetETLeading3", "jetET_cutHoE_1x1_Leading3")
#    plotETRatio("jetETLeading3", "jetET_cutHoE_3x3_Leading3")
#    plotETRatio("jetETLeading4", "jetET_cutHoE_1x1_Leading4")
#    plotETRatio("jetETLeading4", "jetET_cutHoE_3x3_Leading4")

#signaleff_Rate("eff_signal_HTcut", "htSumRates_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)
#signaleff_Rate("eff_signal_NTP_BarrelD4ge1_1_EndcapD4_D7_ORge1_1", "singleJetRates_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)
#signaleff_Rate("eff_signal_NTP_BarrelD4ge1_1_EndcapD4_D7_ORge1_2", "doubleJetRates_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)
#signaleff_Rate("eff_signal_NTP_BarrelD4ge1_1_EndcapD4_D7_ORge1_3", "tripleJetRates_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)
#signaleff_Rate("eff_signal_NTP_BarrelD4ge1_1_EndcapD4_D7_ORge1_4", "quadJetRates_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)
#signaleff_Rate("eff_signal_NTP_BarrelD4ge1_1_EndcapD4_D7_ORge1_4", "htSumRates_Depth_emu", [LLP_m250_500_dict, LLP_m250_1000_dict, LLP_m350_500_dict, LLP_m350_1000_dict], NuGun_dict)



