import ROOT

outpath = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 3000; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
#ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.AddDirectory(0)
ROOT.TProfile.AddDirectory(0)
#ROOT.TProfile.SetDefaultSumw2()
path = "/uscms/home/amercald/nobackup/HCAL/CMSSW_11_0_2/src/HcalTrigger/Validation/result_rates/"
QCD_dict = {"name" : "QCD",
            "color" : ROOT.kBlack,
            "legendlabel" : " QCD"}
NuGun_dict = {"name" : "RelValNuGun",
            "color" : ROOT.kBlack,
            "legendlabel" : "PU Jets"}
LLP500_dict = {"name" : "LLP_MH1000_Ctau500",
               "color" : ROOT.kGreen+2,
               "legendlabel" : "LLP Jets c#tau = 0.5m"}
LLP1000_dict = {"name" : "LLP_MH1000_Ctau1000",
               "color" : ROOT.kBlue,
               "legendlabel" : "LLP c#tau = 1000 mm"}
LLP10000_dict = {"name" : "LLP_10000",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP c#tau = 10000 mm"}
LLP500_350_dict = {"name" : "LLP_MH350_Ctau500",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP Jets, c#tau = 0.5m"}
LLP1000_350_dict = {"name" : "LLP_MH350_Ctau1000",
               "color" : ROOT.kMagenta+2,
               "legendlabel" : "LLP Jets, c#tau = 1m"}
LLP500_250_dict = {"name" : "LLP_MH250_Ctau500",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP Jets, c#tau = 0.5m"}
LLP1000_250_dict = {"name" : "LLP_MH250_Ctau1000",
               "color" : ROOT.kMagenta+2,
               "legendlabel" : "LLP Jets, c#tau = 1m"}
LLP500_350_80_dict = {"name" : "LLP_MH350_80_Ctau500",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP Jets, c#tau = 0.5m"}
LLP1000_350_80_dict = {"name" : "LLP_MH350_80_Ctau1000",
               "color" : ROOT.kMagenta+2,
               "legendlabel" : "LLP Jets, c#tau = 1m"}
LLP500_250_60_dict = {"name" : "LLP_MH250_60_Ctau500",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP Jets, c#tau = 0.5m"}
LLP1000_250_60_dict = {"name" : "LLP_MH250_60_Ctau1000",
               "color" : ROOT.kMagenta+2,
               "legendlabel" : "LLP Jets, c#tau = 1m"}



file_list = [LLP500_dict, LLP1000_dict, QCD_dict]#, LLP10000_dict]
#file_list = [LLP500_dict, QCD_dict]
file_list_gen = [LLP500_dict, LLP1000_dict]
file_list_350 = [LLP500_350_dict, LLP1000_350_dict, QCD_dict]
file_list_250 = [LLP500_250_dict, LLP1000_250_dict, QCD_dict]
file_list_350_NuGun = [LLP500_350_80_dict, LLP1000_350_80_dict, NuGun_dict]
file_list_250_NuGun = [LLP500_250_60_dict, LLP1000_250_60_dict, NuGun_dict]

#file_list_350= [LLP500_350_dict, QCD_dict]

def plotDepthProfile(histname):
    
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
#    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    hist_list = []
    print("plotting "+histname)

    for filename in file_list:
#        print("using rates_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_"+filename["name"]+".root")
        hist = file.Get(histname)
        profilehist = hist.ProfileX()
        profilehist.SetLineColor(filename["color"])
        profilehist.SetMarkerColor(filename["color"])
        profilehist.SetMarkerStyle(20)
        profilehist.SetMarkerSize(1.7)
        profilehist.SetLineWidth(3)
        profilehist.SetStats(0)
        profilehist.GetXaxis().SetTitle("HCAL Depth")
        profilehist.GetYaxis().SetTitle("TP Energy Fraction")
        legend.AddEntry(profilehist, filename["legendlabel"], "l")
        hist_list.append(profilehist)
    hcounter = 0   
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("h")
        else:
            ihist.Draw("h same")
        hcounter += 1
        legend.Draw("same")
    
    c1.SaveAs(outpath+histname+"_profile_All.pdf")
    del c1

def plotGenDepthProfile(currentlist, histname, bghistname_barrel, bghistname_endcap):
    
    if currentlist == file_list: mass = "1000"
    elif currentlist == file_list_350 or currentlist == file_list_350_NuGun: 
        mass = "350"
        mass_LLP = "80"
    elif currentlist == file_list_250 or currentlist == file_list_250_NuGun: 
        mass = "250"
        mass_LLP = "60"
    else: mass = "UNKNOWN_MASS"

    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
 #   ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
#    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    legend = ROOT.TLegend(0.4, 0.55, 0.8, 0.92)
    hist_list = []
    print("plotting "+histname)

    for filename in currentlist:
#        print("using rates_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_"+filename["name"]+".root")
        fill = False
        if filename["name"] == "QCD" or filename["name"] == "RelValNuGun":
            if "Endcap" in histname: adjustedhistname = bghistname_endcap
            elif "Barrel" in histname: adjustedhistname = bghistname_barrel
        else:
            fill = True
            adjustedhistname = histname
        print("file: "+filename["name"]+" histname: "+adjustedhistname)
        hist = file.Get(adjustedhistname)
        profilehist = hist.ProfileX()
        profilehist.SetLineColor(filename["color"])
#        profilehist.SetMarkerColor(filename["color"])
#        profilehist.SetMarkerStyle(21)
#        profilehist.SetMarkerSize(1.7)
        profilehist.SetLineWidth(3)
        if fill:
            profilehist.SetFillStyle(1001)
#            profilehist.SetFillColorAlpha(ROOT.kMagenta-10, 0.35)
        else:
            profilehist.SetFillStyle(1001)
            profilehist.SetFillColorAlpha(16, 0.35)
        profilehist.SetStats(0)
        profilehist.GetXaxis().SetTitle("HCAL Depth")
        profilehist.GetYaxis().SetTitle("TP Energy Fraction")
        profilehist.GetYaxis().SetRangeUser(0, 1)
        profilehist.SetTitle("M_{H} = "+mass+" GeV, M_{X} = "+mass_LLP+" GeV")
        legend.AddEntry(profilehist, filename["legendlabel"], "l")
        hist_list.append(profilehist)
    hcounter = 0   
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("h")
        else:
            ihist.Draw("h same")
        hcounter += 1
        legend.SetTextSize(0.045)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetHeader("#bf{CMS #scale[0.8]{#it{Preliminary}}}")
        legend.Draw("h same")
    
    
    c1.SaveAs(outpath+histname+"_"+mass+"_profile_All.pdf")
    
    fileOUT = ROOT.TFile.Open(histname+"_"+mass+".root", "RECREATE")
    c1.Write()
    del c1

def plotGenDepthEProfile(currentlist, histname, bghistname_barrel, bghistname_endcap):
    
    if currentlist == file_list: mass = "1000"
    elif currentlist == file_list_350 or currentlist == file_list_350_NuGun: mass = "350"
    elif currentlist == file_list_250 or currentlist == file_list_250_NuGun: mass = "250"
    else: mass = "UNKNOWN_MASS"

    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
 #   ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    hist_list = []
    print("plotting "+histname)

    for filename in currentlist:
#        print("using rates_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_"+filename["name"]+".root")
        if filename["name"] == "QCD" or filename["name"] == "RelValNuGun":
            if "Endcap" in histname: adjustedhistname = bghistname_endcap
            elif "Barrel" in histname: adjustedhistname = bghistname_barrel
        else:
            adjustedhistname = histname
        print("file: "+filename["name"]+" histname: "+adjustedhistname)
        hist = file.Get(adjustedhistname)
        profilehist = hist.ProfileX()
        profilehist.SetLineColor(filename["color"])
        profilehist.SetMarkerColor(filename["color"])
        profilehist.SetMarkerStyle(21)
        profilehist.SetMarkerSize(1.7)
        profilehist.SetLineWidth(2)
        profilehist.SetStats(0)
        profilehist.GetXaxis().SetTitle("HCAL Depth")
        profilehist.GetYaxis().SetTitle("Energy (GeV)")
        profilehist.GetYaxis().SetRangeUser(0, 10)
        legend.AddEntry(profilehist, filename["legendlabel"], "l")
        hist_list.append(profilehist)
    hcounter = 0   
    for ihist in hist_list:
        if hcounter == 0:
            ihist.Draw("h")
        else:
            ihist.Draw("h same")
        hcounter += 1
        legend.Draw("h same")
        
    c1.SaveAs(outpath+histname+"_"+mass+"_profile_All.pdf")
    del c1


def plotTPenergy(currentlist, histname):

    if currentlist == file_list: mass = "1000"
    elif currentlist == file_list_350 : mass = "350"
    else: mass = "UNKNOWN_MASS"

    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
 #   ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogy()
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    histcounter = 0
    for filename in currentlist:
        file = ROOT.TFile.Open(path+"rates_"+filename["name"]+".root")
        hist = file.Get(histname)
        hist.SetLineColor(filename["color"])
        hist.SetMarkerColor(filename["color"])
        hist.SetMarkerStyle(21)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(2)
        hist.SetStats(0)
        hist.Rebin(2)
        legend.AddEntry(hist, filename["legendlabel"], "l")
        if histcounter == 0:
            hist.Draw("h")
        else:
            hist.Draw("h same")
        histcounter += 1
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mass+"_All.pdf")

def plotBetaGamma(currentlist, histname):

    if currentlist == file_list: mass = "1000"
    elif currentlist == file_list_350 : mass = "350"
    else: mass = "UNKNOWN_MASS"

    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
 #   ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogy()
    legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
    histcounter = 0
    for filename in currentlist:
        if filename == QCD_dict: continue
        file = ROOT.TFile.Open(path+"rates_"+filename["name"]+".root")
        hist = file.Get(histname)
        hist.SetLineColor(filename["color"])
        hist.SetMarkerColor(filename["color"])
        hist.SetMarkerStyle(21)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(2)
        hist.SetStats(0)
        hist.Rebin(2)
        legend.AddEntry(hist, filename["legendlabel"], "l")
        if histcounter == 0:
            hist.Draw("h")
        else:
            hist.Draw("h same")
        histcounter += 1
    legend.Draw("same")
    c1.SaveAs(outpath+histname+"_"+mass+"_All.pdf")
                    
if __name__ == "__main__":
    
    hist_list = ["energyDepth_Jet1_Barrel", "energyDepth_Jet2_Barrel", "energyDepth_Jet3_Barrel", "energyDepth_Jet4_Barrel", "energyDepth_Jet1_Endcap", "energyDepth_Jet2_Endcap", "energyDepth_Jet3_Endcap", "energyDepth_Jet4_Endcap", "energyDepth_Barrel", "energyDepth_Endcap"]
    hist_list_genmatch = ["energyDepth_genMatchInclusive_Barrel", "energyDepth_genMatchInclusive_Endcap", "energyDepth_genMatchTP_Endcap", "energyDepth_genMatchTP_Barrel"]

    #for h in hist_list:
        #plotGenDepthProfile(file_list, h, "energyDepth_Barrel", "energyDepth_Endcap")
    for g in hist_list_genmatch:
        #plotGenDepthProfile(file_list, g, "energyDepth_Barrel", "energyDepth_Endcap")
        #plotGenDepthProfile(file_list_350, g, "energyDepth_Barrel", "energyDepth_Endcap")
        #plotGenDepthProfile(file_list_250, g, "energyDepth_Barrel", "energyDepth_Endcap")
        plotGenDepthProfile(file_list_350_NuGun, g, "energyDepth_Barrel", "energyDepth_Endcap")
        plotGenDepthProfile(file_list_250_NuGun, g, "energyDepth_Barrel", "energyDepth_Endcap")

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_L1_Endcap", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_L1_Endcap", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_L1_Barrel", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_L1_Barrel", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_L1_Endcap", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_L1_Endcap", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_L1_Barrel", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_L1_Barrel", "energyDepth_L1_Barrel", "energyDepth_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_HT120_L1_Endcap", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_HT120_L1_Endcap", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_LowRatio_L1_Endcap", "energyDepth_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_LowRatio_L1_Endcap", "energyDepth_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Barrel", "energyDepth_LowRatio_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_Endcap", "energyDepth_TPge5_Barrel", "energyDepth_TPge5_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_Endcap", "energyDepth_TPge5_Barrel", "energyDepth_TPge5_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_Barrel", "energyDepth_TPge5_Barrel", "energyDepth_TPge5_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_Barrel", "energyDepth_TPge5_Barrel", "energyDepth_TPge5_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HT120_Endcap", "energyDepth_HT120_Barrel", "energyDepth_HT120_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HT120_Endcap", "energyDepth_HT120_Barrel", "energyDepth_HT120_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HT120_Barrel", "energyDepth_HT120_Barrel", "energyDepth_HT120_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HT120_Barrel", "energyDepth_HT120_Barrel", "energyDepth_HT120_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_L1_Endcap", "energyDepth_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_L1_Endcap", "energyDepth_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Barrel", "energyDepth_TPge5_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HT120_L1_Endcap", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HT120_L1_Endcap", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HT120_L1_Barrel", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HT120_L1_Barrel", "energyDepth_HT120_L1_Barrel", "energyDepth_HT120_L1_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_LowRatio_L1_Endcap", "energyDepth_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_LowRatio_L1_Endcap", "energyDepth_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Barrel", "energyDepth_TPge5_LowRatio_L1_Endcap")#
    
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HoEcut_Endcap", "energyDepth_HoEcut_Barrel", "energyDepth_HoEcut_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HoEcut_Endcap", "energyDepth_HoEcut_Barrel", "energyDepth_HoEcut_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_HoEcut_Barrel", "energyDepth_HoEcut_Barrel", "energyDepth_HoEcut_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_HoEcut_Barrel", "energyDepth_HoEcut_Barrel", "energyDepth_HoEcut_Endcap")#

    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap", "energyDepth_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_HoEcut_Endcap", "energyDepth_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Endcap")
    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Endcap")
    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_genMatchInclusive_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Barrel", "energyDepth_TPge5_HoEcut_Endcap")#

    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_HoEcut_Endcap", "energyDepth_TPE_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_HoEcut_Endcap", "energyDepth_TPE_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Endcap")
    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Barrel", "energyDepth_TPE_HoEcut_Endcap")#

    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_TPge5_Endcap", "energyDepth_TPE_TPge5_Barrel", "energyDepth_TPE_TPge5_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_TPge5_Endcap", "energyDepth_TPE_TPge5_Barrel", "energyDepth_TPE_TPge5_Endcap")
    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_TPge5_Barrel", "energyDepth_TPE_TPge5_Barrel", "energyDepth_TPE_TPge5_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_TPge5_Barrel", "energyDepth_TPE_TPge5_Barrel", "energyDepth_TPE_TPge5_Endcap")#

    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_Endcap", "energyDepth_TPE_Barrel", "energyDepth_TPE_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_Endcap", "energyDepth_TPE_Barrel", "energyDepth_TPE_Endcap")
    plotGenDepthEProfile(file_list_350_NuGun, "energyDepth_TPE_genMatchInclusive_Barrel", "energyDepth_TPE_Barrel", "energyDepth_TPE_Endcap")
    plotGenDepthEProfile(file_list_250_NuGun, "energyDepth_TPE_genMatchInclusive_Barrel", "energyDepth_TPE_Barrel", "energyDepth_TPE_Endcap")#


#    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_ET_50_200_genMatchInclusive_Endcap", "energyDepth_ET_50_200_Barrel", "energyDepth_ET_50_200_Endcap")
#    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_ET_50_200_genMatchInclusive_Endcap", "energyDepth_ET_50_200_Barrel", "energyDepth_ET_50_200_Endcap")
#    plotGenDepthProfile(file_list_350_NuGun, "energyDepth_ET_50_200_genMatchInclusive_Barrel", "energyDepth_ET_50_200_Barrel", "energyDepth_ET_50_200_Endcap")
#    plotGenDepthProfile(file_list_250_NuGun, "energyDepth_ET_50_200_genMatchInclusive_Barrel", "energyDepth_ET_50_200_Barrel", "energyDepth_ET_50_200_Endcap")#

    #plotTPenergy(file_list, "hcalTP_emu")
    #plotTPenergy(file_list_350, "hcalTP_emu")
    #plotBetaGamma(file_list, "betagammaLLP")
    #plotBetaGamma(file_list_350, "betagammaLLP")

