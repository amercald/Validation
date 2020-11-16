import ROOT

outpath = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
#ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.AddDirectory(0)
ROOT.TProfile.AddDirectory(0)
#ROOT.TProfile.SetDefaultSumw2()
path = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/depth_studies/"
QCD_dict = {"name" : "QCD",
            "color" : ROOT.kBlack,
            "legendlabel" : " QCD"}
LLP500_dict = {"name" : "LLP_500",
               "color" : ROOT.kGreen+2,
               "legendlabel" : "LLP c#tau = 500 mm"}
LLP1000_dict = {"name" : "LLP_1000",
               "color" : ROOT.kBlue,
               "legendlabel" : "LLP c#tau = 1000 mm"}
LLP10000_dict = {"name" : "LLP_10000",
               "color" : ROOT.kRed,
               "legendlabel" : "LLP c#tau = 10000 mm"}
LLP500_350_dict = {"name" : "350_LLP_500",
               "color" : ROOT.kGreen+2,
               "legendlabel" : "LLP c#tau = 500 mm"}
LLP1000_350_dict = {"name" : "350_LLP_1000",
               "color" : ROOT.kBlue,
               "legendlabel" : "LLP c#tau = 1000 mm"}

file_list = [LLP500_dict, LLP1000_dict, QCD_dict]#, LLP10000_dict]
#file_list = [LLP500_dict, QCD_dict]
file_list_gen = [LLP500_dict, LLP1000_dict]
file_list_350 = [LLP500_350_dict, LLP1000_350_dict, QCD_dict]
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
#        print("using rates_depth_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_depth_"+filename["name"]+".root")
        hist = file.Get(histname)
        profilehist = hist.ProfileX()
        profilehist.SetLineColor(filename["color"])
        profilehist.SetMarkerColor(filename["color"])
        profilehist.SetMarkerStyle(20)
        profilehist.SetMarkerSize(1.7)
        profilehist.SetLineWidth(2)
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
    elif currentlist == file_list_350 : mass = "350"
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
#        print("using rates_depth_"+filename["nam+".root"+" with color "+str(filename["color"]))
        file = ROOT.TFile.Open(path+"rates_depth_"+filename["name"]+".root")
        if filename["name"] == "QCD":
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
        profilehist.GetYaxis().SetTitle("TP Energy Fraction")
        profilehist.GetYaxis().SetRangeUser(0, 1)
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
        file = ROOT.TFile.Open(path+"rates_depth_"+filename["name"]+".root")
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
        file = ROOT.TFile.Open(path+"rates_depth_"+filename["name"]+".root")
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

    for h in hist_list:
        plotGenDepthProfile(file_list, h, "energyDepth_Barrel", "energyDepth_Endcap")
    for g in hist_list_genmatch:
        plotGenDepthProfile(file_list, g, "energyDepth_Barrel", "energyDepth_Endcap")
        plotGenDepthProfile(file_list_350, g, "energyDepth_Barrel", "energyDepth_Endcap")

    plotTPenergy(file_list, "hcalTP_emu")
    plotTPenergy(file_list_350, "hcalTP_emu")
    plotBetaGamma(file_list, "betagammaLLP")
    plotBetaGamma(file_list_350, "betagammaLLP")

