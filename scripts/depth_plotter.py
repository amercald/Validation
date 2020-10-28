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
file_list = [QCD_dict, LLP500_dict, LLP1000_dict, LLP10000_dict]

def plotDepthProfile(histname):
    
    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);

    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])    
    ROOT.gPad.SetGridy()
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

def basicplotter():

    c1 = ROOT.TCanvas("c1", "c1", XCANVAS, YCANVAS);
    
    file1 = ROOT.TFile.Open(path+"rates_depth_QCD.root")
    hist1 = file1.Get("energyDepth_Jet1_Barrel")
    phist1 = hist1.ProfileX()
    phist1.SetLineColor(ROOT.kRed)
    phist1.Draw("h")
    file2 = ROOT.TFile.Open(path+"rates_depth_LLP_500.root")
#    c1.cd()
    hist2 = file2.Get("energyDepth_Jet1_Barrel")
    phist2 = hist2.ProfileX()
    phist2.SetLineColor(ROOT.kBlue)

    phist2.Draw("h same")
    c1.SaveAs(outpath+"test.pdf")

    
                    
if __name__ == "__main__":
    
    hist_list = ["energyDepth_Jet1_Barrel", "energyDepth_Jet2_Barrel", "energyDepth_Jet3_Barrel", "energyDepth_Jet4_Barrel", "energyDepth_Jet1_Endcap", "energyDepth_Jet2_Endcap", "energyDepth_Jet3_Endcap", "energyDepth_Jet4_Endcap", "energyDepth_Barrel", "energyDepth_Endcap"]

#    basicplotter()

    for h in hist_list:
        plotDepthProfile(h)

