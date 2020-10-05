import ROOT

outpath = "/afs/cern.ch/user/a/amercald/private/HCAL/Validation_10_6_X/CMSSW_10_6_0/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
path = "/afs/cern.ch/user/a/amercald/private/HCAL/Validation_10_6_X/CMSSW_10_6_0/src/HcalTrigger/Validation/HoE_studies/"
variableplots = ["jetEta", "jetET", "njets"]
file = ROOT.TFile.Open(path+"NuGunRates_def.root")
file_list = ["LLP", "QCD", "NuGun"]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]

def plotHEEnergy(histname, filename):
    print("opening "+path+"rates_hoe_"+filename+".root")
    file = ROOT.TFile.Open(path+"rates_hoe_"+filename+".root")
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
    
    hist.Draw("colz")

    c1.SaveAs(outpath+histname+"_2D_"+filename+".pdf")
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
    
    c2.SaveAs(outpath+histname+"_ECAL_"+filename+".pdf")
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
    
    c3.SaveAs(outpath+histname+"_HCAL_"+filename+".pdf")

def plotHoE(histname, filename):
    file = ROOT.TFile.Open(path+"rates_hoe_"+filename+".root")
    hist = file.Get(histname)
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()
    
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.7)
    hist.SetLineWidth(3)
    
    hist.GetXaxis().SetTitle("H/(H+E)")
    hist.SetStats(0)
    hist.Draw("h")    
    c1.SaveAs(outpath+histname+"_"+filename+".pdf")

def plotHoESame(histname):

    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTopMargin(magicMargins["T"])
    ROOT.gPad.SetBottomMargin(magicMargins["B"])
    ROOT.gPad.SetLeftMargin(magicMargins["L"])
    ROOT.gPad.SetRightMargin(magicMargins["R"])
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTicks()

    for f in range(len(file_list)):
        file = ROOT.TFile.Open(path+"rates_hoe_"+file_list[f]+".root")
        hist = file.Get(histname)
    
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1.7)
        hist.SetLineWidth(3)
        hist.SetMarkerColor(colors[f])
        hist.SetLineColor(colors[f])
        hist.GetXaxis().SetTitle("H/(H+E)")
        hist.SetStats(0)
        if f == 0:
            hist.Draw("h")    
        else:
            hist.Draw("h same")
    c1.SaveAs(outpath+histname+"_All.pdf")

    
def plotvariable(filename):
    file = ROOT.TFile.Open(path+"rates_hoe_"+filename+".root")
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

def plotvariablesame():
    file_list = ["LLP", "QCD", "NuGun"]
    for var in variableplots:
        print("plotting variable "+var)
        c5 = ROOT.TCanvas("%s"%(var), "%s"%(var), XCANVAS, YCANVAS);
        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()
        if var != "jetEta":
            ROOT.gPad.SetLogy()
#        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
        legend = ROOT.TLegend(0.7, 0.82, 0.9, 0.92)
        
        for f in range(len(file_list)):
            print("opening "+path+"rates_hoe_"+file_list[f]+".root")
            file = ROOT.TFile.Open(path+"rates_hoe_"+file_list[f]+".root")
            varhist = file.Get(var)
            varhist.SetStats(0)
            xtitle = ""
            if var == "jetEta":
                xtitle = "#eta"
            elif var == "jetET":
                xtitle = "E_{T}"
            elif var == "njets":
                xtitle = "N_{jets}"
            varhist.SetMarkerStyle(20)
            varhist.SetMarkerSize(1.7)
            varhist.SetLineWidth(3)
            varhist.SetMarkerColor(colors[f])
            varhist.SetLineColor(colors[f])
            varhist.GetXaxis().SetTitle(xtitle)
            if f == 0:
                varhist.Draw("h")
            else:
                varhist.Draw("h same")
            legend.AddEntry(varhist, file_list[f], "l")
        legend.Draw("same")
        c5.SaveAs(outpath+var+"_All.pdf")
        del c5
        del legend

    

if __name__ == "__main__":
    for f in file_list:
        plotHEEnergy("HEEnergytotal_1x1_emu_AllJet", f)
        plotHEEnergy("HEEnergytotal_3x3_emu_AllJet", f)
        plotvariable(f)
        plotHoE("HovEtotal_1x1_emu_AllJets", f)
        plotHoE("HovEtotal_3x3_emu_AllJets", f)

    plotvariablesame()
    plotHoESame("HovEtotal_1x1_emu_AllJets")
    plotHoESame("HovEtotal_3x3_emu_AllJets")
