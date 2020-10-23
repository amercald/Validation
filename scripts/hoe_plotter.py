import ROOT

outpath = "/afs/cern.ch/user/a/amercald/private/HCAL/Validation_10_6_X/CMSSW_10_6_0/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TProfile.SetDefaultSumw2()
path = "/afs/cern.ch/user/a/amercald/private/HCAL/Validation_10_6_X/CMSSW_10_6_0/src/HcalTrigger/Validation/HoE_studies/"
variableplots = ["jetEta", "jetET", "jetEtaLeading1", "jetEtaLeading2", "jetEtaLeading3", "jetEtaLeading4", "njets", "jetETLeading1", "jetETLeading2", "jetETLeading3", "jetETLeading4"]
#file = ROOT.TFile.Open(path+"NuGunRates_def.root")
#file_list = ["LLP", "QCD", "NuGun"]
file_list = ["LLP_10000mm", "LLP_1000mm", "LLP_500mm", "QCD"]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kBlack]
hist_dict = {}

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
    if filename == "NuGun":
        hist.GetXaxis().SetRangeUser(0, 100)
        hist.GetYaxis().SetRangeUser(0, 100)
    
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
    #hist.SetStats(0)
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
    legend = ROOT.TLegend(0.15, 0.77, 0.4, 0.92)
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
        legendentr = file_list[f]
        if file_list[f][:3] == "LLP":
            legendentr = "LLP c#tau = "+file_list[f][4:]
        legend.AddEntry(hist, legendentr, "l")
        if f == 0:
            hist.Draw("h")    
        else:
            hist.Draw("h same")
        hist.GetYaxis().SetRangeUser(1, 120)
    legend.Draw("same")
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
    for var in variableplots:
        print("plotting variable "+var)
        c5 = ROOT.TCanvas("%s"%(var), "%s"%(var), XCANVAS, YCANVAS);
        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()
        if not "Eta" in var:
            ROOT.gPad.SetLogy()
#        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
        legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
        
        for f in range(len(file_list)):
            print("opening "+path+"rates_hoe_"+file_list[f]+".root")
            file = ROOT.TFile.Open(path+"rates_hoe_"+file_list[f]+".root")
            varhist = file.Get(var)
            varhist.SetStats(0)
            xtitle = ""
            if "jetEta" in var:
                xtitle = "#eta"
            elif "jetET" in var:
                xtitle = "E_{T}"
            elif "njets" in var:
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
            legendentr = file_list[f]
            if file_list[f][:3] == "LLP":
                legendentr = "LLP c#tau = "+file_list[f][4:]
            legend.AddEntry(varhist, legendentr, "l")
        legend.Draw("same")
        c5.SaveAs(outpath+var+"_All.pdf")
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
        file = ROOT.TFile.Open(path+"rates_hoe_"+file_list[f]+".root")
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
    file = ROOT.TFile.Open(path+"rates_hoe_"+filename+".root")
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
    c1.SaveAs(outpath+histname+"_2D_"+filename+".pdf")
    del c1
    #Draw Profile Histogram on top
    c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
    #ROOT.gPad.SetLogz()
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
    c1.SaveAs(outpath+histname+"_Profile_"+filename+".pdf")

    

    if not histname in hist_dict.keys():
        hist_dict[histname] = {}
    
    for xbin in range(hist.GetNbinsX()):
        if not xbin in hist_dict[histname].keys():
            hist_dict[histname][xbin] = []
        binSize = 50
        histETbin = hist.ProjectionY("",xbin, xbin+1)
        histETbin.SetTitle("H/(H+E) ET from"+str(xbin*binSize)+" - "+str((xbin+1)*binSize)+" "+histname+" "+filename)
        histETbin.GetXaxis().SetTitle("H/(H+E)")
        hist_dict[histname][xbin].append(histETbin)
        
def draw_hist_dict():
    
    for histname in hist_dict.keys():
        for ETbin in hist_dict[histname].keys():
            c1 = ROOT.TCanvas("%s"%(histname), "%s"%(histname), XCANVAS, YCANVAS);
            ROOT.gPad.SetTopMargin(magicMargins["T"])
            ROOT.gPad.SetBottomMargin(magicMargins["B"])
            ROOT.gPad.SetLeftMargin(magicMargins["L"])
            ROOT.gPad.SetRightMargin(magicMargins["R"])
            ROOT.gPad.SetGridy()
            ROOT.gPad.SetTicks()
            ROOT.gPad.SetLogy()
            legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
            legend.SetFillStyle(0)
            binSize = 50

            for histnum in range(4):
                ihist = hist_dict[histname][ETbin][histnum]
                ihist.SetMarkerStyle(20)
                ihist.SetMarkerSize(1.7)
                ihist.SetLineWidth(3)
                ihist.SetMarkerColor(colors[histnum])
                ihist.SetLineColor(colors[histnum])
                ihist.SetStats(0)
                ihist.Draw("h same")
                legend.AddEntry(ihist, file_list[histnum], "l")
            legend.Draw("same")
            c1.SaveAs(outpath+histname+"_ET_"+str(ETbin*binSize)+"_"+str((ETbin+1)*binSize)+"_All.pdf")
            del c1

if __name__ == "__main__":
    for f in file_list:
        #plotHEEnergy("HEEnergytotal_1x1_emu_AllJet", f)
        #plotHEEnergy("HEEnergytotal_3x3_emu_AllJet", f)
        #plotvariable(f)
        #plotHoE("HovEtotal_1x1_emu_AllJets", f)
        #plotHoE("HovEtotal_3x3_emu_AllJets", f)
        #plotHoE("HovEtotal_1x1_emu", f)
        #plotHoE("HovEtotal_3x3_emu", f)
        #plotHEEnergy("HEEnergytotal_1x1_emu_Leading1", f)
        #plotHEEnergy("HEEnergytotal_3x3_emu_Leading1", f)        
        #plotHEEnergy("HEEnergytotal_1x1_emu_Leading2", f)
        #plotHEEnergy("HEEnergytotal_3x3_emu_Leading2", f)        
        #plotHEEnergy("HEEnergytotal_1x1_emu_Leading3", f)
        #plotHEEnergy("HEEnergytotal_3x3_emu_Leading3", f)
        #plotHEEnergy("HEEnergytotal_1x1_emu_Leading4", f)
        #plotHEEnergy("HEEnergytotal_3x3_emu_Leading4", f)                
        
        plotHoEvsET("HovEtotal_1x1_ET_emu_Leading1", f)
        plotHoEvsET("HovEtotal_3x3_ET_emu_Leading1", f)        
        plotHoEvsET("HovEtotal_1x1_ET_emu_Leading2", f)
        plotHoEvsET("HovEtotal_3x3_ET_emu_Leading2", f)        
        plotHoEvsET("HovEtotal_1x1_ET_emu_Leading3", f)
        plotHoEvsET("HovEtotal_3x3_ET_emu_Leading3", f)
        plotHoEvsET("HovEtotal_1x1_ET_emu_Leading4", f)
        plotHoEvsET("HovEtotal_3x3_ET_emu_Leading4", f)
    
    draw_hist_dict()
        
    plotvariablesame()
    plotHoESame("HovEtotal_1x1_emu_AllJets")
    plotHoESame("HovEtotal_3x3_emu_AllJets")
    plotHoESame("HovEtotal_1x1_emu_Leading1")
    plotHoESame("HovEtotal_1x1_emu_Leading2")
    plotHoESame("HovEtotal_1x1_emu_Leading3")
    plotHoESame("HovEtotal_1x1_emu_Leading4")
    plotHoESame("HovEtotal_3x3_emu_Leading1")
    plotHoESame("HovEtotal_3x3_emu_Leading2")
    plotHoESame("HovEtotal_3x3_emu_Leading3")
    plotHoESame("HovEtotal_3x3_emu_Leading4")

    plotETRatio("jetETLeading1", "jetET_cutHoE_1x1_Leading1")
    plotETRatio("jetETLeading1", "jetET_cutHoE_3x3_Leading1")
    plotETRatio("jetETLeading2", "jetET_cutHoE_1x1_Leading2")
    plotETRatio("jetETLeading2", "jetET_cutHoE_3x3_Leading2")
    plotETRatio("jetETLeading3", "jetET_cutHoE_1x1_Leading3")
    plotETRatio("jetETLeading3", "jetET_cutHoE_3x3_Leading3")
    plotETRatio("jetETLeading4", "jetET_cutHoE_1x1_Leading4")
    plotETRatio("jetETLeading4", "jetET_cutHoE_3x3_Leading4")



