import ROOT

outpath = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/plots/"
magicMargins = {"T" : 0.07, "B" : 0.11, "L" : 0.11, "R" : 0.13}
XCANVAS = 2400; YCANVAS = 2400
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TProfile.SetDefaultSumw2()
path = "/afs/cern.ch/user/a/amercald/private/HCAL/test/g14_merge/CMSSW_10_6_0/src/HcalTrigger/Validation/result_rates/"
#variableplots = ["jetEta", "jetET", "jetEtaLeading1", "jetEtaLeading2", "jetEtaLeading3", "jetEtaLeading4", "njets", "jetETLeading1", "jetETLeading2", "jetETLeading3", "jetETLeading4"]
variableplots = ["hJet1x1ov5x5", "DeltaRLLP", "nGenParticles"]

QCD_dict = {"filename" : "QCD", "color" : ROOT.kBlack, "legendlabel" : "QCD"}
NuGun_dict = {"filename" : "NuGun", "color" : ROOT.kBlack+1, "legendlabel" : "NuGun"}

LLP_m1000_500_dict = {"filename" : "LLP_MH1000_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
LLP_m1000_1000_dict = {"filename" : "LLP_MH1000_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}

LLP_m350_500_dict = {"filename" : "LLP_MH350_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
LLP_m350_1000_dict = {"filename" : "LLP_MH350_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}

LLP_m250_500_dict = {"filename" : "LLP_MH250_Ctau500", "color" : ROOT.kRed, "legendlabel" : "LLP c#tau = 500 mm"}
LLP_m250_1000_dict = {"filename" : "LLP_MH250_Ctau1000", "color" : ROOT.kBlue, "legendlabel" : "LLP c#tau = 1000 mm"}

file_list_1000 = { "list" : [QCD_dict, LLP_m1000_500_dict, LLP_m1000_1000_dict], "name" : "MH_1000"}
file_list_350 = { "list" : [QCD_dict, LLP_m350_500_dict, LLP_m350_1000_dict], "name" : "MH_350"}
file_list_250 = { "list" : [QCD_dict, LLP_m250_500_dict, LLP_m250_1000_dict], "name" : "MH_250"}

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
        if "qcd" in f["filename"].lower():
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
        hist.GetYaxis().SetRangeUser(1, 10**4)
        #hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        if hcounter == 0:
            hist.Draw("h")    
        else:
            hist.Draw("h same")
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
        hist.GetYaxis().SetRangeUser(1, 10**4)
       # hist.GetXaxis().SetRangeUser(0, 1)
#        hist.Rebin(2)
        legendentr = f["legendlabel"]
        legend.AddEntry(hist, legendentr, "l")
        if hcounter == 0:
            hist.Draw("h")    
        else:
            hist.Draw("h same")
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

def plotvariablesame(file_list_dict):
    file_list =file_list_dict["list"]
    mh = file_list_dict["name"]

    for var in variableplots:
        print("plotting variable "+var)
        c5 = ROOT.TCanvas("%s"%(var), "%s"%(var), XCANVAS, YCANVAS);
        ROOT.gPad.SetTopMargin(magicMargins["T"])
        ROOT.gPad.SetBottomMargin(magicMargins["B"])
        ROOT.gPad.SetLeftMargin(magicMargins["L"])
        ROOT.gPad.SetRightMargin(magicMargins["R"])
        #ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks()
        if not "Eta" in var:
            ROOT.gPad.SetLogy()
#        legend = ROOT.TLegend(0.12, 0.82, 0.55, 0.92)
        legend = ROOT.TLegend(0.65, 0.77, 0.9, 0.92)
        
        hcounter = 0
        for f in file_list:
            file = ROOT.TFile.Open(path+"rates_"+f["filename"]+".root")
            varhist = file.Get(var)
            varhist.SetStats(0)
            xtitle = ""
            if "jetEta" in var:
                xtitle = "#eta"
            elif "jetET" in var:
                xtitle = "E_{T}"
            elif "njets" in var:
                xtitle = "N_{jets}"
            elif "hJet1x1ov5x5" in var:
                xtitle = "E_{1#times1} / E_{5#times5}"
                varhist.SetTitle("")
            varhist.SetMarkerStyle(20)
            varhist.SetMarkerSize(1.7)
            varhist.SetLineWidth(3)
            varhist.SetMarkerColor(f["color"])
            varhist.SetLineColor(f["color"])
            varhist.GetXaxis().SetTitle(xtitle)
            varhist.Rebin(2)
            if hcounter == 0:
                varhist.Draw("h")
            else:
                varhist.Draw("h same")
            legendentr = f["legendlabel"]
            legend.AddEntry(varhist, legendentr, "l")
            hcounter += 1
        legend.Draw("same")
        c5.SaveAs(outpath+var+"_"+mh+"_All.pdf")
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
    file_list_list = [file_list_1000, file_list_350, file_list_250]

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
            
#    draw_hist_dict()
    for l in file_list_list:
        plotvariablesame(l)
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



