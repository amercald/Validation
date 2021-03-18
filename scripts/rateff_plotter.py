import ROOT
from plotters import rootPlotter as rP

myPlotter = rP()

def_360_dict =    {  "effbin" : 2,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 361 }
def_120_dict =    {  "effbin" : 3,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 121 }
def_HTscan_360_dict =    {  "effbin" : 6,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 361 }



#        HoE_TP05_ONLYTAG_dict = { "effbin" : 5,
#                                  "rateshist" : "htSumRates_HoE_TP05_emu",
#                                  "ratebin": 120}
HoE_TP05_OR360_dict = { "effbin" : 4,
                  "rateshist" : "htSumRates_HoE_TP05_ORHT360_emu",
                  "ratebin" : 121 }

HoE_TP1_OR360_dict = { "effbin" : 5,
                 "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                 "ratebin" : 121 }

HoE_TP2_OR360_dict = { "effbin" : 6,
                 "rateshist" : "htSumRates_HoE_TP2_ORHT360_emu",
                 "ratebin" : 121 }

HoE_TP3_OR360_dict = { "effbin" : 7,
                 "rateshist" : "htSumRates_HoE_TP3_ORHT360_emu",
                 "ratebin" : 121 }

HoE_TP4_OR360_dict = { "effbin" : 8,
                 "rateshist" : "htSumRates_HoE_TP4_ORHT360_emu",
                 "ratebin" : 121 }

HoE_TP5_OR360_dict = { "effbin" : 9,
                 "rateshist" : "htSumRates_HoE_TP5_ORHT360_emu",
                 "ratebin" : 121 }

HoE_TP05_dict = { "effbin" : 10,
                  "rateshist" : "htSumRates_HoE_TP05_emu",
                  "ratebin" : 121 }

HoE_TP1_dict = { "effbin" : 11,
                 "rateshist" : "htSumRates_HoE_TP1_emu",
                 "ratebin" : 121 }

HoE_TP2_dict = { "effbin" : 12,
                 "rateshist" : "htSumRates_HoE_TP2_emu",
                 "ratebin" : 121 }

HoE_TP3_dict = { "effbin" : 13,
                 "rateshist" : "htSumRates_HoE_TP3_emu",
                 "ratebin" : 121 }

HoE_TP4_dict = { "effbin" : 14,
                 "rateshist" : "htSumRates_HoE_TP4_emu",
                 "ratebin" : 121 }

HoE_TP5_dict = { "effbin" : 15,
                 "rateshist" : "htSumRates_HoE_TP5_emu",
                 "ratebin" : 121 }

rateff_TPE_OR360_list = [def_360_dict, def_120_dict, HoE_TP05_OR360_dict, HoE_TP1_OR360_dict, HoE_TP2_OR360_dict, HoE_TP3_OR360_dict, HoE_TP4_OR360_dict, HoE_TP5_OR360_dict]
rateff_TPE_list = [def_360_dict, def_120_dict, HoE_TP05_dict, HoE_TP1_dict, HoE_TP2_dict, HoE_TP3_dict, HoE_TP4_dict, HoE_TP5_dict]

HoE_Ratio02_OR360_dict = { "effbin" : 4,
                     "rateshist" : "htSumRates_HoE_Ratio02_ORHT360_emu",
                     "ratebin" : 121}

HoE_Ratio04_OR360_dict = { "effbin" : 5,
                     "rateshist" : "htSumRates_HoE_Ratio04_ORHT360_emu",
                     "ratebin" : 121}

HoE_Ratio06_OR360_dict = { "effbin" : 6,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 121}
HoE_Ratio08_OR360_dict = { "effbin" : 7,
                     "rateshist" : "htSumRates_HoE_Ratio08_ORHT360_emu",
                     "ratebin" : 121}
HoE_Ratio02_dict = { "effbin" : 8,
                     "rateshist" : "htSumRates_HoE_Ratio02_emu",
                     "ratebin" : 121}

HoE_Ratio04_dict = { "effbin" : 9,
                     "rateshist" : "htSumRates_HoE_Ratio04_emu",
                     "ratebin" : 121}

HoE_Ratio06_dict = { "effbin" : 10,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                     "ratebin" : 121}
HoE_Ratio08_dict = { "effbin" : 11,
                     "rateshist" : "htSumRates_HoE_Ratio08_emu",
                     "ratebin" : 121}

rateff_Ratio_OR360_list = [def_360_dict, def_120_dict, HoE_Ratio02_OR360_dict, HoE_Ratio04_OR360_dict, HoE_Ratio06_OR360_dict, HoE_Ratio08_OR360_dict]
rateff_Ratio_list = [def_360_dict, def_120_dict, HoE_Ratio02_dict, HoE_Ratio04_dict, HoE_Ratio06_dict, HoE_Ratio08_dict]

HoE_HTscan_noLLPtag_120_dict = { "effbin" : 2,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 121}

HoE_HTscan_noLLPtag_180_dict = { "effbin" : 3,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 181}

HoE_HTscan_noLLPtag_240_dict = { "effbin" : 4,
                     "rateshist" : "htSumRates_emu",
                        "ratebin" : 241}

HoE_HTscan_noLLPtag_300_dict = { "effbin" : 5,
                     "rateshist" : "htSumRates_emu",
                        "ratebin" : 301}

HoE_HTscan_noLLPtag_360_dict = { "effbin" : 6,
                     "rateshist" : "htSumRates_emu",
                        "ratebin" : 361}


HoE_HTscan_120_dict = { "effbin" : 12,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                     "ratebin" : 121}

HoE_HTscan_180_dict = { "effbin" : 13,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                     "ratebin" : 181}

HoE_HTscan_240_dict = { "effbin" : 14,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 241}

HoE_HTscan_300_dict = { "effbin" : 15,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 301}

HoE_HTscan_360_dict = { "effbin" : 16,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 361}

HoE_HTscan_ORHT360_120_dict = { "effbin" : 7,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 121}

HoE_HTscan_ORHT360_180_dict = { "effbin" : 8,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 181}

HoE_HTscan_ORHT360_240_dict = { "effbin" : 9,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 241}

HoE_HTscan_ORHT360_300_dict = { "effbin" : 10,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 301}

HoE_HTscan_ORHT360_360_dict = { "effbin" : 11,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 361}

rateff_HTscan_noLLPtag_list = [def_HTscan_360_dict, def_120_dict, HoE_HTscan_noLLPtag_120_dict, HoE_HTscan_noLLPtag_120_dict, HoE_HTscan_noLLPtag_180_dict, HoE_HTscan_noLLPtag_240_dict, HoE_HTscan_noLLPtag_300_dict, HoE_HTscan_noLLPtag_360_dict]

rateff_HTscan_list = [def_HTscan_360_dict, def_120_dict, HoE_HTscan_120_dict, HoE_HTscan_180_dict, HoE_HTscan_240_dict, HoE_HTscan_300_dict, HoE_HTscan_360_dict]

rateff_HTscan_ORHT360_list = [def_HTscan_360_dict, def_120_dict, HoE_HTscan_ORHT360_120_dict, HoE_HTscan_ORHT360_180_dict, HoE_HTscan_ORHT360_240_dict, HoE_HTscan_ORHT360_300_dict, HoE_HTscan_ORHT360_360_dict]


def main():
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_ORHT360_Ratio_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_ORHT360_Ratio_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio", legendentry = " HT > 120 with LLP ID")

    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio_Gen", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio_Gen", legendentry = " HT > 120 with LLP ID")

    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_ORHT360_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_ORHT360_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE", legendentry = " HT > 120 with LLP ID")

    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_Gen", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_Gen", legendentry = " HT > 120 with LLP ID")

#    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_noLLPtag", legendentry = " HT > threshold")
#    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_noLLPtag", legendentry = " HT > threshold")

 #   myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_noLLPtag_Gen", legendentry = " HT > threshold")
 #   myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_noLLPtag_Gen", legendentry = " HT > threshold")


    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan", legendentry = " HT > threshold with LLP ID")
    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan", legendentry = " HT > threshold with LLP ID")

    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_Gen", legendentry = " HT > threshold with LLP ID")
    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_Gen", legendentry = " HT > threshold with LLP ID")

    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_ORHT360", legendentry = " HT > threshold with LLP ID")
    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_ORHT360", legendentry = " HT > threshold with LLP ID")

    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_ORHT360_Gen", legendentry = " HT > threshold with LLP ID")
    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_ORHT360_Gen", legendentry = " HT > threshold with LLP ID")

    myPlotter.plotJetRates("singleJetRates", "Single Jet")
    myPlotter.plotJetRates("doubleJetRates", "Double Jet")
    myPlotter.plotJetRates("tripleJetRates", "Triple Jet")
    myPlotter.plotJetRates("quadJetRates", "Quadruple Jet")
    myPlotter.plotJetRates("htSumRates", "L1 H_{T}")

if __name__ == "__main__":

    main()

    
