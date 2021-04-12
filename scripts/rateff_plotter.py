import ROOT
from plotters import rootPlotter as rP

myPlotter = rP()


#dictionaries for rate vs eff plots

def_360_dict =    {  "effbin" : 2,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 361, "def" : True}
def_120_dict =    {  "effbin" : 3,
                     "rateshist" : "htSumRates_emu",
                     "ratebin" : 121, "def" : True }

HoE_TP05_OR360_dict = { "effbin" : 4,
                  "rateshist" : "htSumRates_HoE_TP05_ORHT360_emu",
                  "ratebin" : 121 , "def" : False}

HoE_TP1_OR360_dict = { "effbin" : 5,
                 "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP2_OR360_dict = { "effbin" : 6,
                 "rateshist" : "htSumRates_HoE_TP2_ORHT360_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP3_OR360_dict = { "effbin" : 7,
                 "rateshist" : "htSumRates_HoE_TP3_ORHT360_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP4_OR360_dict = { "effbin" : 8,
                 "rateshist" : "htSumRates_HoE_TP4_ORHT360_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP5_OR360_dict = { "effbin" : 9,
                 "rateshist" : "htSumRates_HoE_TP5_ORHT360_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP05_dict = { "effbin" : 10,
                  "rateshist" : "htSumRates_HoE_TP05_emu",
                  "ratebin" : 121 , "def" : False}

HoE_TP1_dict = { "effbin" : 11,
                 "rateshist" : "htSumRates_HoE_TP1_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP2_dict = { "effbin" : 12,
                 "rateshist" : "htSumRates_HoE_TP2_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP3_dict = { "effbin" : 13,
                 "rateshist" : "htSumRates_HoE_TP3_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP4_dict = { "effbin" : 14,
                 "rateshist" : "htSumRates_HoE_TP4_emu",
                 "ratebin" : 121 , "def" : False}

HoE_TP5_dict = { "effbin" : 15,
                 "rateshist" : "htSumRates_HoE_TP5_emu",
                 "ratebin" : 121 , "def" : False}

rateff_TPE_OR360_list = [def_360_dict, def_120_dict, HoE_TP05_OR360_dict, HoE_TP1_OR360_dict, HoE_TP2_OR360_dict, HoE_TP5_OR360_dict]
rateff_TPE_list = [def_360_dict, def_120_dict, HoE_TP05_dict, HoE_TP1_dict, HoE_TP2_dict, HoE_TP5_dict]

HoE_Ratio02_OR360_dict = { "effbin" : 4,
                     "rateshist" : "htSumRates_HoE_Ratio02_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_Ratio04_OR360_dict = { "effbin" : 5,
                     "rateshist" : "htSumRates_HoE_Ratio04_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_Ratio06_OR360_dict = { "effbin" : 6,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 121, "def" : False}
HoE_Ratio08_OR360_dict = { "effbin" : 7,
                     "rateshist" : "htSumRates_HoE_Ratio08_ORHT360_emu",
                     "ratebin" : 121, "def" : False}
HoE_Ratio02_dict = { "effbin" : 8,
                     "rateshist" : "htSumRates_HoE_Ratio02_emu",
                     "ratebin" : 121, "def" : False}

HoE_Ratio04_dict = { "effbin" : 9,
                     "rateshist" : "htSumRates_HoE_Ratio04_emu",
                     "ratebin" : 121, "def" : False}

HoE_Ratio06_dict = { "effbin" : 10,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                     "ratebin" : 121, "def" : False}
HoE_Ratio08_dict = { "effbin" : 11,
                     "rateshist" : "htSumRates_HoE_Ratio08_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP2_Ratio02_OR360_dict = { "effbin" : 12,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio02_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP2_Ratio04_OR360_dict = { "effbin" : 13,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio04_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP2_Ratio06_OR360_dict = { "effbin" : 14,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio06_ORHT360_emu",
                     "ratebin" : 121, "def" : False}
HoE_TP2_Ratio08_OR360_dict = { "effbin" : 15,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio08_ORHT360_emu",
                     "ratebin" : 121, "def" : False}
HoE_TP2_Ratio02_dict = { "effbin" : 16,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio02_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP2_Ratio04_dict = { "effbin" : 17,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio04_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP2_Ratio06_dict = { "effbin" : 18,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio06_emu",
                     "ratebin" : 121, "def" : False}
HoE_TP2_Ratio08_dict = { "effbin" : 19,
                     "rateshist" : "htSumRates_HoE_TP2_Ratio08_emu",
                     "ratebin" : 121, "def" : False}


rateff_Ratio_OR360_list = [def_360_dict, def_120_dict, HoE_Ratio02_OR360_dict, HoE_Ratio04_OR360_dict, HoE_Ratio06_OR360_dict, HoE_Ratio08_OR360_dict]
rateff_Ratio_list = [def_360_dict, def_120_dict, HoE_Ratio02_dict, HoE_Ratio04_dict, HoE_Ratio06_dict, HoE_Ratio08_dict]

rateff_TP2_Ratio_OR360_list = [HoE_TP2_Ratio02_OR360_dict, HoE_TP2_Ratio04_OR360_dict, HoE_TP2_Ratio06_OR360_dict, HoE_TP2_Ratio08_OR360_dict]
rateff_TP2_Ratio_list = [HoE_TP2_Ratio02_dict, HoE_TP2_Ratio04_dict, HoE_TP2_Ratio06_dict, HoE_TP2_Ratio08_dict]

HoE_HTscan_noLLPtag_120_dict = { "effbin" : 2,
                     "rateshist" : "htSumRates_emu",
                                 "ratebin" : 121, "def" : True}

HoE_HTscan_noLLPtag_180_dict = { "effbin" : 3,
                     "rateshist" : "htSumRates_emu",
                                 "ratebin" : 181, "def" : True}

HoE_HTscan_noLLPtag_240_dict = { "effbin" : 4,
                     "rateshist" : "htSumRates_emu",
                                 "ratebin" : 241, "def" : True}

HoE_HTscan_noLLPtag_300_dict = { "effbin" : 5,
                     "rateshist" : "htSumRates_emu",
                                 "ratebin" : 301, "def" : True}

HoE_HTscan_noLLPtag_360_dict = { "effbin" : 6,
                     "rateshist" : "htSumRates_emu",
                                 "ratebin" : 361, "def" : True}


HoE_HTscan_120_dict = { "effbin" : 12,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 121, "def" : False}

HoE_HTscan_180_dict = { "effbin" : 13,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                     "ratebin" : 181, "def" : False}

HoE_HTscan_240_dict = { "effbin" : 14,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 241, "def" : False}

HoE_HTscan_300_dict = { "effbin" : 15,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 301, "def" : False}

HoE_HTscan_360_dict = { "effbin" : 16,
                     "rateshist" : "htSumRates_HoE_Ratio06_emu",
                        "ratebin" : 361, "def" : False}

HoE_HTscan_ORHT360_120_dict = { "effbin" : 7,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_HTscan_ORHT360_180_dict = { "effbin" : 8,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                     "ratebin" : 181, "def" : False}

HoE_HTscan_ORHT360_240_dict = { "effbin" : 9,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 241, "def" : False}

HoE_HTscan_ORHT360_300_dict = { "effbin" : 10,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 301, "def" : False}

HoE_HTscan_ORHT360_360_dict = { "effbin" : 11,
                     "rateshist" : "htSumRates_HoE_Ratio06_ORHT360_emu",
                        "ratebin" : 361, "def" : False}

HoE_TP1_HTscan_ORHT360_120_dict = { "effbin" : 17,
                     "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                     "ratebin" : 121, "def" : False}

HoE_TP1_HTscan_ORHT360_180_dict = { "effbin" : 18,
                     "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                     "ratebin" : 181, "def" : False}

HoE_TP1_HTscan_ORHT360_240_dict = { "effbin" : 19,
                     "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                        "ratebin" : 241, "def" : False}

HoE_TP1_HTscan_ORHT360_300_dict = { "effbin" : 20,
                     "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                        "ratebin" : 301, "def" : False}

HoE_TP1_HTscan_ORHT360_360_dict = { "effbin" : 21,
                     "rateshist" : "htSumRates_HoE_TP1_ORHT360_emu",
                        "ratebin" : 361, "def" : False}


rateff_HTscan_noLLPtag_list = [HoE_HTscan_noLLPtag_120_dict, HoE_HTscan_noLLPtag_180_dict, HoE_HTscan_noLLPtag_240_dict, HoE_HTscan_noLLPtag_300_dict, HoE_HTscan_noLLPtag_360_dict]

rateff_HTscan_list = [HoE_HTscan_noLLPtag_120_dict, HoE_HTscan_noLLPtag_360_dict, HoE_HTscan_120_dict, HoE_HTscan_180_dict, HoE_HTscan_240_dict, HoE_HTscan_300_dict, HoE_HTscan_360_dict]

rateff_HTscan_ORHT360_list = [HoE_HTscan_noLLPtag_120_dict, HoE_HTscan_noLLPtag_360_dict, HoE_HTscan_ORHT360_120_dict, HoE_HTscan_ORHT360_180_dict, HoE_HTscan_ORHT360_240_dict, HoE_HTscan_ORHT360_300_dict, HoE_HTscan_ORHT360_360_dict]

rateff_TP1_HTscan_ORHT360_list = [HoE_TP1_HTscan_ORHT360_120_dict, HoE_TP1_HTscan_ORHT360_180_dict, HoE_TP1_HTscan_ORHT360_240_dict, HoE_TP1_HTscan_ORHT360_300_dict, HoE_TP1_HTscan_ORHT360_360_dict]

HoE_Jetscan_noTag_20_triple_dict = { "effbin" : 2,
                               "rateshist" : "tripleJetRates_emu",
                               "ratebin" : 21, "def" : False}

HoE_Jetscan_noTag_30_triple_dict = { "effbin" : 3,
                               "rateshist" : "tripleJetRates_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_noTag_40_triple_dict = { "effbin" : 4,
                               "rateshist" : "tripleJetRates_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_noTag_50_triple_dict = { "effbin" : 5,
                               "rateshist" : "tripleJetRates_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_noTag_60_triple_dict = { "effbin" : 6,
                               "rateshist" : "tripleJetRates_emu",
                               "ratebin" : 61, "def" : False}


HoE_Jetscan_20_triple_dict = { "effbin" : 7,
                               "rateshist" : "tripleJetRates_HoE_TP05_emu",
                               "ratebin" : 21, "def" : False}

HoE_Jetscan_30_triple_dict = { "effbin" : 8,
                               "rateshist" : "tripleJetRates_HoE_TP05_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_40_triple_dict = { "effbin" : 9,
                               "rateshist" : "tripleJetRates_HoE_TP05_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_50_triple_dict = { "effbin" : 10,
                               "rateshist" : "tripleJetRates_HoE_TP05_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_60_triple_dict = { "effbin" : 11,
                               "rateshist" : "tripleJetRates_HoE_TP05_emu",
                               "ratebin" : 61, "def" : False}

HoE_Jetscan_TP1_20_triple_dict = { "effbin" : 12,
                               "rateshist" : "tripleJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 21, "def" : False}

HoE_Jetscan_TP1_30_triple_dict = { "effbin" : 13,
                               "rateshist" : "tripleJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_TP1_40_triple_dict = { "effbin" : 14,
                               "rateshist" : "tripleJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_TP1_50_triple_dict = { "effbin" : 15,
                               "rateshist" : "tripleJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_TP1_60_triple_dict = { "effbin" : 16,
                               "rateshist" : "tripleJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 61, "def" : False}

HoE_Jetscan_noTag_20_quad_dict = { "effbin" : 17,
                               "rateshist" : "quadJetRates_emu",
                               "ratebin" : 21, "def" : False}

HoE_Jetscan_noTag_30_quad_dict = { "effbin" : 18,
                               "rateshist" : "quadJetRates_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_noTag_40_quad_dict = { "effbin" : 19,
                               "rateshist" : "quadJetRates_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_noTag_50_quad_dict = { "effbin" : 20,
                               "rateshist" : "quadJetRates_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_noTag_60_quad_dict = { "effbin" : 21,
                               "rateshist" : "quadJetRates_emu",
                               "ratebin" : 61, "def" : False}

HoE_Jetscan_20_quad_dict = { "effbin" : 22,
                               "rateshist" : "quadJetRates_HoE_TP05_emu",
                               "ratebin" : 21, "def" : False}


HoE_Jetscan_30_quad_dict = { "effbin" : 23,
                               "rateshist" : "quadJetRates_HoE_TP05_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_40_quad_dict = { "effbin" : 24,
                               "rateshist" : "quadJetRates_HoE_TP05_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_50_quad_dict = { "effbin" : 25,
                               "rateshist" : "quadJetRates_HoE_TP05_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_60_quad_dict = { "effbin" : 26,
                               "rateshist" : "quadJetRates_HoE_TP05_emu",
                               "ratebin" : 61, "def" : False}


HoE_Jetscan_TP1_20_quad_dict = { "effbin" : 27,
                               "rateshist" : "quadJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 21, "def" : False}

HoE_Jetscan_TP1_30_quad_dict = { "effbin" : 28,
                               "rateshist" : "quadJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 31, "def" : False}

HoE_Jetscan_TP1_40_quad_dict = { "effbin" : 29,
                               "rateshist" : "quadJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 41, "def" : False}

HoE_Jetscan_TP1_50_quad_dict = { "effbin" : 30,
                               "rateshist" : "quadJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 51, "def" : False}

HoE_Jetscan_TP1_60_quad_dict = { "effbin" : 31,
                               "rateshist" : "quadJetRates_HoE_TP1_Ratio06_emu",
                               "ratebin" : 61, "def" : False}

Jetscan_noTag_triple_list = [HoE_Jetscan_noTag_20_triple_dict, HoE_Jetscan_noTag_30_triple_dict, HoE_Jetscan_noTag_40_triple_dict, HoE_Jetscan_noTag_50_triple_dict, HoE_Jetscan_noTag_60_triple_dict]
Jetscan_triple_list = [HoE_Jetscan_20_triple_dict, HoE_Jetscan_30_triple_dict, HoE_Jetscan_40_triple_dict, HoE_Jetscan_50_triple_dict, HoE_Jetscan_60_triple_dict]
Jetscan_TP1_triple_list = [HoE_Jetscan_TP1_20_triple_dict, HoE_Jetscan_TP1_30_triple_dict, HoE_Jetscan_TP1_40_triple_dict, HoE_Jetscan_TP1_50_triple_dict, HoE_Jetscan_TP1_60_triple_dict]

Jetscan_noTag_quad_list = [HoE_Jetscan_noTag_20_quad_dict, HoE_Jetscan_noTag_30_quad_dict, HoE_Jetscan_noTag_40_quad_dict, HoE_Jetscan_noTag_50_quad_dict, HoE_Jetscan_noTag_60_quad_dict]
Jetscan_quad_list = [HoE_Jetscan_20_quad_dict, HoE_Jetscan_30_quad_dict, HoE_Jetscan_40_quad_dict, HoE_Jetscan_50_quad_dict, HoE_Jetscan_60_quad_dict]
Jetscan_TP1_quad_list = [HoE_Jetscan_TP1_20_quad_dict, HoE_Jetscan_TP1_30_quad_dict, HoE_Jetscan_TP1_40_quad_dict, HoE_Jetscan_TP1_50_quad_dict, HoE_Jetscan_TP1_60_quad_dict]


#dictionaries for rate plots

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

rate_Ratio02_cut = {"histname" : "_HoE_Ratio02_emu",
                    "legendentry" : "LLP ID, Ratios > 0.2",
                    "color" : ROOT.kMagenta }

rate_TP2_Ratio02_cut = {"histname" : "_HoE_TP2_Ratio02_emu",
                    "legendentry" : "LLP ID, Ratios > 0.2 and E_{T,TP} > 2",
                        "color" : ROOT.kRed}


rate_Ratio06_cut = {"histname" : "_HoE_Ratio06_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6",
                    "color" : ROOT.kOrange }

rate_TP2_Ratio06_cut = {"histname" : "_HoE_TP2_Ratio06_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6 and E_{T,TP} > 2",
                    "color" : ROOT.kBlue }

rate_Ratio02_ORHT360_cut = {"histname" : "_HoE_Ratio02_ORHT360_emu",
                    "legendentry" : "LLP ID, Ratios > 0.2",
                    "color" : ROOT.kMagenta }

rate_TP2_Ratio02_ORHT360_cut = {"histname" : "_HoE_TP2_Ratio02_ORHT360_emu",
                    "legendentry" : "LLP ID, Ratios > 0.2 and E_{T,TP} > 2",
                        "color" : ROOT.kRed}

rate_TP1_Ratio06_cut = {"histname" : "_HoE_TP1_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6 and E_{T,TP} > 1",
                        "color" : ROOT.kRed}
rate_TP1_Ratio06_cut_jets = {"histname" : "_HoE_TP1_Ratio06_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6 and E_{T,TP} > 1",
                        "color" : ROOT.kRed}

rate_TP1_Ratio06_ORHT360_cut = {"histname" : "_HoE_TP1_ORHT360_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6 and E_{T,TP} > 1",
                        "color" : ROOT.kRed}


rate_Ratio06_ORHT360_cut = {"histname" : "_HoE_Ratio06_ORHT360_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6",
                    "color" : ROOT.kOrange }

rate_TP2_Ratio06_ORHT360_cut = {"histname" : "_HoE_TP2_Ratio06_ORHT360_emu",
                    "legendentry" : "LLP ID, Ratios > 0.6 and E_{T,TP} > 2",
                    "color" : ROOT.kBlue }

        

def main():
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360", extra_list = rateff_TP2_Ratio_OR360_list, extralabel = " E_{T,TP} > 2")
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360", extra_list = rateff_TP2_Ratio_OR360_list, extralabel = " E_{T,TP} > 2")

    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_ORHT360_Ratio_Gen", legendentry = " HT > 120 with LLP ID or HT > 360", extra_list = rateff_TP2_Ratio_OR360_list, extralabel = " E_{T,TP} > 2")
    myPlotter.plotRatesVsEff(rateff_Ratio_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_ORHT360_Ratio_Gen", legendentry = " HT > 120 with LLP ID or HT > 360", extra_list = rateff_TP2_Ratio_OR360_list, extralabel = " E_{T,TP} > 2")

    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio", legendentry = " HT > 120 with LLP ID", extra_list = rateff_TP2_Ratio_list, extralabel = " E_{T,TP} > 2")
    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio", legendentry = " HT > 120 with LLP ID", extra_list = rateff_TP2_Ratio_list, extralabel = "E_{T,TP} > 2")

    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_Ratio_Gen", legendentry = " HT > 120 with LLP ID", extra_list = rateff_TP2_Ratio_list, extralabel = " E_{T,TP} > 2")
    myPlotter.plotRatesVsEff(rateff_Ratio_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Ratio_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_Ratio_Gen", legendentry = " HT > 120 with LLP ID", extra_list = rateff_TP2_Ratio_list, extralabel = " E_{T,TP} > 2")

    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_ORHT360", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_ORHT360_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")
    myPlotter.plotRatesVsEff(rateff_TPE_OR360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_ORHT360_Gen", legendentry = " HT > 120 with LLP ID or HT > 360")

    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE", legendentry = " HT > 120 with LLP ID")

    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_TPE_Gen", legendentry = " HT > 120 with LLP ID")
    myPlotter.plotRatesVsEff(rateff_TPE_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_TPE_Gen", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_TPE_Gen", legendentry = " HT > 120 with LLP ID")

    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_noLLPtag", legendentry = " HT > threshold", def_rates = False)
    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_noLLPtag", legendentry = " HT > threshold", def_rates = False)

    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_noLLPtag_Gen", legendentry = " HT > threshold", def_rates = False)
    myPlotter.plotRatesVsEff(rateff_HTscan_noLLPtag_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_noLLPtag_Gen", legendentry = " HT > threshold", def_rates = False)


    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan", legendentry = " HT > threshold with LLP ID", extra_list = rateff_HTscan_ORHT360_list, extralabel = " or HT > 360")
    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan", legendentry = " HT > threshold with LLP ID", extra_list = rateff_HTscan_ORHT360_list, extralabel = " or HT > 360")

    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_Gen", legendentry = " HT > threshold with LLP ID", extra_list = rateff_HTscan_ORHT360_list, extralabel = " or HT > 360")
    myPlotter.plotRatesVsEff(rateff_HTscan_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_Gen", legendentry = " HT > threshold with LLP ID", extra_list = rateff_HTscan_ORHT360_list, extralabel = " or HT > 360")

    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_ORHT360", legendentry = " HT > threshold with LLP ID or HT > 360", extra_list = rateff_TP1_HTscan_ORHT360_list, extralabel = "; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_ORHT360", legendentry = " HT > threshold with LLP ID or HT > 360", extra_list = rateff_TP1_HTscan_ORHT360_list, extralabel = "; E_{T,TP} > 1")

    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_HTscan_ORHT360_Gen", legendentry = " HT > threshold with LLP ID or HT > 360", extra_list = rateff_TP1_HTscan_ORHT360_list, extralabel = "; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(rateff_HTscan_ORHT360_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_HTscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_HTscan_ORHT360_Gen", legendentry = " HT > threshold with LLP ID or HT > 360", extra_list = rateff_TP1_HTscan_ORHT360_list, extralabel = "; E_{T,TP} > 1")

    myPlotter.plotRatesVsEff(Jetscan_noTag_triple_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Jetscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_triple_Jetscan", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_triple_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_triple_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(Jetscan_noTag_triple_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Jetscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_triple_Jetscan", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_triple_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_triple_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")

    myPlotter.plotRatesVsEff(Jetscan_noTag_triple_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_Jetscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_triple_Jetscan_Gen", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_triple_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_triple_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(Jetscan_noTag_triple_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_Jetscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_triple_Jetscan_Gen", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_triple_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_triple_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")

    myPlotter.plotRatesVsEff(Jetscan_noTag_quad_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Jetscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_quad_Jetscan", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_quad_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_quad_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(Jetscan_noTag_quad_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Jetscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_quad_Jetscan", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_quad_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_quad_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")

    myPlotter.plotRatesVsEff(Jetscan_noTag_quad_list, myPlotter.NuGun_dict, [myPlotter.LLP_m250_60_500_dict, myPlotter.LLP_m250_60_1000_dict], "effJetID_HoE_DepthFB_Gen_Jetscan", title = "M_{H} = 250 GeV, M_{X} = 60 GeV", appendname = "_250_quad_Jetscan_Gen", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_quad_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_quad_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")
    myPlotter.plotRatesVsEff(Jetscan_noTag_quad_list, myPlotter.NuGun_dict, [myPlotter.LLP_m350_80_500_dict, myPlotter.LLP_m350_80_1000_dict], "effJetID_HoE_DepthFB_Gen_Jetscan", title = "M_{H} = 350 GeV, M_{X} = 80 GeV", appendname = "_350_quad_Jetscan_Gen", legendentry = " Jet E_{T} > threshold", def_rates = False, extra_list = Jetscan_quad_list, extralabel = " and LLP ID", extra_list2 = Jetscan_TP1_quad_list, extralabel2 = " and LLP ID; E_{T,TP} > 1")


#    myPlotter.plotJetRates("singleJetRates", "Single Jet", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut])
#    myPlotter.plotJetRates("doubleJetRates", "Double Jet", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut])
#    myPlotter.plotJetRates("tripleJetRates", "Triple Jet", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut])
#    myPlotter.plotJetRates("quadJetRates", "Quadruple Jet", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut])

    myPlotter.plotJetRates("singleJetRates", "Single Jet", "", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut, rate_TP1_Ratio06_cut_jets])
    myPlotter.plotJetRates("doubleJetRates", "Double Jet", "", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut, rate_TP1_Ratio06_cut_jets])
    myPlotter.plotJetRates("tripleJetRates", "Triple Jet", "", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut, rate_TP1_Ratio06_cut_jets])
    myPlotter.plotJetRates("quadJetRates", "Quadruple Jet", "", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut, rate_TP1_Ratio06_cut_jets])
    myPlotter.plotJetRates("htSumRates", "L1 H_{T}; H_{T}", "", [rate_no_cut, rate_HoE_cut, rate_Ratio06_cut, rate_TP1_Ratio06_cut])
    myPlotter.plotJetRates("htSumRates", "L1 H_{T}, LLP ID includes OR H_{T} > 360; H_{T}", "ORHT360", [rate_no_cut, rate_HoE_cut, rate_Ratio06_ORHT360_cut, rate_TP1_Ratio06_ORHT360_cut])

if __name__ == "__main__":

    main()

    
