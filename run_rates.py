#!/usr/bin/python

import argparse
import os
import sys

rates_dict = {
"NuGun" : "RelValNuGun_PU_step1_ALL.root",
"QCD"   : "SMP-PhaseIITDRFall17DR-00002_step1_1-10.root",
"LLP_MH250_Ctau500" : "MH-250_MFF-120_CTau-500mm_step1.root",
"LLP_MH250_Ctau1000" : "MH-250_MFF-120_CTau-1000mm_step1.root",
"LLP_MH350_Ctau500" : "MH-350_MFF-160_CTau-500mm_step1.root",
"LLP_MH350_Ctau1000" : "MH-350_MFF-160_CTau-1000mm_step1.root",
"LLP_MH1000_Ctau500" : "MH-1000_MFF-450_CTau-500mm_step1.root",
"LLP_MH1000_Ctau1000" : "MH-1000_MFF-450_CTau-1000mm_step1.root",
}

def run_rates(name):
    print("Running rates for "+name)
    os.system("rates.exe "+name+" dir_ntuple/L1Ntuple_"+rates_dict[name])


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", action="store_true", default=False, help="Run over all samples")
    parser.add_argument("-nugun", action="store_true", default=False, help="Run over Neutrino Gun")
    parser.add_argument("-qcd", action="store_true", default=False, help="Run over QCD.")
    parser.add_argument("-m1000", action="store_true", default=False, help="Run over LLP with MH=1000")
    parser.add_argument("-m350", action="store_true", default=False, help="Run over LLP with MH=350")
    parser.add_argument("-m250", action="store_true", default=False, help="Run over LLP with MH=350")
    args = parser.parse_args()

    if args.a:
        for name in rates_dict.keys():
            run_rates(name)

    if args.nugun:
        run_rates("NuGun")

    if args.qcd:
        run_rates("QCD")

    if args.m1000:
        run_rates("LLP_MH1000_Ctau500")
        run_rates("LLP_MH1000_Ctau1000")

    if args.m350:
        run_rates("LLP_MH350_Ctau500")
        run_rates("LLP_MH350_Ctau1000")

    if args.m250:
        run_rates("LLP_MH250_Ctau500")
        run_rates("LLP_MH250_Ctau1000")

    os.system("mkdir -p result_rates")
    os.system("mv rates*.root result_rates")

if __name__ == "__main__":
    main()




        
    

