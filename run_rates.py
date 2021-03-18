#!/usr/bin/env python

import argparse
import os
import sys

from dict_ntuples import dict as rates_dict

 
def run_rates(name):
    print("\nRunning rates for "+name+"\n")
    os.system("rates.exe "+name+" "+rates_dict[name])
 #   os.system("mv rates_"+name+".root  result_rates/")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", nargs="+", help="List of filenames")
    args = parser.parse_args()

#    os.system("mkdir -p result_rates")
    for f in args.l:
        if f not in rates_dict.keys():
            sys.exit(f+" was not found in dictionary")
        else:
            run_rates(f)

if __name__ == "__main__":
    main()
