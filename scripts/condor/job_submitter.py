#!/usr/bin/env python

from dict_ntuples import dict

import argparse
import os
import sys

def write_jdl(name):
    os.system("rm -f "+name+".jdl")
    with open(name+".jdl", "w") as file:
        file.write("universe = vanilla\n")
        file.write("Executable = submit_run_rates.sh\n")
        file.write("Transfer_Input_Files = CMSSW_11_0_2.tgz\n")
        file.write("should_transfer_files = YES\n")
        file.write("when_to_transfer_output = ON_EXIT\n")
        file.write("Output = "+name+"_$(Cluster)_$(Process).stdout\n")
        file.write("Error = "+name+"_$(Cluster)_$(Process).stderr\n")
        file.write("Log = "+name+"_$(Cluster)_$(Process).log\n")
        file.write("Arguments = "+name+"\n")
        file.write("Queue\n")

#    os.system("condor_submit "+name+".jdl")

def main():

    parser = argparse.ArgumentParser()
#    parser.add_argument("-o", help="Output path for files")
    parser.add_argument("-l", nargs="+", help ="List of filetypes to run")
    args = parser.parse_args()

#    outpath = args.o
#    os.system("mkdir -p "+outpath)

    #os.system("rm -f CMSSW_11_0_2.tgz")
    #os.system("cd $CMSSW_BASE/..")
    #os.system("tar -zcvf CMSSW_11_0_2.tgz CMSSW_11_0_2")
    #os.system("mv CMSSW_11_0_2.tgz CMSSW_11_0_2/src/HcalTrigger/Validation/submission_dir")
    os.system("./make_tar.sh")


    for name in args.l:
        if name not in dict.keys():
            sys.exit(name+" was not found in the dictionary.")
        else: 
            write_jdl(name)

if __name__ == "__main__":
    main()
