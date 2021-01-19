#!/bin/bash

printf "Beginning\n"

ls -l

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820

base_dir=`pwd`

printf "Unpacking tar\n"

tar -xf CMSSW_11_0_2.tgz

printf "Moving into base dir\n"

cd CMSSW_11_0_2

mkdir -p src
cd src

printf "Compiling...\n"

scram b ProjectRename
eval `scramv1 runtime -sh`

cd HcalTrigger/Validation/

printf "Here's where I am and what I see\n"
echo `pwd`
ls -l

printf "Running\n"

./run_rates.py -l $1

printf "Done running, moving files" 
ls -l

mv rates*.root ${base_dir}
cd ${base_dir}
printf "Here is output"
ls -l
