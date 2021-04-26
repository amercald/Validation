rm -f CMSSW_11_0_2.tgz

cd $CMSSW_BASE

cd ..

echo Making tarball of $CMSSW_BASE

tar --exclude-caches-all --exclude-vcs -zcf CMSSW_11_0_2.tgz CMSSW_11_0_2 --exclude=tmp --exclude="*.root" --exclude="*.pdf" --exclude="*.png"

mv $CMSSW_BASE.tgz $CMSSW_BASE/src/HcalTrigger/Validation/scripts/condor
