rm -f CMSSW_11_0_2.tgz

cd $CMSSW_BASE

cd ..

echo Making tarball of $CMSSW_BASE

#tar -zcf CMSSW_11_0_2.tgz CMSSW_11_0_2 --exclude="Filename*.root"

tar --exclude-caches-all --exclude-vcs -zcf $CMSSW_BASE.tgz $CMSSW_BASE --exclude=src --exclude=tmp

mv $CMSSW_BASE.tgz $CMSSW_BASE/src/HcalTrigger/Validation/scripts/condor
