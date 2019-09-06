    #!/bin/bash
# . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/makeCards.sh baseline 0 0 2 /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/
catName=$1
runSR=$2
signal=$3
nLeps=$4
macroLoc=$5


dirName="${catName}"
if [ "${runSR}" = "1" ] ; then
    dirName="${dirName}_TopCR"
fi
if [ "${runSR}" = "2" ] && [ "${nLeps}" = "1" ] ; then
    dirName="${dirName}_QGCR"
fi
if [ "${runSR}" = "2" ] && [ "${nLeps}" = "2" ] ; then
    dirName="${dirName}_NonTopCR"
fi
if [ "${signal}" = "0" ] ; then
    dirName="${dirName}_rad"
fi
if [ "${signal}" = "1" ] ; then
    dirName="${dirName}_blk"
fi
if [ "${signal}" = "2" ] ; then
    dirName="${dirName}_rad"
fi

mkdir ${dirName}
cd ${dirName}
mkdir plots

if [ "${nLeps}" = "1" ] ; then
    RCMD="root -b -q '${macroLoc}/makeCard.C+(${runSR},${signal})'"
fi
if [ "${nLeps}" = "2" ] ; then
    RCMD="root -b -q '${macroLoc}/makeDileptonCard.C+(${runSR},${signal})'"
fi

eval $RCMD
. comp.sh
cd ..
RCMD="scp -r ${dirName} bwstone@cmslpc26.fnal.gov:/uscms/home/bwstone/nobackup/TreeAnalyzer/CMSSW_8_3_0/src/combine/"
if [ "${runSR}" = "0" ]; then
    eval $RCMD
fi
cd ${dirName}
RCMD="text2workspace.py combinedCard.txt -o combined.root"
eval $RCMD
