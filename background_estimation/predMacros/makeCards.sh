#!/bin/bash
# . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/makeCards.sh nonCond_replaceDRforPTRatio 0 0 /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/
catName=$1
runSR=$2
condSignal=$3
macroLoc=$4

mkdir ${catName}
cd ${catName}
mkdir plots
RCMD="root -b -q '${macroLoc}/makeCard.C+(${runSR},${condSignal})'"
eval $RCMD
. comp.sh
cd ..
RCMD="scp -r ${catName} cmslpc26.fnal.gov:/uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/work/combineWork/"
if [ "${runSR}" = "0" ]; then
    eval $RCMD
fi
cd ${catName}
RCMD="text2workspace.py combinedCard.txt -o combined.root"
eval $RCMD
