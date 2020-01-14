    #!/bin/bash
# . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/makeCards.sh baseline 0 0 2 /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/
# . /Users/brentstone/Dropbox/Physics/GitRepos/TreeAnalyzer/background_estimation/predMacros/makeCards.sh baseline 0 0 2 /Users/brentstone/Dropbox/Physics/GitRepos/TreeAnalyzer/background_estimation/predMacros/

catName=$1
runSR=$2
signal=$3
channel=$4      # 0 for 1l + 2l, 1 for 1l, 2 for 2l
macroLoc=$5


dirName="${catName}"
if [ "${runSR}" = "1" ] ; then
    dirName="${dirName}_TopCR"
fi
if [ "${runSR}" = "2" ] ; then
    dirName="${dirName}_NonTopCR"
fi
if [ "${signal}" = "0" ] ; then
    dirName="${dirName}_rad"
fi
if [ "${signal}" = "1" ] ; then
    dirName="${dirName}_blk"
fi

mkdir ${dirName}
cd ${dirName}
mkdir plots

RCMD="root -b -q '${macroLoc}/makeCard.C+(${runSR},${signal},${channel})'"
eval $RCMD

. comp.sh
cd ..
RCMD="scp -r ${dirName} bwstone@cmslpc26.fnal.gov:/uscms/home/bwstone/nobackup/TreeAnalyzer/combine/"
if [ "${runSR}" = "0" ]; then
    eval $RCMD
fi
cd ${dirName}
RCMD="text2workspace.py combinedCard.txt -o combined.root"
eval $RCMD
