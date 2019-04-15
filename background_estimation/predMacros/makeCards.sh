    #!/bin/bash
# . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/makeCards.sh baseline 0 0 /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/
catName=$1
runSR=$2
signal=$3
macroLoc=$4


dirName="${catName}"
if [ "${runSR}" = "1" ] ; then
    dirName="${dirName}_TopCR"
fi
if [ "${runSR}" = "2" ] ; then
    dirName="${dirName}_QGCR"
fi
if [ "${signal}" = "0" ] ; then
    dirName="${dirName}_rad"
fi
if [ "${signal}" = "1" ] ; then
    dirName="${dirName}_blk"
fi
if [ "${condSignal}" = "0" ] ; then
    dirName="${dirName}_nc"
fi

mkdir ${dirName}
cd ${dirName}
mkdir plots
RCMD="root -b -q '${macroLoc}/makeCard.C+(${runSR},${signal})'"
eval $RCMD
. comp.sh
cd ..
RCMD="scp -r ${dirName} cmslpc26.fnal.gov:/uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/work/combineWork/"
if [ "${runSR}" = "0" ]; then
    eval $RCMD
fi
cd ${dirName}
RCMD="text2workspace.py combinedCard.txt -o combined.root"
eval $RCMD
