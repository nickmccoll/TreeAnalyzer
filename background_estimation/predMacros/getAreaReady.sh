#!/bin/bash
# ./getAreaReady.sh /uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/work/analyzer_running/compiled/ trees/ /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/skimTree.C
inputdir=$1
outputdir=$2
skimLoc=$3

OUTSRV="cmslpc26.fnal.gov:"

mkdir ${outputdir}
mkdir ${outputdir}/../bkgInputs
mkdir ${outputdir}/../bkgInputsTopCR
mkdir ${outputdir}/../bkgInputsQGCR
mkdir ${outputdir}/../signalInputs
mkdir ${outputdir}/../signalInputsNoCond
scp ${OUTSRV}${inputdir}/* ${outputdir}/
mkdir ${outputdir}/mcPieces
mkdir ${outputdir}/dataPieces
mv ${outputdir}/data*.root ${outputdir}/dataPieces/
mv ${outputdir}/*.root ${outputdir}/mcPieces/
mv ${outputdir}/mcPieces/out*radion*.root ${outputdir}/
hadd ${outputdir}/betrees_data.root ${outputdir}/dataPieces/*.root
hadd ${outputdir}/betrees_mc.root ${outputdir}/mcPieces/*.root
RCMD="root -b -q '${skimLoc}(\"${outputdir}/betrees_mc.root\",\"${outputdir}/betrees\")'"
eval $RCMD

mkdir ${outputdir}/bkgCompLMT
mkdir ${outputdir}/bkgCompAB
RCMD="root -b -q '${skimLoc}(\"${outputdir}/betrees_mc.root\",\"${outputdir}/bkgCompLMT/betrees\",\"hbbCSVCat>=4\",true)'"
eval $RCMD
RCMD="root -b -q '${skimLoc}(\"${outputdir}/betrees_mc.root\",\"${outputdir}/bkgCompAB/betrees\",\"hbbCSVCat==1\",true)'"
eval $RCMD