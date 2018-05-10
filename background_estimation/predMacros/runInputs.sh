#!/bin/bash
# . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/runInputs.sh . /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/

startDir=$1
macroLoc=$2

cd ${startDir}/signalInputs
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(3)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(0)' &"
eval $RCMD
cd ../signalInputsNoCond
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(0)' &"
eval $RCMD
wait
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(2)' &"
eval $RCMD
cd ../signalInputs
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(1)' &"
eval $RCMD
cd ../

cd ${startDir}/bkgInputs
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(5)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(0)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(1)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(2)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(3)' &"
eval $RCMD
cd .. 

cd ${startDir}/bkgInputsTopCR
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,5)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,0)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,1)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,2)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,3)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(true,-1)' &"
eval $RCMD
cd .. 

wait
cd ${startDir}/bkgInputsQGCR
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,5)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,0)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,1)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,2)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,3)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(false,-1)' &"
eval $RCMD
cd .. 

wait
cd ${startDir}/bkgInputs
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(4)' &"
eval $RCMD
cd .. 

