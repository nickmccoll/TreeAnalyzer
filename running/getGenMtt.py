#!/usr/bin/python
#import os
#import re
import argparse
import glob
from ROOT import gROOT, TFile
gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Prepare and submit ntupling jobs')
parser.add_argument("-d", "--dataDir", dest="dataDir", default="/eos/uscms/store/user/lpchh/HHWWbb_trees/", help="directory with ntuples")
parser.add_argument("-o", "--outFile", dest="outFile", default="mtt.root", help="Output file")

args = parser.parse_args()

filenames = glob.glob(args.dataDir+'/*/*/*/*.root')

f0 = TFile.Open(filenames[0])
hist0 = f0.Get('genMttFiller/genMtt0')
hist1 = f0.Get('genMttFiller/genMtt1')
hist2 = f0.Get('genMttFiller/genMtt2')

isFile1 = True
for fname in filenames:
    file = TFile.Open(fname)
    if not file:
        continue
    if isFile1:
        isFile1 = False
        continue
    
    h0 = file.Get('genMttFiller/genMtt0')
    h1 = file.Get('genMttFiller/genMtt1')
    h2 = file.Get('genMttFiller/genMtt2')

    hist0.Add(h0,1)
    hist1.Add(h1,1)
    hist2.Add(h2,1)

    file.Close()

outF = TFile(args.outFile,'RECREATE')
outF.cd()
hist0.Write()
hist1.Write()
hist2.Write()
outF.Close()

exit()
