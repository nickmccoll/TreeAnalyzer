#!/usr/bin/python
import os
import sys
import re
import argparse
import subprocess
from ROOT import gROOT, TFile
gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Prepare and submit ntupling jobs')
parser.add_argument("-i", "--inputData", dest="inputData", default="datasets.conf", help="Input dataset. [Default: datasets.conf]")
parser.add_argument("-o", "--outputData", dest="outputData", default="procdatasets.conf", help="output dataset. [Default: procdatasets.conf]")
parser.add_argument("-d", "--dataDir", dest="dataDir", default="/eos/uscms/store/user/${USER}/13TeV/processed", help="Location of data [Default: \"/eos/uscms/store/user/${USER}/13TeV/ntuples\"]")
parser.add_argument("-t", "--histName", dest="histName", default="eventCounter/eventWeights", help="data tree [Default: \"treeMaker/Events\"]")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def get_num_mc_events(filelist, prefix=''):
    totposentries = 0
    totnegentries = 0
    for filename in filelist :
        filepath = '/'.join(['%s' % prefix, '%s' % filename])
        file = TFile.Open(filepath)
        if not file:
            continue
        hist = file.Get(args.histName)
        totnegentries += hist.GetBinContent(3)
        totposentries += hist.GetBinContent(2)
    return totposentries, totnegentries
def get_num_data_events(filelist, prefix=''):
    totposentries = 0
    for filename in filelist :
        filepath = '/'.join(['%s' % prefix, '%s' % filename])
        file = TFile.Open(filepath)
        if not file:
            continue
        hist = file.Get(args.histName)
        totposentries += hist.GetBinContent(2)
    return totposentries

def addSample(name,datarun,cross,dasName,configs) :
    if(re.match(".+-ext\d*", name)) : return
    	
    files = []
    	
    if args.dataDir.startswith("/eos/uscms/store/user") or args.dataDir.startswith("/store/user") :
    	cmd = (" eos root://cmseos.fnal.gov find -f %s | egrep '.*%s(-ext[0-9]*|)_[0-9]*.root'" % ( args.dataDir, name))
    	prefix = "root://cmseos:1094/"
    else:
    	cmd = ("find %s | egrep '.*%s(-ext[0-9]*|)_[0-9]*.root'" % (args.dataDir, name))
    	prefix = "/"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = ps.communicate()
    filelist = result[0].rstrip('\n').split('\n')
    totNumE = 0
    if datarun.startswith("MC"):
    	nposevents, nnegevents = get_num_mc_events(filelist, prefix)
    	print "Sample " + name + " has " + str(nposevents) + " positive and " + str(nnegevents) + " negative weight events"
    	totNumE = nposevents - nnegevents
    else :
    	nposevents = get_num_data_events(filelist, prefix)
    	print "Sample " + name + " has " + str(nposevents) + " number of events"
    	totNumE = nposevents
    
    if totNumE > 0 :
       cfgLine = ("%s\t%s\t%s\t%s\t5\t%s" % (name,str(1 if datarun.startswith("MC") else 0 ),cross,str(totNumE),args.dataDir))
       configs.append(cfgLine)

outputLines = []
inputData = open(args.inputData, "r")
for line in inputData:
	if re.match("\s*#.*", line) : 
		continue
	match = re.match("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s*", line)
	if not match : 
		print "Do not understand:"
		print line
		continue
	addSample(match.group(1),match.group(3),match.group(5),
			match.group(7),outputLines)
print "Creating new config file: " + args.outputData
with open(args.outputData, "w") as script:
	for line in outputLines:
		script.write("%s\n" % line)

print "Done!"
exit()
