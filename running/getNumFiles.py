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
parser.add_argument("-d", "--dataDir", dest="dataDir", default="/eos/uscms/store/user/${USER}/13TeV/processed", help="Location of data [Default: \"/eos/uscms/store/user/${USER}/13TeV/ntuples\"]")
if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()


def addSample(name,sample,datarun,cross,numE,numF,cmdLine,dasName,configs) :
# 	if(re.match(".+-ext\d*", name)) : return

	files = []
	
	if args.dataDir.startswith("/eos/uscms/store/user") or args.dataDir.startswith("/store/user") :
		cmd = (" eos root://cmseos.fnal.gov find -f %s | egrep '.*%s(-ext[0-9]*|)_[0-9]*.root'  | wc -l" % ( args.dataDir, name))
		prefix = "root://cmseos:1094/"
	else:
		cmd = ("find %s -f | egrep '.*%s(-ext[0-9]*|)_[0-9]*.root' | wc -l" % (args.dataDir, name))
		prefix = ""
	
	ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
	result = ps.communicate()
	
	dascmd =  ("dasgoclient --query=\"summary dataset=%s\"  |  sed \'s/.*nevents\\\":\\([0-9]*\\),\\\".*/\\1/\'" % dasName)
	dasResult = subprocess.Popen(dascmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	cfgLine = ("%s\t%s\t%s" % (name,result[0].rstrip('\n'),dasResult[0].rstrip('\n')))
	configs.append(cfgLine)

outputLines = []
inputData = open(args.inputData, "r")
for line in inputData:
	if re.match("\s*#.*", line) : 
		continue
	match = re.match("^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
	if not match : 
		print "Do not understand:"
		print line
		continue
	addSample(match.group(1),match.group(2),match.group(3),match.group(4),
			match.group(5),match.group(6),match.group(7),match.group(8),outputLines)
print "Number of files per dataset: "
for line in outputLines:
	print line

print "Done!"
exit()
