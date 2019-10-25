#!/usr/bin/python
import os
import sys
import re
import argparse
import subprocess
import cmd

parser = argparse.ArgumentParser(description='Copy eos to eos')
parser.add_argument("-i", "--inputDir",      dest="inputDir",  help="file to be run. [Default: runSomething.C]")
parser.add_argument("-o", "--outputDir",     dest="outputDir",  help="Should we setup a condor job? [Default: False]")
parser.add_argument("-d", "--dryRun",        dest="dryRun", action='store_true', default=False, help="Should we include the parameters to calculate weight? [Default: False]")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


cmd = (" eos root://cmseos.fnal.gov find -f --xurl %s" % ( args.inputDir))
ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
result = ps.communicate()
fileList = result[0].rstrip('\n').split('\n')          

for f in fileList:
    mtch = re.match(r".*{INDIR}(.*)".format(INDIR=args.inputDir),f)
    cpcmd = ("xrdcp %s root://cmseos.fnal.gov/%s/%s" % (f,args.outputDir,mtch.group(1)))
    if args.dryRun is True:
        print cpcmd
    else:
        print ("copy %s" %(f))
        ps2 = subprocess.Popen(cpcmd, shell=True)
        ps2.wait()
            