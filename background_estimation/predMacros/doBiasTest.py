#!/usr/bin/env python
import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import os, sys, re, optparse,pickle,shutil,json,random
import subprocess

parser = optparse.OptionParser()

parser.add_option("-i","--inputFile",dest="input",help="Input file",default='')
parser.add_option("-N","--nToys",dest="toys",type=int,help="number of tens of toys",default=1)
parser.add_option("-r","--r",dest="r",type=float,help="Generate Signal",default=0)
parser.add_option("-m","--mass",dest="mass",type=int,help="Mass",default=0)
parser.add_option("-l","--label",dest="label",help="label")
parser.add_option("-c","--combine",dest="combine",help="combine executable")
parser.add_option("-f","--freq",dest="freq",type=int,help="frequentistFit",default=0)
parser.add_option("-s","--skipInitialFit",dest="skipInitialFit",type=int,help="skipInitialFit",default=0)
parser.add_option("-o","--onlyInitialFit",dest="onlyInitialFit",type=int,help="onlyInitialFit",default=0)

(options,args) = parser.parse_args()

if options.freq!=0 and options.skipInitialFit==0:
    rcmd = "/bin/bash -l -c '{combine} -m {mass} -M MultiDimFit --saveWorkspace   --freezeParameters r --setParameters r=0 -n {tag}_{mass} {card}'".format(combine=options.combine,mass=options.mass,card=options.input,tag=options.label)
    ps = subprocess.Popen(rcmd, shell=True)
    ps.wait()
    fitName="higgsCombine{tag}_{mass}.MultiDimFit.mH{mass}.root".format(tag=options.label,mass=options.mass)
    fitNameOut="higgsCombine{tag}_{mass}_forBias.MultiDimFit.mH{mass}.root".format(tag=options.label,mass=options.mass)
    fin = ROOT.TFile(fitName,"read");
    w=fin.Get("w")
    w.loadSnapshot("MultiDimFit")    ;
    w.var("r").setConstant(0);
    w.var("r").setVal(options.r);
    w.var("r").setMin(0);
    w.var("r").setMax(1);
    args = ROOT.RooArgSet(w.allVars());
    args.add(w.allCats());
    w.saveSnapshot("MultiDimFit2",args);
    fout = ROOT.TFile(fitNameOut,"recreate");
    fout.WriteTObject(w,"w");
    fout.Close();
    fin.Close();

if options.onlyInitialFit == 0:
    for i in range(0,options.toys):
        seed=int(201606+random.random()*10101982)
    
        if options.freq==0:
            rcmd = "/bin/bash -l -c '{combine} -m {mass} -M MaxLikelihoodFit --expectSignal={r} --bypassFrequentistFit -t 10 --seed {seed}  -n {tag}_{i}_{mass} --rMin=-1 --rMax=1 --skipBOnlyFit  {card}'".format(combine=options.combine,mass=options.mass,seed=seed,card=options.input,r=options.r,tag=options.label,i=i)
        else:
            fitName="higgsCombine{tag}_{mass}_forBias.MultiDimFit.mH{mass}.root".format(tag=options.label,mass=options.mass)
            rcmd = "/bin/bash -l -c '{combine} -m {mass} -M MaxLikelihoodFit --expectSignal={r} --bypassFrequentistFit -t 10 --seed {seed} --snapshotName MultiDimFit2 -n {tag}_{i}_{mass} --rMin=-1 --rMax=1 --skipBOnlyFit  {card}'".format(combine=options.combine,mass=options.mass,seed=seed,card=fitName,r=options.r,tag=options.label,i=i)
            ps = subprocess.Popen(rcmd, shell=True)
            ps.wait()
