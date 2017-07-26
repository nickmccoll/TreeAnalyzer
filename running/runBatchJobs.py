#!/usr/bin/python
import os
import sys
import re
import argparse
import subprocess
import cmd

parser = argparse.ArgumentParser(description='Prepare and submit ntupling jobs')
parser.add_argument("-m", "--macro",         dest="macro", default="runSomething.C", help="file to be run. [Default: runSomething.C]")
parser.add_argument("-b", "--runBatch",      dest="runBatch", action='store_true', default=False, help="Should we setup a condor job? [Default: False]")
parser.add_argument("-w", "--computeWeight", dest="compW", action='store_true', default=False, help="Should we include the parameters to calculate weight? [Default: False]")
parser.add_argument("-i", "--input",         dest="input", default="procdatasets.conf", help="input config or directory [Default: procdatasets.conf]")
parser.add_argument("-o", "--outputDir",     dest="outdir", default="out", help="Output directory for ntuples. [Default: \"out\"]")
parser.add_argument("-j", "--jobdir"       , dest="jobdir", default="jobs", help="Directory for job files  [Default: jobs]")
parser.add_argument("-r", "--runningDir"    , dest="runningDir", default="running/", help="Where to find helper files for running  [Default: running/]")
parser.add_argument("-s", "--runScript"    , dest="runScript", default="running/batchScript.sh", help="File that tells condor how to run  [Default: running/batchScript.sh]")
parser.add_argument("-c", "--compCommand"  , dest="compCommand", default="running/saComp.C", help="File with to make the macro compilable for batch (use none if not needed) [Default: running/saComp.C]")
parser.add_argument("-l", "--libDir"       , dest="libDir", default="$CMSSW_BASE/../TreeAnalyzer/framework/", help="Include location for batch (use none if not needed) [Default: $CMSSW_BASE/../TreeAnalyzer/framework/]")
parser.add_argument("-n", "--numFiles",      dest="numFiles", default=5, help="Number of files per job if no config [Default: 5]")
parser.add_argument("-t", "--treeInt",       dest="treeInt" , default=1, help="treeInt if no config [Default: 1]")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()





def compileSAMacro() :
    compM = args.macro
    bareM = os.path.splitext(os.path.basename(args.macro))[0]
    exeM = os.path.normpath(os.path.join(os.path.join(os.getcwd(), args.jobdir),"comp.exe"))
    if args.compCommand != "" :
        compM = os.path.normpath(os.path.join(args.jobdir,bareM+"_copy.C"))
        ps = subprocess.Popen(("cp %s %s" % (args.macro,compM)), shell=True)
        ps.wait()  
        with open(compM, "a") as of :
            with open(args.compCommand,"r") as inf :
                for line in inf:
                    of.write(line.replace("__MACRO__NAME__", bareM))
    envCMD = "root-config --cflags --libs"
    ps = subprocess.Popen(envCMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    envStr = ps.communicate()[0]    
    compCmd = ("g++ {MCR} {ENV} -o {OBJ} -lGenVector" .format (MCR=compM,ENV=envStr.rstrip(), OBJ=exeM))
    if args.libDir != "" :
        compCmd += (" -I{INC} {INC}/*.a "  .format( INC=os.path.normpath(args.libDir)))
    ps = subprocess.Popen(compCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print "Compiling macro:"
    output = ps.communicate()            
    if not os.path.isfile(exeM) :
        print "could not find " + exeM
        print output[0]
        print output[1]
        exit()            
    else:
        print "Compiled!"
    if args.compCommand != "" :
        ps = subprocess.Popen("rm %s" % compM,shell=True)
        ps.wait()
    
    return exeM

def compileLOCMacro() :
    compCmd = ("root -l -b -q \'%s/compileMacro.C(\"%s\")\'" % (args.runningDir, args.macro))
    ps = subprocess.Popen(compCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print "Compiling macro:"
    output = ps.communicate()            
    libName = re.sub(r'(.+)\.(\w+)', r'\1_\2.so', args.macro)
    if not os.path.isfile(libName) :
        print "could not find " + libName
        print output[0]
        print output[1]
        exit()            
    else:
        print "Compiled!"
    return libName
            

def prepareSampleJob(libName,outList, name, filelist, nFilesPerJob, treeInt, weightJob = False,  xsec = 1, numE = 1):
    
    pcmName = re.sub(r'(.+)\.so', r'\1_ACLiC_dict_rdict.pcm', libName)
    dName   = re.sub(r'(.+)\.so', r'\1.d', libName)
    
    nFiles = len(filelist)
    iF = 0
    iF2 = 0
    iJ = 0
    while iF < nFiles:
        if(iF2 == 0):
            jobfiles = open("{0}/files_{1}_{2}.txt".format(args.jobdir,name, iJ), "w")    
        if filelist[iF].startswith("/eos/uscms/store/user") :
            jobfiles.write("root://cmseos:1094/%s" % (re.match(r'/eos/uscms(.*)',filelist[iF]).group(1) ) )        
        else :
            jobfiles.write(filelist[iF])        
        jobfiles.write("\n")
        iF2 += 1
        iF += 1
        if(iF2 == int(nFilesPerJob)  or  iF ==  nFiles):
            jobfiles.close()
            
            inputF = "files_{0}_{1}.txt".format(name,iJ)
            outputF = "out_{0}_{1}.root".format(name,iJ)
            if not args.runBatch: 
                inputF  = os.path.normpath(os.path.join(os.path.join(os.getcwd(), args.jobdir),inputF))
                outputF = os.path.normpath(os.path.join(os.path.join(os.getcwd(), args.outdir),outputF))
                if weightJob  :
                    CMD = "root -l -b -q \'{cfg}+(\"{INF}\",{TreeInt},\"{OUTF}\",{xs},{nE})\'".format( cfg= args.macro,INF=inputF,TreeInt=treeInt,OUTF=outputF,xs=xsec,nE=numE)
                else :
                    CMD = "root -l -b -q \'{cfg}+(\"{INF}\",{TreeInt},\"{OUTF}\")\'".format( cfg=args.macro,INF=inputF,TreeInt=treeInt,OUTF=outputF)
                outList.append(CMD + " &")
            
            else :
                if weightJob  :
                    CMD = "./{MCR} {INF} {TreeInt} {OUTF} {xs} {nE}".format( MCR=os.path.basename(libName),INF=inputF,TreeInt=treeInt,OUTF=outputF,xs=xsec,nE=numE)
                else :
                    CMD = "./{MCR} {INF} {TreeInt} {OUTF}".format( MCR=os.path.basename(libName),INF=inputF,TreeInt=treeInt,OUTF=outputF)

                jobscript = open("{0}/submit_{1}_{2}".format(args.jobdir,name,iJ), "w")
                jobscript.write("""
universe                = vanilla
Requirements            = (Arch == "X86_64") && (OpSys == "LINUX")
Executable              = {runscript}
Arguments               = "'{cmd}' {outputdir} {CMSSWVERS}"
Output                  = logs/job_{samp}_{num}.out
Error                   = logs/job_{samp}_{num}.err
Log                     = logs/job_{samp}_{num}.log
use_x509userproxy       = true
initialdir              = {jobdir}
Should_Transfer_Files   = YES
transfer_input_files    = {workdir}/{libdir},{workdir}/files_{samp}_{num}.txt
WhenToTransferOutput    = ON_EXIT
Queue
""".format(
                runscript=args.runScript, 
                cmd=CMD, outputdir=args.outdir,CMSSWVERS=os.path.expandvars("$CMSSW_VERSION"),
                samp=name, num=iJ,
                jobdir=args.jobdir,workdir=os.path.normpath(os.path.join(os.getcwd(), args.jobdir)),                 
                libdir=os.path.basename(libName)            
                ))
                jobscript.close()
                outList.append("condor_submit {jobdir}/submit_{samp}_{num}".format(jobdir=args.jobdir,samp=name,num=iJ))

            iF2 = 0
            iJ += 1

def getFileList(inputDir,name="") :    
    if name == "":
        grepCmd = "egrep '.*root'"
    else :
        grepCmd = "egrep '.*%s(-ext[0-9]*|())(_[0-9]*|()).root'" % ( name)
        
    if inputDir.startswith("/eos/uscms/store/user") or inputDir.startswith("/store/user") :
        cmd = (" eos root://cmseos.fnal.gov find -f %s | %s" % ( inputDir,grepCmd))
        prefix = "root://cmseos:1094/"
    else:
        cmd = ("find %s | %s" % (inputDir,grepCmd))
        prefix = ""
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = ps.communicate()
    return result[0].rstrip('\n').split('\n')                



os.system("mkdir -p %s" % args.jobdir)
os.system("mkdir -p %s/logs" % args.jobdir)
if args.outdir.startswith("/eos/cms/store/user") or args.outdir.startswith("/store/user") :
    os.system("eosmkdir -p %s" % (args.outdir))
else :
    os.system("mkdir -p %s" % args.outdir)
    
if args.runBatch:
    libName = compileSAMacro()
else :
    libName = compileLOCMacro()
    
outList = []
if os.path.isdir(args.input) :
    fileList = getFileList(args.input,"")
    prepareSampleJob(libName,outList,"all", fileList, args.numFiles,args.treeInt)
else :
    inputData = open(args.input, "r")
    for line in inputData:
        if re.match("\s*#.*", line) : 
            continue
        match = re.match("^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
        if not match : 
            print "Do not understand:"
            print line
            continue
        fileList = getFileList(match.group(6),match.group(1))
        prepareSampleJob(libName,outList,match.group(1), fileList, match.group(5), match.group(2), 
                         args.compW, match.group(3),match.group(4))
        
if args.runBatch:
    subscript = open("submitall.sh", "w")
    subscript.write("#!/bin/bash\n")
    for line in outList :
        subscript.write("%s\n" % line)
else:
    for line in outList :
        print line

print "Done!"
exit()
