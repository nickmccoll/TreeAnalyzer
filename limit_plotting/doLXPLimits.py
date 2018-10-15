#!/usr/bin/env python

import ROOT
import os, sys, re, optparse,pickle,shutil,json,random




parser = optparse.OptionParser()

parser.add_option("-s","--step",dest="step",type=float,help="step for mass points",default=100.0)
parser.add_option("-m","--min",dest="min",type=float,help="minimum Mass point",default=800.0)
parser.add_option("-M","--max",dest="max",type=float,help="maximum Mass point",default=3500.0)
parser.add_option("-o","--options",dest="options",help="Combine Options",default='-M AsymptoticLimits')
parser.add_option("-q","--queue",dest="queue",help="Batch Queue",default='1nw')
parser.add_option("-r","--randomSeeds",dest="randomize",type=int, help="randomize seeds",default=0)
(options,args) = parser.parse_args()


STEPS = int((options.max-options.min)/options.step)

massPoints=[]

for i in range(0,STEPS+1):
    massPoints.append(options.min+i*options.step)




for i,m in enumerate(massPoints):

    if options.randomize:
        suffixOpts=" -s {rndm}".format(rndm = int(random.random()*950000))
    else:
        suffixOpts=" "
        
    opt = " --rMax 0.5 "
    if m > 850 :
        opt = " --rMax 0.3 "
    if m > 950 :
        opt = " --rMax 0.2 "
    if m > 1150 :
        opt = " --rMax 0.13 "
    if m > 1550 :
        opt = " --rMax 0.08 "
    if m > 1950 :
        opt = " --rMax 0.06 "
    


    f=open("submit_{i}.sh".format(i=i),'w')
    execScript = 'cd {cwd} \n'.format(cwd=os.getcwd())
    execScript += 'eval `scramv1 runtime -sh` \n'
    execScript += "combine -m {mass} {options}  {file}\n".format(mass=m,options=options.options+suffixOpts+opt,file=args[0])
    f.write(execScript)
    f.close()
    os.system('chmod +x submit_{i}.sh'.format(i=i))

    if options.queue!="local":
        os.system('bsub -q {queue} submit_{i}.sh '.format(queue=options.queue,i=i))
    else:    
        os.system('sh submit_{i}.sh '.format(i=i))


