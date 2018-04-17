import os
import ROOT
from ROOT import gStyle,gROOT,gPad
from CMS_lumi import *
from tdrstyle import *
from RooPlotter import *
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

import optparse
parser = optparse.OptionParser()
parser.add_option("-i","--input",dest="inputFile",default='combined.root',help="input root datacard")
parser.add_option("-r","--fixR",dest="fixR",type=float,help="fix r in the fit")
parser.add_option("-s","--signalType",dest="signalType",default='',help="XWW or XWZ")
parser.add_option("-u","--doUncBand",dest="doUncBand",type=int,default=0,help="do uncertainty band")
(options,args) = parser.parse_args()



def saveCanvas(canvas,name):
  canvas.SaveAs(name+".root")
  canvas.SaveAs(name+".C")
  canvas.SaveAs(name+".pdf")
  #canvas.SaveAs(name+".png")
  canvas.SaveAs(name+".eps")
  os.system("convert -density 150 -quality 100 "+name+".eps "+name+".png")
  os.system("rm "+name+".eps")

def cmsLabel(canvas):
  #cmslabel_prelim(canvas,'2016',11)
  cmslabel_final(canvas,'2016',11)



directory='PlotsPostFit_'+options.signalType
os.system("mkdir -p "+directory)

plotter=RooPlotter(options.inputFile)
# plotter.w.Print()
s = options.signalType
if s=='radHH':
    plotter.fix("MH",1500)
elif s=='XWW':
    plotter.fix("MH",1350)
elif s=='XWZ':
    plotter.fix("MH",1450)
if options.fixR is not None:
    plotter.fix("r",options.fixR)

plotter.prefit()

if s=='XWW':
    plotter.addContribution("XWW",True,"X #rightarrow WW",3,1,ROOT.kOrange+10,0,ROOT.kWhite)
elif s=='XWZ':
    plotter.addContribution("XWZ",True,"X #rightarrow WZ",3,1,ROOT.kMagenta,0,ROOT.kWhite)
elif s=='radHH':
    plotter.addContribution("radHH",True,"X #rightarrow HH",3,1,ROOT.kMagenta,0,ROOT.kWhite)
# plotter.addContribution("resW",False,"W+V",1,1,ROOT.TColor.GetColor("#0F5500"),1001,ROOT.TColor.GetColor("#60B037"))
# plotter.addContribution("nonRes",False,"W+jets",1,1,ROOT.TColor.GetColor("#0041AA"),1001,ROOT.TColor.GetColor("#A5D2FF"),"_opt")
plotter.addContribution("mt",False,"#it{m}_{t} bkg.",1,1,ROOT.TColor.GetColor("#000000"),1001,ROOT.TColor.GetColor("#993366"))
plotter.addContribution("mw",False,"#it{m}_{W} bkg.",1,1,ROOT.TColor.GetColor("#000000"),1001,ROOT.TColor.GetColor("#ff66cc"))
plotter.addContribution("losttw",False,"lost t/W bkg.",1,1,ROOT.TColor.GetColor("#000000"),1001,ROOT.TColor.GetColor("#ffffcc"),"_opt")
plotter.addContribution("qg",False,"q/g bkg.",1,1,ROOT.TColor.GetColor("#000000"),1001,ROOT.TColor.GetColor("#9999ff"),"_opt")







for c in ['L','M','T']:
#for c in ['Wplus','Wminus']:
    if c=='vbf':
        pur=['NP']
    else:
        pur=['LP','HP']
    for p in pur:
        for l in ['mu','e']:
#            continue

#             label="W #rightarrow "+(("e","#mu")[l=='mu'])+"#nu"
            label=" "
            
            binLabel = "std_"+ l +"_"+c+"_"+p +"_full_13TeV"

            plotter.drawBinned("MR","#it{m}_{HH} [GeV]",label,binLabel,[0,0],options.doUncBand,1,"")
            cmsLabel(plotter.canvas)
            saveCanvas(plotter.canvas,directory+"/postFitMVV_"+s+"_"+binLabel)

            plotter.drawBinned("MJ","#it{m}_{H#rightarrowbb} [GeV]",label,binLabel,[0,0],options.doUncBand,0,"")
            cmsLabel(plotter.canvas)
            saveCanvas(plotter.canvas,directory+"/postFitMJJ_"+s+"_"+binLabel)

            plotter.drawBinned("MR","#it{m}_{HH} [GeV]",label,binLabel,[0,0],options.doUncBand,0,"MJ:sig:110:140")
            cmsLabel(plotter.canvas)
            saveCanvas(plotter.canvas,directory+"/postFitMVVW_"+s+"_"+binLabel)

            plotter.frame.GetXaxis().SetRangeUser(1000.,2000.)
            plotter.frame.GetYaxis().SetRangeUser(0.,1.1*(plotter.frame.getHist("datapoints").GetY()[9]+plotter.frame.getHist("datapoints").GetEYhigh()[9]))
            plotter.frame2.GetXaxis().SetRangeUser(1000.,2000.)
            plotter.line.SetX1(1000.)
            plotter.line.SetX2(2000.)
            saveCanvas(plotter.canvas,directory+"/postFitMVVWZoom_"+s+"_"+binLabel)

#            plotter.drawBinned("MLNuJ","m_{WV} (GeV)",label,c+"_"+l+"_"+p+"_13TeV",[0,10000],options.doUncBand,c!='vbf',"")
#            cmsLabel(plotter.canvas)
#            saveCanvas(plotter.canvas,directory+"/postFitMVV_"+s+"_"+c+"_"+l+"_"+p+"_13TeV")

#            plotter.drawBinned("MJ","m_{jet} (GeV)",label,c+"_"+l+"_"+p+"_13TeV",[64,106],options.doUncBand,0,"")
#            cmsLabel(plotter.canvas)
#            saveCanvas(plotter.canvas,directory+"/postFit_"+s+"_"+c+"_"+l+"_"+p+"_13TeV")

            plotter.drawBinned("MR","#it{m}_{HH} [GeV]",label,binLabel,[0,0],options.doUncBand,c!='vbf',"MJ:low:30:80")
            cmsLabel(plotter.canvas)
            saveCanvas(plotter.canvas,directory+"/postFitMVVLo_"+s+"_"+binLabel)

            plotter.drawBinned("MR","#it{m}_{HH} [GeV]",label,binLabel,[0,0],options.doUncBand,c!='vbf',"MJ:high:160:210")
            cmsLabel(plotter.canvas)
            saveCanvas(plotter.canvas,directory+"/postFitMVVHi_"+s+"_"+binLabel)



#for l in ['mu','e']:
#    for p in ['both']:
#        for c in ['vbf']:
#            plotter.drawProjection("MJ","m_{jet} [GeV]",c+"_"+l+"_"+p+"_13TeV",1,0)
#            saveCanvas(plotter.canvas,directory+"/postfitMJJ"+c+"_"+l+"_"+p)
#            plotter.drawProjection("MLNuJ","m_{WV} [GeV]",c+"_"+l+"_"+p+"_13TeV",1,0)
#            saveCanvas(plotter.canvas,directory+"/postfitMVV"+c+"_"+l+"_"+p)



#plotter=RooPlotter("LNuJJ_topPreFit_HP.root")    
#plotter.prefit()
#plotter.addContribution("topRes",True,"t#bar{t}",1,1,ROOT.kRed,0,ROOT.kWhite)
#plotter.addContribution("topNonRes",False,"non-resonant t#bar{t}",1,1,ROOT.kBlack,1001,ROOT.kGreen-5)
#plotter.drawStack("MJ","m_{jet} [GeV]","top_mu_HP_13TeV","top_mu_HP_13TeV")







