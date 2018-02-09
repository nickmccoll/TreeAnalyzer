#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#include<vector>
#include<utility>
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
namespace CutConstants{

class CutStr : public std::string{
public:
    CutStr(std::string name,std::string cut) : std::string(name), cut(cut){}
    CutStr(std::string name,std::string cut, std::string title) : std::string(name), cut(cut), title(title){}
    std::string cut;
    std::string title;
};


std::string hhFilename = "HHlnujj";

CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N*lep_N*btag_N");

CutStr bkgTS("bkgTS" , "(process==2||process==5||process==7)");
CutStr bkgWS("bkgWS" , "(process==3||process==4||process==6||process==8)");
CutStr aQCD ("aQCD"  , "process!=8");

CutStr wjjBC("wjjBC" , "wjjTau2o1<0.55&&wjjMass>10");
CutStr exA  ("exA"   , "wlnuDR<3.2&&wwDM<2");
CutStr bV   ("bV"    , "nAK4Btags==0");


CutStr hbbMCS("hbbMass","hbbMass","#it{m}_{H#rightarrowbb} [GeV]");
CutStr hhMCS ("hhMass" ,"hhMass","#it{m}_{HH} [GeV]");

unsigned int nHbbMassBins   =90;
double minHbbMass = 30  ;
double maxHbbMass = 210 ;
unsigned int nHHMassBins   =168;
double minHHMass  = 800 ;
double maxHHMass  = 5000;

unsigned int nInclHHMassBins   =280;
double minInclHHMass  = 0   ;
double maxInclHHMass  = 7000;

unsigned int nInclHbbMassBins   =125;
double minInclHbbMass  = 0   ;
double maxInclHbbMass  = 250;


CutStr hhRange  ("hhRange" , hhMCS.cut+">"+ASTypes::flt2Str(minHHMass)+"&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxHHMass));
CutStr hhInclRange  ("hhInclRange" , hhMCS.cut+">"+ASTypes::flt2Str(minInclHHMass)+"&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxInclHHMass));
CutStr hbbRange ("hbbRange", hbbMCS.cut+">"+ASTypes::flt2Str(minHbbMass)+"&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxHbbMass));
CutStr hbbInclRange  ("hbbInclRange" , hbbMCS.cut+">"+ASTypes::flt2Str(minInclHbbMass)+"&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxInclHbbMass));

CutStr hhBinning ("hhBinning" ,ASTypes::int2Str(nHHMassBins)+","+ASTypes::flt2Str(minHHMass)+","+ASTypes::flt2Str(maxHHMass));
CutStr hbbBinning ("hbbBinning" ,ASTypes::int2Str(nHbbMassBins)+","+ASTypes::flt2Str(minHbbMass)+","+ASTypes::flt2Str(maxHbbMass));

CutStr hhInclBinning ("hhInclBinning" ,ASTypes::int2Str(nInclHHMassBins)+","+ASTypes::flt2Str(minInclHHMass)+","+ASTypes::flt2Str(maxInclHHMass));
CutStr hbbInclBinning ("hbbInclBinning" ,ASTypes::int2Str(nInclHbbMassBins)+","+ASTypes::flt2Str(minInclHbbMass)+","+ASTypes::flt2Str(maxInclHbbMass));

enum BKGModels  {BKG_QG, BKG_LOSTTW, BKG_MW, BKG_MT};
std::vector<CutStr > bkgSels = {
        CutStr("qg"    ,"hbbWQuark==0","q/g bkg."),
        CutStr("losttw","hbbWQuark>0&&hbbWQuark<=3","lost t/W bkg."),
        CutStr("mw"     ,"hbbWQuark==4","m_{W} bkg."),
        CutStr("mt"     ,"hbbWQuark==5","m_{t} bkg.")
};

enum LEPCuts  {LEP_EMU, LEP_E, LEP_MU};
std::vector<CutStr> lepSels = {
        CutStr("emu","isMuon>=0"),
        CutStr("e"  ,"isMuon==0"),
        CutStr("mu" ,"isMuon==1")
};
enum PURCuts  {PUR_I, PUR_LMT, PUR_L, PUR_M,PUR_T};
std::vector<CutStr > purSels = {
        CutStr("I"  ,"1.0"),
        CutStr("LMT","hbbCSVCat>=4"),
        CutStr("L"  ,"hbbCSVCat==4"),
        CutStr("M"  ,"hbbCSVCat==5"),
        CutStr("T"  ,"hbbCSVCat==6")
};
//enum HADCuts  {HAD_FULL, HAD_LWW, HAD_LB,HAD_LTMB,HAD_NONE};
enum HADCuts  {HAD_FULL, HAD_LWW, HAD_LB,HAD_LTMB,HAD_NONE};
std::vector<CutStr > hadSels = {
        CutStr("full",exA.cut+"&&"+wjjBC.cut+"&&"+bV.cut),
        CutStr("lWW" ,wjjBC.cut+"&&"+bV.cut),
        CutStr("lb"  ,exA.cut+"&&"+wjjBC.cut),
        CutStr("ltmb",exA.cut),
        CutStr("none","1.0")
};

std::vector<double> resPTBins = {600,700,750,800,850,900,1000,1100,1250,1500,1750,2000,2500,3000,3500,4000};

CutStr radionSig("radHH","radHH");
std::vector<int> signalMassBins = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};


}


#endif

