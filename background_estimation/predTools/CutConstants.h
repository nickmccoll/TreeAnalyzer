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


enum PROC  {TTBAR,WJETS,QCD,OTHER};
std::vector<CutStr > processes = {
        CutStr("ttbar"     ,"process==2","t#bar{t}"),
        CutStr("wjets"     ,"process==3","W+jets"),
        CutStr("qcd"       ,"process==8","Multijet"),
        CutStr("other"     ,"(process>1&&!(process==2||process==3||process==8))","Other SM")
};

enum REGION  {REG_SR, REG_TOPCR, REG_QGCR};

CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N*lep_N*btag_N");

CutStr aQCD ("aQCD"  , "process!=8");

CutStr wjjBC("wjjBC" , "wjjTau2o1<0.75");
CutStr exA  ("exA"   , "(hwwPT/hhMass>0.3)&&wwDM<125.0");
CutStr bV   ("bV"    , "nAK4Btags==0");
CutStr abV  ("abV"   , "nAK4Btags!=0");
CutStr preSel("preSel"  , "passPre==1");


CutStr hbbMCS("hbbMass","hbbMass","#it{m}_{b#bar{b}} [GeV]");
CutStr hhMCS ("hhMass" ,"hhMass","#it{m}_{HH} [GeV]");

unsigned int nHbbMassBins   =90;
double minHbbMass = 30  ;
double maxHbbMass = 210 ;
unsigned int nHHMassBins   =132;
double minHHMass  = 700;
double maxHHMass  = 4000;

unsigned int nInclHHMassBins   =200;
double minInclHHMass  = 0   ;
double maxInclHHMass  = 5000;

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
        CutStr("qg"    ,"hbbWQuark==0","Q/G bkg."),
        CutStr("losttw","hbbWQuark>0&&hbbWQuark<=3","Lost t/W bkg."),
        CutStr("mw"     ,"hbbWQuark==4","m_{W} bkg."),
        CutStr("mt"     ,"hbbWQuark==5","m_{t} bkg.")
};

enum LEPCats  {LEP_EMU, LEP_E, LEP_MU};
std::vector<CutStr> lepCats = {
        CutStr("emu","isMuon>=0","e#mu"),
        CutStr("e"  ,"isMuon==0","e"),
        CutStr("mu" ,"isMuon==1","#mu")
};

enum BTAGCats  {BTAG_LMT, BTAG_L, BTAG_M, BTAG_T};
std::vector<CutStr > btagCats = {
        CutStr("LMT","hbbCSVCat>=4","bLMT"),
        CutStr("L"  ,"hbbCSVCat==4","bL"),
        CutStr("M"  ,"hbbCSVCat==5","bM"),
        CutStr("T"  ,"hbbCSVCat==6","bT")
};

std::vector<CutStr > qgBtagCats = {
        CutStr("LMT","hbbCSVCat==1",""),
        CutStr("L"  ,"hbbCSVCat==1","")
};
CutStr inclBtagCat("I","hbbCSVCat>=0");

enum   PURCats {PURE_I, PURE_LP, PURE_HP};
std::vector<CutStr > purCats = {
        CutStr("I","1.0","LHP"),
        CutStr("LP" ,"wjjTau2o1>=0.55","LP"),
        CutStr("HP"  ,"wjjTau2o1<0.55","HP")
};

enum HADCuts  {HAD_NONE,HAD_LB,HAD_LT,HAD_LTMB,HAD_FULL};
std::vector<CutStr > hadCuts = {
        CutStr("none",preSel.cut,"-ExB -#it{m}_{D} -#it{p}_{T}/#it{m} -#tau_{0.75}"),
        CutStr("lb"  ,preSel.cut+"&&"+exA.cut+"&&"+wjjBC.cut,"-ExB"),
        CutStr("lt"  ,preSel.cut+"&&"+exA.cut+"&&"+bV.cut,"-#tau_{0.75}"),
        CutStr("ltmb",preSel.cut+"&&"+exA.cut,"-ExB -#tau_{0.75}"),
        CutStr("full",preSel.cut+"&&"+exA.cut+"&&"+wjjBC.cut+"&&"+bV.cut,"")

};

std::string getCategoryLabel(const LEPCats lep, const BTAGCats btag, const PURCats pur ){
    return lepCats[lep].title+", "+btagCats[btag].title+", "+purCats[pur].title;
}



std::string getCategoryLabel(const std::string& inStr){

    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(inStr);
    while (std::getline(tokenStream, token, '_'))
    {tokens.push_back(token);}
    std::string title = "";
    auto getTitle = [&](const std::string in, const std::vector<CutStr >& opts ) -> std::string{
        for(const auto& c : opts) if(in == c ) return c.title;
        return "";
    };
    std::string lep,btag,pur,ex;
    if(tokens.size()){lep = getTitle(tokens[0],lepCats);}
    if(tokens.size()>1){btag = getTitle(tokens[1],btagCats);}
    if(tokens.size()>2){pur = getTitle(tokens[2],purCats);}
    if(tokens.size()>3){ex = getTitle(tokens[3],hadCuts);}
    if(lep.size()){
        title += lep;
        if(btag.size() || pur.size()) title+=", ";
    }
    if(btag.size()){
        title += btag;
        if(pur.size())
            title += ", ";
    }
    if(pur.size()){
        title += pur;
        if(title == lepCats[LEP_EMU].title + ", "+btagCats[BTAG_LMT].title + ", "+purCats[PURE_I].title )
            title = "All categories";
        if(ex.size())
            title += ", ";
    }


    title += ex;
    return title;
}

std::vector<double> resPTBins = {600,700,750,800,850,900,1000,1100,1250,1500,1750,2000,2500,3000,3500,4000};

enum SIGNALS  {RADION,BLKGRAV};
std::vector<CutStr > signals = {
        CutStr("radHH"     ,"radion_hh_bbinc","radion"),
        CutStr("blkHH"     ,"blkgrv_hh_bbinc","bulk graviton")
};
std::vector<std::vector<int> > signalMassBins = {
        {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500},
        {600 ,650 ,700 ,800 ,900 ,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500}
};

//Constants for models when building limits
std::string MOD_MJ("MJ");
std::string MOD_MR("MR");
std::string MOD_MS("MH");

CutStr sigMCS("mx","mx","#it{m}_{#it{X}} [GeV]");


}


#endif

