#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#include<vector>
#include<utility>
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
namespace CutConstants{

class CutStr : public std::string{
public:
	CutStr() {}
    CutStr(std::string name,std::string cut) : std::string(name), cut(cut){}
    CutStr(std::string name,std::string cut, std::string title) : std::string(name),
            cut(cut), title(title){}

    std::string cut;
    std::string title;
};
std::string hhFilename    = "HHbb1o2l";

enum PROC  {TTBAR,WJETS,QCD,OTHER};
std::vector<CutStr > processes = {
        CutStr("ttbar"     ,"process==2","t#bar{t}"),
        CutStr("wjets"     ,"process==3","W+jets"),
        CutStr("qcd"       ,"process==8","Multijet"),
        CutStr("other"     ,"(process>1&&!(process==2||process==3||process==8))","Other SM")
};

enum REGION  {REG_SR, REG_TOPCR, REG_QGCR};

//CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N*lep_N*btag_N");
CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N");

CutStr aQCD ("aQCD"  , "process!=8");

CutStr wjjBC("wjjBC" , "wjjTau2o1<0.75");
CutStr exA  ("exA"   , "(hwwPT/hhMass>0.3)&&(hwwChi2<=11)");
CutStr bV   ("bV"    , "nAK4Btags==0");
CutStr abV  ("abV"   , "nAK4Btags!=0");

CutStr dPhiC ("dPhiC"  , "abs(dPhi_metll)<(3.14159/2)");
CutStr mllV ("mllV"  , "(dilepMass>6)&&(dilepMass<75)");
CutStr metC ("metC"  , "met/hhMass>0.1");
CutStr dRC   ("dRC"    , "dilepDR<1.0");
CutStr drCrC  ("drCrC"   , "dilepDR>0.4");

CutStr preSel("preSel"  , "passPre==1");

CutStr hbbMCS("hbbMass","hbbMass","#it{m}_{b#bar{b}} [GeV]");
CutStr hhMCS ("hhMass" ,"hhMass","#it{m}_{HH} [GeV]");

unsigned int nHbbMassBins   =30;
double minHbbMass = 30  ;
double maxHbbMass = 210 ;
unsigned int nHHMassBins   =86;
double minHHMass  = 700;
double maxHHMass  = 4000;
std::vector<double> hhMassBins = {
        700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
        1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,
        1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,
        2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,
        3100,3175,3250,3325,3400,3475,3550,3625,3700,3775,3850,3925,4000
};


unsigned int nInclHHMassBins   =128;
double minInclHHMass  = 0   ;
double maxInclHHMass  = 5050;
std::vector<double> inclHHMassBins = {0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,
        400,425,450,475,500,525,550,575,600,625,650,675,
        700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
        1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,
        1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,
        2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,
        3100,3175,3250,3325,3400,3475,3550,3625,3700,3775,3850,3925,4000,
        4075,4150,4225,4300,4375,4450,4525,4600,4675,4750,4825,4900,4975,5050
};
unsigned int nInclHbbMassBins   =42;
double minInclHbbMass  = 0   ;
double maxInclHbbMass  = 252;


CutStr hhRange  ("hhRange", hhMCS.cut+">"+ASTypes::flt2Str(minHHMass)+
        "&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxHHMass));
CutStr hhInclRange  ("hhInclRange" , hhMCS.cut+">"+ASTypes::flt2Str(minInclHHMass)+
        "&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxInclHHMass));
CutStr hbbRange ("hbbRange", hbbMCS.cut+">"+ASTypes::flt2Str(minHbbMass)+
        "&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxHbbMass));
CutStr hbbInclRange  ("hbbInclRange", hbbMCS.cut+">"+ASTypes::flt2Str(minInclHbbMass)+
        "&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxInclHbbMass));

std::string getHbbBinningString(bool inclusive){
    using namespace ASTypes;
    if(inclusive)
        return int2Str(nInclHbbMassBins)+","+flt2Str(minInclHbbMass)+","+flt2Str(maxInclHbbMass);
    else
        return int2Str(nHbbMassBins)+","+flt2Str(minHbbMass)+","+flt2Str(maxHbbMass);
}
std::string getHHBinningString(bool inclusive){
    using namespace ASTypes;
    const std::vector<double> * bins = &(inclusive ? inclHHMassBins : hhMassBins);
    std::string outString="";
    for(unsigned int iB = 0; iB < bins->size(); ++iB){
        if(iB) outString += ",";
        outString += flt2Str((*bins)[iB]);
    }
    return outString;
}

enum BKGModels  {BKG_QG, BKG_LOSTTW, BKG_MW, BKG_MT};
std::vector<CutStr > bkgSels = {
        CutStr("qg"    ,"hbbWQuark==0","q/g bkg."),
        CutStr("losttw","(hbbWQuark>0)&&(hbbWQuark<=3)","Lost t/W bkg."),
        CutStr("mw"     ,"hbbWQuark==4","#it{m}_{W} bkg."),
        CutStr("mt"     ,"hbbWQuark==5","#it{m}_{t} bkg.")
};

enum BKGModels2L  {BKG_NONTOP,BKG_TOP};
std::vector<CutStr > llBkgSels = {
		CutStr("nontop","hbbDecayTypeMC==0","NonTop bkg"),
		CutStr("top","hbbDecayTypeMC>0","Top bkg")

};

enum LEPCats  {LEP_EMU, LEP_E, LEP_MU};
std::vector<CutStr> lepCats = {
        CutStr("emu","isMuon>=0","e#mu"),
        CutStr("e"  ,"isMuon==0","e"),
        CutStr("mu" ,"isMuon==1","#mu")
};

enum DILEPCats  {LEP_INCL, LEP_SF, LEP_OF};
std::vector<CutStr> dilepCats = {
        CutStr("IF","((isMuon1>=0)&&(isMuon2>=0))","incl flavor"),
        CutStr("SF","((isMuon1==0)==(isMuon2==0))","same flavor"),
        CutStr("OF"  ,"((isMuon1==0)!=(isMuon2==0))","opposite flavor"),
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

enum SELCuts  {SEL_NONE,SEL_RPhiB,SEL_FULL};
std::vector<CutStr > selCuts = {
        CutStr("none",preSel.cut,"-ExB -#it{#DeltaR}_{ll} -#it{M}_{ll} -#it{MET} -#it{#Delta#Phi}_{met,ll}"),
		CutStr("R_phi_b",preSel.cut+"&&"+dRC.cut+"&&"+mllV.cut+"&&"+metC.cut,"full relax B, phi"),
        CutStr("full",preSel.cut+"&&"+bV.cut+"&&"+dRC.cut+"&&"+dPhiC.cut+"&&"+mllV.cut+"&&"+metC.cut,"")

};

struct CatIterator{
    bool firstBin = true;
    LEPCats  l =LEP_EMU;
    BTAGCats b =BTAG_LMT;
    PURCats  p =PURE_I;
    HADCuts  h =HAD_NONE;
    bool is(const LEPCats  cl ) const {return l == cl;}
    bool is(const BTAGCats cb ) const {return b == cb;}
    bool is(const PURCats  cp ) const {return p == cp;}
    bool is(const HADCuts  ch ) const {return h == ch;}
    std::string name() const {
        return lepCats[l] +"_"+btagCats[b]+"_"+purCats[p] +"_"+hadCuts[h];
    }
    std::string cut() const {
        return "("+lepCats[l].cut+"&&"+btagCats[b].cut+"&&"+purCats[p].cut +"&&"+hadCuts[h].cut+")";
    }
    bool getBin() {
        if(firstBin){
            firstBin = false;
            return true;
        }
        if(h < hadCuts.size()-1){
            h = HADCuts(h+1);
            return true;
        }
        else if(p < purCats.size()-1){
            p= PURCats(p+1);
            h=HAD_NONE;
            return true;
        }
        else if(b < btagCats.size()-1){
            b= BTAGCats(b+1);
            p=PURE_I;
            h=HAD_NONE;
            return true;
        }
        else if(l < lepCats.size()-1){
            l= LEPCats(l+1);
            b=BTAG_LMT;
            p=PURE_I;
            h=HAD_NONE;
            return true;
        }
        return false;
    }
    void reset() {
        firstBin = true;
        l =LEP_EMU;
        b =BTAG_LMT;
        p =PURE_I;
        h =HAD_NONE;
    }
};

struct DilepCatIterator{
    bool firstBin = true;
    DILEPCats  l =LEP_INCL;
    BTAGCats b =BTAG_LMT;
    SELCuts  s =SEL_NONE;
    bool is(const DILEPCats  cl ) const {return l == cl;}
    bool is(const BTAGCats cb ) const {return b == cb;}
    bool is(const SELCuts  cs ) const {return s == cs;}
    std::string name() const {
        return dilepCats[l] +"_"+btagCats[b] +"_"+selCuts[s];
    }
    std::string cut() const {
        return "("+dilepCats[l].cut+"&&"+btagCats[b].cut +"&&"+selCuts[s].cut+")";
    }
    bool getBin() {
        if(firstBin){
            firstBin = false;
            return true;
        }
        if(s < selCuts.size()-1){
            s = SELCuts(s+1);
            return true;
        }
        else if(b < btagCats.size()-1){
            b= BTAGCats(b+1);
            s=SEL_NONE;
            return true;
        }
        else if(l < dilepCats.size()-1){
            l= DILEPCats(l+1);
            b=BTAG_LMT;
            s=SEL_NONE;
            return true;
        }
        return false;
    }
    void reset() {
        firstBin = true;
        l =LEP_INCL;
        b =BTAG_LMT;
        s =SEL_NONE;
    }
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
        if(title == lepCats[LEP_EMU].title + ", "+btagCats[BTAG_LMT].title +
                ", "+purCats[PURE_I].title )
            title = "All categories";
        if(ex.size())
            title += ", ";
    }


    title += ex;
    return title;
}



std::vector<double> resPTBins = {600,700,750,800,850,900,1000,1100,1250,1500,1750,
                                 2000,2500,3000,3500,4000};

enum SIGNALS  {RADION,BLKGRAV,DILEP_RADION};
std::vector<CutStr > signals = {
        CutStr("radHH"     ,"Radion_hh_hVVhbb_inclusive_narrow","radion"),
        CutStr("blkHH"     ,"BulkGravTohhTohVVhbb_narrow","bulk graviton"),
		CutStr("raddHHcomb","spin0")
};
std::vector<std::vector<int> > signalMassBins = {
        {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500},
        {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500},
		{800,900,1000,2000,2500,3000}
};

//Constants for models when building limits
std::string MOD_MJ("MJ");
std::string MOD_MR("MR");
std::string MOD_MS("MH");

CutStr sigMCS("mx","mx","#it{m}_{#it{X}} [GeV]");
//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG#Higgs_cross_sections_and_decay_b
const double HtoBBBF = 0.5824;
const double HtoWWBF = 0.2137;
const double HtoZZBF = 0.02619;
const double HHtobbVVBF = 2*HtoBBBF*(HtoWWBF+HtoZZBF);


}


#endif

