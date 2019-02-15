#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_DILEPTONCUTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_DILEPTONCUTCONSTANTS_H
#include<vector>
#include<utility>
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
namespace DileptonCutConstants{

class CutStr : public std::string{
public:
    CutStr(std::string name,std::string cut) : std::string(name), cut(cut){}
    CutStr(std::string name,std::string cut, std::string title) : std::string(name), cut(cut), title(title){}
    std::string cut;
    std::string title;
};
std::string hhFilename = "HHlnulnu";

enum PROC  {TTBAR,WJETS,ZJETS,OTHER};
std::vector<CutStr > processes = {
        CutStr("ttbar"     ,"process==2","t#bar{t}"),
        CutStr("wjets"     ,"process==3","W+jets"),
		CutStr("zjets"     ,"process==4","Z+jets"),
        CutStr("other"     ,"(process>1&&!(process==2||process==3||process==4))","Other SM")
};

enum REGION  {REG_SR, REG_TOPCR, REG_QGCR};

CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N*lep_N*btag_N");

CutStr aQCD ("aQCD"  , "process!=8");

CutStr bV   ("bV"    , "nAK4Btags==0");
CutStr abV  ("abV"   , "nAK4Btags!=0");

CutStr dR   ("dR"    , "dilepDR<1.6");
CutStr dPhi ("dPhi"  , "(dPhi_met_dilep<(3.14159/2))&&(dPhi_met_dilep>(-3.14159/2))");
CutStr mllV ("mllV"  , "(dilepMass>12)&&(dilepMass<75)");
CutStr metC ("metC"  , "met>40");
CutStr preSel("preSel"  , "passPre==1");


CutStr hbbMCS("hbbMass","hbbMass","#it{m}_{b#bar{b}} [GeV]");
CutStr hhMCS ("hhMass" ,"hhMass","#it{m}_{HH} [GeV]");

unsigned int nHbbMassBins   =90; // 90 originally
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

// for now, only one bkg model is being considered
enum BKGModels  {BKG_NOM};
std::vector<CutStr > bkgSels = {
		CutStr("nomBKG","0==0","Nom bkg.")
};

enum LEPCats  {LEP_INCL, LEP_SF, LEP_OF};
std::vector<CutStr> lepCats = {
        CutStr("inclFlav","((isMuon1>=0)&&(isMuon2>=0))","incl flavor"),
        CutStr("sameFlav","((isMuon1==0)==(isMuon2==0))","same flavor"),
        CutStr("oppFlav"  ,"((isMuon1==0)!=(isMuon2==0))","opposite flavor"),
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

enum HADCuts  {HAD_NONE,HAD_RPhiB,HAD_FULL};
std::vector<CutStr > hadCuts = {
        CutStr("none",preSel.cut,"-ExB -#it{m}_{D} -#it{p}_{T}/#it{m} -#tau_{0.75}"),
		CutStr("R_phi_b",preSel.cut+"&&"+dR.cut+"&&"+mllV.cut+"&&"+metC.cut,"full relax B, phi"),
        CutStr("full",preSel.cut+"&&"+bV.cut+"&&"+dR.cut+"&&"+dPhi.cut+"&&"+mllV.cut+"&&"+metC.cut,"")

};
std::string getCategoryLabel(const LEPCats lep, const BTAGCats btag){
    return lepCats[lep].title+", "+btagCats[btag].title;
}


std::string getCategoryLabel(const std::string& inStr){

    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(inStr);
    while (std::getline(tokenStream, token, '_')) {
    	tokens.push_back(token);
    }
    std::string title = "";
    auto getTitle = [&](const std::string in, const std::vector<CutStr >& opts ) -> std::string{
        for(const auto& c : opts) if(in == c ) return c.title;
        return "";
    };
    std::string lep,btag,ex;
    if(tokens.size()){lep = getTitle(tokens[0],lepCats);}
    if(tokens.size()>1){btag = getTitle(tokens[1],btagCats);}
    if(tokens.size()>2){ex = getTitle(tokens[2],hadCuts);}
    if(lep.size()){
        title += lep;
        if(btag.size()) title+=", ";
    }
    if(btag.size()){
        title += btag;
    }

    title += ex;
    return title;
}

std::vector<double> resPTBins = {600,700,750,800,850,900,1000,1100,1250,1500,1750,2000,2500,3000,3500,4000};

enum SIGNALS  {RADION /*,BLKGRAV*/};
std::vector<CutStr > signals = {
        CutStr("radHH"     ,"radion_hh_bbinc","radion"),
//        CutStr("blkHH"     ,"blkgrv_hh_bbinc","bulk graviton")
};
std::vector<std::vector<int> > signalMassBins = {
        {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500}/*,
        {600 ,650 ,700 ,800 ,900 ,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500}*/
};

//Constants for models when building limits
std::string MOD_MJ("MJ");
std::string MOD_MR("MR");
std::string MOD_MS("MH");

CutStr sigMCS("mx","mx","#it{m}_{#it{X}} [GeV]");


}


#endif

