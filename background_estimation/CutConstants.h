#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H

#include <string>
namespace CutConstants{
std::string nomW = "xsec*trig_N*pu_N*lep_N*btag_N";
std::string resS  = "hbbWQuark!=0";
std::string nresS = "hbbWQuark==0";
std::string hbbBC = "hbbCSVCat>=4&&hbbMass>30";
std::string wjjBC = "wjjTau2o1<0.55&&wjjMass>10";
std::string exA   = "wlnuDR<3.2&&wwDM<2";
std::string bV    = "nAK4Btags==0";
std::string aQCD  = "process!=8";

std::vector<std::string > leps  = {"emu","e","mu"};
std::vector<std::string > lepsS = {"isMuon>=0","isMuon==0","isMuon==1"};
std::vector<std::string > purs  = {"LMT","L","M","T"};
std::vector<std::string > pursS = {"hbbCSVCat>=4", "hbbCSVCat==4","hbbCSVCat==5","hbbCSVCat==6"};

}


#endif

