
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../predTools/InputsHelper.h"
#include "../predTools/makePlots.C"


void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doSignal = false, bool isData = false){
    std::vector<PlotVar> vars;
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );

    std::vector<PlotSamp> samps;
    samps.emplace_back(processes[QCD],processes[QCD].cut);
    samps.emplace_back(processes[WJETS],processes[WJETS].cut);

    std::vector<PlotSel> sels;

    std::vector<CutStr > srPCrBtagCats = btagCats;
    for(const auto& b : qgBtagCats) srPCrBtagCats.push_back(b);

    for(const auto& l : lepCats)
        for(const auto& b : srPCrBtagCats)
            for(const auto& p : purCats)
                for(const auto& h : hadCuts){
                    sels.emplace_back(l+"_"+b+"_"+p+"_"+h         ,l.cut + "&&"+b.cut+"&&"+p.cut+"&&"+h.cut+"&&"+hhRange.cut+"&&"+hbbRange.cut);
                }

    std::string outFileName=filename+"_"+name+ "_getQCDRatio.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut, isData ? std::string("1.0") : nomW.cut);
}

void go(){
    makeFittingDistributions("bkg",hhFilename,"../trees/betrees_mc.root");
}

#endif
void getQCDRatio(){
    go();

}
