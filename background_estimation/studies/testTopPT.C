
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../predTools/InputsHelper.h"
#include "../predTools/makePlots.C"


void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0"){
    std::vector<PlotVar> vars;

    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
//    vars.emplace_back("avgTopPT" ,";average top #it{p}_{T} [GeV]","avgTopPT",60,0,3000);
//    vars.emplace_back("avgTopPT" ,";average top #it{p}_{T} [GeV]","avgTopPT",60,0,3000,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
//    vars.emplace_back("avgTopPT" ,";average top #it{p}_{T} [GeV]","avgTopPT",60,0,3000,"hbbPT",std::string("#it{p}_{T,H#rightarrowbb} [GeV]"),"hbbPT",60,0,3000);

    std::vector<PlotSamp> samps;
    samps.emplace_back(bkgSels[BKG_LOSTTW],bkgSels[BKG_LOSTTW].cut);
    samps.emplace_back(bkgSels[BKG_MW]    ,bkgSels[BKG_MW].cut);
    samps.emplace_back(bkgSels[BKG_MT]    ,bkgSels[BKG_MT].cut);
    samps.emplace_back("mt_o_losttw"      ,std::string("(!(")+bkgSels[BKG_MT].cut+"))");
    samps.emplace_back("all"              ,"1.0");

    std::vector<PlotSel> sels;
    for(const auto& l : lepCats)
        for(const auto& b : btagCats)
            for(const auto& p : purCats)
                for(const auto& h : hadCuts){
                    if(h != hadCuts[HAD_FULL]) continue;
                    sels.emplace_back(l+"_"+b+"_"+p+"_"+h         ,l.cut + "&&"+b.cut+"&&"+p.cut+"&&"+h.cut+"&&"+hhRange.cut+"&&"+hbbRange.cut);
//                    sels.emplace_back(l+"_"+b+"_"+p+"_"+h +"_p50"         ,"(1+0.5*avgTopPT/2500)*("+l.cut + "&&"+b.cut+"&&"+p.cut+"&&"+h.cut+"&&"+hhRange.cut+"&&"+hbbRange.cut+")");
//                    sels.emplace_back(l+"_"+b+"_"+p+"_"+h +"_m50"         ,"(1-0.5*avgTopPT/2500)*("+l.cut + "&&"+b.cut+"&&"+p.cut+"&&"+h.cut+"&&"+hhRange.cut+"&&"+hbbRange.cut+")");
                }

    std::string outFileName=filename+"_"+name+ "_testTopPT.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,nomW.cut);
}

void go(){
    makeFittingDistributions("ttbar",hhFilename,"../trees/betrees_mc.root",
            processes[TTBAR].cut+"&&hbbWQuark!=0"+"&&"+hbbRange.cut+"&&"+hhRange.cut+"&&"+hadCuts[HAD_FULL].cut);
}

#endif
void testTopPT(){
    go();

}
