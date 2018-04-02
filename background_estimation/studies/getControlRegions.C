
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../predTools/InputsHelper.h"
#include "../predTools/makePlots.C"


void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doSignal = false){
    std::vector<PlotVar> vars;
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );

    std::vector<PlotSamp> samps;

    if(doSignal){
        samps.emplace_back(name,"1.0");
    } else {
        samps.emplace_back("bkg",aQCD.cut);
        samps.emplace_back("bkg_wQCD","1.0");
        samps.emplace_back(bkgSels[BKG_QG],bkgSels[BKG_QG].cut+"&&"+ aQCD.cut);
        samps.emplace_back(bkgSels[BKG_LOSTTW],bkgSels[BKG_LOSTTW].cut);
        samps.emplace_back(bkgSels[BKG_MW],bkgSels[BKG_MW].cut);
        samps.emplace_back(bkgSels[BKG_MT],bkgSels[BKG_MT].cut);
        samps.emplace_back("QCD",bkgSels[BKG_QG].cut+"&&!("+ aQCD.cut+")");
    }


    std::vector<PlotSel> sels;
    sels.emplace_back("sr",exA.cut+"&&"+wjjBC.cut+"&&"+bV.cut);
    sels.emplace_back("ab",exA.cut+"&&"+wjjBC.cut+"&&"+"nAK4Btags>0");
    sels.emplace_back("ab_noM",exA.cut+"&&"+"nAK4Btags>0");
    sels.emplace_back("ab_noEx",wjjBC.cut + "&&nAK4Btags>0");
    sels.emplace_back("ab_noM_noEx","nAK4Btags>0");

    std::string outFileName=filename+"_"+name+ "_findCRPlots.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,nomW.cut);
}


void go(){
    std::string genSel = hhRange.cut+"&&"+hbbRange.cut;
    genSel+= "&&"+lepCats[LEP_EMU].cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+purCats[PURE_I].cut+"&&"+ hadCuts[HAD_NONE].cut;
    makeFittingDistributions("bkg",hhFilename,"../trees/betrees_LMT_mc.root",genSel,false);
    makeFittingDistributions("m1000",hhFilename,"../trees/out_radion_hh_bbinc_m1000_0.root",genSel,true);
    makeFittingDistributions("m2000",hhFilename,"../trees/out_radion_hh_bbinc_m2000_0.root",genSel,true);
    makeFittingDistributions("m3000",hhFilename,"../trees/out_radion_hh_bbinc_m3000_0.root",genSel,true);
}

#endif
void getControlRegions(){
    go();

}
