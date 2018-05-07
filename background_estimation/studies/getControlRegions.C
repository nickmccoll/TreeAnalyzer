
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

    if(doSignal){
        samps.emplace_back(name,"1.0");
    } else if(isData){
        samps.emplace_back("data","1.0");
    } else {
        samps.emplace_back("bkg",aQCD.cut);
        samps.emplace_back("bkg_wQCD","1.0");
        samps.emplace_back(bkgSels[BKG_QG],bkgSels[BKG_QG].cut+"&&"+ aQCD.cut);
        samps.emplace_back(bkgSels[BKG_QG]+"_wQCD",bkgSels[BKG_QG].cut);

        samps.emplace_back(bkgSels[BKG_LOSTTW],bkgSels[BKG_LOSTTW].cut);
        samps.emplace_back(bkgSels[BKG_MW],bkgSels[BKG_MW].cut);
        samps.emplace_back(bkgSels[BKG_MT],bkgSels[BKG_MT].cut);

        samps.emplace_back("QCD",bkgSels[BKG_QG].cut+"&&(process==8)");
        samps.emplace_back("wjets",bkgSels[BKG_QG].cut+"&&(process==3)");
        samps.emplace_back("qg_ttbar",bkgSels[BKG_QG].cut+"&&(process==2)");
        samps.emplace_back("qg_other",bkgSels[BKG_QG].cut+"&&!(process==2||process==3||process==8)");

        samps.emplace_back("ttbar","process==2");
        samps.emplace_back("singlet","process==5");
        samps.emplace_back("other","!(process==2||process==3||process==5||process==8)");
    }
    std::vector<PlotSel> sels;

    auto addSels = [&](const std::string& prevcut,const std::string& prefix) {
        if(!isData){
            sels.emplace_back(prefix+"sr",prevcut + "&&"+btagCats[BTAG_I].cut+"&&"+nSJs.cut+"&&"+exA.cut+"&&"+wjjBC.cut+"&&"+bV.cut);
        }
        //For ttbar cr
        sels.emplace_back(prefix+"ab"         ,prevcut + "&&"+btagCats[BTAG_I].cut+"&&"+nSJs.cut+"&&"+exA.cut+"&&"+wjjBC.cut+"&&"+"nAK4Btags>0");
        sels.emplace_back(prefix+"ab_noM"     ,prevcut + "&&"+btagCats[BTAG_I].cut+"&&"+nSJs.cut+"&&"+exA.cut+"&&"+"nAK4Btags>0");
        sels.emplace_back(prefix+"ab_noEx"    ,prevcut + "&&"+btagCats[BTAG_I].cut+"&&"+nSJs.cut+"&&"+wjjBC.cut + "&&nAK4Btags>0");
        sels.emplace_back(prefix+"ab_noM_noEx",prevcut + "&&"+btagCats[BTAG_I].cut+"&&"+nSJs.cut+"&&"+"nAK4Btags>0");
        //for wjets cr

        sels.emplace_back(prefix+"hbb1",  prevcut + "&&"+exA.cut+"&&"+nSJs.cut+"&&"+wjjBC.cut+"&&"+bV.cut+"&&hbbCSVCat==1");
        sels.emplace_back(prefix+"hbb1_noM_noEx"  ,prevcut + "&&"+bV.cut+"&&"+nSJs.cut+"&&hbbCSVCat==1");
        sels.emplace_back(prefix+"hbb1_noM_noEx_noNJ"  ,prevcut + "&&"+bV.cut+"&&hbbCSVCat==1");
    };
    addSels("1.0","");
    addSels(lepCats[LEP_E].cut,"e_");
    addSels(lepCats[LEP_MU].cut,"mu_");

    std::string outFileName=filename+"_"+name+ "_findCRPlots.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut, isData ? std::string("1.0") : nomW.cut);
}


void go(){
    std::string genSel = hhRange.cut+"&&"+hbbRange.cut;
    genSel+= "&&"+lepCats[LEP_EMU].cut+"&&"+purCats[PURE_I].cut;
    makeFittingDistributions("data",hhFilename,"../trees/betrees_data.root",genSel,false,true);
    makeFittingDistributions("bkg",hhFilename,"../trees/betrees_mc.root",genSel,false);
    makeFittingDistributions("m1000",hhFilename,"../trees/out_radion_hh_bbinc_m1000_0.root",genSel,true);
    makeFittingDistributions("m2000",hhFilename,"../trees/out_radion_hh_bbinc_m2000_0.root",genSel,true);
    makeFittingDistributions("m3000",hhFilename,"../trees/out_radion_hh_bbinc_m3000_0.root",genSel,true);
}

#endif
void getControlRegions(){
    go();

}
