
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"

void makeDataDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    std::vector<PlotVar> vars;
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    }
    std::vector<PlotSamp> samps = { {name,"1.0"}};
    std::vector<PlotSel> sels;
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(h != hadCuts[HAD_FULL] ) continue;
        sels.emplace_back(l +"_"+b+"_"+p +"_"+h,
                l.cut +"&&"+b.cut+"&&"+p.cut+"&&"+h.cut);
    }
    std::string outFileName=filename+"_"+name+ (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,"1.0");
}

#endif

void makeBKGCRInputs(bool doTopRegion = true, int bkgToDo = BKG_QG, std::string treeDir = "../trees/"){
    if(doTopRegion){
        hadCuts[HAD_NONE].cut = nSJs.cut;
        hadCuts[HAD_LB].cut   = nSJs.cut+"&&"+wjjBC.cut;
        hadCuts[HAD_LT].cut   = nSJs.cut+ "&&"+abV.cut;
        hadCuts[HAD_LTMB].cut = nSJs.cut ;
        hadCuts[HAD_FULL].cut = nSJs.cut + "&&"+abV.cut+"&&"+wjjBC.cut;


//        hadCuts[HAD_NONE].cut = nSJs.cut;
//        hadCuts[HAD_LB].cut   = nSJs.cut+"&&"+wjjBC.cut+"&&"+exA.cut;
//        hadCuts[HAD_LT].cut   = nSJs.cut+ "&&"+abV.cut+"&&"+exA.cut;
//        hadCuts[HAD_LTMB].cut = nSJs.cut +"&&"+exA.cut;
//        hadCuts[HAD_FULL].cut = nSJs.cut + "&&"+abV.cut+"&&"+wjjBC.cut+"&&"+exA.cut;


        hhFilename +="_TopCR";
        go(static_cast<BKGModels>(bkgToDo),treeDir+"/bkgCompLMT/");
        if(bkgToDo < 0) makeDataDistributions("data",hhFilename,treeDir+"betrees_data.root","1.0",false);
    } else {
        btagCats = qgBtagCats;
        hhFilename +="_QGCR";
        go(static_cast<BKGModels>(bkgToDo),treeDir+"/bkgCompAB/");
        if(bkgToDo < 0) makeDataDistributions("data",hhFilename,treeDir+"betrees_data.root","1.0",false);
    }

}
