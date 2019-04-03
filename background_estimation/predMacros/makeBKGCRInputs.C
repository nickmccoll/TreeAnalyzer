
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"
#endif

void makeBKGCRInputs(bool doTopRegion = true, int bkgToDo = BKG_QG, std::string treeDir = "/Users/brentstone/Dropbox/Physics/HHbbWW/BEtrees/SingleLepton17/"){
    if(doTopRegion){
        hadCuts[HAD_NONE].cut = preSel.cut;
        hadCuts[HAD_LB].cut   = preSel.cut+"&&"+wjjBC.cut;
        hadCuts[HAD_LT].cut   = preSel.cut+ "&&"+abV.cut;
        hadCuts[HAD_LTMB].cut = preSel.cut ;
        hadCuts[HAD_FULL].cut = preSel.cut + "&&"+abV.cut+"&&"+wjjBC.cut;


//        hadCuts[HAD_NONE].cut = preSel.cut;
//        hadCuts[HAD_LB].cut   = preSel.cut+"&&"+wjjBC.cut+"&&"+exA.cut;
//        hadCuts[HAD_LT].cut   = preSel.cut+ "&&"+abV.cut+"&&"+exA.cut;
//        hadCuts[HAD_LTMB].cut = preSel.cut +"&&"+exA.cut;
//        hadCuts[HAD_FULL].cut = preSel.cut + "&&"+abV.cut+"&&"+wjjBC.cut+"&&"+exA.cut;


        hhFilename +="_TopCR";
        go(bkgToDo,treeDir+"/bkgCompLMT/");
    } else {
        btagCats = qgBtagCats;
        hhFilename +="_QGCR";
        go(bkgToDo,treeDir+"/bkgCompAB/");
    }

}
