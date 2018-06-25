
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"
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
        go(bkgToDo,treeDir+"/bkgCompLMT/");
    } else {
        btagCats = qgBtagCats;
        hhFilename +="_QGCR";
        go(bkgToDo,treeDir+"/bkgCompAB/");
    }

}
