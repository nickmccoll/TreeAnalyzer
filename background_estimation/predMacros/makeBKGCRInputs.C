
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"
#endif

void makeBKGCRInputs(bool doTopRegion = true, int bkgToDo = BKG_QG, int channel = 0, std::string treeDir = "../trees/"){
	if(doTopRegion){
	    hadCuts[HAD_NONE].cut = preSel1.cut;
	    hadCuts[HAD_LB].cut   = preSel1.cut+ "&&"+wjjBC.cut;
	    hadCuts[HAD_LT].cut   = preSel1.cut+ "&&"+abV.cut;
	    hadCuts[HAD_LTMB].cut = preSel1.cut;
	    hadCuts[HAD_FULL].cut = preSel1.cut + "&&"+abV.cut+"&&"+wjjBC.cut;

        selCuts[SEL_NONE].cut  = preSel2.cut;
        selCuts[SEL_RPhiB].cut = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut;
        selCuts[SEL_FULL].cut  = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut+"&&"+dPhiC.cut;

	    hhFilename +="_TopCR";
	    go(bkgToDo,channel,treeDir+"/bkgCompLMT/");
	} else {
	    btagCats = qgBtagCats;
	    hhFilename +="_NonTopCR";
	    go(bkgToDo,channel,treeDir+"/bkgCompAB/");
	}

}
