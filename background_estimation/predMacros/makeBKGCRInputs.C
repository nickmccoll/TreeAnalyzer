
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"
#endif

void makeBKGCRInputs(int year, bool doTopRegion = true, int bkgToDo = BKG_QG, int channel = 0, std::string treeDir = "../trees/"){
	if(doTopRegion){
	    selCuts1[SEL1_NONE].cut = preSel1.cut;
	    selCuts1[SEL1_LB].cut   = preSel1.cut+ "&&"+wjjBC.cut;
	    selCuts1[SEL1_LT].cut   = preSel1.cut+ "&&"+abV.cut;
	    selCuts1[SEL1_LTMB].cut = preSel1.cut;
	    selCuts1[SEL1_FULL].cut = preSel1.cut + "&&"+abV.cut+"&&"+wjjBC.cut;

        selCuts2[SEL2_NONE].cut  = preSel2.cut;
        selCuts2[SEL2_RPhiB].cut = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut;
        selCuts2[SEL2_FULL].cut  = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut+"&&"+dPhiC.cut;

	    hhFilename +="_TopCR";
	    go(year,bkgToDo,channel,treeDir+"/bkgCompLMT/");
	} else {
	    btagCats = qgBtagCats;
	    hhFilename +="_NonTopCR";
	    go(year,bkgToDo,channel,treeDir+"/bkgCompAB/");
	}

}
