#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "HistoPlotting/include/Plotter.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Configuration/interface/FillerConstants.h"

using namespace std;

void getDilepBEVars(int trigPreSel=0) {
	TString fArea = "/Users/brentstone/Dropbox/Physics/HHbbWW/BEtrees/Dilepton17/";
    HistGetter plotter;

    TFile *f1 = new TFile(fArea+"betrees_data.root");
    TFile *f2 = new TFile(fArea+"betrees_mc.root");
//    TFile *f2 = new TFile(fArea+"OldSelTrees/bkg.root");

    TString MX = "1000";
    TFile *fw0 = new TFile(fArea+"out_Radion_hh_hVVhbb_inclusive_narrow_M-"+MX+"_TuneCP5_13TeV-madgraph-pythia8_0.root");
    TFile *ft0 = new TFile(fArea+"out_GluGluToRadionToHHTo2B2Tau_M-"+MX+"_narrow_13TeV-madgraph_0.root");
    TFile *fw2 = new TFile(fArea+"out_BulkGravTohhTohVVhbb_narrow_M-"+MX+"_TuneCP5_13TeV-madgraph-pythia8_0.root");

    vector<TFile*> files = {f1,f2,fw0,ft0,fw2};
    vector<TString> fileNames = {"","","bbWWspin0","bbttspin0","bbWWspin2"};

    float dilepMass, hbbMass, hhMass, dilepDR, xsec, dPhi_metll, met, trigN, puN, weight, ptww;
    UChar_t hbbCSVCat, isMuon1, isMuon2, nAK4Btags, numBinHbb, passPre, process, dataset, nLepsTT, hbbDecayTypeMC;
    bool isSignal = true;

    map<UChar_t,TString> dataMap = { {7,"JetHT"}, {8,"MET"}, {11,"SingleElectron"}, {12,"SingleMuon"}, {13,"SinglePhoton"} };

    auto passBaseline = [&]() {
    	if (!passPre) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (dilepMass < 6 || dilepMass > 75) return false;
    	if (dilepDR > 1.0) return false;
    	if (met / hhMass < 0.1) return false;
    	if (fabs(dPhi_metll) > TMath::PiOver2()) return false;

    	return true;
    };

    auto passSR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat < 4) return false;
    	return true;
    };

    auto passNonTopCR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags > 0) return false;
    	if (hbbCSVCat != 1) return false;
    	return true;
    };

    auto passTopCR = [&]() {
    	if (!passPre) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (met / hhMass < 0.1) return false;
    	if (fabs(dPhi_metll) > TMath::PiOver2()) return false;
    	if (hbbCSVCat < 4) return false;
    	if (dilepMass < 6 || dilepMass > 75) return false;

    	if (nAK4Btags == 0) return false;
    	if (dilepDR < 0.4) return false;

    	return true;
    };

    auto passTestReg = [&](int step) {
    	if (!passPre) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (hbbCSVCat < 4) return false;
    	if (nAK4Btags != 0) return false;
    	if (fabs(dPhi_metll) > TMath::PiOver2()) return false;

    	if (dilepMass < 6) return false;

    	if (step == 1) {
        	if (met / hhMass < 0.1) return false;
    	} else if (step == 2) {
        	if (dilepDR > 1.0) return false;
    	} else if (step > 2) {
        	if (dilepDR > 1.0) return false;
        	if (met / hhMass < 0.1) return false;
    	}

    	return true;
    };

    auto pltVars = [&](TString pref, TString bCat, TString fN) {
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);
        TString hhS = "mhh2000toInf";
        if (hhMass > 700 && hhMass < 1000) hhS = "mhh700to1000";
        else if (hhMass > 1000 && hhMass < 1500) hhS = "mhh1000to1500";
        else if (hhMass > 1500 && hhMass < 2000) hhS = "mhh1500to2000";
        plotter.getOrMake1DPre(pref+"_"+bCat+"_"+hhS,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);

        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mhh",";M_{HH}",132,700,4000)->Fill(hhMass,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mll",";M_{ll}",500,0,500)->Fill(dilepMass,weight);

        TString lch = "emu";
        if (isMuon1 && isMuon2) lch = "mumu";
        else if (!isMuon1 && !isMuon2) lch = "ee";
        plotter.getOrMake1DPre(pref+"_"+lch+"_"+bCat,fN+"_mll",";M_{ll}",500,0,500)->Fill(dilepMass,weight);
    };

    auto mkPlots = [&](TString procName, bool isData, TString fN, TString region) {

//        TString bCat;
//        if (hbbCSVCat==3) bCat = "LL";
//        else if (hbbCSVCat==4) bCat = "bL";
//        else if (hbbCSVCat==5) bCat = "bM";
//        else if (hbbCSVCat==6) bCat = "bT";

    	TString ttS = TString::Format("%d",int(nLepsTT));
    	TString numbS = TString::Format("nb%d_",numBinHbb);

    	if (hbbDecayTypeMC==0) numbS = "bkgNonTop_";
    	else                   numbS = "bkgTop_";

    	pltVars(procName,region,fN);
    	if (procName=="ttbar") pltVars(procName+ttS,region,fN);

    	if (isData) pltVars("data",region,fN);
    	else if (!isSignal) pltVars("bkg" ,region,fN);

        pltVars(procName,region,numbS+fN);
    	if (procName=="ttbar") pltVars(procName+ttS,region,numbS+fN);

        if (isData) pltVars("data",region,numbS+fN);
        else if (!isSignal) pltVars("bkg" ,region,numbS+fN);

    	TString cat1l = "tw";
    	if (hbbDecayTypeMC == 0) return;
    	if (hbbDecayTypeMC == 4) cat1l = "mw";
    	if (hbbDecayTypeMC == 5) cat1l = "mt";
    	if (hbbDecayTypeMC > 5) cat1l = TString::Format("decaytype%d",hbbDecayTypeMC);

    	if (process==2 && nLepsTT >= 1 && passSR()) {
        	pltVars(procName+ttS+"_"+cat1l,region,fN);
        }
    	pltVars("bkg_"+cat1l,region,fN);

    };

    bool isData;
    for (unsigned int k=0; k<files.size(); k++) {
        TTree *t = (TTree*)files[k]->Get("treeMaker/Events");

        if (k==0) isData = true;
        else isData = false;
        if (k>1) isSignal = true;
        else isSignal = false;

        t->SetBranchAddress("dilepDR",&dilepDR);
        t->SetBranchAddress("hbbMass",&hbbMass);
        t->SetBranchAddress("hhMass",&hhMass);
        t->SetBranchAddress("isMuon1",&isMuon1);
        t->SetBranchAddress("isMuon2",&isMuon2);
        t->SetBranchAddress("dPhi_metll",&dPhi_metll);
        t->SetBranchAddress("met",&met);
        t->SetBranchAddress("hbbCSVCat",&hbbCSVCat);
        t->SetBranchAddress("nAK4Btags",&nAK4Btags);
        t->SetBranchAddress("passPre",&passPre);
        t->SetBranchAddress("dilepMass",&dilepMass);
        t->SetBranchAddress("hwwPT",&ptww);

        if(!isData) {
        	t->SetBranchAddress("process",&process);
        	t->SetBranchAddress("xsec",&xsec);
        	t->SetBranchAddress("pu_N",&puN);
        	t->SetBranchAddress("trig_N",&trigN);
            t->SetBranchAddress("numBinHbb",&numBinHbb);
            t->SetBranchAddress("nLepsTT",&nLepsTT);
            t->SetBranchAddress("hbbDecayTypeMC",&hbbDecayTypeMC);
        } else {
        	t->SetBranchAddress("dataset",&dataset);
        	weight = 1;
        }

        for (unsigned int i=0; i<t->GetEntries(); i++) {
            if (i%100000 == 0) printf("processing evt %d\n",i);
            t->GetEntry(i);

            TString pref;
            if (isData) pref = dataMap[dataset];
            else if (isSignal) pref = "m"+MX;
            else pref = FillerConstants::MCProcessNames[process];

            if(!isData) {
            	weight = xsec*trigN*puN;
//            	weight = xsec;
                if(process == 2) weight *= 0.836984; // ttbar scale factor

            	if(passSR()) mkPlots(pref,isData,fileNames[k],"SR");

            	bool passBT = (hbbCSVCat >= 4);
            	bool passBV = (nAK4Btags == 0);
            	bool passMll = dilepMass > 6 && dilepMass < 75;
            	bool passMet = met/hhMass > 0.1;
            	bool passDR = dilepDR < 1.0;
            	bool passDP = fabs(dPhi_metll) < TMath::PiOver2();

            	if (hbbMass > 30 && hbbMass < 210 && hhMass > 700 && hhMass < 4000) {
            		if (passBV && passMll && passMet && passDR && passDP) mkPlots(pref,isData,fileNames[k],"relaxBT");
            		if (passBT && passBV && passDP)                       mkPlots(pref,isData,fileNames[k],"relaxMll_MET_DR");
            		if (passBT && passMll && passMet && passDR)           mkPlots(pref,isData,fileNames[k],"relaxB_DP");
            		if (passBT && passMll && passDP && passDR)            mkPlots(pref,isData,fileNames[k],"relaxB_MET");
            		if (passBT && passMet && passDR)                      mkPlots(pref,isData,fileNames[k],"relaxB_DP_Mll");
            		if (passBT && passMet && passMll)                     mkPlots(pref,isData,fileNames[k],"relaxB_DP_DR");
            		if (passBT)                       mkPlots(pref,isData,fileNames[k],"wBT");
            		if (passBT && passMet)                      mkPlots(pref,isData,fileNames[k],"relaxB_DP_Mll_DR");

            	}

            	if(passTestReg(0)) mkPlots(pref,isData,fileNames[k],"baseTest");
            	if(passTestReg(1)) mkPlots(pref,isData,fileNames[k],"baseTest_wMet");
            	if(passTestReg(2)) mkPlots(pref,isData,fileNames[k],"baseTest_wDR");
            	if(passTestReg(3)) mkPlots(pref,isData,fileNames[k],"baseTest_wMetDR");

            }

            if (!isSignal) {
                if(passNonTopCR())  mkPlots(pref,isData,fileNames[k],"NonTopCR");
                if(passTopCR()) mkPlots(pref,isData,fileNames[k],"TopCR");
            }

        }
    }
    plotter.write("debug2l.root");

}
