#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "HistoPlotting/include/Plotter.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Configuration/interface/FillerConstants.h"

using namespace std;

void getBEVars() {
	TString fArea = "/Users/brentstone/Dropbox/Physics/HHbbWW/BEtrees/2017/";
    HistGetter plotter;

    TFile *f1 = new TFile(fArea+"betrees_data.root");
    TFile *f2 = new TFile(fArea+"betrees_mc.root");
    vector<TFile*>  files = {f1,f2};
    vector<TString> fileNames = {"",""};
    vector<bool>    isDatas = {true,false};

    float hwwChi, hbbMass, hhMass, xsec, ht, lepETA, trigN, puN, weight;
    float dilepMass, dilepDR, llMetDphi, met, tau21, hwwPt, hwwLi;
    UChar_t hbbCSVCat, hbbDecayType, nAK4Btags, lepChan,
		isMuon, isMuon1, isMuon2, process, dataset, nLepsTT;

    map<UChar_t,TString> dataMap = { {7,"JetHT"}, {8,"MET"}, {11,"SingleElectron"}, {12,"SingleMuon"}, {13,"SinglePhoton"} };

    auto pltVars = [&](TString pref, TString bCat, TString fN) {
        TString hhS = "mhh2000toInf";
        if (hhMass > 700 && hhMass < 1000) hhS = "mhh700to1000";
        else if (hhMass > 1000 && hhMass < 1500) hhS = "mhh1000to1500";
        else if (hhMass > 1500 && hhMass < 2000) hhS = "mhh1500to2000";

        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat+"_"+hhS,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);

        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mhh",";M_{HH}",132,700,4000)->Fill(hhMass,weight);
    };

    auto mkPlots = [&](TString procName, bool isData, TString fN, TString region, int nQuarks) {
    	if (!isData) {
        	TString bName = "tw";
        	if(nQuarks == 0) bName = "qg";
        	else if (nQuarks==4) bName = "mw";
        	else if (nQuarks==5) bName = "mt";


//        	TString bCat;
//        	else if (hbbCSVCat==4) bCat = "bL";
//        	else if (hbbCSVCat==5) bCat = "bM";
//        	else if (hbbCSVCat==6) bCat = "bT";

        	pltVars(bName,region,fN);
        	if(nQuarks==0 && process!=8) pltVars(bName+"_noQCD",region,fN);

        	TString ttS = TString::Format("%d",int(nLepsTT));
        	if(process==2) {
        		pltVars(procName+ttS,region,fN);
            	pltVars(procName+"_"+bName,region,fN);
            	pltVars(procName+ttS+"_"+bName,region,fN);
        	}

    	}
    	pltVars(procName,region,fN);

    	if (isData) pltVars("data",region,fN);
    	else        pltVars("bkg" ,region,fN);

    };

    auto passSR1 = [&]() {
    	if (lepChan != 1) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat < 4) return false;
    	if (tau21 > 0.75) return false;
    	if (hwwPt / hhMass < 0.3) return false;
    	if (hwwChi > 11) return false;

    	return true;
    };

    auto passSR2 = [&]() {
    	if (lepChan != 2) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat < 4) return false;
    	if (dilepMass < 6 || dilepMass > 75) return false;
    	if (dilepDR > 1.0) return false;
    	if (llMetDphi > TMath::PiOver2()) return false;
    	if (met/hhMass < 0.1) return false;

    	return true;
    };

    auto passTopCR1 = [&]() {
    	if (lepChan != 1) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags == 0) return false;
    	if (hbbCSVCat < 4) return false;

    	return true;
    };

    auto passTopCR2 = [&]() {
    	if (lepChan != 2) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags == 0) return false;
    	if (hbbCSVCat < 4) return false;
    	if (dilepDR < 0.4) return false;

    	return true;
    };

    auto passQgCR = [&]() {
    	if (lepChan != 1) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat != 1) return false;

    	return true;
    };

    auto passNonTopCR = [&]() {
    	if (lepChan != 2) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat != 1) return false;

    	return true;
    };

    auto pass1lQCDTest = [&]() {
    	if (lepChan != 1) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700) return false;
    	if (nAK4Btags != 0) return false;
    	if (hbbCSVCat < 4) return false;
    	if (tau21 > 0.55) return false;
    	if (hwwPt / hhMass < 0.3) return false;
    	if (hwwChi > 11) return false;
    	if (isMuon==1) return false;

    	return true;
    };

    bool isData = false;
    for (unsigned int k=0; k<files.size(); k++) {
        TTree *t = (TTree*)files[k]->Get("treeMaker/Events");

        isData = isDatas[k];

        t->SetBranchAddress("hwwChi",&hwwChi);
        t->SetBranchAddress("hwwLi",&hwwLi);
        t->SetBranchAddress("hbbMass",&hbbMass);
        t->SetBranchAddress("hhMass",&hhMass);
        t->SetBranchAddress("hwwPT",&hwwPt);
        t->SetBranchAddress("wjjTau2o1",&tau21);
        t->SetBranchAddress("ht",&ht);
        t->SetBranchAddress("hbbCSVCat",&hbbCSVCat);
        t->SetBranchAddress("nAK4Btags",&nAK4Btags);
        t->SetBranchAddress("isMuon",&isMuon);
        t->SetBranchAddress("lepETA",&lepETA);

        t->SetBranchAddress("dilepMass",&dilepMass);
        t->SetBranchAddress("llMetDphi",&llMetDphi);
        t->SetBranchAddress("met",&met);

        t->SetBranchAddress("dilepDR",&dilepDR);

        t->SetBranchAddress("lepChan",&lepChan);
        t->SetBranchAddress("isMuon1",&isMuon1);
        t->SetBranchAddress("isMuon2",&isMuon2);

        if(!isData) {
        	t->SetBranchAddress("process",&process);
        	t->SetBranchAddress("xsec",&xsec);
        	t->SetBranchAddress("pu_N",&puN);
        	t->SetBranchAddress("trig_N",&trigN);
        	t->SetBranchAddress("hbbDecayType",&hbbDecayType);
            t->SetBranchAddress("nLepsTT",&nLepsTT);

        } else {
        	t->SetBranchAddress("dataset",&dataset);
        	weight = 1;
        }

        for (unsigned int i=0; i<t->GetEntries(); i++) {
            if (i%100000 == 0) printf("processing evt %d\n",i);
            t->GetEntry(i);
            TString pref = isData ? dataMap[dataset] : FillerConstants::MCProcessNames[process];

            if(!isData) {
            	weight = xsec*trigN*puN;
                if(process == 2) weight *= 0.836984; // ttbar scale factor

            	if(passSR1()) mkPlots(pref,isData,fileNames[k],"SR1",hbbDecayType);
            	if(passSR2()) mkPlots(pref,isData,fileNames[k],"SR2",hbbDecayType);
            }

            if(passQgCR())      mkPlots(pref,isData,fileNames[k],"QgCR",hbbDecayType);
            if(passNonTopCR())  mkPlots(pref,isData,fileNames[k],"NonTopCR",hbbDecayType);
            if(passTopCR1())     mkPlots(pref,isData,fileNames[k],"TopCR1",hbbDecayType);
            if(passTopCR2())     mkPlots(pref,isData,fileNames[k],"TopCR2",hbbDecayType);

            if(pass1lQCDTest())     mkPlots(pref,isData,fileNames[k],"qcdTest",hbbDecayType);

        }
    }
    plotter.write("debugNew.root");

}
