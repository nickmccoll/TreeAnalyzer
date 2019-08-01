#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "HistoPlotting/include/Plotter.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Configuration/interface/FillerConstants.h"

using namespace std;

void getBEVars(int trigPreSel=0) {
	TString fArea = "/Users/brentstone/Dropbox/Physics/HHbbWW/BEtrees/SingleLepton17/";
    HistGetter plotter;

    TFile *f1 = new TFile(fArea+"betrees_data.root");
    TFile *f2 = new TFile(fArea+"betrees_mc.root");
    vector<TFile*> files = {f1,f2};
    vector<TString> fileNames = {"",""};

    float hwwChi, hbbMass, hhMass, tau21, xsec, ht, hwwPt, ht_chs, ht_pup, lepETA, trigN, puN, weight;
    UChar_t hbbCSVCat, hbbWQuark, nAK4Btags, passPre, isMuon, passTrigCHS,
	        passTrigPUP_loose, passTrigPUP_tight, process, dataset;

    map<UChar_t,TString> dataMap = { {7,"JetHT"}, {8,"MET"}, {11,"SingleElectron"}, {12,"SingleMuon"}, {13,"SinglePhoton"} };

    auto passBaseline = [&]() {
    	if (!passPre) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (hwwChi > 11) return false;
    	if (tau21 > 0.75) return false;
    	if (hwwPt/hhMass < 0.3) return false;
    	return true;
    };

    auto passSR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags > 0) return false;
    	if (hbbCSVCat < 4) return false;
    	return true;
    };

    auto passTopCR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags == 0) return false;
    	if (hbbCSVCat < 4) return false;
    	return true;
    };

    auto passQgCR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags > 0) return false;
    	if (hbbCSVCat != 1) return false;
    	return true;
    };

    auto pltVars = [&](TString pref, TString bCat, TString fN) {
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mhh",";M_{HH}",132,700,4000)->Fill(hhMass,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_ht",";H_{T}",100,400,3000)->Fill(ht,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_eta",";#eta",40,-2.4,2.4)->Fill(lepETA,weight);
    };

    auto mkPlots = [&](TString procName, bool isData, TString fN, TString region, int nQuarks) {
    	if (!isData) {
        	TString bName = "tw";
        	if(nQuarks == 0) bName = "qg";
        	else if (nQuarks==4) bName = "mw";
        	else if (nQuarks==5) bName = "mt";
        	else if (nQuarks > 5) printf("some shit is wrong\n");

//        	TString bCat;
//        	if (hbbCSVCat==3) bCat = "LL";
//        	else if (hbbCSVCat==4) bCat = "bL";
//        	else if (hbbCSVCat==5) bCat = "bM";
//        	else if (hbbCSVCat==6) bCat = "bT";

        	pltVars(bName,region,fN);
    	}

    	pltVars(procName,region,fN);
    	if (isData) pltVars("data",region,fN);
    	else        pltVars("bkg" ,region,fN);

    };


//    double htBins[] = {0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1700,2000,2500,3000};
//    int nBins = 18;
//    TH1F *hht = new TH1F("hht","hht",nBins,htBins);

    bool isData;
    for (unsigned int k=0; k<files.size(); k++) {
        TTree *t = (TTree*)files[k]->Get("treeMaker/Events");

        if (k==0) isData = true;
        else isData = false;

        t->SetBranchAddress("hwwChi2",&hwwChi);
        t->SetBranchAddress("hbbMass",&hbbMass);
        t->SetBranchAddress("hhMass",&hhMass);
        t->SetBranchAddress("hwwPT",&hwwPt);
        t->SetBranchAddress("wjjTau2o1",&tau21);
        t->SetBranchAddress("ht",&ht);
        t->SetBranchAddress("ht_pup_30",&ht_pup);
        t->SetBranchAddress("ht_chs_20",&ht_chs);
        t->SetBranchAddress("hbbCSVCat",&hbbCSVCat);
        t->SetBranchAddress("nAK4Btags",&nAK4Btags);
        t->SetBranchAddress("passPre",&passPre);
        t->SetBranchAddress("isMuon",&isMuon);
        t->SetBranchAddress("lepETA",&lepETA);

        t->SetBranchAddress("passTrigCHS_30",&passTrigCHS);
        t->SetBranchAddress("passTrigPUP_30",&passTrigPUP_loose);
        t->SetBranchAddress("passTrigPUP_tight",&passTrigPUP_tight);

        if(!isData) {
        	t->SetBranchAddress("process",&process);
        	t->SetBranchAddress("xsec",&xsec);
        	t->SetBranchAddress("pu_N",&puN);
        	t->SetBranchAddress("trig_N",&trigN);
        	t->SetBranchAddress("hbbWQuark",&hbbWQuark);
        } else {
        	t->SetBranchAddress("dataset",&dataset);
        	weight = 1;
        }

        for (unsigned int i=0; i<t->GetEntries(); i++) {
            if (i%100000 == 0) printf("processing evt %d\n",i);
            t->GetEntry(i);
            TString pref = isData ? dataMap[dataset] : FillerConstants::MCProcessNames[process];
            if(!passTrigPUP_tight) continue;

            if(!isData) {
            	weight = xsec*trigN*puN;
            	if(passSR()) mkPlots(pref,isData,fileNames[k],"SR",hbbWQuark);
            }
            if(passQgCR())  mkPlots(pref,isData,fileNames[k],"qgCR",hbbWQuark);
            if(passTopCR()) mkPlots(pref,isData,fileNames[k],"topCR",hbbWQuark);

            if (passQgCR() && passTopCR()) printf("some shit is wronggg\n");

        }
    }
    plotter.write("outDebug.root");

}
