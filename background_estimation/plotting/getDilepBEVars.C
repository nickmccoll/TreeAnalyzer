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

    map<UChar_t,TString> dataMap = { {7,"JetHT"}, {8,"MET"}, {11,"SingleElectron"}, {12,"SingleMuon"}, {13,"SinglePhoton"} };

    auto passBaseline = [&]() {
    	if (!passPre) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (dilepMass < 6 || dilepMass > 75) return false;
    	if (dilepDR > 1.0) return false;
    	if (met / hhMass < 0.1) return false;
//    	if (met < 40) return false;
    	if (fabs(dPhi_metll) > TMath::PiOver2()) return false;
//    	if (ptww / hhMass < 0.25) return false;
    	return true;
    };

    auto passSR = [&]() {
    	if (!passBaseline()) return false;
    	if (nAK4Btags != 0) return false;
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

    auto passTopDrCR = [&]() {
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

    auto pltVars = [&](TString pref, TString bCat, TString fN) {
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mhh",";M_{HH}",132,700,4000)->Fill(hhMass,weight);

//        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mll",";M_{ll}",200,0,200)->Fill(dilepMass,weight);
//        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_dRll",";#DeltaR_{ll}",100,0,5)->Fill(dilepDR,weight);
//        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_ptom",";MET",200,0,1)->Fill(ptww/hhMass,weight);
//        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_dPhi_metll",";dPhi",100,-3.14,3.14)->Fill(dPhi_metll,weight);
//        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_hbbDecayType",";",8,-0.5,7.5)->Fill(hbbDecayTypeMC,weight);

    };

    auto mkPlots = [&](TString procName, bool isData, TString fN, TString region) {

//        TString bCat;
//        if (hbbCSVCat==3) bCat = "LL";
//        else if (hbbCSVCat==4) bCat = "bL";
//        else if (hbbCSVCat==5) bCat = "bM";
//        else if (hbbCSVCat==6) bCat = "bT";

    	TString ttS = TString::Format("%d",int(nLepsTT));
    	TString numbS = TString::Format("nb%d_",numBinHbb);

    	if (hbbDecayTypeMC==0) numbS = "bkgRad_";
    	else                   numbS = "bkgQk_";

    	pltVars(procName,region,fN);
    	if (procName=="ttbar") pltVars(procName+ttS,region,fN);

    	if (isData) pltVars("data",region,fN);
    	else        pltVars("bkg" ,region,fN);

        pltVars(procName,region,numbS+fN);
    	if (procName=="ttbar") pltVars(procName+ttS,region,numbS+fN);

        if (isData) pltVars("data",region,numbS+fN);
        else        pltVars("bkg" ,region,numbS+fN);

    };

    bool isData, isSignal;
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
                if(process == 2) weight *= 0.831493; // ttbar scale factor

            	if(passSR()) mkPlots(pref,isData,fileNames[k],"SR");
            }

            if(passQgCR())  mkPlots(pref,isData,fileNames[k],"qgCR");
            if(passTopCR()) mkPlots(pref,isData,fileNames[k],"topCR");
            if(passTopDrCR()) mkPlots(pref,isData,fileNames[k],"TopDrCR");

//            if(passStudyR()) mkPlots(pref,isData,fileNames[k],"specR");
        }
    }
    plotter.write("debug.root");

}
