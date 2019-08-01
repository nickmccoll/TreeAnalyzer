#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

TH1* getTotBkgDist(vector<pair<TString,TFile*>> bkgs, TString hist) {

//	if (bkgs.size() < 2) {cout<<"why u doin this? less than 2 samps"<<endl; return 0;}
	TH1 *h1 = (TH1F*)bkgs[0].second->Get(bkgs[0].first+"_"+hist);
	if (!h1) {
		cout<<"bad hist: "<<bkgs[0].first+"_"+hist<<endl;
		return 0;
	}
	TH1 *hb = (TH1F*)h1->Clone();

	for (unsigned int i = 1; i < bkgs.size(); i++) {
		TH1 *h = (TH1*)bkgs[i].second->Get(bkgs[i].first+"_"+hist);
		if (!h) {
			cout<<"bad hist: "<<bkgs[i].first+"_"+hist<<endl;
			continue;
		}
		hb->Add(h,1);
	}
	return hb;
}

TH1* getTotRadDist(TFile *fW, TFile *fT, TString hist) {

	TH1 *hW = (TH1*)fW->Get(hist);
	TH1 *hT = (TH1*)fT->Get(hist);
	if (!hW) cout<<"no bbWW hist: "<<hist<<endl;
	if (!hT) cout<<"no bbTT hist: "<<hist<<endl;
	if (!hW || !hT) return 0;

	hW->Add(hT,0.26);
	return hW;
}

void drawSelVars() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/Vars/";
	TFile *fww = new TFile(prePath+"bbWW_spin0.root");
	TFile *ftt = new TFile(prePath+"bbtt_spin0.root");
	TFile *fbkg = new TFile(prePath+"bkg.root");

	TFile *fwwO = new TFile(prePath+"bbWW_spin0_old.root");
	TFile *fttO = new TFile(prePath+"bbtt_spin0_old.root");
	TFile *fbkgO = new TFile(prePath+"bkg_old.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double htBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000};
	int nhtBins = 17;
	double hhBins[] = {700,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2500,3000};
	int nhhBins = 12;

	vector<pair<TString,TFile*>> FbkgSamps = { {"ttbar1",fbkg}, {"wjets",fbkg} };
	vector<pair<TString,TFile*>> bkgSamps = { {"ttbar1",fbkg}, {"wjets",fbkg}, {"ttbar2",fbkg}, {"zjets",fbkg}, {"singlet",fbkg}, {"diboson",fbkg} };
	vector<pair<TString,TFile*>> otherBkgs = { {"wjets",fbkg}, {"singlet",fbkg}, {"diboson",fbkg}, {"ttx",fbkg}, {"hx",fbkg} };

	vector<TString> vars = {"met","metOht","metOrht","dRll","dEtall","dPhill","metOmhh","Mww","ptwwomhh"};
//	vector<TString> vars = {"Mbb","Mhh"};
	TString histS = "SR_fullSel_";

	for (const auto& var : vars) {
		Plotter *p = new Plotter();

		TH1 *hww = (TH1*)fww->Get("m1000_"+histS+var);
		TH1 *htt = (TH1*)ftt->Get("m1000_"+histS+var);
		TH1 *hww3 = (TH1*)fww->Get("m3000_"+histS+var);
		TH1 *htt3 = (TH1*)ftt->Get("m3000_"+histS+var);

		TH1 *ht1 = (TH1*)fbkg->Get("ttbar1_"+histS+var);
		TH1 *ht2 = (TH1*)fbkg->Get("ttbar2_"+histS+var);
		TH1 *hz  = (TH1*)fbkg->Get("zjets_"+histS+var);
		TH1 *ho  = getTotBkgDist(otherBkgs,histS+var);

		p->addStackHist(ht1,"ttbar 1");
		p->addStackHist(ht2,"ttbar 2");
		p->addStackHist(hz,"zjets");
		p->addStackHist(ho,"other");

		p->addHistLine(hww,"bbWW 1000 GeV");
		p->addHistLine(htt,"bb#tau#tau 1000 GeV");
		p->addHistLine(hww3,"bbWW 3000 GeV");
		p->addHistLine(htt3,"bb#tau#tau 3000 GeV");

		p->normalize();
		p->draw(false,histS+var);

		ho->Add(ht1,1); ho->Add(ht2,1); ho->Add(hz,1);
		if (var=="Mww"||var=="yww"||var=="Mll") continue;
		bool cutGt = true;
		if (var == "dPhi_metll" || var == "dRll"|| var == "dEtall"|| var == "dPhill") cutGt = false;
		TGraph *gw = PlotTools::getRocCurve(hww,ho,cutGt,"sig","bkg");
		TGraph *gt = PlotTools::getRocCurve(htt,ho,cutGt,"sig","bkg");
		Plotter *pR = new Plotter();
		pR->addGraph(gw,"bbWW 1 TeV");
		pR->addGraph(gt,"bb#tau#tau 1 TeV");

		pR->draw(false,"ROC "+var);

	}

//	TString hhS = "SR_fullSel_Mhh";
//	vector<pair<TFile*,TFile*>> files = { {fww,fwwO}, {ftt,fttO}, {fbkg,fbkgO} };
//	vector<TString> ids = {"bbWW","bbtt","bkg"};
//	vector<pair<TString,TFile*>> allBkgs = { {"ttbar1",fbkg},{"ttbar2",fbkg},{"zjets",fbkg},{"wjets",fbkg}, {"singlet",fbkg}, {"diboson",fbkg}, {"ttx",fbkg}, {"hx",fbkg} };
//	vector<pair<TString,TFile*>> allBkgsO = { {"ttbar1",fbkgO},{"ttbar2",fbkgO},{"zjets",fbkgO},{"wjets",fbkgO}, {"singlet",fbkgO}, {"diboson",fbkgO}, {"ttx",fbkgO}, {"hx",fbkgO} };
//
//	for (unsigned int i=0; i<ids.size(); i++) {
//		Plotter *pp = new Plotter();
//		TH1 *hold=0;
//		TH1 *hnew=0;
//
//		if (i==(files.size()-1)) {
//			hold = getTotBkgDist(allBkgsO,hhS);
//			hnew = getTotBkgDist(allBkgs,hhS);
//		} else {
//			hold = (TH1*)files[i].first->Get("m1000_"+hhS);
//			hnew = (TH1*)files[i].second->Get("m1000_"+hhS);
//		}
//
//		pp->addHistLine(hold,"MET > 40");
//		pp->addHistLine(hnew,"MET / Mhh > 0.05");
//		pp->drawSplitRatio(0,"Met comparison "+ids[i],false,false,"Met comparison "+ids[i]);
//	}

	return;
}
