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

void makeSignalPlots(TFile *f, TString preN, vector<TString> selsN, TString postN,
		TString denS, TString canName="slurm") {

	vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500};
	if (canName.Contains("bbtt")) masses = {800,900,1000,2000,2500,3000};

	Plotter *p = new Plotter();

	auto getLegStr=[&](TString sel) {
		TString s = sel;
		if (preN.BeginsWith("ee")) {
			s.ReplaceAll("_m1_L_m2_L","");
			s.ReplaceAll("_m0p20","");
			s.ReplaceAll("ID_","");
		} else if (preN.BeginsWith("mumu")) {
			s.ReplaceAll("ID_e1_MVA_e2_MVA_","");
			s.ReplaceAll("_e0p20","");
		} else {
			cout<<"emu getLegStr"<<endl;
		}
		return s;
	};

	auto getEffHistByMass=[&](TString num, TString den, vector<int> masses) {
		TH1F *hs = new TH1F(num+canName,num+canName,5000,0,5000);
		for (const auto& mass : masses) {
			TString m = TString::Format("m%i_",mass);
			TH1 *hn = (TH1*)f->Get(m+num);
			TH1 *hd = (TH1*)f->Get(m+den);

//			TH1 *hn = getTotRadDist(fWW,ftt,m+num);
//			TH1 *hd = getTotRadDist(fWW,ftt,m+den);
			if (!hn || !hd) continue;
			double numF = hn->Integral(0,hn->GetNbinsX()+1);
			double denF = hd->Integral(0,hd->GetNbinsX()+1);
//			cout<<"num = "<<numF<<", den = "<<denF<<endl;
			delete hn; delete hd;
			hs->Fill(mass,numF/denF);
		}
		return hs;
	};

	for (const auto& sel : selsN) {
		TH1F *h = getEffHistByMass(preN+sel+postN,denS,masses);
		TString legS = getLegStr(sel);
		p->addHist(h,legS,-1,1,4,20,1,true,false);
		delete h;
	}
	p->setXTitle("M_{X}");
	p->draw(false,canName);

}

void makeBkgPlots(vector<pair<TString,TFile*>> bkgs, TString pre, vector<TString> sels, TString post,
		TString canName="slurm",float rebin=0, int nR = 0, double * rebins = 0) {

	Plotter *p = new Plotter();

	auto getLegStr=[&](TString sel) {
		TString s = sel;
		if (pre.BeginsWith("ee")) {
			s.ReplaceAll("_m1_L_m2_L","");
			s.ReplaceAll("_m0p20","");
			s.ReplaceAll("ID_","");
		} else if (pre.BeginsWith("mumu")) {
			s.ReplaceAll("ID_e1_MVA_e2_MVA_","");
			s.ReplaceAll("_e0p20","");
		} else {
			cout<<"emu getLegStr"<<endl;
		}
		return s;
	};

	for (const auto& sel : sels) {
		TH1 *hn = getTotBkgDist(bkgs,pre+sel+post);
    	if (hn==0) {cout<<"bad num: "<<pre+sel+post<<endl; continue;}

    	hn = (TH1*)hn->Clone();
    	PlotTools::toOverflow(hn);
    	PlotTools::toUnderflow(hn);

        if(rebin > 0) PlotTools::rebin(hn,rebin);
        else if(rebins) hn = PlotTools::rebin(hn,nR,rebins);

        TString legS = getLegStr(sel);
        p->addHist(hn,legS);
	}
    p->draw(false,canName);

}

void drawDileptonSelection() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/lepSel17_lnulnu/";
	TFile *frad = new TFile(prePath+"bbWW_spin0.root");
	TFile *ftt = new TFile(prePath+"bbtt_spin0.root");
	TFile *fbkg = new TFile(prePath+"llSel_bkg.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double htBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000};
	int nhtBins = 17;
	double hhBins[] = {700,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2500,3000};
	int nhhBins = 12;

	vector<pair<TString,TFile*>> sigSamps = { {"m1000",frad}, {"m3000",frad} };
	vector<pair<TString,TFile*>> FbkgSamps = { {"ttbar1",fbkg}, {"wjets",fbkg} };
	vector<pair<TString,TFile*>> bkgSamps = { {"ttbar1",fbkg}, {"wjets",fbkg}, {"ttbar2",fbkg}, {"zjets",fbkg}, {"singlet",fbkg}, {"diboson",fbkg} };

	vector<TString> elSels = {};
	vector<TString> muSels = {};

// ----------------------------------------------------------------------------------

	elSels = {"ID_e1_L_e2_L_m1_L_m2_L","ID_e1_M_e2_M_m1_L_m2_L",
			"ID_e1_T_e2_T_m1_L_m2_L","ID_e1_MVA_e2_MVA_m1_L_m2_L"};

	makeSignalPlots(frad,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbWW eID 1");

	makeSignalPlots(ftt,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt eID 1");

	makeBkgPlots(FbkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","fake bkg eID 1",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","tot bkg eID 1",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	elSels = {"ID_e1_MVA_e2_MVA_m1_L_m2_L","ID_e1_MVA_e2_M_m1_L_m2_L",
			"ID_e1_M_e2_MVA_m1_L_m2_L","ID_e1_M_e2_M_m1_L_m2_L"};

	makeSignalPlots(frad,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbWW eID 2");

	makeSignalPlots(ftt,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt eID 2");

	makeBkgPlots(FbkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","fake bkg eID 2",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","tot bkg eID 2",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	muSels = {"ID_e1_MVA_e2_MVA_m1_L_m2_L","ID_e1_MVA_e2_MVA_m1_M_m2_M",
			"ID_e1_MVA_e2_MVA_m1_L_m2_M","ID_e1_MVA_e2_MVA_m1_M_m2_L","ID_e1_MVA_e2_MVA_m1_S_m2_S"};

	makeSignalPlots(frad,"mumu_passIP_",muSels,"_miniIso0p2_pass2lSel_sip1",
			"mumu_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbWW muID 1");

	makeSignalPlots(ftt,"mumu_passIP_",muSels,"_miniIso0p2_pass2lSel_sip1",
			"mumu_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt muID 1");

	makeBkgPlots(FbkgSamps,"mumu_passIP_",muSels,"_miniIso0p2_pass2lSel_mhh","fake bkg muID 1",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"mumu_passIP_",muSels,"_miniIso0p2_pass2lSel_mhh","tot bkg muID 1",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	muSels = {"miniIso_e0p20_m0p20","pfIso_e0p20_m0p20","miniIso_e0p20_m0p15","miniIso_e0p20_m0p10"};

	makeSignalPlots(frad,"mumu_passIP_ID_eMVA_mL_",muSels,"_pass2lSel_sip1",
			"mumu_passIP_ID_eMVA_mL_inclIso_pass2lSel_sip1","bbWW muISO 1");

	makeSignalPlots(ftt,"mumu_passIP_ID_eMVA_mL_",muSels,"_pass2lSel_sip1",
			"mumu_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt muISO 1");

	makeBkgPlots(FbkgSamps,"mumu_passIP_ID_eMVA_mL_",muSels,"_pass2lSel_mhh","fake bkg muISO 1",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"mumu_passIP_ID_eMVA_mL_",muSels,"_pass2lSel_mhh","tot bkg muISO 1",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	elSels = {"miniIso_e0p20_m0p20","pfIso_e0p20_m0p20","miniIso_e0p15_m0p20","miniIso_e0p10_m0p20",
			"pfIso_e0p25_m0p20"};

	makeSignalPlots(frad,"ee_passIP_ID_eMVA_mL_",elSels,"_pass2lSel_sip1",
			"ee_passIP_ID_eMVA_mL_inclIso_pass2lSel_sip1","bbWW eISO 1");

	makeSignalPlots(ftt,"ee_passIP_ID_eMVA_mL_",elSels,"_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt eISO 1");

	makeBkgPlots(FbkgSamps,"ee_passIP_ID_eMVA_mL_",elSels,"_pass2lSel_mhh","fake bkg eISO 1",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"ee_passIP_ID_eMVA_mL_",elSels,"_pass2lSel_mhh","tot bkg eISO 1",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	elSels = {"ID_e1_MVA_e2_MVA_m1_L_m2_L","ID_e1_MVA_e2_M_m1_L_m2_L",
			"ID_e1_M_e2_MVA_m1_L_m2_L","ID_e1_M_e2_M_m1_L_m2_L"};

	makeSignalPlots(frad,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbWW eID 2");

	makeSignalPlots(ftt,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_sip1",
			"ee_passIP_ID_incl_miniIso0p2_pass2lSel_sip1","bbtt eID 2");

	makeBkgPlots(FbkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","fake bkg eID 2",0,nhhBins,hhBins);
	makeBkgPlots(bkgSamps,"ee_passIP_",elSels,"_miniIso0p2_pass2lSel_mhh","tot bkg eID 2",0,nhhBins,hhBins);

// ----------------------------------------------------------------------------------

	return;
}
