#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

TH1* getTotBkgDist(vector<pair<TString,TFile*>> bkgs, TString hist) {

	if (bkgs.size() < 2) {cout<<"why u doin this? less than 2 samps"<<endl; return 0;}
	TH1 *h1 = (TH1F*)bkgs[0].second->Get(bkgs[0].first+"_"+hist);
	TH1 *hb = (TH1F*)h1->Clone();

	for (unsigned int i = 1; i < bkgs.size(); i++) {
		TH1 *h = (TH1*)bkgs[i].second->Get(bkgs[i].first+"_"+hist);
		hb->Add(h,1);
	}

	return hb;
}

TH1* getTotRadDist(TString hist) {

	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/lepSel17_lnulnu/";
	TFile *fW = new TFile(prePath+"llSel_rad.root");
	TFile *fT = new TFile(prePath+"llSel_bbtautau.root");

	TH1 *hW = (TH1*)fW->Get(hist);
	TH1 *hT = (TH1*)fT->Get(hist);
	if (!hW) cout<<"no W hist: "<<hist<<endl;
	if (!hT) cout<<"no T hist: "<<hist<<endl;
	if (!hW || !hT) return 0;

	hW->Add(hT,0.29);
	return hW;
}

void drawEffs(TFile *f, bool isSig, TString pre, vector<TString> sels, TString numS, TString denS,
		TString canName="slurm", float rebin=0, int nR = 0, double * rebins = 0) {

	vector<int> masses = {800,900,1000,2000,2500,3000};
//	vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	Plotter *p = new Plotter();

	auto getLegStr=[&](TString sel) {
		TString s = sel;
		if (pre.BeginsWith("ee")) {
			s.ReplaceAll("_m1_M_m2_M","");
			s.ReplaceAll("_m0p20","");
			s.ReplaceAll("ID_","");
		} else if (pre.BeginsWith("mumu")) {
			s.ReplaceAll("ID_e1_M_e2_M_","");
			s.ReplaceAll("_e0p20","");
		} else {
			cout<<"emu getLegStr"<<endl;
		}
		return s;
	};

	auto getEffHistByMass=[&](TString num, TString den, TFile *f, vector<int> masses) {
		TH1 *hs = new TH1F(num+den,num+den,5000,0,5000);
		for (const auto& mass : masses) {
			TString m = TString::Format("m%i_",mass);
			TH1 *hn = getTotRadDist(m+num);
			TH1 *hd = getTotRadDist(m+num);

//			TH1 *hn = (TH1*)f->Get(m+num);
//			TH1 *hd = (TH1*)f->Get(m+den);
			if (!hn || !hd) continue;
			double numF = hn->Integral();
			double denF = hd->Integral();
			hs->Fill(mass,numF/denF);
		}
		return hs;
	};

	if (isSig) {
		for (const auto& sel : sels) {
			TH1 *h = getEffHistByMass(pre+sel+numS,denS,f,masses);
			TString legS = getLegStr(sel);
			p->addHist(h,legS,-1,1,4,20,1,true,false);
			delete h;
		}
		p->setXTitle("M_{X}");
		p->draw(false,canName);
	} else {
		vector<pair<TString,TFile*>> bkgs = { {"ttbar",f}, {"wjets",f}/*, {"zjets",f}*/ };
//		TH1 *hd = getTotBkgDist(bkgs,denS);
//		if (hd==0) {cout<<"bad den: "<<denS<<endl; return;}

//		hd = (TH1*)hd->Clone();
//		PlotTools::toOverflow(hd);
//		PlotTools::toUnderflow(hd);
//	    if(rebin > 0) PlotTools::rebin(hd,rebin);
//	    else if(rebins) hd = PlotTools::rebin(hd,nR,rebins);

		for (const auto& sel : sels) {
			TH1 *hn = getTotBkgDist(bkgs,pre+sel+numS);
	    	if (hn==0) {cout<<"bad num: "<<pre+sel+numS<<endl; continue;}

	    	hn = (TH1*)hn->Clone();
	    	PlotTools::toOverflow(hn);
	    	PlotTools::toUnderflow(hn);

	        if(rebin > 0) PlotTools::rebin(hn,rebin);
	        else if(rebins) hn = PlotTools::rebin(hn,nR,rebins);

//	        hn->Divide(hn,hd,1,1,"b");
	        TString legS = getLegStr(sel);
	        p->addHist(hn,legS);
		}
	    p->draw(false,canName);
	}
}

void drawDileptonSelection() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/lepSel17_lnulnu/";
	TFile *fbgv = new TFile(prePath+"llSel_blk.root");
	TFile *frad = new TFile(prePath+"llSel_rad.root");
	TFile *ftautau = new TFile(prePath+"llSel_bbtautau.root");
	TFile *fcomb = new TFile(prePath+"llSel_radComb.root");

	TFile *fbkg = new TFile(prePath+"llSel_fake.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double htBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000};
	int nhtBins = 17;
	double hhBins[] = {500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2500,3000};
	int nhhBins = 14;

	vector<pair<TString,TFile*>> sigSamps = { {"m1000",fbgv}, {"m3000",fbgv} };
	vector<pair<TString,TFile*>> bkgSamps = { {"ttbar",fbkg}, {"wjets",fbkg} };

	vector<TString> isos = {"miniIso","pfIso","trkIso"};
	vector<TString> elIds = {"MVA","L","M","T"};
	vector<TString> muIds = {"L","M","T"};

	vector<TString> elSels1 = {"ID_e1_M_e2_L_m1_M_m2_M","ID_e1_M_e2_M_m1_M_m2_M",
			"ID_e1_T_e2_L_m1_M_m2_M","ID_e1_M_e2_MVA_m1_M_m2_M"};

	vector<TString> elSels2 = {"ID_e1_L_e2_L_m1_M_m2_M","ID_e1_M_e2_M_m1_M_m2_M",
			"ID_e1_T_e2_T_m1_M_m2_M","ID_e1_MVA_e2_MVA_m1_M_m2_M"};

	vector<TString> elSels3 = {"ID_e1_MVA_e2_MVA_m1_M_m2_M","ID_e1_MVA_e2_M_m1_M_m2_M",
			"ID_e1_M_e2_MVA_m1_M_m2_M","ID_e1_M_e2_M_m1_M_m2_M"};

	vector<TString> muSels1 = {"ID_e1_MVA_e2_MVA_m1_L_m2_L","ID_e1_MVA_e2_MVA_m1_M_m2_M",
			"ID_e1_MVA_e2_MVA_m1_L_m2_M","ID_e1_MVA_e2_MVA_m1_M_m2_L","ID_e1_MVA_e2_MVA_m1_S_m2_S"};

	vector<vector<TString>> elSels = {elSels1,elSels2,elSels3};
	vector<vector<TString>> muSels = {muSels1};

//	for (unsigned int i = 0; i < elSels.size(); i++) {
//		TString iS = TString::Format("%d",i);
//		drawEffs(fbgv, true,"ee_passIP_",elSels[i],"_miniIso0p2_pass2lSel_ht","ee_passIP_ID_incl_miniIso0p2_pass2lSel_ht","el id sig "+iS);
//		drawEffs(ftautau, true,"ee_passIP_",elSels[i],"_miniIso0p2_pass2lSel_ht","ee_passIP_ID_incl_miniIso0p2_pass2lSel_ht","el id sig bbtautau "+iS);
//
//		drawEffs(fbkg,false,"ee_passIP_",elSels[i],"_miniIso0p2_pass2lSel_mhh","","el id bkg "+iS,0,nhhBins,hhBins);
//	}
//	for (unsigned int i = 0; i < muSels.size(); i++) {
//		TString iS = TString::Format("%d",i);
//		drawEffs(fbgv, true,"mumu_passIP_",muSels[i],"_miniIso0p2_pass2lSel_ht","mumu_passIP_ID_incl_miniIso0p2_pass2lSel_ht","mu id sig "+iS);
//		drawEffs(ftautau, true,"mumu_passIP_",muSels[i],"_miniIso0p2_pass2lSel_ht","mumu_passIP_ID_incl_miniIso0p2_pass2lSel_ht","mu id sig bbtautau "+iS);
//		drawEffs(fbkg,false,"mumu_passIP_",muSels[i],"_miniIso0p2_pass2lSel_mhh","","mu id bkg "+iS,0,nhhBins,hhBins);
//	}

//	elSels1 = {"miniIso_e0p10_m0p20","miniIso_e0p15_m0p20","miniIso_e0p20_m0p20","miniIso_e0p25_m0p20"};
//	muSels1 = {"miniIso_e0p20_m0p10","miniIso_e0p20_m0p15","miniIso_e0p20_m0p20","miniIso_e0p20_m0p25"};
//	vector<TString> muSels2 = {"miniIso_e0p20_m0p20","pfIso_e0p20_m0p40","miniIso_e0p20_m0p40","pfIso_e0p20_m0p20"};
//
//	elSels = {elSels1};
//	muSels = {muSels1,muSels2};
//
//	for (unsigned int i = 0; i < elSels.size(); i++) {
//		TString iS = TString::Format("%d",i);
//		drawEffs(fbgv, true,"ee_passIP_ID_eMVA_mM_",elSels[i],"_pass2lSel_ht","ee_passIP_ID_eMVA_mM_inclIso_pass2lSel_ht","el iso sig "+iS);
//		drawEffs(fbkg,false,"ee_passIP_ID_eMVA_mM_",elSels[i],"_pass2lSel_mhh","","el iso bkg "+iS,0,nhhBins,hhBins);
//	}
//	for (unsigned int i = 0; i < muSels.size(); i++) {
//		TString iS = TString::Format("%d",i);
//		drawEffs(fbgv, true,"mumu_passIP_ID_eMVA_mM_",muSels[i],"_pass2lSel_ht","mumu_passIP_ID_eMVA_mM_inclIso_pass2lSel_ht","mu iso sig "+iS);
//		drawEffs(fbkg,false,"mumu_passIP_ID_eMVA_mM_",muSels[i],"_pass2lSel_mhh","","mu iso bkg "+iS,0,nhhBins,hhBins);
//	}


//	drawEffs(fbgv, true,"mumu_passIP_",sels,"_miniIso0p2_ht","mumu_passIP_ID_incl_miniIso0p2_ht","el id sig");
//	drawEffs(fbkg,false,"mumu_passIP_",sels,"_miniIso0p2_ht","mumu_passIP_ID_incl_miniIso0p2_ht","el id bkg",0,nhtBins,htBins);

//	drawEffs(fbgv, true,"mumu_passIP_ID_eMVA_mM_",sels,"_ht", "mumu_passIP_ID_eMVA_mM_inclIso_ht","muon iso sig");
//	drawEffs(fbkg,false,"mumu_passIP_ID_eMVA_mM_",sels,"_ht","baseline_evts","muon iso bkg",0,nhtBins,htBins);

	vector<TString> vars = {"sip","met"};
	TString hist = "sip_9999p00_ID_eMVA_mM_miniIso0p2_pass2lSel_";
	TString hist_old = "sip_4p00_ID_eMVA_mM_miniIso0p2_pass2lSel_";
	for (const auto& var : vars) {
		Plotter *p= new Plotter;
		Plotter *pb = new Plotter();
		Plotter *ps = new Plotter();

		TH1F *hb = (TH1F*)getTotBkgDist(bkgSamps,hist+var);
		pb->addHistLine(hb,"bkg");
		TH1F *hb_old = (TH1F*)getTotBkgDist(bkgSamps,hist_old+var);
		pb->addHistLine(hb_old,"bkg SIP < 4");

		if (var=="met") {
			cout<<"new bkg"<<endl<<hist<<endl<<hb->Integral()<<endl<<endl;
			cout<<"old bkg"<<endl<<hist_old<<endl<<hb_old->Integral()<<endl<<endl;
		}

		vector<TString> masses = {"m1000","m3000"};
		for (const auto& mass : masses) {
			TH1F *h1 = (TH1F*)frad->Get(mass+"_"+hist+var);
			TH1F *h2 = (TH1F*)ftautau->Get(mass+"_"+hist+var);

			TH1F *hh1 = (TH1F*)h1->Clone();
			h1->Add(h2,0.29); // applying ratio of B(H->tautau) / B(H->WW)

			if (var=="met") {
				cout<<"new "<<mass<<endl<<hist<<endl<<h1->Integral()<<endl<<endl;
				cout<<"old "<<mass<<endl<<hist_old<<endl<<hh1->Integral()<<endl<<endl;
			}

			TGraph *g = PlotTools::getRocCurve(h1,hb,false,"sig","bkg");
			p->addGraph(g,mass);
			ps->addHistLine(h1,mass+" w TauTau");
			ps->addHistLine(hh1,mass+" bbWW");

		}
//		p->draw(false,"ROC "+var);
		pb->draw(false,"bkg "+var);
		ps->draw(false,"sig "+var);
	}

//	TH1 *h1 = (TH1*)frad->Get("m1000_"+hist+"sip");
//	TH1 *h1t = (TH1*)ftautau->Get("m1000_"+hist+"sip");
//	h1->Add(h1t,0.29);
//
//	TH1 *h3 = (TH1*)frad->Get("m3000_"+hist+"sip");
//	TH1 *h3t = (TH1*)ftautau->Get("m3000_"+hist+"sip");
//	h3->Add(h3t,0.29);
//
//	TH1 *hb = (TH1*)getTotBkgDist(bkgSamps,hist+"sip");
//
//	h1 = (TH1*)h1->Clone();
//	h3 = (TH1*)h3->Clone();
//	hb = (TH1*)hb->Clone();
//
//	PlotTools::toOverflow(h1);
//	PlotTools::toOverflow(h3);
//	PlotTools::toOverflow(hb);
//
//	TH1* i1 = PlotTools::getIntegral(h1,false,true);
//	TH1* i3 = PlotTools::getIntegral(h3,false,true);
//	TH1* ib = PlotTools::getIntegral(hb,false,true);
//
//	Plotter *p1 = new Plotter();
//	p1->addHist(i1,"");
//	p1->draw(false,"m1000 integral");
//	Plotter *p3 = new Plotter();
//	p3->addHist(i3,"");
//	p3->draw(false,"m3000 integral");
//	Plotter *pb = new Plotter();
//	pb->addHist(ib,"");
//	pb->draw(false,"fake bkg integral");

	return;
}
