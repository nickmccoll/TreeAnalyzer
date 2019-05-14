#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

void drawEffs(TFile *f, bool isSig, TString pre, vector<TString> sels, TString numS, TString denS,
		TString canName="slurm", float rebin=0, int nR = 0, double * rebins = 0) {

	vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	Plotter *p = new Plotter();

	auto getTotBkgDist=[&](vector<pair<TString,TFile*>> bkgs, TString hist)->TH1* {

		if (bkgs.size() < 2) {cout<<"why u doin this? less than 2 samps"<<endl; return 0;}
		TH1 *h1 = (TH1F*)bkgs[0].second->Get(bkgs[0].first+"_"+hist);
		TH1 *hb = (TH1F*)h1->Clone();

		for (unsigned int i = 1; i < bkgs.size(); i++) {
			TH1 *h = (TH1*)bkgs[i].second->Get(bkgs[i].first+"_"+hist);
			hb->Add(h,1);
		}

		return hb;
	};

	auto getEffHistByMass=[&](TString num, TString den, TFile *f, vector<int> masses) {
		TH1 *hs = new TH1F(num+den,num+den,5000,0,5000);
		for (const auto& mass : masses) {
			TString m = TString::Format("m%i_",mass);
			double numF = ((TH1*)f->Get(m+num))->Integral();
			double denF = ((TH1*)f->Get(m+den))->Integral();
			hs->Fill(mass,numF/denF);
		}
		return hs;
	};

	if (isSig) {
		for (const auto& sel : sels) {
			TH1 *h = getEffHistByMass(pre+sel+numS,denS,f,masses);

			p->addHist(h,sel,-1,1,4,20,1,true,false);
		}
		p->setXTitle("M_{X}");
		p->draw(false,canName);
	} else {
		vector<pair<TString,TFile*>> bkgs = { {"ttbar",f}, {"wjets",f}/*, {"zjets",f}*/ };
		TH1 *hd = getTotBkgDist(bkgs,denS);
		if (hd==0) {cout<<"bad den: "<<denS<<endl; return;}

		hd = (TH1*)hd->Clone();
		PlotTools::toOverflow(hd);
		PlotTools::toUnderflow(hd);
	    if(rebin > 0) PlotTools::rebin(hd,rebin);
	    else if(rebins) hd = PlotTools::rebin(hd,nR,rebins);

		for (const auto& sel : sels) {
			TH1 *hn = getTotBkgDist(bkgs,pre+sel+numS);
	    	if (hn==0) {cout<<"bad num: "<<pre+sel+numS<<endl; continue;}

	    	hn = (TH1*)hn->Clone();
	    	PlotTools::toOverflow(hn);
	    	PlotTools::toUnderflow(hn);

	        if(rebin > 0) PlotTools::rebin(hn,rebin);
	        else if(rebins) hn = PlotTools::rebin(hn,nR,rebins);

	        hn->Divide(hn,hd,1,1,"b");
	        p->addHist(hn,sel);
		}
	    p->draw(false,canName);
	}
}

void drawDileptonSelection() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/";
	TFile *fbgv = new TFile(prePath+"llSel_blk.root");
	TFile *frad = new TFile(prePath+"llSel_rad.root");
	TFile *fbkg = new TFile(prePath+"llSel_bkg.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double hhBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000};
	int nhhBins = 17;
	double htBins[] = {400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000};
	int nhtBins = 14;

	vector<pair<TString,TFile*>> sigSamps = { {"m1000",fbgv}, {"m3000",fbgv} };
	vector<pair<TString,TFile*>> bkgSamps = { {"ttbar",fbkg}/*, {"qcd",fqcd}, {"wjets",fwj}, {"singlet",fst} */};

	vector<TString> isos = {"miniIso","pfIso","trkIso"};
	vector<TString> elIds = {"MVA","L","M","T"};
	vector<TString> muIds = {"L","M","T"};

	vector<TString> sels = {"ID_e1_H_e2_L_m1_M_m2_M","ID_e1_M_e2_M_m1_M_m2_M",
			"ID_e1_T_e2_L_m1_M_m2_M","ID_e1_T_e2_T_m1_M_m2_M"};

//	drawEffs(fbgv, true,"ee_passIP_",sels,"_miniIso0p2_ht","ee_passIP_ID_incl_miniIso0p2_ht","el id sig");
//	drawEffs(fbkg,false,"ee_passIP_",sels,"_miniIso0p2_ht","baseline_evts","el id bkg",0,nhtBins,htBins);
	drawEffs(fbgv, true,"mumu_passIP_",sels,"_miniIso0p2_ht","mumu_passIP_ID_incl_miniIso0p2_ht","el id sig");
	drawEffs(fbkg,false,"mumu_passIP_",sels,"_miniIso0p2_ht","baseline_evts","el id bkg",0,nhtBins,htBins);

	sels = {"miniIso_e0p10_m0p25","miniIso_e0p10_m0p20","pfIso_e0p10_m0p20","pfIso_e0p10_m0p25"};
//	drawEffs(fbgv, true,"mumu_passIP_ID_eMVA_mM_",sels,"_ht", "mumu_passIP_ID_eMVA_mM_inclIso_ht","muon iso sig");
//	drawEffs(fbkg,false,"mumu_passIP_ID_eMVA_mM_",sels,"_ht","baseline_evts","muon iso bkg",0,nhtBins,htBins);



	return;
}
