#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "Configuration/interface/FillerConstants.h"
#include "TFile.h"
#include "TH1.h"

using namespace std;

void plotNm1Eff(TFile *f, TString smp, TString denS, vector<TString> numS, TString var, TString canName="slurm", float rebin=0, int nR = 0, double * rebins = 0) {
	Plotter *p = new Plotter();

	TString denHist = smp+"_"+denS+"_"+var;
	TH1 *hd = (TH1*)f->Get(denHist);
	if (hd==0) {cout<<"bad den: "<<denHist<<endl; return;}

	hd = (TH1*)hd->Clone();
	PlotTools::toOverflow(hd);
	PlotTools::toUnderflow(hd);

    if(rebin > 0) PlotTools::rebin(hd,rebin);
    else if(rebins) hd = PlotTools::rebin(hd,nR,rebins);

//    p->addHist(hd,"den");

    for (const auto& num : numS) {
    	TString numHist = smp+"_"+denS+"_"+num+"_"+var;
    	TH1 *hn = (TH1F*)f->Get(numHist);
    	if (hn==0) {cout<<"bad num: "<<numHist<<endl; continue;}

    	hn = (TH1*)hn->Clone();
    	PlotTools::toOverflow(hn);
    	PlotTools::toUnderflow(hn);

        if(rebin > 0) PlotTools::rebin(hn,rebin);
        else if(rebins) hn = PlotTools::rebin(hn,nR,rebins);

        hn->Divide(hn,hd,1,1,"b");
        p->addHist(hn,num);
    }
    p->draw(false,canName);
}

void plotDistribution(TFile *fs, map<TFile*,TString> fbs, TString mass, TString hist, TString canName="slurm", float rebin=0, int nR=0, double *rebins = 0) {
	Plotter *p = new Plotter();
	TH1 *hs = (TH1*)fs->Get(mass+"_"+hist);
	if (hs==0) {cout<<"bad sig: "<<mass+"_"+hist<<endl; return;}

	hs = (TH1*)hs->Clone();
	PlotTools::toOverflow(hs);
	PlotTools::toUnderflow(hs);

    if(rebin > 0) PlotTools::rebin(hs,rebin);
    else if(rebins) hs = PlotTools::rebin(hs,nR,rebins);

    p->addHistLine(hs,mass);

    for (const auto& f : fbs) {
    	TH1 *h = (TH1*)f.first->Get("bkg_"+hist);
    	if (h==0) {cout<<"bad bkg: "<<"bkg_"+hist<<endl; return;}
    	h = (TH1*)h->Clone();
    	PlotTools::toOverflow(h);
    	PlotTools::toUnderflow(h);

        if(rebin > 0) PlotTools::rebin(h,rebin);
        else if(rebins) h = PlotTools::rebin(h,nR,rebins);
        p->addStackHist(h,f.second);
    }
    p->normalize();
    p->draw(false,canName);

}

TH1* getHistByMass(TString pre, TFile *f, vector<int> masses, TString sel, TString post) {
	TH1F *hs = new TH1F(pre+sel,pre+sel,100,0,5000);
	for (const auto& mass : masses) {
		TString m = TString::Format("m%i_",mass);
//		cout<<m+pre+"_"+sel+post+"_pt"<<endl;
		double num = ((TH1F*)f->Get(m+pre+"_"+sel+post+"_pt"))->Integral();
		double den;
//		if (post=="") den = ((TH1F*)f->Get(m+pre+"_pt"))->Integral();
//		else          den = ((TH1F*)f->Get(m+pre+"_"+post+"_pt"))->Integral();
//		cout<<m+pre+"_maxEta2p5_pt"<<endl;
		den = ((TH1F*)f->Get(m+pre+"_maxEta2p5_pt"))->Integral();
		hs->Fill(mass,num/den);
	}
	return hs;
}

void drawElectronIDvars() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/";
//	TFile *frad = new TFile(prePath+"eID_radion.root");
	TFile *fbgv = new TFile(prePath+"eID_bulkgrav.root");
	TFile *ftt = new TFile(prePath+"eID_ttbar.root");
	TFile *fwj = new TFile(prePath+"eID_wjets.root");
	TFile *fqcd = new TFile(prePath+"eID_qcd.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double hhBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3500,5000};
	int nhhBins = 18;

	vector<TString> selVars = {"relax_hoe","relax_dEtaSeed","relax_dPhiIn","relax_sigmaIetaIeta","relax_1oEm1op","TightID","MVAID"};
//	plotNm1Eff(fbgv,"m1200","passIPandISO",selVars,"pt","BG 1.2 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);
//	plotNm1Eff(fbgv,"m2000","passIPandISO",selVars,"pt","BG 2 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);
//	plotNm1Eff(frad,"m1400","passIPandISO",selVars,"pt","RAD 1.4 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);
//	plotNm1Eff(fbgv,"m4500","passIPandISO",selVars,"pt","BG 4.5 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);

	map<TFile*,TString> sigSamps = { /*{frad,"radion"},*/ {fbgv,"graviton"} };
	map<TFile*,TString> bkgSamps = { {fqcd,"QCD"} };
//	plotDistribution(fbgv,bkgs,"m4500","passIPandISO_hoe");

	vector<TFile*> files = {fqcd};
	vector<TString> name = {"QCD"};
	vector<TString> isos = {"miniIso0p1_","miniIso0p2_","miniIsoFP0p1_","miniIsoFP0p2_"};
	vector<TString> ids = {"MedID","TightID","MVAID"};
	vector<TString> etas = {"","maxEta1p5_","maxEta2p1_"};
	vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};

	// ID --------------------------------------- //
//	for (const auto& f : sigSamps) {
//		Plotter *p = new Plotter();
//		for (const auto& id : ids) {
//			TH1 *h = getHistByMass("passIP_maxEta2p1",f.first,masses,"miniIso0p1_"+id,"");
//			p->addHist(h,id+"+miniIso 0.1",-1,1,4,20,1,true,false);
//		}
//		p->setXTitle("M_{X}");
//		p->draw(false,f.second+" ID efficiencies: absEta < 2.1, miniIso<0.1");
//	}
//	for (unsigned int i=0; i<files.size(); i++) {
//		Plotter *p = new Plotter();
//		TH1 *hd = (TH1*)(files[i])->Get("bkg_passIP_maxEta2p1_hhmass");
//		cout<<endl<<hd->Integral()<<endl;
//		hd = (TH1*)hd->Clone();
//		PlotTools::toOverflow(hd);
//		hd = PlotTools::rebin(hd,nhhBins,hhBins);
//
//		p->addHist(hd,"den");
//		delete hd;
//		for (const auto& id : ids) {
//			TH1 *hn = (TH1*)(files[i])->Get("bkg_passIP_maxEta2p1_miniIso0p1_"+id+"_hhmass");
//			if (hn==0) cout<<"bad num: "<<"bkg_passIP_maxEta2p1_miniIso0p1_"+id+"_hhmass"<<endl;
//			cout<<hn->Integral()<<endl;
//			hn = (TH1*)hn->Clone();
//			PlotTools::toOverflow(hn);
//			hn = PlotTools::rebin(hn,nhhBins,hhBins);
//			p->addHist(hn,id+"+miniIso 0.1");
//			delete hn;
//		}
//		p->setXTitle("M_{HH}");
//		p->drawRatio(0,"stack",true,false,name[i]+" ID efficiencies: absEta < 2.1, miniIso<0.1");
//	}

	// ISO ---------------------------------------- //
//	for (const auto& f : sigSamps) {
//		Plotter *p = new Plotter();
//		for (const auto& iso : isos) {
//			TH1 *h = getHistByMass("passIP_maxEta2p1",f.first,masses,iso+"MVAID","");
//			p->addHist(h,iso+"+ MVA90 ID",-1,1,4,20,1,true,false);
//		}
//		p->setXTitle("M_{X}");
//		p->draw(false,f.second+" Iso efficiencies");
//	}
//	for (const auto& b : bkgSamps) {
//		Plotter *p = new Plotter();
//		TH1 *hd = (TH1*)b.first->Get("bkg_passIP_maxEta2p1_hhmass");
//		hd = (TH1*)hd->Clone();
//		PlotTools::toOverflow(hd);
//		hd = PlotTools::rebin(hd,nhhBins,hhBins);
//
//		p->addHist(hd,"den");
//
//		for (const auto& iso : isos) {
//			TH1 *hn = (TH1*)b.first->Get("bkg_passIP_maxEta2p1_"+iso+"MVAID_hhmass");
//			if (hn==0) cout<<"bad num: "<<"bkg_passIP_maxEta2p1_"+iso+"MVAID_hhmass"<<endl;
//			hn = (TH1*)hn->Clone();
//			PlotTools::toOverflow(hn);
//			hn = PlotTools::rebin(hn,nhhBins,hhBins);
//			p->addHist(hn,iso+" + MVA90 ID");
//		}
//		p->drawRatio(0,"stack",true,false,b.second+" Iso efficiencies");
//	}

	// Any combination of numerators
	vector<TString> combs = {"maxEta2p5_miniIso0p2_TightID","maxEta2p1_miniIsoFP0p1_MVAID","maxEta1p5_miniIso0p2_MVAID","maxEta2p1_miniIso0p1_MVAID"};


	for (const auto& f : sigSamps) {
		Plotter *p = new Plotter();
		for (const auto& comb : combs) {
			TH1 *h = getHistByMass("passIP",f.first,masses,comb,"");
			p->addHist(h,comb,-1,1,4,20,1,true,false);
		}
		p->setXTitle("M_{X}");
		p->draw(false,f.second+" Isalurm");
	}

	for (const auto& b : bkgSamps) {
		Plotter *p = new Plotter();
		TH1 *hd = (TH1*)b.first->Get("bkg_passIP_maxEta2p5_hhmass");
		hd = (TH1*)hd->Clone();
		PlotTools::toOverflow(hd);
		hd = PlotTools::rebin(hd,nhhBins,hhBins);

		p->addHist(hd,"den");

		for (const auto& comb : combs) {
			TH1 *hn = (TH1*)b.first->Get("bkg_passIP_"+comb+"_hhmass");
			if (hn==0) cout<<"bad num: "<<"bkg_passIP_"+comb+"_hhmass"<<endl;
			hn = (TH1*)hn->Clone();
			PlotTools::toOverflow(hn);
			hn = PlotTools::rebin(hn,nhhBins,hhBins);
			p->addHist(hn,comb);
		}
		p->drawRatio(0,"stack",true,false,b.second+" Iso-Eta efficiencies 0.1");
	}

	// ETA ---------------------------------------- //


	return;
}
