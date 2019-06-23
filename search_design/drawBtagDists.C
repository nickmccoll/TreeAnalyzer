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
    p->drawRatio(false,canName);
}

void plotDistribution(vector<pair<TString,TFile*>> fs, vector<pair<TString,TFile*>> fbs,
		TString hist, TString canName="slurm", float rebin=0, int nR=0, double *rebins = 0) {

	auto tfHist = [&](TString n, TH1* h) {
		if (h->GetNbinsX() != 10) return h;
		TH1F *hn = new TH1F(n,n,6,-0.5,5.5);
		for (unsigned int i=1; i < 10; i++) {
			if (i==6 || i==9 || i==10) hn->Fill(5,h->GetBinContent(i));
			else if (i==5 || i==8) hn->Fill(4,h->GetBinContent(i));
			else if (i==4 || i==7) hn->Fill(3,h->GetBinContent(i));
			else hn->Fill((i-1),h->GetBinContent(i));
		}
		return (TH1*)hn;
	};

	Plotter *p = new Plotter();
	for (const auto& f : fs) {
		TH1 *hs = (TH1*)f.second->Get(f.first+"_"+hist);
		if (hs==0) {cout<<"bad sig: "<<f.first+"_"+hist<<endl; return;}

		hs = (TH1*)hs->Clone();
		hs = tfHist(f.first+hist,hs);
		PlotTools::toOverflow(hs);
		PlotTools::toUnderflow(hs);

		if(rebin > 0) PlotTools::rebin(hs,rebin);
		else if(rebins) hs = PlotTools::rebin(hs,nR,rebins);
		p->addHistLine(hs,f.first);
	}

    for (const auto& f : fbs) {
    	TH1 *h = (TH1*)f.second->Get(f.first+"_"+hist);
    	if (h==0) {cout<<"bad bkg: "<<f.first+"_"+hist<<endl; return;}
    	h = (TH1*)h->Clone();
    	h = tfHist(f.first+hist,h);
    	PlotTools::toOverflow(h);
    	PlotTools::toUnderflow(h);

        if(rebin > 0) PlotTools::rebin(h,rebin);
        else if(rebins) h = PlotTools::rebin(h,nR,rebins);
        p->addStackHist(h,f.first);
    }

    p->normalize();
    p->drawSplitRatio(-1,canName,false,false,canName);

}

TH1F *getTotBkgDist(vector<pair<TString,TFile*>> bkgs, TString hist) {

	TH1F *h1 = (TH1F*)bkgs[0].second->Get(bkgs[0].first+"_"+hist);
	TH1F *hb = (TH1F*)h1->Clone();

	for (unsigned int i = 1; i < bkgs.size(); i++) {
		TH1F *h = (TH1F*)bkgs[i].second->Get(bkgs[i].first+"_"+hist);
		hb->Add(h,1);
	}

	return hb;
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

void printEfficiencies(vector<pair<TString,TFile*>> sigSamps,vector<pair<TString,TFile*>> bkgSamps) {
	vector<TString> bCats = {"Deep_","CSV_"};
	vector<TString> lCats = {"","mu_","el_"};
	vector<TString> jpts = {"jetPt20_","jetPt25_","jetPt30_"};
	vector<TString> flvs = {"Fb_","Fc_","Fl_"};
	vector<TString> wps = {"bLoose_","bMed_","bTight_"};

	auto printeffs = [&](pair<TString,TFile*> smp) {
		cout<<"Sample "<<smp.first<<" ///////////////////////////"<<endl<<endl;
		for(const auto& fl : flvs) for(const auto& l : lCats) for(const auto& pt : jpts) {
			cout<<"Den = "<<l+fl+pt<<" ("<<smp.first<<")"<<endl;
			double den = ((TH1F*)smp.second->Get(smp.first+"_"+l+fl+pt+"bIncl_pt_eta"))->Integral();
			for (const auto& b : bCats) for(const auto& wp : wps) {
				TH1F *h = (TH1F*)smp.second->Get(smp.first+"_"+l+fl+pt+wp+b+"pt_eta");
				double num = 0;
				if(h) num = ((TH1F*)smp.second->Get(smp.first+"_"+l+fl+pt+wp+b+"pt_eta"))->Integral();
				cout<<b+wp<<" eff = "<<(num/den)<<endl;
			}
			cout<<endl;
		}
	};

	for (const auto& ss : sigSamps) printeffs(ss);
	for (const auto& bs : bkgSamps) printeffs(bs);
}

void drawBtagDists() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/btag17_lnuqq/";
	TFile *fbgv = new TFile(prePath+"btag_blk.root");
	TFile *ftt = new TFile(prePath+"btag_ttbar.root");
	TFile *fwj = new TFile(prePath+"btag_WJets.root");
	TFile *fqcd = new TFile(prePath+"btag_QCD.root");
	TFile *fst = new TFile(prePath+"btag_ST.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
	double hhBins[] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3500,5000};
	int nhhBins = 18;
	vector<int> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};

	vector<pair<TString,TFile*>> sigSamps = { {"m1000",fbgv}, {"m3000",fbgv} };
	vector<pair<TString,TFile*>> bkgSamps = { {"ttbar",ftt}, {"qcd",fqcd}, {"wjets",fwj}};
//	printEfficiencies(sigSamps,bkgSamps);

	vector<TString> bcats = {"csv","deep"};
	vector<TString> lcats = {""/*,"el_","mu_"*/};
	for (const auto& b : bcats) for (const auto& l : lcats) {
		plotDistribution(sigSamps,bkgSamps,l+"sj_"+b+"_bCats",l+b);
	}
	for (const auto& sig : sigSamps) {
		Plotter *p = new Plotter();
		for (const auto& b : bcats) {
			TH1 *h = (TH1*)sig.second->Get(sig.first+"_sj_"+b+"_bCats");
			TH1 *hb = (TH1*)getTotBkgDist(bkgSamps,"sj_"+b+"_bCats");

			h->Scale(1/h->Integral());
			hb->Scale(1/hb->Integral());
			h->Divide(h,hb,1,1,"b");
			p->addHistLine(h,b);
		}
		p->draw(false,"S/B "+sig.first);
	}

//	TString hist1 = "deep_bCats";
//	TString hist2 = "csv_bCats";
//	Plotter *p = new Plotter();
//	TH1F *hx1 = (TH1F*)fbgv->Get("m2000_"+hist1);
//	TH1F *hy1  = getTotBkgDist(bkgSamps,hist1);
//	TGraph *g1 = PlotTools::getRocCurve(hx1,hy1,true,"signal","bkg");
//	TH1F *hx2 = (TH1F*)fbgv->Get("m2000_"+hist2);
//	TH1F *hy2 = getTotBkgDist(bkgSamps,hist2);
//	TGraph *g2 = PlotTools::getRocCurve(hx2,hy2,true,"signal","bkg");
//	p->addGraph(g2,"CSV");
//	p->addGraph(g1,"Deep");
//	p->draw(false,"fck u");

	vector<TString> WPs = {"bLoose_","bMed_","bTight_"};
	vector<TString> flvs = {"Fl_","Fb_"};
	vector<TString> bcats2 = {"Deep_","CSV_"};
	vector<TString> jpts = {"jetPt20_"};

//	for (const auto& b : bcats2) for(const auto& wp : WPs) for(const auto& l : lcats) for(const auto& j : jpts){
//		Plotter *p = new Plotter();
//		for (const auto& sig : sigSamps) {
//			TH1 *h = (TH1*)sig.second->Get(sig.first+"_"+j+l+wp+b+"nSepAK4");
//			p->addHistLine(h,sig.first);
//		}
//		for (const auto& bkg : bkgSamps) {
//			TH1 *h = (TH1*)bkg.second->Get(bkg.first+"_"+j+l+wp+b+"nSepAK4");
//			p->addStackHist(h,bkg.first);
//		}
//		p->normalize();
//		p->draw(false,j+l+wp+b);
//	}

	vector<pair<TString,TFile*>> samps = { {"ttbar",ftt},  {"wjets",fwj} };
	for (const auto& b : bcats2) for (const auto& flv : flvs) for (const auto& smp : samps){
		Plotter *p = new Plotter();
		cout<<smp.first+"_sj_"+flv+"bIncl_"+b+"mbb"<<endl;
		TH1 *den = (TH1*)smp.second->Get(smp.first+"_SR_sj_"+flv+"bIncl_mbb");
		p->addHist(den,"incl");
		double count_den = den->Integral();
		double count_num = 0;
		for (const auto& wp : WPs) {
			TH1 *num = (TH1*)smp.second->Get(smp.first+"_SR_sj_"+flv+wp+b+"mbb");
			p->addHist(num,wp);
			if (wp=="bLoose_") count_num += num->Integral();
		}
		printf("eff = %f\n",(count_num/count_den));
		p->drawRatio(0,b+flv+smp.first,false,false,b+flv+smp.first);
		cout<<endl;
	}



	return;
}
