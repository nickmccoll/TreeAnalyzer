#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "Configuration/interface/FillerConstants.h"
#include "TFile.h"
#include "TH1.h"

using namespace std;

void plotNm1Eff(TFile *f, TString smp, TString denS, vector<TString> numS, TString var, TString canName, float rebin=0, int nR = 0, double * rebins = 0) {
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

void drawElectronIDvars() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/";
	TFile *frad = new TFile(prePath+"eID_rad.root");
	TFile *fbgv = new TFile(prePath+"eID_bg.root");

	int nLepBins = 17;
	double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};

	vector<TString> selVars = {"relax_hoe","relax_dEtaSeed","relax_dPhiIn","relax_sigmaIetaIeta","relax_1oEm1op","TightID","MVAID"};
	plotNm1Eff(fbgv,"m1200","passIPandISO",selVars,"pt","BG 1.2 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);
	plotNm1Eff(fbgv,"m2000","passIPandISO",selVars,"pt","BG 2 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);
	plotNm1Eff(frad,"m1400","passIPandISO",selVars,"pt","RAD 1.4 TeV N-1 Electron Tight ID eff",0,nLepBins,lepBins);


	return;
}
