#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "HistoPlotting/include/Plotter.h"
//#include "HistoPlotting/include/PlotTools.h"
using namespace std;

// run using: rr -b -q 'makePUSF.C(year)'

void makePUSF(int year, bool showPlots=false) {
    TString path = "/Users/brentstone/Dropbox/Physics/HHbbWW/pileup/";
    vector<TString> wps = {"nom","up","down"};
    
    TString outS = TString::Format("/Users/brentstone/Dropbox/Physics/HHbbWW/pileup/puSF_%d.root",year);
    TFile *fout = new TFile(outS,"RECREATE");
    TFile *ft   = new TFile(path+TString::Format("mcFiles/puDistMC_%d.root",year));

    Plotter *p = new Plotter();
    Plotter *pp = new Plotter();
    TH1F *htt = (TH1F*)ft->Get("ttbar_numTruePUInt");
    htt->Scale(1/htt->Integral());
    if(showPlots) p->addHist(htt,"ttbar MC");

    for (const auto& wp : wps) {
        TFile *f = new TFile(path+TString::Format("dataFiles/pileup_%d_",year)+wp+".root");
        TH1F *h = (TH1F*)f->Get("pileup");
        h = (TH1F*)h->Clone("puSF_"+wp);
        h->Scale(1/h->Integral());

	if(showPlots) {
	    pp->addHist(h,wp);
	    if(wp == "nom") p->addHist(h,"data");
	}

        h->Divide(htt);
        fout->cd();
        h->Write();
    }
    fout->Close();
    cout<<"Done!"<<endl;

    if(showPlots) {
	cout<<"slurm"<<endl;
/*
	TFile *fb   = new TFile(path+"puDist_bkg2017.root");
	TH1 *h = (TH1*)fb->Get("ttbar_numTruePUInt");
	vector<TString> smps = {"qcd","singlet","zjets","wjets","ttx","hx","diboson"};
	for(const auto& s : smps) h->Add((TH1*)fb->Get(s+"_numTruePUInt"),1);
	h->Scale(1./h->Integral());
	TFile *ff = new TFile(path+"dataFiles/pileup_2017_nom.root");
	TH1 *h1 = (TH1*)ff->Get("pileup"); h1 = (TH1*)h1->Clone("pu");
	TH1 *h2 = (TH1*)h1->Clone("pu2");
	h1->Scale(1./h1->Integral());
	h2->Scale(1./h2->Integral());
	h1->Divide(h);
	h2->Divide(htt);
//	ff->Close();

	TFile *ffu = new TFile(path+"dataFiles/pileup_2017_up.root");
	TFile *ffd = new TFile(path+"dataFiles/pileup_2017_down.root");
        TH1 *hu = (TH1*)ffu->Get("pileup"); hu = (TH1*)hu->Clone("hup");
        TH1 *hd = (TH1*)ffd->Get("pileup"); hd = (TH1*)hd->Clone("hdn");
	hu->Scale(1./hu->Integral()); hu->Divide(htt);
	hd->Scale(1./hd->Integral()); hd->Divide(htt);

	Plotter *pi = new Plotter();
	pi->addHist(h2,"nom data / ttbar");
	pi->addHist(h1,"nom data / all bkg");
	pi->addHist(hu,"up data / ttbar");
	pi->addHist(hd,"dn data / ttbar");
	pi->drawSplitRatio(0,"stack",0,0,"slurm");
*/
        p->drawSplitRatio(0,"stack",0,0,"data / MC");
	pp->drawSplitRatio(0,"stack",0,0,"data PU scenario");
    }
}
