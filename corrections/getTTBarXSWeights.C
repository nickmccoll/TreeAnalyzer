#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

TH1* getMttSubset(TH1* h, double low, double high) {

	TH1 *hc = (TH1*)h->Clone();
	for (unsigned int i=1; i < h->GetNbinsX()+1; i++) {

		if ( i >= h->FindFixBin(low) && i < h->FindFixBin(high) ) continue;
		hc->SetBinContent(i, 0.0);
		hc->SetBinError(i, 0.0);
	}
	return hc;
}

void getTTBarXSWeights(int year) {
	TString yS = TString::Format("%d",year);
	TFile *fl = new TFile("genmtt_lep_"+yS+".root");
	TFile *f1 = new TFile("genmtt_1000_"+yS+".root");
	TFile *f7 = new TFile("genmtt_700_"+yS+".root");

	double lumi = 1000;
	if (year == 2016) lumi *= 35.9;
	else if (year == 2017) lumi *= 41.5;
	else if (year == 2018) lumi *= 59.7;
	else throw std::invalid_argument("invalid year - 2016, 2017, or 2018 only");

	double xs0 = 380.094;
	double xs1 = 364.3508;
	double xs2 = 83.71;

//	TH1 *hl_0 = (TH1*)fl->Get("genMtt0");
//	TH1 *hl_1 = (TH1*)fl->Get("genMtt1");
//	TH1 *hl_2 = (TH1*)fl->Get("genMtt2");
//
//	TH1 *h1_0 = (TH1*)f1->Get("genMtt0");
//	TH1 *h1_1 = (TH1*)f1->Get("genMtt1");
//	TH1 *h1_2 = (TH1*)f1->Get("genMtt2");
//
//	TH1 *h7_0 = (TH1*)f7->Get("genMtt0");
//	TH1 *h7_1 = (TH1*)f7->Get("genMtt1");
//	TH1 *h7_2 = (TH1*)f7->Get("genMtt2");
//
//	TH1 *hl_0_700to1000 = getMttSubset(hl_0,700,1000);
//	TH1 *h7_0_700to1000 = getMttSubset(h7_0,700,1000);

	vector<TString> nleps = {"0","1","2"};
	vector<double> xsecs = {xs0,xs1,xs2};

	for (int i = 0; i < nleps.size(); i++) {

	    TH1 *hl = (TH1*)fl->Get("genMtt"+nleps[i]);
	    TH1 *h1 = (TH1*)f1->Get("genMtt"+nleps[i]);
	    TH1 *h7 = (TH1*)f7->Get("genMtt"+nleps[i]);

	    double defwt = lumi * xsecs[i] / hl->Integral();

	    TH1 *hl_700to1000 = getMttSubset(hl,700,1000);
	    TH1 *hl_1000toInf = getMttSubset(hl,1000,9999);

	    TH1 *hl_0to700 = getMttSubset(hl,0,700);
//	    TH1 *hl_960to1010 = getMttSubset(hl,960,1010);

	    TH1 *h7_700to1000 = getMttSubset(h7,700,1000);
	    TH1 *h1_1000toInf = getMttSubset(h1,1000,9999);

	    double w7 = hl_700to1000->Integral() / ( hl_700to1000->Integral() + h7_700to1000->Integral() );
	    double w1 = hl_1000toInf->Integral() / ( hl_1000toInf->Integral() + h1_1000toInf->Integral() );

	    cout<<year<<": nleps = "<<nleps[i]<<endl;
	    cout<<"w7_"<<i<<" = "<<w7*defwt*1000/lumi<<endl;
	    cout<<"w1_"<<i<<" = "<<w1*defwt*1000/lumi<<endl;
	    cout<<endl;

	    hl_700to1000->Add(h7_700to1000,1); hl_700to1000->Scale(defwt*w7);
	    hl_1000toInf->Add(h1_1000toInf,1); hl_1000toInf->Scale(defwt*w1);

	    TH1 *hreco = (TH1*)hl->Clone(); hreco->Reset();
	    hl_0to700->Scale(defwt);

	    hreco->Add(hl_0to700,1);
	    hreco->Add(hl_700to1000,1);
	    hreco->Add(hl_1000toInf,1);

	    Plotter *pt = new Plotter();
	    Plotter *ptr = new Plotter();

	    hl->Scale(defwt);
	    pt->addHist(hl,"original lep-binned");
	    pt->addHist(hreco,"new stitched");
	    pt->draw(false,"mtt "+nleps[i]);

	    ptr->addHist(hl,"original lep-binned");
	    ptr->addHist(hreco,"new stitched");
	    ptr->drawRatio(0,"rat "+nleps[i],false,false,"rat "+nleps[i]);
	}

}
