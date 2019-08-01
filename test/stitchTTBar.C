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

		if ( i >= h->FindFixBin(low) && i <= h->FindFixBin(high) ) continue;
		hc->SetBinContent(i, 0.0);
		hc->SetBinError(i, 0.0);
	}
	return hc;
}

void stitchTTBar() {
	TString pre = "~/Dropbox/Physics/HHbbWW/plots/MttXsec/";
	TFile *fl = new TFile(pre+"ttbar_lep.root");
	TFile *f1 = new TFile(pre+"ttbar_1000toInf.root");
	TFile *f7 = new TFile(pre+"ttbar_700to1000.root");

	TString wtS = "defwt";
	double lumi = 41.53*1000;

	double xs0 = 380.094;
	double xs1 = 364.3508;
	double xs2 = 83.71;


	TH1 *hl = (TH1*)fl->Get("mtt_hand");
	TH1 *h1 = (TH1*)f1->Get("mtt_hand");
	TH1 *h7 = (TH1*)f1->Get("mtt_hand");


	TH1 *hl_0 = (TH1*)fl->Get("mtt_sgnwt_hand_0lep");
	TH1 *hl_1 = (TH1*)fl->Get("mtt_sgnwt_hand_1lep");
	TH1 *hl_2 = (TH1*)fl->Get("mtt_sgnwt_hand_2lep");

	TH1 *h1_0 = (TH1*)f1->Get("mtt_sgnwt_hand_0lep");
	TH1 *h1_1 = (TH1*)f1->Get("mtt_sgnwt_hand_1lep");
	TH1 *h1_2 = (TH1*)f1->Get("mtt_sgnwt_hand_2lep");

	TH1 *h7_0 = (TH1*)f7->Get("mtt_sgnwt_hand_0lep");
	TH1 *h7_1 = (TH1*)f7->Get("mtt_sgnwt_hand_1lep");
	TH1 *h7_2 = (TH1*)f7->Get("mtt_sgnwt_hand_2lep");

	TH1 *hl_0_710to960 = getMttSubset(hl_0,710,960);
	TH1 *h7_0_710to960 = getMttSubset(h7_0,710,960);

	vector<TString> nleps = {"0lep","1lep","2lep"};
	vector<double> xsecs = {xs0,xs1,xs2};

	for (int i = 0; i < nleps.size(); i++) {

		TH1 *hl = (TH1*)fl->Get("mtt_sgnwt_hand_"+nleps[i]);
		TH1 *h1 = (TH1*)f1->Get("mtt_sgnwt_hand_"+nleps[i]);
		TH1 *h7 = (TH1*)f7->Get("mtt_sgnwt_hand_"+nleps[i]);

		double defwt = lumi * xsecs[i] / hl->Integral();

		TH1 *hl_710to960 = getMttSubset(hl,710,960);
		TH1 *hl_1010toIn = getMttSubset(hl,1010,99999);

		TH1 *hl_0to710 = getMttSubset(hl,0,710);
		TH1 *hl_960to1010 = getMttSubset(hl,960,1010);

		TH1 *h7_710to960 = getMttSubset(h7,710,960);
		TH1 *h1_1010toIn = getMttSubset(h1,1010,99999);

		double w7 = hl_710to960->Integral() / ( hl_710to960->Integral() + h7_710to960->Integral() );
		double w1 = hl_1010toIn->Integral() / ( hl_1010toIn->Integral() + h1_1010toIn->Integral() );

		cout<<nleps[i]<<endl;
		cout<<"w7_"<<i<<" = "<<w7*defwt*1000/lumi<<endl;
		cout<<"w1_"<<i<<" = "<<w1*defwt*1000/lumi<<endl;

		hl_710to960->Add(h7_710to960,1); hl_710to960->Scale(defwt*w7);
		hl_1010toIn->Add(h1_1010toIn,1); hl_1010toIn->Scale(defwt*w1);

		TH1 *hlwt = (TH1*)fl->Get("mtt_hand_"+nleps[i]);

		TH1 *hrec = (TH1*)hlwt->Clone(); hrec->Reset();
		hl_0to710->Scale(defwt);
		hl_960to1010->Scale(defwt);
		hrec->Add(hl_0to710,1);
		hrec->Add(hl_710to960,1);
		hrec->Add(hl_960to1010,1);
		hrec->Add(hl_1010toIn,1);

		// mhh shit
		TH1 *mlwt = (TH1*)fl->Get("mhh_defwt_"+nleps[i]);
		TH1 *mrec = (TH1*)mlwt->Clone(); mrec->Reset();

		TH1 *ml_0to710    = (TH1*)fl->Get("mhh_sgnwt_"+nleps[i]+"_0to710");
		TH1 *ml_710to960  = (TH1*)fl->Get("mhh_sgnwt_"+nleps[i]+"_710to960");
		TH1 *ml_960to1010 = (TH1*)fl->Get("mhh_sgnwt_"+nleps[i]+"_960to1010");
		TH1 *ml_1010toInf = (TH1*)fl->Get("mhh_sgnwt_"+nleps[i]+"_1010toInf");

		TH1 *m1 = (TH1*)f1->Get("mhh_sgnwt_"+nleps[i]+"_1010toInf");
		TH1 *m7 = (TH1*)f7->Get("mhh_sgnwt_"+nleps[i]+"_710to960");

		ml_0to710->Scale(defwt);
		ml_710to960->Add(m7,1); ml_710to960->Scale(w7*defwt);
		ml_960to1010->Scale(defwt);
		ml_1010toInf->Add(m1,1); ml_1010toInf->Scale(w1*defwt);

		mrec->Add(ml_0to710,1);
		mrec->Add(ml_710to960,1);
		mrec->Add(ml_960to1010,1);
		mrec->Add(ml_1010toInf,1);

		Plotter *pt = new Plotter();
		Plotter *ph = new Plotter();

		pt->addHist(hlwt,"original");
		pt->addHist(hrec,"stitched");

		cout<<endl<<"Bin errors:"<<endl;
		cout<<"og: "<<hlwt->GetBinError(hlwt->FindFixBin(1100))<<endl;
		cout<<"st: "<<hrec->GetBinError(hrec->FindFixBin(1100))<<endl;
		cout<<endl;

		ph->addHist(mlwt,"original");
		ph->addHist(mrec,"stitched");

		pt->draw(false,"mtt "+nleps[i]);
//		pt->drawSplitRatio(0,"mtt "+nleps[i],false,false,"mtt "+nleps[i]);
		ph->drawSplitRatio(0,"mhh "+nleps[i],false,false,"mhh "+nleps[i]);

	}


}
