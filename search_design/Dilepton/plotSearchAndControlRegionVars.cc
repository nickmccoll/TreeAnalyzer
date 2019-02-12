#include "TROOT.h"
#include "TH1.h"
#include "TSystem.h"
#include "TH2.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "/Users/brentstone/Dropbox/Physics/GitRepos/HistoPlotting/include/Plotter.h"
#include "/Users/brentstone/Dropbox/Physics/GitRepos/HistoPlotting/include/PlotTools.h"

using namespace std;

void plotSearchAndControlRegionVars() {
    vector<TString> chans = {"","ee_","mumu_","emu_"};
    vector<TString> vars = {"Mww","Mbb","Mhh"};
    vector<TString> regions = {"SR_","TopCR_","QgCR_"};
    vector<TString> btags = {"","btagLMT_","btagL_","btagM_","btagT_"};

    TFile *fdy = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/dyll.root");
    TFile *ftt = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/ttbar_2l.root");
    TFile *fdata = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/data.root");
    TFile *f1l = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/fake_1l.root");
    TFile *fother = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/other_bkg.root");

    for (const TString& reg : regions) {
    	for (const TString& var : vars) {
    		vector<TObject*> canvases;
    		gROOT->SetBatch(true);

    		for (TString& ch : chans) {
    			for (const auto& btag : btags) {

    				if (reg=="QgCR_" && btag.Contains("b")) continue;
    				if (reg!="QgCR_" && btag == "") continue;

    				TString suff;
    				//if (var=="mll" || var=="Mbb" || var=="ht" || var=="met" || var=="dr_ll") suff = ch+reg+btag+"fullSel_"+var;
    				suff = ch+reg+btag+"fullSel_"+var;

    				cout<<suff<<endl;
    				Plotter *p = new Plotter();

    				TH1F* htt = (TH1F*)ftt->Get("ttbar_"+suff);

    				TH1F* hdy = (TH1F*)fdy->Get("zjets_"+suff);

    				TH1F *h1l = new TH1F("1l"+suff,"1l"+suff,htt->GetNbinsX(),htt->GetBinLowEdge(1),htt->GetBinLowEdge(htt->GetNbinsX()+1));
    				TH1F* h1l_w = (TH1F*)f1l->Get("wjets_"+suff);
    				TH1F* h1l_t = (TH1F*)f1l->Get("ttbar_"+suff);
    				if (h1l_w) h1l->Add(h1l_w,1);
    				if (h1l_t) h1l->Add(h1l_t,1);

    				TH1F *h0 = new TH1F("other"+suff,"other"+suff,htt->GetNbinsX(),htt->GetBinLowEdge(1),htt->GetBinLowEdge(htt->GetNbinsX()+1));
    				TH1F* h0_db = (TH1F*)fother->Get("diboson_"+suff);
    				TH1F* h0_t = (TH1F*)fother->Get("singlet_"+suff);
    				TH1F* h0_qc = (TH1F*)fother->Get("qcd_"+suff);
    				TH1F* h0_ttx = (TH1F*)fother->Get("ttX_"+suff);
    				TH1F* h0_hx = (TH1F*)fother->Get("hx_"+suff);
    				if (h0_db) h0->Add(h0_db,1);
    				if (h0_t)  h0->Add(h0_t,1);
    				if (h0_qc) h0->Add(h0_qc,1);
    				if (h0_ttx) h0->Add(h0_ttx,1);
    				if (h0_hx) h0->Add(h0_hx,1);

    				p->addStackHist(htt,"ttbar 2l");
    				p->addStackHist(hdy,"DY 2l");
    				p->addStackHist(h1l,"1l fake");
    				p->addStackHist(h0,"other");

    				if (reg !="SR_") {
    					TH1F* data = (TH1F*)fdata->Get("data_"+suff);
    					if (data) p->addHist(data,"data");
    					else cout<<endl<<"NO DATA IN "<<suff<<endl<<endl;
    				}

    				p->addText(ch+btag,0.2,0.8);
    				p->setBotMinMax(0,2);
    				if (var=="Mww") p->rebin(2);
    				if (!reg.Contains("Qg") && ((ch=="" && btag.Contains("LMT")) || (ch=="" && !btag.Contains("LMT")) || (ch.Contains("_") && btag.Contains("LMT"))) ) {
    					canvases.push_back(p->draw(false,suff));
    				} else if (reg.Contains("Qg")) canvases.push_back(p->drawSplitRatio(-1,suff,false,false,suff));
    			}
    		}
    		gROOT->SetBatch(false);
    		Drawing::drawAll(canvases,reg+var,"");
    	}
    }

/*    vector<TString> vars2 = {"mll","ht","dPhi_metLL","dr_ll","numB"};
    vector<TString> cuts = {"passDR","passDphi","passB","passDR_Dphi","passDR_B","passDphi_B","passDR_Dphi_B"};
    vector<TString> runs = {"B","C","D","E","F","G","H"};

    for (const auto& cut :cuts) {
    	for (const auto& var : vars2) {
    		for (const auto& run : runs) {
    			cout<<"Run "<<run<<endl;
    		vector<TObject*> canvases;

    		for (const auto& ch:chans) {
    			TString suff = ch+"QgCR_"+cut+"_"+var;
    			cout<<suff<<endl;
                Plotter *p = new Plotter();

                TH1F* htt = (TH1F*)ftt->Get("ttbar_"+suff);

                TH1F* hdy = (TH1F*)fdy->Get("zjets_"+suff);

                TH1F *h1l = new TH1F("1l"+suff,"1l"+suff,htt->GetNbinsX(),htt->GetBinLowEdge(1),htt->GetBinLowEdge(htt->GetNbinsX()+1));
                TH1F* h1l_w = (TH1F*)f1l->Get("wjets_"+suff);
                TH1F* h1l_t = (TH1F*)f1l->Get("ttbar_"+suff);
                if (h1l_w) h1l->Add(h1l_w,1);
                if (h1l_t) h1l->Add(h1l_t,1);

                TH1F *h0 = new TH1F("other"+suff,"other"+suff,htt->GetNbinsX(),htt->GetBinLowEdge(1),htt->GetBinLowEdge(htt->GetNbinsX()+1));
                TH1F* h0_db = (TH1F*)fother->Get("diboson_"+suff);
                TH1F* h0_t = (TH1F*)fother->Get("singlet_"+suff);
                TH1F* h0_qc = (TH1F*)fother->Get("qcd_"+suff);
                TH1F* h0_ttx = (TH1F*)fother->Get("ttX_"+suff);
                TH1F* h0_hx = (TH1F*)fother->Get("hx_"+suff);
                if (h0_db) h0->Add(h0_db,1);
                if (h0_t)  h0->Add(h0_t,1);
                if (h0_qc) h0->Add(h0_qc,1);
                if (h0_ttx) h0->Add(h0_ttx,1);
                if (h0_hx) h0->Add(h0_hx,1);

                p->addStackHist(htt,"ttbar 2l");
                p->addStackHist(hdy,"DY 2l");
                p->addStackHist(h1l,"1l fake");
                p->addStackHist(h0,"other");

                TFile *f = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/SR_CR/data_"+run+".root");
                TH1F* data = (TH1F*)f->Get("data_"+suff);
                p->addHist(data,"data");

                p->addText(suff+"_"+run,0.2,0.8);
                canvases.push_back(p->drawSplitRatio(-1,"bkg",false,false,suff+run));
    		}
    		Drawing::drawAll(canvases,cut+var+"_run"+run,"");
    		}
    	}
    } */
}
