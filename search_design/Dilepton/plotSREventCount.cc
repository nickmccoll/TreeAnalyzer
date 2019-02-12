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

void plotSREventCount() {
//	gROOT->SetBatch(true);

    vector<TString> tableRows = {"SR_btagLMT_","SR_btagL_","SR_btagM_","SR_btagT_","ee_SR_btagLMT_","mumu_SR_btagLMT_","emu_SR_btagLMT_"};
    vector<TString> bkgCats = {"incl_","Mbb100to150_","MhhGt2000_","Mbb100to150_MhhGt2000_"};
    vector<TString> vars = {"Mbb"};
    vector<TString> cutCats = {"","relax_dPhi_","relax_Bveto_","relax_Btag_"};
    vector<TString> btagCats = {"","BtagCat_1_","BtagCat_2_","BtagCat_3_"};
    map<TString,TString> panes = {{"SR_btagLMT_","Incl"},{"SR_btagL_","Loose Btag"},{"SR_btagM_","Medium Btag"},{"SR_btagT_","Tight Btag"},{"mumu_SR_btagLMT_","Same Flavor"},{"emu_SR_btagLMT_","Opp Flavor"}};
    map<TString,int> SRs = {{"mumu_SR_btagL_",1},{"emu_SR_btagL_",2},{"mumu_SR_btagM",3},{"emu_SR_btagM_",4},{"mumu_SR_btagT_",5},{"emu_SR_btagT_",6}};


    TFile *fdy = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/dyll_evtCount.root");
    TFile *ftt = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/ttbar_2l_evtCount.root");
    TFile *f1l = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/fake_1l_evtCount.root");
    TFile *fother = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/other_evtCount.root");
    TFile *fsig = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/radion_evtCount.root");
    TFile *fbkg = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/bkgTot_evtCount.root");

    for (const auto& var : vars) {
    	vector<TObject*> canvases;
    	gROOT->SetBatch(true);
    	Plotter *pp = new Plotter();
    	for (const auto& pane : panes) {
    		Plotter *p = new Plotter();

    		TString histID = pane.first+"incl_"+var;
//    		cout<<histID<<endl;
    		TH1F *hdy = (TH1F*)fdy->Get("bkg_"+histID);
    		TH1F *htt = (TH1F*)ftt->Get("bkg_"+histID);
    		TH1F *h1l = (TH1F*)f1l->Get("bkg_"+histID);
    		TH1F *hot = (TH1F*)fother->Get("bkg_"+histID);
    		TH1F *hTO = (TH1F*)fbkg->Get("bkg_"+histID);

    		if (pane.first.BeginsWith("mu")) { // same flavor leptons
    			TString eeID = "ee_SR_btagLMT_incl_"+var;
    			hdy->Add((TH1F*)fdy->Get("zjets_"+eeID));
    			htt->Add((TH1F*)ftt->Get("ttbar_"+eeID));
    			h1l->Add((TH1F*)f1l->Get("bkg_"+eeID));
    			hot->Add((TH1F*)fother->Get("bkg_"+eeID));
    			hTO->Add((TH1F*)fbkg->Get("bkg_"+eeID));
    		}

    		p->addHistLine(hTO,"total bkg");
    		p->addHistLine(hdy,"Drell-Yan");
    		p->addHistLine(htt,"tbar 2l");
    		p->addHistLine(h1l,"1l fakes");
    		p->addHistLine(hot,"others");

    		if (var=="Mbb") {
        		if (pane.second == "Incl") {
            		pp->addStackHist(htt,"tbar 2l");
            		pp->addStackHist(hdy,"Drell-Yan");
            		pp->addStackHist(h1l,"1l fakes");
            		pp->addStackHist(hot,"others");
        		}
    		}

    		p->normalize();
		    p->setBotMinMax(0,2);
            p->addText(pane.second,0.2,0.775,0.116);
            if (var=="Mbb") p->rebin(4);
            else            p->rebin(5);
    		canvases.push_back(p->drawSplitRatio(0,pane.second+var,false,false,pane.second+var));
    	}
    	gROOT->SetBatch(false);

    	pp->draw(false,"inclusive stats");

    	Plotter *pALL = new Plotter();
    	Plotter *pB = new Plotter();
    	for (const auto& cut : cutCats) {
        	TString histID = "SR_btagLMT_"+cut+"incl_"+var;
        	TH1F *hTO = (TH1F*)fbkg->Get("bkg_"+histID);
        	pALL->addHistLine(hTO,cut);
        }
        pALL->normalize();
    	pALL->setBotMinMax(0,2);
    	pALL->rebin(4);
        pALL->drawSplitRatio(0,var+"tot bkg",false,false,var+"tot bkg");
        for (const auto& cat : btagCats) {
        	TString histID = "SR_btagLMT_"+cat+"incl_"+var;
        	TH1F *hTO = (TH1F*)fbkg->Get("bkg_"+histID);
        	pB->addHistLine(hTO,cat);
        }
        pB->normalize();
    	pB->setBotMinMax(0,2);
    	pB->rebin(4);
        pB->drawSplitRatio(0,var+" btag Comp",false,false,var+" btag Comp");
    	Drawing::drawAll(canvases,var+"","");
    }

//    cout<<"new_section"<<endl;
    for (const auto& row : tableRows) {
//    	cout<<row<<endl;
    	vector<double> yields;
    	vector<double> dy_yields;
    	vector<double> tt_2l_yields;
    	vector<double> tt_1l_yields;
    	vector<double> wjets_yields;
    	vector<double> other_yields;

    	double n_800  =  ((TH1F*)fsig->Get("m800_"+row+"incl_ht"))->Integral();
    	double n_1000 = ((TH1F*)fsig->Get("m1000_"+row+"incl_ht"))->Integral();
    	double n_2000 = ((TH1F*)fsig->Get("m2000_"+row+"incl_ht"))->Integral();
    	double n_3000 = ((TH1F*)fsig->Get("m3000_"+row+"incl_ht"))->Integral();

//    	cout<<row<<endl<<endl;
    	yields.push_back(n_800);
    	yields.push_back(n_1000);
    	yields.push_back(n_2000);
    	yields.push_back(n_3000);

//    	cout<<"n_800 = "<<n_800<<endl;
//    	cout<<"n_1000 = "<<n_1000<<endl;
//    	cout<<"n_2000 = "<<n_2000<<endl;
//    	cout<<"n_3000 = "<<n_3000<<endl<<endl;

    	for (const auto& cat : bkgCats) {
//    		cout<<cat<<endl;
    	    TH1F *hdy = (TH1F*)fdy->Get("zjets_"+row+cat+"ht");
    	    TH1F *htt = (TH1F*)ftt->Get("ttbar_"+row+cat+"ht");
    	    TH1F *h1_w = (TH1F*)f1l->Get("wjets_"+row+cat+"ht");
    	    TH1F *h1_t = (TH1F*)f1l->Get("ttbar_"+row+cat+"ht");
            TH1F* h0_db = (TH1F*)fother->Get("diboson_"+row+cat+"ht");
            TH1F* h0_t = (TH1F*)fother->Get("singlet_"+row+cat+"ht");
            TH1F* h0_qc = (TH1F*)fother->Get("qcd_"+row+cat+"ht");
            TH1F* h0_ttx = (TH1F*)fother->Get("ttX_"+row+cat+"ht");
            TH1F* h0_hx = (TH1F*)fother->Get("hx_"+row+cat+"ht");

            double n_bkg = 0;
	        double n_other = 0;

            if (hdy)    n_bkg += hdy->Integral();
            if (htt)    n_bkg += htt->Integral();
            if (h1_w)   n_bkg += h1_w->Integral();
            if (h1_t)   n_bkg += h1_t->Integral();;
            if (h0_db)  {n_bkg += h0_db->Integral(); n_other += h0_db->Integral();}
            if (h0_t)   {n_bkg += h0_t->Integral();  n_other += h0_t->Integral();}
//            if (h0_qc) n_bkg += h0_qc->Integral();
            if (h0_ttx) {n_bkg += h0_ttx->Integral(); n_other += h0_ttx->Integral();}
            if (h0_hx)  {n_bkg += h0_hx->Integral();  n_other += h0_hx->Integral();}

//             cout<<cat<<endl;
//             cout<<"n_bkg = "<<n_bkg<<endl<<endl;
            yields.push_back(n_bkg);
	    
	        if (hdy) dy_yields.push_back(hdy->Integral());
	        else     dy_yields.push_back(0);

	        if (htt) tt_2l_yields.push_back(htt->Integral());
            else     tt_2l_yields.push_back(0);

	        if (h1_t) tt_1l_yields.push_back(h1_t->Integral());
            else      tt_1l_yields.push_back(0);

	        if (h1_w) wjets_yields.push_back(h1_w->Integral());
            else      wjets_yields.push_back(0);

	        other_yields.push_back(n_other);
//	        if (h0_qc) other_yields.push_back(h0_qc->Integral());
//	        else      other_yields.push_back(0);
    	}
    	// fill vector with different S/B
    	yields.push_back(yields[0]/yields[5]); // 800 GeV div Hbb 100-150 bkg
    	yields.push_back(yields[1]/yields[5]); // 1 TeV div Hbb 100-150 bkg
    	yields.push_back(yields[2]/yields[7]); // 2 TeV div Hbb 100-150 && HH > 2000 bkg
    	yields.push_back(yields[3]/yields[7]); // 3 TeV div Hbb 100-150 && HH > 2000 bkg

//	cout<<"Tot BKG"<<endl;
//    	for (auto& k : yields) cout<<k<<", ";
//    	cout<<endl;

//	cout<<"Drell-Yan"<<endl;
//	for (auto& k : dy_yields) cout<<k<<", ";
//        cout<<endl;
//
//	cout<<"ttbar 2l"<<endl;
//        for (auto& k : tt_2l_yields) cout<<k<<", ";
//        cout<<endl;
//
//	cout<<"ttbar 1l"<<endl;
//        for (auto& k : tt_1l_yields) cout<<k<<"	";
//        cout<<endl;
//
//	cout<<"W+jets"<<endl;
//        for (auto& k : wjets_yields) cout<<k<<", ";
//        cout<<endl;
//
//	cout<<"Other"<<endl;
//        for (auto& k : other_yields) cout<<k<<", ";
//        cout<<endl;
    }

}
