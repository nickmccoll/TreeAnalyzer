#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TFile.h"
#include <iostream>
#include "TString.h"
#include "/Users/brentstone/Dropbox/Physics/GitRepos/HistoPlotting/include/Plotter.h"
using namespace std;

// run using: rr -b -q 'makePUSF.C(year)'

void makePUSF(int year) {
    TString path = "/Users/brentstone/Dropbox/Physics/HHbbWW/pileup/";
    vector<TString> wps = {"nom","up","down"};
    
    TString outS = TString::Format("/Users/brentstone/Dropbox/Physics/HHbbWW/pileup/puSF_%d.root",year);
    TFile *fout = new TFile(outS,"RECREATE");
    TFile *ft   = new TFile(path+TString::Format("mcFiles/ttbar_%d.root",year));
    
    vector<TString> fpustrs = {"nom","up","down"};

    TH1F *htt = (TH1F*)ft->Get("ttbar_numTruePUInt");
    htt->Scale(1/htt->Integral());

    for (const auto& wp : wps) {
	TFile *f = new TFile(path+TString::Format("dataFiles/pileup_%d_",year)+wp+".root");
        TH1F *h = (TH1F*)f->Get("pileup");
        h = (TH1F*)h->Clone("puSF_"+wp);
        h->Scale(1/h->Integral());
        h->Divide(htt);
	fout->cd();
        h->Write();
    }
    fout->Close();

}
