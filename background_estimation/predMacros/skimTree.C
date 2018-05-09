#include "TFile.h"
#include "TTree.h"
#include "../predTools/CutConstants.h"

void skimTree(std::string inFile,std::string outFile, std::string skimString){
    TFile * fin = new TFile(inFile.c_str(),"read");
    TTree * tin = 0;
    fin->GetObject("treeMaker/Events",tin);
    if(tin == 0) return;
    TFile * fout = new TFile(outFile.c_str(),"recreate");
    TDirectory * cdtof = fout->mkdir("treeMaker");
    cdtof->cd();
    TTree *tout = tin->CopyTree(skimString.c_str());
    tout->Write();
    delete fin;
    delete fout;
}

void skimTree(std::string inFile,std::string outPrefix){
    TFile * fin = new TFile(inFile.c_str(),"read");
    TTree * tin = 0;
    fin->GetObject("treeMaker/Events",tin);
    if(tin == 0) return;

    for(const auto& b : CutConstants::bkgSels){
        TFile * fout = new TFile((outPrefix+"_" + b+".root").c_str(),"recreate");
        TDirectory * cdtof = fout->mkdir("treeMaker");
        cdtof->cd();
        TTree *tout = tin->CopyTree(b.cut.c_str());
        tout->Write();
        delete fout;
    }
}
