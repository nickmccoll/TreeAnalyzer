#include "TFile.h"
#include "TTree.h"
#include "../predTools/CutConstants.h"
using namespace CutConstants;
//void skimTree(std::string inFile,std::string outFile, std::string skimString){
//    TFile * fin = new TFile(inFile.c_str(),"read");
//    TTree * tin = 0;
//    fin->GetObject("treeMaker/Events",tin);
//    if(tin == 0) return;
//    TFile * fout = new TFile(outFile.c_str(),"recreate");
//    TDirectory * cdtof = fout->mkdir("treeMaker");
//    cdtof->cd();
//    TTree *tout = tin->CopyTree(skimString.c_str());
//    tout->Write();
//    delete fin;
//    delete fout;
//}

void skimTree(std::string inFile,std::string outPrefix,std::string skimString, bool cutByBackground){
    TFile * fin = new TFile(inFile.c_str(),"read");
    TTree * tin = 0;
    fin->GetObject("treeMaker/Events",tin);
    if(tin == 0) return;


    if(cutByBackground){
        for(const auto& b : bkgSels){
            TFile * fout = new TFile((outPrefix+"_" + b+".root").c_str(),"recreate");
            TDirectory * cdtof = fout->mkdir("treeMaker");
            cdtof->cd();
            std::string tss = skimString.size()?  b.cut+"&&("+skimString+")" : b.cut;
            TTree *tout = tin->CopyTree(tss.c_str());
            tout->Write();
            delete fout;
        }
    } else {
        TFile * fout = new TFile((outPrefix+".root").c_str(),"recreate");
        TDirectory * cdtof = fout->mkdir("treeMaker");
        cdtof->cd();
        TTree *tout = tin->CopyTree(skimString.c_str());
        tout->Write();
        delete fout;
    }

    delete fin;
}
