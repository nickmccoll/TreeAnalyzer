//#include "TFile.h"
//#include "TTree.h"
//#include "string"
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
