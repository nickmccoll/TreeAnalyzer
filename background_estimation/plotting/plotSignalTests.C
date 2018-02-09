#include "../CutConstants.h"
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
using namespace CutConstants;
using namespace ASTypes;

void testMJJFits(std::string name, std::string filename) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....

    std::vector<std::string> sels = {"emu_L_ltmb","emu_M_ltmb","emu_T_ltmb"};

    for(const auto& s : sels){
        TFile * fo = new TFile((filename+"_"+name+"_"+s+"_MJJ_fit.root").c_str(),"read");
        TH1 * hof = 0;
        fo->GetObject("histo",hof);

        auto addCan = [&](const std::string& name,std::vector<TObject*>& list){
            TCanvas * can= 0;
            fo->GetObject(name.c_str(),can);
            if(can != 0) list.push_back(can);
        };
        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
            TGraphErrors * can= 0;
            fo->GetObject(name.c_str(),can);
            if(can == 0) return;
            can->GetYaxis()->SetTitle(name.c_str());
            list.push_back(can);
        };

        std::vector<TObject*> fitPads;
        for(const auto& sB : signalMassBins){ addCan(std::string("can_m") +int2Str(sB), fitPads);}
        std::vector<TObject*> paramPads;
        addGraph("meanSMJJ", paramPads);
        addGraph("sigmaSMJJ", paramPads);
        addGraph("alphaSMJJ", paramPads);
        addGraph("alpha2SMJJ", paramPads);
        addGraph("nSMJJ", paramPads);
        addGraph("n2SMJJ", paramPads);
        addGraph("slopeSMJJ", paramPads);
        addGraph("fESMJJ", paramPads);
        addGraph("chi2", paramPads);
        Drawing::drawAll(fitPads, (s + ": fits").c_str());
        Drawing::drawAll(paramPads, (s + ": params").c_str());

    }

}



void plotSignalTests(){
    std::string filename = hhFilename;

//    test2DCondTemplate(bkgSels[BKG_LOSTTW],filename);
//    testMJJKern(bkgSels[BKG_LOSTTW],filename);
//    test2DTemplate(bkgSels[BKG_LOSTTW],filename);
//    test2DFits(bkgSels[BKG_LOSTTW],filename);

//        test2DCondTemplate(bkgSels[BKG_QG],filename);
//        testMJJKern(bkgSels[BKG_QG],filename);
//        test2DTemplate(bkgSels[BKG_QG],filename);
    testMJJFits(radionSig,filename);
}
