#include "../CutConstants.h"
#include <vector>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"
using namespace CutConstants;
using namespace ASTypes;


void test1DFits(std::string name, std::string filename, std::string varName, std::string fitName, const std::vector<std::string>& sels, const std::vector<std::string>& canNames) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....

    for(const auto& s : sels){
        TFile * fo = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".root").c_str(),"read");
        if(fo == 0) continue;

        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");

        auto addCan = [&](const std::string& name,std::vector<TObject*>& list){
            TCanvas * can= 0;
            fo->GetObject(name.c_str(),can);
            if(can != 0) list.push_back(can);
        };
        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
            TGraphErrors * can= 0;
            fo->GetObject(name.c_str(),can);

            if(ff){
                ff->GetObject((name).c_str(),can);
            }
            //            if(!can) fo->GetObject(name.c_str(),can);
            if(can == 0) return;
            can->GetYaxis()->SetTitle(name.c_str());
            list.push_back(can);
        };

        std::vector<TObject*> fitPads;
        for(const auto& sB : canNames){ addCan(sB, fitPads);}
        std::vector<TObject*> paramPads;
        auto vN=[&](const std::string& v )->std::string{ return v +varName;};
        addGraph(vN("mean"  ), paramPads);
        addGraph(vN("sigma" ), paramPads);
        addGraph(vN("alpha" ), paramPads);
        addGraph(vN("alpha2"), paramPads);
        addGraph(vN("n"     ), paramPads);
        addGraph(vN("n2"    ), paramPads);
        addGraph(vN("slope" ), paramPads);
        addGraph(vN("fE"    ), paramPads);
        addGraph("chi2", paramPads);
        Drawing::drawAll(fitPads, (s + ": fits").c_str());
        Drawing::drawAll(paramPads, (s + ": params").c_str());

    }

}

void compFitParams(const std::string& name, const std::string& filename, const std::string& varName, const std::string& fitName,
        const std::string& canPre, const std::vector<std::string>& compSels) {
    gROOT->SetBatch(true);

    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<std::pair<std::string,TFile *>> files;
    for(const auto& s : compSels){
        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");
        files.emplace_back(s,ff);
    }


    auto addGraph = [&](const std::string& pname,std::vector<TObject*>& list){
        TCanvas * c = 0;
        int nLines = 0;
        for(auto& f : files){
            if(c == 0){
                TGraphErrors * graph= 0;
                f.second->GetObject((pname).c_str(),graph);
                if(!graph) continue;
                c = new TCanvas((canPre+"_"+f.first+"_"+pname).c_str());
                graph->GetYaxis()->SetTitle(pname.c_str());
                graph->Draw("AP");
                nLines++;
            }
            if(true){
                TF1 * fit= 0;
                f.second->GetObject((pname +"_func").c_str(),fit);
                if(!fit) continue;
                fit->SetLineColor(nLines);
                nLines++;
                fit->SetLineStyle(7);
                fit->Draw("SAME");
            }
        }
        if(c == 0) return;
        list.push_back(c);
    };

    std::vector<TObject*> paramPads;
    auto vN=[&](const std::string& v )->std::string{ return v +varName;};
    addGraph(vN("mean"  ), paramPads);
    addGraph(vN("sigma" ), paramPads);
    addGraph(vN("alpha" ), paramPads);
    addGraph(vN("alpha2"), paramPads);
    addGraph(vN("n"     ), paramPads);
    addGraph(vN("n2"    ), paramPads);
    addGraph(vN("slope" ), paramPads);
    addGraph(vN("fE"    ), paramPads);
    addGraph("chi2", paramPads);

    gROOT->SetBatch(false);
    Drawing::drawAll(paramPads, (canPre+ " : params").c_str());

}
