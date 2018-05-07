#include "../predTools/CutConstants.h"
#include <vector>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TROOT.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"
using namespace CutConstants;
using namespace ASTypes;


std::vector<TObject*> test1DFits(std::string name, std::string filename, std::string varName, std::string fitName, const std::vector<std::string>& sels, const std::vector<std::string>& canNames) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<TObject*> writeables;
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
            if(ff){
                ff->GetObject((name).c_str(),can);
            }
            if(!can) fo->GetObject(name.c_str(),can);
            if(!can) return;
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
        auto *c   =Drawing::drawAll(fitPads, (s + "_fits").c_str());
        auto * c1 =Drawing::drawAll(paramPads, (s + "_params").c_str());
        writeables.push_back(c);
        writeables.push_back(c1);
    }
    return writeables;
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


std::vector<TObject*> test2DModel(std::vector<CutStr> types, std::string filename, const std::vector<std::string>& sels,const std::vector<double>& bins, bool binInY = true, bool addData = false ) {
    std::vector<TObject*> writeables;
    std::vector<TFile*> distFiles;
    std::vector<TFile*> tempFiles;
    for(const auto& t: types){
        TFile * fY = new TFile((filename+"_"+t+"_distributions.root").c_str(),"read");
        distFiles.push_back(fY);
        TFile * fF = new TFile((filename+"_"+t+"_2D_template_debug.root").c_str(),"read");
        tempFiles.push_back(fF);
    }

    TFile * fD = 0;
    if(addData) {
        fD  = new TFile((filename+"_data_distributions.root").c_str(),"read");
    }

  for(const auto& s :sels){
      std::vector<TH2*> hs;
      std::vector<std::string> hNs;
      TH2 * dataHist = 0;
      if(addData){
          fD->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dataHist);
      }

      TH2 * dh = 0;
      for(unsigned int iT = 0; iT < types.size(); ++iT){
          const auto& t = types[iT];
          TFile * fY =distFiles[iT];
          TH2 * h = 0;
          fY->GetObject((t+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
          if(h==0) continue;
          if(dh==0)dh= (TH2*)h->Clone();
          else dh->Add(h);

          TFile * fF =tempFiles[iT];
          TH2 * hF = 0;
          fF->GetObject((t+"_"+s).c_str(),hF);
          if(hF==0) continue;
          hF->Scale(h->Integral()/hF->Integral());
          hs.push_back(hF);
          hNs.push_back(t.title);
      }

      const TAxis * ax = hs[0]->GetXaxis();
      if(binInY) ax = hs[0]->GetYaxis();

      for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
        int binL = ax->FindFixBin(bins[iB]);
        int binH = ax->FindFixBin(bins[iB+1]) -1;
        auto proj =[&](TH2* h, const std::string& title) ->TH1*{
            return binInY ? h->ProjectionX(  (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH) :  h->ProjectionY( (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH);
        };
        Plotter * p = new Plotter();
        if(dataHist){
            auto dataHist1 = proj(dataHist,"data");
            p->addHist(dataHist1,"data");
            double max = 0;
            for(int iB = 1; iB <= dataHist1->GetNbinsX();++iB) if(dataHist1->GetBinContent(iB) > max) max =dataHist1->GetBinContent(iB);
            p->setMinMax(.0001,3*max);

        } else {
            auto dh1 = proj(dh,"MC");
            p->addHist(dh1,"MC");
        }

        for(unsigned int iH = 0; iH < hs.size(); ++iH){
            TH1 * h = proj(hs[iH],hNs[iH]);
            for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
            p->addStackHist(h,hNs[iH].c_str());
        }
        p->setUnderflow(false);
        p->setOverflow(false);
        p->rebin(2);

        p->setXTitle( (binInY ? hbbMCS : hhMCS) .title.c_str());
        p->setYTitle("N. of events");
        p->setCMSLumi();
//        auto * c = p->draw(false,(s + ": "+flt2Str(bins[iB]) +"-"+flt2Str(bins[iB+1])).c_str());
//           c->SetLogy();
//           c->Update();

         p->setBotMinMax(0,2);
         p->setYTitleBot("N/N(template)");
         auto * c = p->drawSplitRatio(-1,"stack",false,false,(s + "_"+flt2Str(bins[iB]) +"_"+flt2Str(bins[iB+1])).c_str());
         writeables.push_back(c);

        // c->GetPad(1)->SetLogy();
        // c->GetPad(1)->Update();
  }
  }
  return writeables;
}

