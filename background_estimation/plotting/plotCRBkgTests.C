#include "../predTools/CutConstants.h"
#include <vector>
#include "TFile.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;
std::vector<TObject*> writeables;

class Dummy {
public:
    Dummy(const std::string& outName = "") : outName(outName) {};
    ~Dummy() {
        if(outName.size()){
            TFile * f = new TFile(outName.c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
            }
            f->Close();
        }
    }
    std::string outName;
};





std::vector<TObject*> comparePreAndPostFits(std::vector<CutStr> types, std::string preFitFilenames, std::string postFitFilenamePrefix, const std::vector<std::string>& sels,const std::vector<double>& bins, bool binInY = true, bool addData = false, bool doIntegral =false ) {
    std::vector<TObject*> writeables;
    std::vector<TFile*> distFiles;
    std::vector<TFile*> tempFiles;
    for(const auto& t: types){
        TFile * fY = new TFile((preFitFilenames+"_"+t+"_distributions.root").c_str(),"read");
        distFiles.push_back(fY);
        TFile * fF = new TFile((preFitFilenames+"_"+t+"_2D_template_debug.root").c_str(),"read");
        tempFiles.push_back(fF);
    }

    TFile * fD = 0;
    if(addData) {
        fD  = new TFile((preFitFilenames+"_data_distributions.root").c_str(),"read");
    }

  for(const auto& s :sels){
      //start with postfit
      TFile * fpf = new TFile((postFitFilenamePrefix+"_"+s+"_13TeV_comp.root").c_str(),"read");

      std::vector<TH2*> pfHs;
      std::vector<std::string> pfHNs;

      for(unsigned int iT = 0; iT < types.size(); ++iT){
          const auto& t = types[iT];
          TH2 * h = 0;
          fpf->GetObject((t).c_str(),h);
          if(h==0) continue;
          pfHs.push_back(h);
          pfHNs.push_back(t.title);
      }
      //Now data
      TH2 * dataHist = 0;
      if(addData){
          fD->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dataHist);
      }

      //and prefit
      std::vector<TH2*> hs;
      std::vector<std::string> hNs;
      TH2 * totalPrefit=0;
      TH2 * totalMC=0;
      for(unsigned int iT = 0; iT < types.size(); ++iT){
          const auto& t = types[iT];
          TFile * fY =distFiles[iT];
          TH2 * h = 0;
          fY->GetObject((t+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
          if(h==0) continue;
          if(totalMC ==0) totalMC= (TH2*)h->Clone();
          else totalMC->Add(h);

          TFile * fF =tempFiles[iT];
          TH2 * hF = 0;
          fF->GetObject((t+"_"+s).c_str(),hF);
          if(hF==0) continue;
          hF->Scale(h->Integral()/hF->Integral());
          hs.push_back(hF);
          hNs.push_back(t.title);
          if(totalPrefit ==0) totalPrefit= (TH2*)hF->Clone();
          else totalPrefit->Add(hF);
      }

      const TAxis * ax = hs[0]->GetXaxis();
      if(binInY) ax = hs[0]->GetYaxis();

      for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
        int binL = ax->FindFixBin(bins[iB]);
        int binH = ax->FindFixBin(bins[iB+1]) -1;
        if(binH<=binL) continue;
        auto proj =[&](TH2* h, const std::string& title) ->TH1*{
            auto h1 = binInY ? h->ProjectionX(  (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH) :  h->ProjectionY( (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH);
            if(doIntegral) return PlotTools::getIntegral(h1,true,false);
            return h1;
        };
        Plotter * p = new Plotter();

        if(dataHist){
            auto dataHist1 = proj(dataHist,"data");
            p->addHist(dataHist1,"data");
            double max = 0;
            for(int iB = 1; iB <= dataHist1->GetNbinsX();++iB) if(dataHist1->GetBinContent(iB) > max) max =dataHist1->GetBinContent(iB);
            p->setMinMax(.0001,3*max);

        }
        if(totalPrefit){
            auto prefit1 = proj(totalPrefit,"prefit");
            p->addHistLine(prefit1,"prefit");
        }

        for(unsigned int iH = 0; iH < pfHs.size(); ++iH){
            TH1 * h = proj(pfHs[iH],pfHNs[iH]);
            for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
            p->addStackHist(h,pfHNs[iH].c_str());
        }
        p->setUnderflow(false);
        p->setOverflow(false);
//        p->rebin(2);

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
         c->GetPad(1)->Update();
  }
  }
  return writeables;
}



void plotCRBkgTests(int step = 0, bool topCR = true,  std::string outName = ""){
    std:: string inName =  "bkgInputs" ;
    std::vector<std::string> srList = {"emu_LMT_I_full","e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"};
    if(topCR == true){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
    }
    else{
        inName =  "bkgInputsQGCR";
        hhFilename +="_QGCR";
        srList = {"emu_L_I_full","e_L_LP_full","mu_L_LP_full","e_L_HP_full","mu_L_HP_full"};
    }

    if(outName.size()){
        outName += topCR ? std::string("/") +  "TopCR" : std::string("/") +  "QGCR";
    }

    std::string filename = inName +"/"+hhFilename;

    if(step== 0){
        if(outName.size())         outName += "_dataComp.root";
        writeables =  test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,srList,{700,4000},true,true);
        writeables =  test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,srList,{30,210,100,150},false,true);
    }

    if(step== 1){
        if(outName.size())         outName += "_postPreComp.root";
        writeables =  comparePreAndPostFits({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,"limits/baseline_topSyst_reallybigptyunc_QGCR/PlotsPostFit_radHH/postFitMJJMVV_radHH_std",srList,{700,4000},true,true);
        writeables =  comparePreAndPostFits({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,"limits/baseline_topSyst_reallybigptyunc_QGCR/PlotsPostFit_radHH/postFitMJJMVV_radHH_std",srList,{30,210,100,150},false,true);
    }

    Dummy d(outName);


}







