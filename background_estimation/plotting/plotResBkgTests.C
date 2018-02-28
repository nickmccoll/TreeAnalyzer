#include "plotTestHelper.h"
#include "TH1.h"
#include "TH2.h"

  void testHHKern(std::string name, std::string filename) {
      bool withRatio = false;
    TFile *f = new TFile((filename + "_"+name +"_incl_template.root").c_str(),"read");
    std::vector<TH1*> hs;
    std::vector<std::string> hNs;
    TH1* h = 0;
    f->GetObject("histo_data",h);hs.push_back(h);hNs.push_back("MC");
    f->GetObject("histo",h);     hs.push_back(h);hNs.push_back("KDE");
//    f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back("KDE w/o expo. tail smoothing");

    int binL = hs[0]->FindFixBin(minHHMass);
    int binH = hs[0]->FindFixBin(maxHHMass);
     for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral(binL,binH)/hs[iH]->Integral(binL,binH));

    Plotter * p = new Plotter();
    Plotter * pf = new Plotter();
    for(unsigned int iH = 0; iH < hs.size(); ++iH){
        if(iH == 0){
            p->addHist(hs[iH],hNs[iH].c_str());
            pf->addHist(hs[iH],hNs[iH].c_str());
        }
        else  {
            p->addHistLine(hs[iH],hNs[iH].c_str());
            pf->addHistLine(hs[iH],hNs[iH].c_str());
        }
    }
    auto setupPlotter = [&](Plotter * p, std::string name, double rebin =-1){
        p->setMinMax(.0001,(rebin < 0 ? 1.0 : rebin) * hs[0]->Integral()/4);
        p->setUnderflow(false);
        p->setOverflow(false);
        p->setBotMinMax(0,2);
        p->setYTitle("N. of events");
        p->setXTitle(hhMCS.title.c_str());
        if(rebin > 0) p->rebin(rebin);
        if(withRatio){
            auto * c = p->drawSplitRatio(1,"stack",false,false,name);
            c->GetPad(1)->SetLogy();
            c->GetPad(1)->Update();
        } else {
            auto * c = p->draw(false,name);
            c->SetLogy();
            c->Update();
        }
    };

    setupPlotter(p,name + "_HHKDE_C",5);
    setupPlotter(pf,name + "_HHKDE_F");
  }

  void testHHPDFFits(std::string name, std::string filename) {
      std::vector<std::string> sels = {"emu_LMT_lb","mu_L_full","mu_M_full","mu_T_full","e_L_full","e_M_full","e_T_full"};

      TFile * fd = new TFile((filename+"_"+name+"_distributions.root").c_str(),"read");
      TFile * fo = new TFile((filename+"_"+name+"_template.root").c_str(),"read");
      TH1 * hof = 0;
      fo->GetObject("histo",hof);

      for(const auto& s : sels){
        TH1 * hd = 0;
        fd->GetObject(  (name+"_"+s+"_hhMass").c_str(),hd);
        if(hd==0) continue;
        std::vector<TH1*> hs;
        std::vector<TString> hNs;
        hs.push_back((TH1*)hof->Clone());
        hNs.push_back("Baseline KDE before MC fit");

        TFile * ff = new TFile((filename+"_"+name+"_"+s+"_fitTemplate.root").c_str(),"read");
        TH1* h = 0;
        ff->GetObject("histo",h);hs.push_back(h);hNs.push_back("Nominal SR PDF");

        for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(hd->Integral()/hs[iH]->Integral());

        Plotter * p = new Plotter();
        p->addHist(hd,"SR MC");
        for(unsigned int iH = 0; iH < hs.size(); ++iH){
          TH1 * h1D = hs[iH];
          for(int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
          p->addHist(h1D,hNs[iH],-1,1,4,20,1,false,true,false,"E");
        }
        p->setMinMax(.0001,hs[0]->Integral());
        p->setUnderflow(false);
        p->setOverflow(false);
        p->rebin(8);
        p->setXTitle("HH mass [GeV]");
        p->setYTitle("N. of events");
        p->setBotMinMax(0,2);
        // auto * c = p->drawSplitRatio(0,"stack",false,false,s);
        // c->GetPad(1)->SetLogy();
        // c->GetPad(1)->Update();
        auto * c = p->draw(false,(name + "_"+s+"_HHFIT").c_str());
        c->SetLogy();
        c->Update();

  }

  }

  void test2DFits(std::vector<CutStr> types, std::string filename) {
//      std::vector<std::string> sels = {"emu_LMT_ltmb","e_L_full","e_M_full","e_T_full","mu_L_full","mu_M_full","mu_T_full"};
      std::vector<std::string> sels = {"emu_LMT_lb","emu_M_lb","emu_LMT_full","emu_L_full","emu_M_full","emu_T_full"};
//      std::vector<std::string> sels = {"emu_LMT_ltmb"};
      // std::vector<double> bins = {30,50,100,150,210};
      // bool binInY = false;
//      std::vector<double> bins = {800,900,1000,1500,2000,3000,4000,5000};
//      std::vector<double> bins = {800,1000,2000,5000};
            std::vector<double> bins = {800,5000};
        bool binInY = true;




    for(const auto& s :sels){
        std::vector<TH2*> hs;
        std::vector<std::string> hNs;
        TH2 * dh = 0;
        for(const auto& t: types){
          TFile * fY = new TFile((filename+"_"+t+"_distributions.root").c_str(),"read");
          TH2 * h = 0;
          fY->GetObject((t+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
          if(h==0) continue;
          if(dh==0)dh= (TH2*)h->Clone();
          else dh->Add(h);

          TFile * fF = new TFile((filename+"_"+t+"_2D_template_debug.root").c_str(),"read");
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
          auto dh1 = proj(dh,"MC");
          p->addHist(dh1,"MC");
          for(unsigned int iH = 0; iH < hs.size(); ++iH){
              TH1 * h = proj(hs[iH],hNs[iH]);
              for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
              p->addStackHist(h,hNs[iH].c_str());
          }
          p->setUnderflow(false);
          p->setOverflow(false);
          p->rebin(2);
//           p->setMinMax(.0001,dh1->Integral());
          p->setXTitle( (binInY ? hbbMCS : hhMCS) .title.c_str());
          p->setYTitle("N. of events");
          auto * c = p->draw(false,(s + ": "+flt2Str(bins[iB]) +"-"+flt2Str(bins[iB+1])).c_str());
//           c->SetLogy();
//           c->Update();

          // p->setBotMinMax(0,2);
          // auto * c = p->drawSplitRatio(-1,"stack",false,false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
          // c->GetPad(1)->SetLogy();
          // c->GetPad(1)->Update();
    }
    }
  }



  void testBKG1DFits(std::string name, std::string filename, std::string varName, std::string fitName, const std::vector<std::string>& sels){
      std::vector<std::string> canNames;
      for(unsigned int iM = 0; iM+1 < resPTBins.size(); ++iM){
          canNames.push_back(std::string("can_") +flt2Str(resPTBins[iM]) + "to"+flt2Str(resPTBins[iM+1]));
      }
      test1DFits(name,filename,varName,fitName,sels,canNames);
  }




void plotResBkgTests(int step = 0){
    std::string filename = hhFilename;

    if(step == 0)testHHKern(bkgSels[BKG_MW],filename);
    if(step == 1)testHHPDFFits(bkgSels[BKG_MW],filename);
    if(step == 2)testBKG1DFits(bkgSels[BKG_MW],filename,"W","fit1stIt",{"emu_LMT_none"});
    if(step == 3)testBKG1DFits(bkgSels[BKG_MW],filename,"W","fit",{"emu_LMT_none"});
    if(step == 4)test2DFits({bkgSels[BKG_MW] },filename);

    if(step == 5)testHHKern(bkgSels[BKG_MT],filename);
    if(step == 6)testHHPDFFits(bkgSels[BKG_MT],filename);
    if(step == 7)testBKG1DFits(bkgSels[BKG_MT],filename,"T","fit1stIt",{"emu_L_none","emu_M_none","emu_T_none"});
    if(step == 8)testBKG1DFits(bkgSels[BKG_MT],filename,"T","fit",{"emu_L_none","emu_M_none","emu_T_none"});
    if(step == 9)test2DFits({bkgSels[BKG_MT] },filename);
    if(step == 10)test2DFits({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename);
}
