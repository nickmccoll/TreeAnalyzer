{
  vector<unsigned int> fullSigMasses = {
    600,
    800,
    1000,
    1200,
    1400,
    1600,
    1800,
    2000,
    2500,
    3000,
    3500,
    4000,
    4500

  };
    vector<TString> bkgs = {
      "ttbar",
      "wjets",
      "rare"
    };
    vector<TString> bkgNamess = {
      "t#bar{t}",
      "w+jets",
      "other"
    };
  vector<unsigned int> sigMasses = {
    600,
    // 800,
    1000,
    // 1200,
    // 1400,
    // 1600,
    // 1800,
    2000,
    // 2500,
    // 3000,
    // 3500,
    // 4000,
    4500

  };

  TFile * f = new TFile("getWjjJetSel_plots.root","read");
  



  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres, bool doNorm){
    for(auto v : vars)     for(auto p : pres){
      Plotter * plots = new Plotter;      
      
      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;        
        f->GetObject(TString::Format("%s_%s%s",bkgs[iS].Data(),p.Data(),v.Data()),h);        
        if(h == 0) continue;
        plots->addStackHist(h,bkgNamess[iS]);        
      }      
      
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s%s",sigMasses[iM],p.Data(),v.Data()),h);
        if(h ==0)continue;
        if(doNorm){
          float stackNorm = plots->getStackIntegral();
          if(stackNorm) h->Scale(stackNorm/h->Integral(0,-1)); else PlotTools::normalize(h);
        } else if(plots->getTotStack()){
           h->Add(plots->getTotStack());
         }

        plots->addHistLine(h,TString::Format("#it{m}(X) %u GeV",sigMasses[iM]));
    }
    plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    // plots->yAxis()->SetTitleOffset(1.5);
  }
};

auto makeRocs  = [&](std::vector<TString> vars,TString prefix, TString name,  bool doForward = true ){
  Plotter * p = new Plotter();
  Plotter * pSoB = new Plotter();
  for(const auto& v : vars) {
    TH1 *hb = 0;
    Plotter * pI = new Plotter();
    
    for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
      TH1 * h = 0;        
      f->GetObject(TString::Format("%s_%s%s",bkgs[iS].Data(),prefix.Data(),v.Data()),h);        
      if(h == 0) continue;
      if(hb == 0) hb = (TH1*)h->Clone();
      else hb->Add(h);
      TH1 *hI = PlotTools::getIntegral(h,doForward,false);                              
      pI->addStackHist(hI,bkgNamess[iS].Data());
    }
    
    for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
      TH1 * h = 0;
      f->GetObject(TString::Format("m%u_%s%s",sigMasses[iM],prefix.Data(),v.Data()),h);
      if(h ==0)continue;
      
      TGraph * roc = PlotTools::getRocCurve(h,hb,doForward,"signal","bkg");
      p->addGraphLine(roc,TString::Format("#it{m}(X) %u GeV, %s",sigMasses[iM],v.Data()));        
      TH1 *hI = PlotTools::getIntegral(h,doForward,true);                              
      float stackNorm = hb->Integral(0,-1);
      hI->Scale(stackNorm);
      pI->addHistLine(hI,TString::Format("#it{m}(X) %u GeV",sigMasses[iM]));      
  }      
          
  pI->draw(false,TString::Format("%s_int_%s",name.Data(), v.Data()));    
  }
  p->draw(false,name.Data());
};

  auto effPlots = [&](TString name, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%u",name.Data(),iN),";#it{m}(X) [TeV]",40,0.550,4.550));
      }

      for(auto m : fullSigMasses){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_selection",m),h);
        if(h==0){continue;}
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(float(m)/1000.),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHist(hists[iN],cutNames[cuts[iN]],-1,1,4,20,1,true,false);
      }
      p->setMinMax(0,1.0);
      p->drawRatio(false,"stack",false,false,name);

            // p->draw(false);
};

std::vector<TString> vars = {"wjj_fj_mass","wjj_fj_sd_mass","wjj_fj_rawsd_mass","wjj_fj_tau2otau1","wjj_fj_csv","wjj_fj_minsdcsv","wjj_fj_maxsdcsv","wjj_fj_oM_csv","wjj_fj_oM_minsdcsv","wjj_fj_oM_maxsdcsv","wjj_fj_oM_tau2otau1","wjj_fj_oM_mass","wjj_fj_oM_sd_mass","wjj_fj_oM_rawsd_mass"};
std::vector<TString> pres = {""};
distPlots("plots",vars,pres,true);

std::vector<TString> rocvars = {"hbb_fj_mass","hbb_fj_sd_mass","hbb_fj_rawsd_mass"};
// std::vector<TString> rocvars = {"hbb_fj_oM_bbcsv","hbb_fj_oM_minsdcsv","hbb_fj_oM_maxMed_minsdcsv"};
// makeRocs(rocvars,"","roc",true);

std::vector<unsigned int> cuts = {0,1,2,3,4};
// std::vector<unsigned int> cuts = {0,1,5,6,7,8};

vector<TString> cutNames = {
  "inclusive",
  "fatjet candidate",
  "+ softdrop mass  > 10 GeV",
  "+ #tau_{3/2} < 0.55",
  "+ no medium CSV subjet",
  "wjj fj cand truth matched",
  "+ softdrop mass  > 10 GeV",
  "+ #tau_{3/2} < 0.55",
  "+ no medium CSV subjet"
};
effPlots("cutflow",cuts,cutNames);


}