rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getCutflow.C+("../testTrees/out_radion_hh_bbinc_m3000_0.root",2,"test_m3000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getCutflow.C+("../testTrees/out_radion_hh_bbinc_m800_0.root",2,"test_m800.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getCutflow.C+("../testTrees/out_radion_hh_bbinc_m1000_0.root",2,"test_m1000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getCutflow.C+("../testTrees/out_t-tbar1l-madgraph_18.root",1,"test_mttbar.root")' &
hadd -f test_all.root test_m*.root

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
    vector<TString> bkgs 
    = {
      "ttbar",
      "wjets",
      // "qcd",
      "other"
    };
    vector<TString> bkgNamess = {
      "t#bar{t}",
      "w+jets",
      // "QCD",
      "other"
    };


    TFile * f = new TFile("getCutflow.root","read");
    // TFile * f = new TFile("test_all.root","read");
  



  auto distPlots = [&](TString name, const std::vector<TString>& vars,const std::vector<unsigned int>& sigMasses, const std::vector<TString>& pres, bool doNorm, float rebin = -1, bool addLumi=false){
    for(auto v : vars)     for(auto p : pres){
      Plotter * plots = new Plotter;      
      
      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s_%s",bkgs[iS].Data(),p.Data(),v.Data()),h);
        if(h == 0) continue;
        plots->addStackHist(h,bkgNamess[iS]);
      }
      
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s_%s",sigMasses[iM],p.Data(),v.Data()),h);
        if(h ==0)continue;
        if(doNorm){

          float stackNorm =  plots->getTotStack() ? plots->getStackIntegral() : 1.0;
          if(stackNorm) h->Scale(stackNorm/h->Integral(0,-1)); else PlotTools::normalize(h);
        } else if(plots->getTotStack()){
           h->Scale(20./1000.);
           h->Add(plots->getTotStack());
         }

        plots->addHistLine(h,TString::Format("#it{m}(X) %.1f TeV",float(sigMasses[iM])/1000.));
    }
    if(rebin > 0) plots->rebin(rebin);
    if(addLumi){
    plots->setCMSLumi(33, "36 fb^{-1} (13 TeV)", "Simulation preliminary" );
    plots->setYTitle("Events / 10 GeV");
    }
    plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    // plots->yAxis()->SetTitleOffset(1.5);
  }
  };
  
      
  std::vector<TString> vars = {"hh_mass","hbbTM_hh_mass"};
  std::vector<TString> pres = {"emu_L","emu_T","emu_LP","emu_VL"};

  vector<unsigned int> sigMasses = {800,1000,1600,2000,3000};
  distPlots("plots",vars,sigMasses,pres,false,20,false);
    
  // std::vector<TString> vars = {"hh900to1100_hbb_mass","hh1400to1600_hbb_mass","hh2500to3000_hbb_mass"};
  // std::vector<TString> pres = {"mu_L","mu_T","e_L","e_T","emu_T","emu_L","emu_LT"};
  // vector<unsigned int> sigMasses = {1000,1600,3000};
  //   for(unsigned int i = 0; i < sigMasses.size(); ++i){
  //     std::vector<TString> v1s = {vars[i]};
  //     std::vector<unsigned int> s1s = {sigMasses[i]};
  //     distPlots(TString::Format("plots_%s",vars[i].Data()),v1s,s1s,pres,false,2,false);
  //
  //   }


  
}


//Eff plots
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
    vector<TString> bkgs 
    = {
      "ttbar",
      "wjets",
      "zjets",
      "singlet",
      "diboson",
      "ttX",
      "qcd"
      // "other"
    };
    vector<TString> bkgNames 
    = {
      "t#bar{t}",
      "W+jets",
      "Z+jets",
      "t",
      "VV",
      "t#bar{t}X",
      "QCD"
      // "other"
    };
    
    // TString bkgNames = ";t#bar{t} / W+jets / Z+jets / t / VV / t#bar{t}X / QCD";

    TFile * f = new TFile("test_all.root","read");
    // TFile * f = new TFile("getCutflow.root","read");
    
    auto effPlots = [&](TString name, const TString& pre, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
        std::vector<TH1*> hists;
        for(unsigned int iN = 0; iN < cuts.size(); ++iN){
          hists.push_back(new TH1F(TString::Format("%s_%s_%u",name.Data(),pre.Data(),iN),";#it{m}(X) [TeV]",40,0.550,4.550));
        }

        for(auto m : fullSigMasses){
          TH1 * h = 0;
          f->GetObject(TString::Format("m%u_%s_evtCounts",m,pre.Data()),h);
          if(h==0){std:: cout <<TString::Format("m%u_%s_evtCounts",m,pre.Data())<<std::endl; continue;}
          for(unsigned int iC = 0; iC < cuts.size(); ++iC){
            hists[iC]->SetBinContent(hists[iC]->FindFixBin(float(m)/1000.),h->GetBinContent(cuts[iC]+1));
          }
        }
        Plotter * p = new Plotter;
        for(unsigned int iN = 0; iN < cuts.size(); ++iN){
          // hists[iN]->Scale(26.18/1000.);
          p->addHist(hists[iN],cutNames[cuts[iN]],-1,1,4,20,1,true,false);
        }
        p->setMinMax(0,1.0);
        p->drawRatio(false,"stack",false,false,name);

              // p->draw(false,name);
  };
  
  auto effBKGPlots = [&](TString name, const TString& pre, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
    const unsigned int iB = bkgs.size();
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%s_%u",name.Data(),pre.Data(),iN),";SM bkg",iB,-0.5,float(iB) -.5));
        auto * ax = hists.back()->GetXaxis();
        for(unsigned int b = 1; b <= iB; ++b) ax->SetBinLabel(b,bkgNames[b -1]);
      }


      for(unsigned int b = 0; b < bkgs.size(); ++b){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s_evtCounts",bkgs[b].Data(),pre.Data()),h);
        if(h==0){std:: cout <<TString::Format("%s_%s_evtCounts",bkgs[b].Data(),pre.Data())<<std::endl; continue;}
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(b),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHistLine(hists[iN],cutNames[cuts[iN]]);
      }
      // p->setMinMax(0,1.0);
      // p->drawRatio(false,"stack",false,false,name);

            p->draw(false,name);
};
  
  // std::vector<unsigned int> cuts = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};  
  // std::vector<unsigned int> cuts = {0,1,2,3,4,5,6,7};
    // std::vector<unsigned int> cuts = {0,8,9,10,11,12,13,14};
        // std::vector<unsigned int> cuts = {0,3,5,7,8,9,10,11,12};  
// std::vector<unsigned int> cuts = {0,3,5,7,8,14,15,16,17};
// std::vector<unsigned int> cuts = {0,20,26,27,28,29,30,31};
std::vector<unsigned int> cuts = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  vector<TString> cutNames = {
    "inclusive",  //emu
    "trigger presel",
    "cleaning",
    "lepton veto",
    "Wjj cand",
    "Wjj SD mass > 10 GeV",
    "Wjj #tau_{2}/#tau_{1} < 0.55" ,
    "Wjj anti b-tag",
    "Hbb cand",   //8
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag",
    "Hbb mass window", //15
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag",
    "HH mass 9-1.1 TeV", //22
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag",
    "HH mass 9-1.1 TeV, Hbb mass window", //29
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag",
    "HH mass 1.4-1.8 TeV", //36
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag",
    "HH mass 1.4-1.8 TeV, Hbb mass window", //43
    "Hbb LT tag",
    "WW #DeltaR < 2",
    "Wl#nu #DeltaR < #pi/2",
    "Hbb loose b-tag",
    "Hbb med b-tag",
    "Hbb tight b-tag", //43
  };
  
  std::vector<TString> presS = {"emu","e","mu","genemu_emu"};
    std::vector<TString> presB = {"emu","e","mu"};
  for(const auto& p : presB){
    effBKGPlots(p+"_bkg",p,cuts,cutNames);
  }
  for(const auto& p : presS){
    effPlots(p +"_sig",p,cuts,cutNames);
  }
  
}
