
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getAntiBTag.C+("../testTrees/out_radion_hh_bbinc_m3000_0.root",2,"test_m3000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getAntiBTag.C+("../testTrees/out_radion_hh_bbinc_m800_0.root",2,"test_m800.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getAntiBTag.C+("../testTrees/out_radion_hh_bbinc_m1000_0.root",2,"test_m1000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getAntiBTag.C+("../testTrees/out_t-tbar1l-madgraph_18.root",1,"test_mttbar.root")' &
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
      // "zjets",
      // "singlet",
      // "diboson",
      // "ttX",
      "qcd",
      "other"
    };
    vector<TString> bkgNames 
    = {
      "t#bar{t}",
      "W+jets",
      // "Z+jets",
      // "t",
      // "VV",
      // "t#bar{t}X",
      "QCD",
      "other"
    };
    
    // TString bkgNames = ";t#bar{t} / W+jets / Z+jets / t / VV / t#bar{t}X / QCD";


    TFile * f = new TFile("getAntiBTag.root","read");
        // TFile * f = new TFile("test_all.root","read");
    
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
        p->setMinMax(0.8,1.0);
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
      p->setMinMax(0,1.0);
      p->drawRatio(false,"stack",false,false,name);

            // p->draw(false,name);
};
  

// std::vector<unsigned int> cuts = {0,1,2,3,4,5,6,7,8};
// std::vector<unsigned int> cuts = {27,28,29,30,31,32,33,34,35};
// std::vector<unsigned int> cuts = {36,37,38,39,40,41,42,43,44};

std::vector<unsigned int> cuts = {2,3,4,5,6,7,8};
// std::vector<unsigned int> cuts = {0,8};

// std::vector<unsigned int> cuts = {29,30,31,32,33,34,35};
  vector<TString> cutNames = {
    "inclusive", //0 HBBI LI
    "SJL",
    "SJM",
    "SJL JL",
    "SJL JM",
    "SJM JL",
    "SJM JM",
    "JL",
    "JM",    //8
    "inclusive", //9 MF
    "SJL",
    "SJM",
    "SJL JL",
    "SJL JM",
    "SJM JL",
    "SJM JM",
    "JL",
    "JM",    //17
        "inclusive", //18 ML
      "SJL",
      "SJM",
      "SJL JL",
      "SJL JM",
      "SJM JL",
      "SJM JM",
      "JL",
      "JM",    //26        
      "inclusive", //27 MM
      "SJL",
      "SJM",
      "SJL JL",
      "SJL JM",
      "SJM JL",
      "SJM JM",
      "JL",
      "JM",    //35
        
        "inclusive", //36 MM
        "SJL",
        "SJM",
        "SJL JL",
        "SJL JM",
        "SJM JL",
        "SJM JM",
        "JL",
        "JM"    //44
  };
  
  std::vector<TString> presS = {"pwl_30_tID_emu","pwl_20_tID_emu"};
    // std::vector<TString> presB = {"emu_pwlep30","emu_pwlep"};
  for(const auto& p : presS){
    effBKGPlots(p+"_bkg",p,cuts,cutNames);
  }
  for(const auto& p : presS){
    effPlots(p +"_sig",p,cuts,cutNames);
  }
  
}

{

  vector<TString> bkgs 
  = {
    "ttbar",
    // "wjets",
    // "zjets",
    // "singlet",
    // "diboson",
    // "ttX",
    // "qcd",
    "other"
  };
  
  TFile * f = new TFile("getAntiBTag.root","read");
  std::vector<TString> abT = {"sjM","jM"};
  std::vector<TString> pres = {"L","M","T","L_hbbTM","M_hbbTM","T_hbbTM"};
  
  for(const auto& b : bkgs){
  for(const auto& p : pres){
    Plotter * p2 = new Plotter();
      for(const auto& a : abT){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_pwl_30_tID_emu_%s_%s_hh_mass",b.Data(),p.Data(),a.Data()),h);
        if(h == 0) continue;
        p2->addHist(h,a);
      }
      p2->rebin(5);
      p2->drawSplitRatio(0,"stack",false,false,TString::Format("%s_%s",b.Data(),p.Data()));
  }
  }
  
  

}
