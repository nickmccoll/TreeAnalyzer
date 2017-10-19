rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getRecoHHCuts.C+("../testTrees/out_radion_hh_bbinc_m3000_0.root",2,"test_m3000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getRecoHHCuts.C+("../testTrees/out_radion_hh_bbinc_m800_0.root",2,"test_m800.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getRecoHHCuts.C+("../testTrees/out_radion_hh_bbinc_m1000_0.root",2,"test_m1000.root")' &
rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/getRecoHHCuts.C+("../testTrees/out_t-tbar1l-madgraph_18.root",1,"test_mttbar.root")' &
hadd -f test_all.root test_m*.root
  

    


{
  TFile * f = new TFile("test_all.root","read");
  // std::vector<float> masses = {800,1000,3000};
    // std::vector<float> masses = {1000,3000};
      std::vector<float> masses = {1000,3000};
  // std::vector<TString> types = {"noCorr","corrHbb","corrHbbMET"};
  std::vector<TString> types = {"nuH","nuH2","nuW2","nuWsdW"};
  // std::vector<TString> sels = {"hbbL","hbbT","hbbI"}; 
    std::vector<TString> sels = {"hbbI"}; 
  std::vector<TString> vars = {"wlnu_mass","wjj_mass","hWW_mass","hh_mass","wjj_wlnu_dR","hWWv_mass","hWWvm_mass"}; 
  
  for(const auto& v : vars){
    for(unsigned int iS = 0; iS < sels.size(); ++iS){
      Plotter * p = new Plotter();
      for(unsigned int iM = 0; iM < masses.size(); ++iM)
        for(unsigned int iT = 0; iT < types.size(); ++iT){
          TH1 * h = 0; f->GetObject(TString::Format("m%.0f_%s_%s_%s",masses[iM],types[iT].Data(),sels[iS].Data(),v.Data()),h);
          if(h == 0) continue;
          if(iM > 1)
            p->addHistLine(h,TString::Format("M(%.0f), %s",masses[iM],types[iT].Data()),StyleInfo::getLineColor(iT), 2);
          else
            p->addHistLine(h,TString::Format("M(%.0f), %s",masses[iM],types[iT].Data()),-1 );
        }
        p->normalize();
        p->rebin(2);
        p->draw(false,TString::Format("%s_%s",sels[iS].Data(),v.Data()));
    }
  }
  
  
}



{
  Plotter * p = new Plotter();
  p->addHistLine(m1000_hbbI_h900to1100_hWWv_mass,"m1000");
  p->addHistLine(ttbar_hbbI_h900to1100_hWWv_mass,"ttbar");
  p->addHistLine(m1000_hbbI_h900to1100_hWWvsd_mass,"m1000_2");
  p->addHistLine(ttbar_hbbI_h900to1100_hWWvsd_mass,"ttbar_2");
  p->normalize();
  p->draw();
}

{
  Plotter * p = new Plotter();
  p->addHistLine(m1000_hbbI_h900to1100_wjj_wlnu_dR,"m1000");
  p->addHistLine(ttbar_hbbI_h900to1100_wjj_wlnu_dR,"ttbar");
  p->normalize();
  p->draw();
}



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
      // "ttbar"
      // "wjets",
      "qcd"
      // "other"
    };
    vector<TString> bkgNamess = {
      // "t#bar{t}"
      // "w+jets",
      "QCD"
      // "other"
    };
  vector<unsigned int> sigMasses = {
    // 600,
    // 800,
    1000,
    // 1200,
    // 1400,
    1600
    // 1800,
    // 2000,
    // 2500
    // 3000
    // 3500,
    // 4000
    // 4500

  };

    TFile * f = new TFile("getRecoHHCuts.root","read");
  



  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres, bool doNorm, float rebin = -1, bool addLumi=false){
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
  
  auto makeRocs  = [&](std::vector<TString> vars, const std::vector<int>& sigMasses, TString prefix, TString name,  bool doForward = true ){
    Plotter * p = new Plotter();
    Plotter * pSoB = new Plotter();
    for(const auto& v : vars) {
      TH1 *hb = 0;
      Plotter * pI = new Plotter();
    
      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;        
        f->GetObject(TString::Format("%s_%s_%s",bkgs[iS].Data(),prefix.Data(),v.Data()),h);        
        if(h == 0) continue;
        if(hb == 0) hb = (TH1*)h->Clone();
        else hb->Add(h);
      }
      float stackNorm = hb->Integral(0,-1);
      
      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;        
        f->GetObject(TString::Format("%s_%s_%s",bkgs[iS].Data(),prefix.Data(),v.Data()),h);        
        if(h == 0) continue;
        TH1 *hI = PlotTools::getIntegral(h,doForward,false);                              
        hI->Scale(1./stackNorm);
        pI->addStackHist(hI,bkgNamess[iS].Data());
      }
      
    
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%i_%s_%s",sigMasses[iM],prefix.Data(),v.Data()),h);
        if(h ==0)continue;
      
        TGraph * roc = PlotTools::getRocCurve(h,hb,doForward,"signal","bkg");
        p->addGraphLine(roc,TString::Format("#it{m}(X) %i GeV, %s",sigMasses[iM],v.Data()));        
        TH1 *hI = PlotTools::getIntegral(h,doForward,true);                              
        // hI->Scale(stackNorm);
        pI->addHistLine(hI,TString::Format("#it{m}(X) %i GeV",sigMasses[iM]));      
    }      
          
    pI->draw(false,TString::Format("%s_int_%s",name.Data(), v.Data()));    
    }
    p->draw(false,name.Data());
  };
  

  // std::vector<TString> vars = {"hWWv_mass","wjj_wlnu_dR","hWWTM2_mass","WTM_mass","hWWTM_mass",""};
    // std::vector<TString> vars = {"wlnu_pt","hbb_nu_dPhi","lep_nu_dPhi","lep_nu_dR","wjj_wlnu_dRM","wjj_wlnu_dR",""};
        // std::vector<TString> vars = {"wlnu_pt","hbb_nu_dPhi","lep_nu_dPhi","lep_nu_dR",""};
        
                // std::vector<TString> vars = {"hbb_nu_dPhi","lep_nu_dPhi","lep_nu_dR"};
                std::vector<TString> vars = {"hbb_nu_dPhi","hbb_nu_dR","lep_nu_dPhi","lep_nu_dR","hWW_nu_dPhi","hWW_nu_dR","wlnu_nu_dPhi","wlnu_nu_dR"};

                std::vector<TString> pres = {"hbbI_hh900to1100","hbbI_hh1400to1800"};
                // std::vector<TString> pres = {"hbbI_hh900to1100","hbbI_hh1400to1800"};
                                // std::vector<TString> pres = {"hbbI_hh900to1100","hbbT_hh900to1100","hbbI_hh1400to1800"};
                                                                // std::vector<TString> pres = {"hbbI_hh900to1100","hbbT_hh900to1100","hbbI_hh1400to1800"};
                                                                // std::vector<TString> pres = {"hbbI_hhInc"};
  
  // distPlots("plots",vars,pres,true,2,false);
                std::vector<int> allMass = {1000,1600};
                std::vector<int> lowMass = {1000};
                std::vector<int> highMass = {1600};
        for(unsigned int iP = 0; iP < pres.size(); ++iP){
          if(iP <= 0)
            makeRocs(vars,lowMass,pres[iP],pres[iP],false);
          else if (iP >= 1) makeRocs(vars,highMass,pres[iP],pres[iP],false);
          else  makeRocs(vars,allMass,pres[iP],pres[iP],false);
        }
  
}

///do everything
{
  TFile * f = new TFile("getRecoHHCuts.root","read");
  std::vector<int> allMass = {800,1000,2000,3000,4000};
  // std::vector<TString> vars = {"wjj_wlnu_dRM","lep_nu_dR"};
  std::vector<TString> vars = {"hbb_nu_dPhi","lep_nu_dR","hWW_nu_dR","wlnu_nu_dPhi","wlnu_nu_dR"};
  
  TString prefix = "hbbT_hhInc";
  
for(const auto& v: vars ){
  Plotter * p = new Plotter();
    for(const auto& m: allMass ){
      TH1 * h = 0;
      f->GetObject(TString::Format("m%i_%s_%s",m,prefix.Data(),v.Data()),h);
      if(h ==0)continue;
      TH1 *hI = PlotTools::getIntegral(h,false,true);                              
      p->addHistLine(hI,TString::Format("#it{m}(X) %i GeV, %s",m,v.Data()));        
    }
    p->draw(false,v);
  }
  
}

///Try out dists
{
  TFile * f = new TFile("getRecoHHCuts.root","read");
  
  std::vector<TString> samps = {"signal","ttbar","qcd"};
  TString exP = "hbbI_";
  std::vector<TString> pres = {"hWW_pt_lt500","hWW_pt_500to750","hWW_pt_750to1000","hWW_pt_1000to1500","hWW_pt_geq1500","hWWInc"};
  // std::vector<TString> vars = {"wlnu_dR","wlnu_minPToTotPT","wlnu_deltaRPToMH","wlnu_deltaRPToMT","wlnu_dPhi","wlnu_pt"};
    // std::vector<TString> vars = {"hWW_dR","hWW_minPToTotPT","hWW_deltaRPToMH","hWW_deltaRPToMT","hWW_dPhi","hWW_pt"};
    std::vector<TString> vars = {"hWWSD_dR","hWWSD_minPToTotPT","hWWSD_deltaRPToMH","hWWSD_deltaRPToMT","hWWSD_dPhi","hWWSD_pt"};
  
  for(unsigned int iV = 0; iV < vars.size(); ++iV){
    for(unsigned int iS = 0; iS < samps.size(); ++iS){
      Plotter * p = new Plotter();
      for(unsigned int iP = 0; iP < pres.size(); ++iP){
        if(iS && iP+2 ==pres.size() ) continue;
        else if(!iS && iP+1 ==pres.size()) continue;
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s%s_%s",samps[iS].Data(),exP.Data(),pres[iP].Data(),vars[iV].Data()),h);
        if(h ==0)continue;
        p->addHistLine(h,pres[iP]);
      }
      p->normalize();
      TString name = TString::Format("%s_%s",samps[iS].Data(),vars[iV].Data());
      p->draw(false,name);
    }
  }
  
  
}