{
  TFile * f = new TFile("testLSF_all.root");
  std::vector<TString> sigs = {"m1000","m2000"};
  std::vector<TString> bkgCs = {"","htgeq1500_"};
  
  std::vector<TString> sigNs = {"#it{m}_{X} 1 TeV","#it{m}_{X} 2 TeV"};
  
  // std::vector<TString> bkgs = {"ttbar","wjets","qcd"};
  std::vector<TString> bkgs = {"qcd"};
  
  TString sel = "mu_inclPT";
  
  

  // std::vector<TString> vars = {"miniIso","pfIso","lsf"};
    std::vector<TString> vars = {"miniIso","pfIso"};
    std::vector<TString> varNs = {"mini-iso","std. iso"};
    // std::vector<TString> vars = {"likeliBHWW","likeliHWW"};
  int cutGT = 2;
  for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
    Plotter * p = new Plotter();
    for(unsigned int iV = 0; iV < vars.size(); ++iV){
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
        TH1* hs = 0;
        f->GetObject(sigs[iS]+"_"+sel+"_"+vars[iV],hs);     
        if(hs==0)continue;
        hs=(TH1*)hs->Clone();
        TH1* hb = 0;
        f->GetObject(bkgs[iB]+"_"+bkgCs[iS]+sel+"_"+vars[iV],hb);     
        if(hb==0)continue;
        hb=(TH1*)hb->Clone();
        
        
        PlotTools::toUnderflow(hb);
        PlotTools::toOverflow(hb);
        PlotTools::toUnderflow(hs);
        PlotTools::toOverflow(hs);
        
        hs->SetBinContent(0,0);
        hb->SetBinContent(0,0);
        
        auto * roc = PlotTools::getRocCurve(hs,hb,iV==cutGT, "signal eff","bkg. eff");
        // p->addGraph(roc,sigs[iS]+"_"+vars[iV]);        
        p->addGraph(roc,sigNs[iS]+", "+varNs[iV]);

      }                  
    }
    p->setMinMax(0.001,1);
    p->draw(false,bkgs[iB])      ;
  }
}
