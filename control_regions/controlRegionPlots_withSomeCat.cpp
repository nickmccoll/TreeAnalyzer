///Categorization
{
// TH2 * h2 = bkg_categories;
TH2 * h2 = signal_categories;
new TCanvas();
h2->Draw("COLZ");
h2 = (TH2*)h2->Clone();
for(unsigned int iX =1 ; iX <= h2->GetNbinsX(); ++iX){
 double tot = h2->Integral(iX,iX,0,-1);
for(unsigned int iY =1 ; iY <= h2->GetNbinsY(); ++iY){
h2->SetBinContent(iX,iY,h2->GetBinContent(iX,iY)/tot);
} 
}
new TCanvas();

h2->Draw("COLZ");

}

{
TH2 * h2 = bkg_categories;
new TCanvas();
h2->Draw("COLZ");
h2 = (TH2*)h2->Clone();
for(unsigned int iX =1 ; iX <= h2->GetNbinsX(); ++iX){
 double tot = h2->Integral(iX,iX,0,-1);
for(unsigned int iY =1 ; iY <= h2->GetNbinsY(); ++iY){
h2->SetBinContent(iX,iY,h2->GetBinContent(iX,iY)/tot);
} 
}
new TCanvas();

h2->Draw("COLZ");

}

{ TFile * f = new TFile("plots_all.root");
  std::vector<TString> sels = {"m1000","m2500"};
  std::vector<TString> bC = {"bL","bM","bT"};
  std::vector<TString> hC = {"wML_wTH","wMH_wTH","wML_wTL","wMH_wTL"};

  
    for(unsigned int iS = 0; iS < sels.size(); ++iS){
      Plotter * sPlot = new Plotter();
      Plotter * bPlot = new Plotter();
      Plotter * rPlot = new Plotter();
      
  for(unsigned int iB = 0; iB < bC.size(); ++iB){
    int lS = 1;
    if(iB == 1)lS = 7;
    if(iB == 2)lS = 9;
      for(unsigned int iH = 0; iH < hC.size(); ++iH){
        TH1 * hs = 0;
        f->GetObject("signal_"+bC[iB]+"_"+hC[iH]+"_"+sels[iS]+"_tauEff",hs);
        TH1 * hb = 0;
        f->GetObject("bkg_"+bC[iB]+"_"+hC[iH]+"_"+sels[iS]+"_tauEff",hb);
        if(hs==0 || hb==0) continue;
        sPlot->addHistLine(hs,bC[iB]+"_"+hC[iH],StyleInfo::getLineColor(iH),lS);
        bPlot->addHistLine(hb,bC[iB]+"_"+hC[iH],StyleInfo::getLineColor(iH),lS);
        
        hs = (TH1*)hs->Clone();
        hs->Divide(hb);
        rPlot->addHistLine(hs,bC[iB]+"_"+hC[iH],StyleInfo::getLineColor(iH),lS);        
      }
  }
  
  sPlot->draw(false,sels[iS] +"_sig");
  bPlot->draw(false,sels[iS] +"_bkg");
  rPlot->draw(false,sels[iS] +"_ratio");
}
}

///CONTROL REGION PLOTS
{
  TFile * f = new TFile("HHlnujj_findCRPlots.root");
  std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  std::vector<TString> selNs = {"search region","invert extra b veto","invert extra b veto, relax extra cuts"};
  std::vector<TString> samps = {"bkg","QCD","m1000","m2000","m3000"};

  for(const auto& samp : samps ){
    Plotter * p = new Plotter();
    for(unsigned int iS = 0; iS < sels.size(); ++iS){
      const auto& sel = sels[iS];
      TH1 * h = 0;
      f->GetObject(samp+"_"+sel+"_hbbMass",h);
      if(h==0)continue;
      p->addHistLine(h,selNs[iS]);
    }
    p->rebin(5);
    p->drawRatio(0,"stack",false,false,samp);
  }
  
}


{
  TFile * f = new TFile("HHlnujj_findCRPlots.root");
  // std::vector<TString> sels = {"sr","ab","ab_noM","ab_noEx","ab_noM_noEx"};
    std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  std::vector<TString> samps = {"qg","losttw","mw","mt"};
  std::vector<TString> sampNs = {"q/g bkg.","lost t/W bkg.","#it{m}_{W} bkg.","#it{m}_{t} bkg."};
    for(const auto& sel : sels ){
          Plotter * p = new Plotter();
          for(unsigned int iS = 0; iS < samps.size(); ++iS){
    const auto& samp = samps[iS];
      TH1 * h = 0;
      f->GetObject(samp+"_"+sel+"_hbbMass",h);
      if(h==0)continue;
      p->addStackHist(h,sampNs[iS]);
    }
    p->draw(false,sel);
  }
  
}