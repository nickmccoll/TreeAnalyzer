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
  //for ttbar cr
  // std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  // std::vector<TString> selNs = {"search region","invert extra b veto","invert extra b veto, relax extra cuts"};
  // std::vector<TString> samps = {"bkg","QCD","m1000","m2000","m3000"};
  
  //For WJets
  std::vector<TString> sels = {"sr","hbb1","hbb1_noM_noEx"};
  std::vector<TString> selNs = {"search region","hbb=0","hbb=0, relax extra cuts"};
  std::vector<TString> samps = {"bkg","qg","m1000","m2000","m3000"};
  
  

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
  bool addData = true;
  TString var = "hhMass";
  //ttbar cr
  // std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  // std::vector<TString> samps = {"qg","losttw","mw","mt"};
  // std::vector<TString> sampNs = {"q/g bkg.","lost t/W bkg.","#it{m}_{W} bkg.","#it{m}_{t} bkg."};
  
  //wjets cr
  std::vector<TString> sels = {"hbb1","e_hbb1","e_hbb1_noM_noEx","mu_hbb1","mu_hbb1_noM_noEx"};
  // std::vector<TString> sels = {"hbbLT4_noM_noEx","hbb3_noM_noEx","hbb2_noM_noEx","hbb1_noM_noEx"};
  // std::vector<TString> sels = {"hbbLT4_noM_noEx_noNJ","hbb3_noM_noEx_noNJ","hbb2_noM_noEx_noNJ","hbb1_noM_noEx_noNJ"};
  std::vector<TString> samps = {"qg_wQCD","losttw","mw","mt"};
  std::vector<TString> sampNs = {"q/g bkg.","lost t/W bkg.","#it{m}_{W} bkg.","#it{m}_{t} bkg."};
  
  // std::vector<TString> samps = {"QCD","wjets","qg_ttbar","qg_singlet","qg_other"};
  // std::vector<TString> sampNs = {"QCD","wjets","ttbar","singlet","other"};
  
    for(const auto& sel : sels ){
          Plotter * p = new Plotter();
          for(unsigned int iS = 0; iS < samps.size(); ++iS){
    const auto& samp = samps[iS];
      TH1 * h = 0;
      f->GetObject(samp+"_"+sel+"_"+var,h);
      if(h==0)continue;
      p->addStackHist(h,sampNs[iS]);
    }
    
    if(addData){
      TH1 * h = 0;
      f->GetObject("data_"+sel+"_"+var,h);
      if(h!=0) p->addHist(h,"data");
    }
    
    p->rebin(5);
    // p->draw(false,sel);
    p->drawSplitRatio(-1,"stack",false,false,sel);
  }
  
}


{
  TFile * f = new TFile("HHlnujj_findCRPlots.root");
  //ttbar cr
  // std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  // std::vector<TString> samps = {"bkg","qg","losttw","mw","mt"};
  // std::vector<TString> sampNs = {"all","q/g bkg.","lost t/W bkg.","#it{m}_{W} bkg.","#it{m}_{t} bkg."};
  
  // //wjets cr
    std::vector<TString> sels = {"e_hbb1","mu_hbb1"};
  // std::vector<TString> sels = {"sr","ab_noM_noEx","hbb1","hbb1_noM_noEx"};
  
  std::vector<TString> samps = {"qg_wQCD","QCD","wjets","qg_ttbar","qg_other"};
  std::vector<TString> sampNs = {"all q/g bkg.","QCD","W+jets","t#bar{t}","other"};
  
    for(const auto& sel : sels ){
          Plotter * p = new Plotter();
          for(unsigned int iS = 0; iS < samps.size(); ++iS){
    const auto& samp = samps[iS];
      TH1 * h = 0;
      f->GetObject(samp+"_"+sel+"_hhMass",h);
      if(h==0)continue;
      p->addHist(h,sampNs[iS]);
    }
    p->rebin(10);
    p->setBotMinMax(0,1);
    p->drawSplitRatio(0,"stack",false,false,sel);
  }
  
}


///QCD Scaling
{
  TFile * f = new TFile("HHlnujj_QCDAll_distributions.root");
  //ttbar cr
  // std::vector<TString> sels = {"sr","ab","ab_noM_noEx"};
  // std::vector<TString> samps = {"bkg","qg","losttw","mw","mt"};
  // std::vector<TString> sampNs = {"all","q/g bkg.","lost t/W bkg.","#it{m}_{W} bkg.","#it{m}_{t} bkg."};
  
  // //wjets cr
    // std::vector<TString> sels = {"emu_LMT_I_none","emu_LMT_I_ltmb","emu_LMT_I_full","e_LMT_I_ltmb","mu_LMT_I_ltmb","emu_L_I_ltmb","emu_M_I_ltmb","emu_T_I_ltmb","emu_LMT_LP_ltmb","emu_LMT_HP_ltmb"};
  // std::vector<TString> sels = {"L_LP_full","L_HP_full","M_LP_full","M_HP_full","T_LP_full","T_HP_full"};
    std::vector<TString> sels = {"e_L_LP_full","e_L_HP_full","e_M_LP_full","e_M_HP_full","e_T_LP_full","e_T_HP_full","mu_L_LP_full","mu_L_HP_full","mu_M_LP_full","mu_M_HP_full","mu_T_LP_full","mu_T_HP_full"};

  
  // std::vector<TString> samps = {"QCDNoIso","QCD"};
  // std::vector<TString> sampNs = {"QCDNoIso","QCD"};
    
    std::vector<TString> samps = {"QCDWOthers","QCDWOnlyOthers","QCD"};
    std::vector<TString> sampNs = {"QCDWOthers","QCDWOnlyOthers","QCD"};
  
    for(const auto& sel : sels ){
          Plotter * p = new Plotter();
          for(unsigned int iS = 0; iS < samps.size(); ++iS){
    const auto& samp = samps[iS];
      TH1 * h = 0;
      f->GetObject(samp+"_"+sel+"_hhMass",h);
      if(h==0)continue;
      p->addHist(h,sampNs[iS]);
    }
    p->rebin(10);
    p->setBotMinMax(0,1);
    p->setMinMax(0.001,1000);
    auto* c= p->drawSplitRatio(0,"stack",false,false,sel);
        c->GetPad(1)->SetLogy();
        c->GetPad(1)->Update();
  }
  
}



// Get QCD Ratio
{
  TFile * f = new TFile("HHlnujj_bkg_getQCDRatio.root");
  TString denSel = "wjets";
  TString numSel = "qcd";
  TString var = "hhMass";
  int rebin = 5;
  
  std::vector<TString> leps = {"e","mu"};
  // std::vector<TString> purs = {"I","LP","HP"};
    std::vector<TString> purs = {"I"};
  std::vector<TString> hads = {"ltmb","full"};
  // std::vector<TString> btags = {"AB","LMT","L","M","T"};
  std::vector<TString> btags = {"AB","LMT"};
  
  
  for(const auto& l : leps ){
    for(const auto& p : purs ){
      for(const auto& h : hads ){
        Plotter * plot = new Plotter();
        for(const auto& b : btags ){
          TH1 * hn = 0;
          f->GetObject(numSel+"_"+l+"_"+b+"_"+p+"_"+h+"_"+var ,hn);          
          TH1 * hd = 0;
          f->GetObject(denSel+"_"+l+"_"+b+"_"+p+"_"+h+"_"+var ,hd);          
          if(hn == 0 || hd==0) continue;
          PlotTools::rebin(hn,rebin);
          PlotTools::rebin(hd,rebin);
          hn->Divide(hd);
          plot->addHist(hn,b);        
        }      
        // p->setMinMax(0,1);
        plot->setUnderflow(false);
        auto* c= plot->draw(false,l+"_"+p+"_"+h);        
      }
    }
  }
}