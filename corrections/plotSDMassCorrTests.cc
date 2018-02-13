{
  TFile * f = new TFile("plots_allM.root","read");
  TString prefix = "goodHBB";
  std::vector<TString> vars = {"rawMass","corrMass"};
  std::vector<TString> varNs = {"#it{m}_{H#rightarrowbb}","C_{SD}*#it{m}_{H#rightarrowbb}"};
  std::vector<float> masses = {1000,2000,3000,4000};
  
  Plotter * p = new Plotter();
  
  for(unsigned int iV = 0; iV < vars.size();++iV){
        for(unsigned int iM = 0; iM < masses.size();++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%.0f_%s_%s",masses[iM],prefix.Data(),vars[iV].Data()),h);
        if(h == 0) continue;
        p->addHistLine(h,TString::Format("#it{m}_{X} %.0f, %s",masses[iM],varNs[iV].Data()),StyleInfo::getLineColor(iM),iV+1);
        
      }
      
  }
  p->setYTitle("arbitrary units");
  p->normalize();
  p->rebin(2);
  p->draw();
  
  
}



{

    // vector<TString> bkgs
    // = {
    //   "ttbar",
    //   "wjets",
    //   "zjets",
    //   "qcd",
    //   "other"
    // };
  // vector<TString> bkgNamess = {
  //   "t#bar{t}",
  //   "w+jets",
  //   "z+jets",
  //   "QCD",
  //   "other"
  // };
    vector<TString> bkgs 
    = {
      "other",
      "zjets",
      "wjets",
      "qcd",
      "ttbar",
    };
    vector<TString> bkgNamess = {
      "other",
      "z+jets",
      "w+jets",
      "QCD",
      "t#bar{t}",




    };
  vector<unsigned int> sigMasses = {
    // 600,
    800,
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
  const float sigNorm = 20.0/1000.0;

    TFile * f = new TFile("testSDMassCorrectionsInData.root","read");
  
    int nBins = 40;
    double bins[] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200};

    auto cutHistogram= [](const TH1* inH,  double nMinX, double nMaxX)->TH1*{
        double binWX = inH->GetXaxis()->GetBinWidth(1);
        std::string name = inH->GetName();
        name +="_cut";
        
      TH1 * outH = new TH1F(name.c_str(),TString(";")+inH->GetXaxis()->GetTitle(),
              (nMaxX-nMinX)/binWX,nMinX,nMaxX);

      for(int iX =1; iX <= inH->GetNbinsX(); ++iX){
        const int outIX =outH->GetXaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(iX));
        if(outIX < 1 || outIX > outH->GetNbinsX() ) continue;
        outH->SetBinContent(outIX,inH->GetBinContent(iX));
        outH->SetBinError(outIX,inH->GetBinError(iX));
      }
      return outH;
    };



  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres,bool doSig = false, float rebin = 0){
    for(auto v : vars) for(auto p : pres){
      Plotter * plots = new Plotter;      
      
      TH1 *ht = 0;
      f->GetObject(TString::Format("bkg_%s_%s",p.Data(),v.Data()),ht);

      TH1 * hd = 0;
      f->GetObject(TString::Format("data_%s_%s",p.Data(),v.Data()),hd);
      // double min =120;
      // double max =210;
      // int lb = hd->FindFixBin(160);
      // int hb = hd->FindFixBin(210);
      double min =30;
      double max =210;      
      int lb = hd->FindFixBin(80);
      int hb = hd->FindFixBin(100);
      float SF = hd->Integral(lb,hb)/ht->Integral(lb,hb);
      
      hd = cutHistogram(hd,min,max);
      ht = cutHistogram(ht,min,max);

      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s_%s",bkgs[iS].Data(),p.Data(),v.Data()),h);
        if(h == 0) continue;
        h = cutHistogram(h,min,max);
        h->SetXTitle("H(bb) corrected SD mass [GeV]"); 
        h->Scale(SF);
        plots->addStackHist(h,bkgNamess[iS]);
      }
        
      TH1 *hs = 0;
      TString t = p;
      f->GetObject(TString::Format("bkg_%s_%s",t.ReplaceAll("corr","scorr").Data(),v.Data()),hs);
      if(hs != 0 && ht != 0){
        hs = cutHistogram(hs,min,max);
        hs->Scale(SF);
        ht->Scale(SF);
        plots->addHistLine(hs,"MC +JER unc.",StyleInfo::getLineColor(4));
        hs = (TH1*)hs->Clone();
        hs->Add(ht,-1);
        TH1 *hu = (TH1*)ht->Clone();
        hu->Add(hs,-1);
        plots->addHistLine(hu,"MC -JER unc.",StyleInfo::getLineColor(3));

      }

      TH1 * scU = 0;
      TH1 * scD = 0;
      TString t2 = p;
      TString t3 = p;
      f->GetObject(TString::Format("bkg_%s_%s",t2.ReplaceAll("corr","sUcorr").Data(),v.Data()),scU);
      f->GetObject(TString::Format("bkg_%s_%s",t3.ReplaceAll("corr","sDcorr").Data(),v.Data()),scD);
      if(scU != 0 && scD != 0){
        scU = cutHistogram(scU,min,max);
        scD = cutHistogram(scD,min,max);
        scU->Scale(SF);
        scD->Scale(SF);
        plots->addHistLine(scU,"MC +JES unc.",StyleInfo::getLineColor(2));
        plots->addHistLine(scD,"MC -JES unc.",StyleInfo::getLineColor(1));
      }
      
      

      

      
      if(doSig)
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s_%s",sigMasses[iM],p.Data(),v.Data()),h);
        if(h ==0)continue;
        h = cutHistogram(h,min,max);
        h->Scale(sigNorm);
        plots->addHistLine(h,TString::Format("#it{m}(X) %.1f TeV",float(sigMasses[iM])/1000.));
    }

        
            if(hd != 0) plots->addHist(hd,"data",StyleInfo::getLineColor(0));
    
    if(rebin > 0){
      plots->rebin(rebin);
    }
    plots->setYTitle("number of events / 5 GeV");
    plots->setBotMinMax(0.0001,1.999);
    plots->setYTitleBot("N/N(MC)");
    plots->setCMSLumi();
    plots->setCMSLumiPosition(0,1.1);
    plots->setCMSLumiExtraText("Preliminary");
    plots->setMinMax(0,1200);
    plots->setBotMinMax(0.5001,1.499);
    plots->setLegendNColumns(2);
    // plots->setLegendPos(0.6,0.6,0.85,0.95);
    // plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    plots->drawSplitRatio(-1,"stack",false,true,TString::Format("%s_%s_%s.pdf",name.Data(),p.Data(),v.Data()));
  }
};

  // vector<TString> vars = { "njets","njetswlep","ht_wlep","ht_wlep_plusmet"};
  vector<TString> vars = { "hbbMass","hh800to1000_hbbMass","hh1500to2000_hbbMass"};
    // vector<TString> pres = { "1m_L_corr"};
    
  vector<TString> pres = { "1m_L_corr","1m_AL_corr"};
  distPlots("plots",vars,pres,false,5);

}

//compare corrections by sample
{
  vector<TString> samps 
  = {
    "wjets",
    "ttbar",
    "m1000",
    "m1600"
  };
  vector<TString> vars 
  = {
    "raw",
    "wCorr",
    "corr"
  };
  vector<TString> sels 
  = {
    "hbbI_hh900to1100",
    "hbbI_hh1400to1800"
  };

  TFile * f = new TFile("testSDMassCorrectionsInData.root","read");
  for(unsigned int iS = 0; iS < samps.size(); ++iS){
    Plotter * p = new Plotter();
    for(unsigned int iSel = 0; iSel < sels.size(); ++iSel){
      for(unsigned int iV = 0; iV < vars.size(); ++iV){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_1l_%s_%s_hbb_mass",samps[iS].Data(),vars[iV].Data(),sels[iSel].Data()),h);        
        if(h == 0) continue;
        p->addHistLine(h,TString::Format("%s, %s",sels[iSel].Data(),vars[iV].Data()),StyleInfo::getLineColor(iV),iSel+1);
      }
      
    }
    p->normalize();
    p->rebin(2);
    p->draw(false,samps[iS]);
  }
  
}


//compare norm
{
  vector<TString> samps 
  = {
    "wjets",
    "ttbar"
  };
  vector<TString> vars 
  = {
    "raw",
    "wCorr",
    "corr"
  };
  
  vector<TString> signals
    = {
      "m1000",
      "m1600"
    };
  
  vector<TString> sels 
  = {
    "hbbI_hh900to1100",
    "hbbI_hh1400to1800"
  };

  TFile * f = new TFile("testSDMassCorrectionsInData.root","read");
  for(unsigned int iS = 0; iS < samps.size(); ++iS){
    Plotter * p = new Plotter();
    for(unsigned int iSel = 0; iSel < sels.size(); ++iSel){
      for(unsigned int iV = 0; iV < vars.size(); ++iV){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_1l_%s_%s_hbb_mass",samps[iS].Data(),vars[iV].Data(),sels[iSel].Data()),h);        
        if(h == 0) continue;
        h = (TH1*)h->Clone();
        TH1 * hs = 0;
        f->GetObject(TString::Format("%s_1l_%s_%s_hbb_mass",signals[iSel].Data(),vars[iV].Data(),sels[iSel].Data()),hs);        
        if(hs == 0) continue;
        hs = (TH1*)hs->Clone();
        // h->Rebin(2);
        // hs->Rebin(2);
        hs->Scale(1./hs->Integral(0,-1));
        h->Scale(1./h->Integral(0,-1));
        hs->Divide(h);
        
        p->addHistLine(hs,TString::Format("%s, %s",signals[iSel].Data(),vars[iV].Data()),StyleInfo::getLineColor(iV),iSel+1);
      }
      
    }
    p->draw(false,samps[iS]);
  }
  
}



//checkBosonMom....used to verify jet cuts

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
  vector<unsigned int> sigMasses = {
    // 600,
    800,
    1000,
    // 1200,
    // 1400,
    1600,
    // 1800,
    // 2000,
    // 2500
    3000
    // 3500,
    // 4000
    // 4500

  };

    TFile * f = new TFile("checkBosonMom.root","read");
  



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
  

  std::vector<TString> vars = {"hh_mass"};
  // std::vector<TString> vars = {"hbb_pt"};
  // std::vector<TString> pres = {"fj_hbbL","fj_hbbT","sdCorr_hbbL","sdCorr_hbbT","sdRaw_hbbL","sdRaw_hbbT"};
// std::vector<TString> pres = {"fj_hbbL_hbb90to140","fj_hbbT_hbb90to140","sdCorr_hbbL_hbb90to140","sdCorr_hbbT_hbb90to140","sdRaw_hbbL_hbb90to140","sdRaw_hbbT_hbb90to140"};
    // std::vector<TString> pres = {"fj_hbbT_hbb90to140_hh900to1100","sdCorr_hbbT_hbb90to140_hh900to1100","sdRaw_hbbT_hbb90to140_hh900to1100"};
    // std::vector<TString> pres = {"fj_hbbL_hbb90to140_hh900to1100","sdCorr_hbbL_hbb90to140_hh900to1100","sdRaw_hbbL_hbb90to140_hh900to1100"};
  // std::vector<TString> pres = {"fj_hbbT_hbb90to140_hh1400to1800","sdCorr_hbbT_hbb90to140_hh1400to1800","sdRaw_hbbT_hbb90to140_hh1400to1800"};
    // std::vector<TString> pres = {"fj_hbbL_hbb90to140_hh1400to1800","sdCorr_hbbL_hbb90to140_hh1400to1800","sdRaw_hbbL_hbb90to140_hh1400to1800"};
    // std::vector<TString> pres = {"fj_hbbT_hbb90to140_hh2500to3500","sdCorr_hbbT_hbb90to140_hh2500to3500","sdRaw_hbbT_hbb90to140_hh2500to3500"};
  
      // std::vector<TString> pres = {"fj_hbbL_hbb90to140_hh900to1100","fj_hbbT_hbb90to140_hh900to1100","fj_hbbL_hbb90to140","fj_hbbT_hbb90to140","fj_hbbL_hh900to1100","fj_hbbT_hh900to1100","fj_hbbL","fj_hbbT"};
  
        std::vector<TString> pres = {"fj_hbbL_hbb90to140","fj_hbbT_hbb90to140","fj_hbb90to140"};
  
  distPlots("plots",vars,pres,true,5,false);
}
