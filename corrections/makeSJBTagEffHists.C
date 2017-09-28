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
      "qcd",
      "rare"
    };
    vector<TString> bkgNamess = {
      "t#bar{t}",
      "w+jets",
      "QCD",
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

  // TFile * f = new TFile("getHHSel_plots.root","read");
    TFile * f = new TFile("sjBTagSFEffHists.root","read");
  



  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres, bool doNorm, float rebin = -1, bool addLumi=false){
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

          float stackNorm =  plots->getTotStack() ? plots->getStackIntegral() : 1.0;
          if(stackNorm) h->Scale(stackNorm/h->Integral(0,-1)); else PlotTools::normalize(h);
        } else if(plots->getTotStack()){
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




std::vector<TString> pres = {"wjj_m1000_","wjj_"};
// std::vector<TString> pres = {"wjjB_m1000_","wjjB_"};
std::vector<TString> vars = {"maxSJETA","minSJPT"};


distPlots("plots",vars,pres,true);




}



{
  TString inFN = "sjBTagSFEffHists.root";
  TString outFN = "sj_csvEff.root";
  TFile * inF = new TFile(inFN,"read");
  TFile * outF = new TFile(outFN,"recreate");
  
  TString denN = "incl";
  std::vector<TString> numNs {"loose","med","tight"};
  std::vector<TString> flvs {"l","c","b"};
  
  double ptBins[] = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,700,800,900,1000};
  int nPTBins = 22;
  
  double ptBins2[] = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,800,1000};
  int nPTBins2 = 20;
  
  float maxX = 900;
  
  double etaBins[] = {0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,2.4};
  int nETABins = 9;
  
  for(unsigned int iF = 0; iF < flvs.size(); ++iF){
    TH2 * hd = 0; 
    inF->GetObject(TString::Format("%s_%s_noxsec",flvs[iF].Data(),denN.Data()),hd);
    if(hd == 0) continue;
    
    for(unsigned int iN = 0; iN < numNs.size(); ++iN){
      TH2 * hn = 0; 
      TString name = TString::Format("%s_%s",flvs[iF].Data(),numNs[iN].Data());
      TString name2 = TString::Format("%s_%s_noxsec",flvs[iF].Data(),numNs[iN].Data());
      inF->GetObject(name2,hn);
      if(hn == 0) continue;
      
      auto doRebin =[&](const TH2* oh, const TString& newName, int nptbins, double * ptbins) ->TH2* {
        TH2* nh = new TH2F(newName,";jet p_{T}[GeV];jet |#eta|",nptbins,ptbins,nETABins,etaBins);
        for(unsigned int iX = 1; iX <= oh->GetNbinsX()+1; ++iX){
          for(unsigned int iY = 1; iY <= oh->GetNbinsY(); ++iY){
            float xVal = oh->GetXaxis()->GetBinCenter( iX == oh->GetNbinsX()+1 ? oh->GetNbinsX() : iX );
            float yVal = oh->GetYaxis()->GetBinCenter( iY );
            float val = oh->GetBinContent(iX,iY);
            nh->Fill(xVal,yVal,val);
          }
        }
        
        for(unsigned int iX = 0; iX <= nh->GetNbinsX()+1; ++iX ){
          nh->SetBinContent(iX,nh->GetNbinsY()+1,0);          
        }
        for(unsigned int iY = 0; iY <= nh->GetNbinsY()+1; ++iY ){
          nh->SetBinContent(0,iY,0);          
        }
        
        return nh;
      };

      hn = doRebin(hn,name + "new_N", iF !=1 ? nPTBins: nPTBins2 , iF !=1 ? ptBins:ptBins2    );
      hn->Divide(doRebin(hd,name + "new_D", iF !=1 ? nPTBins: nPTBins2 , iF !=1 ? ptBins:  ptBins2   ));
            
      // hn->Divide(hd);
      outF->cd();
      hn->Write(name);
      
    }
    
    
  }
  
  outF->Close();
  
}

{
  TString inFN = "sj_csvEff.root";
  TFile * inF = new TFile(inFN,"read");
  
  std::vector<TString> numNs {"loose","med"};
  std::vector<TString> flvs {"l","c","b"};
  
  for(unsigned int iF = 0; iF < flvs.size(); ++iF){

    for(unsigned int iN = 0; iN < numNs.size(); ++iN){
      TString name = TString::Format("%s_%s",flvs[iF].Data(),numNs[iN].Data());
      TH2 * hn = 0; 
      inF->GetObject(name,hn);
      if(hn == 0) continue;
      new TCanvas(name,name); hn->Draw("COLZ"); 
    }
    
    
  }

  
}
