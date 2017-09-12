

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
    vector<pair<TString,TString> > samps 
    = {
      std::pair<TString,TString>("ttbar","t#bar{t}"         ),
      std::pair<TString,TString>("m600" ,"#it{m}(X) 0.6 TeV"),
      // std::pair<TString,TString>("m800" ,"#it{m}(X) 0.8 TeV"),
      std::pair<TString,TString>("m1000","#it{m}(X) 1.0 TeV"),
      // std::pair<TString,TString>("m1200","#it{m}(X) 1.2 TeV"),
      // std::pair<TString,TString>("m1400","#it{m}(X) 1.4 TeV"),
      std::pair<TString,TString>("m1600","#it{m}(X) 1.6 TeV"),
      // std::pair<TString,TString>("m1800","#it{m}(X) 1.8 TeV"),
      // std::pair<TString,TString>("m2000","#it{m}(X) 2.0 TeV"),
      // std::pair<TString,TString>("m2500","#it{m}(X) 2.5 TeV"),
      std::pair<TString,TString>("m3000","#it{m}(X) 3.0 TeV"),
      // std::pair<TString,TString>("m3500","#it{m}(X) 3.5 TeV"),
      // std::pair<TString,TString>("m4000","#it{m}(X) 4.0 TeV"),
      std::pair<TString,TString>("m4500","#it{m}(X) 4.5 TeV")
    };
    // vector<TString> sampNamess = {
    //   "t#bar{t}"         ,
    //   "#it{m}(X) 0.6 TeV" ,
    //   "#it{m}(X) 0.8 TeV" ,
    //   "#it{m}(X) 1.0 TeV",
    //   "#it{m}(X) 1.2 TeV",
    //   "#it{m}(X) 1.4 TeV",
    //   "#it{m}(X) 1.6 TeV",
    //   "#it{m}(X) 1.8 TeV",
    //   "#it{m}(X) 2.0 TeV",
    //   "#it{m}(X) 2.5 TeV",
    //   "#it{m}(X) 3.0 TeV",
    //   "#it{m}(X) 3.5 TeV",
    //   "#it{m}(X) 4.0 TeV",
    //   "#it{m}(X) 4.5 TeV"
    // };
    
    TFile * f = new TFile("trueLepParam.root","read");
  



  auto distPlots = [&](TString name,const TString& num, const TString& den, float rebin = 0, int nR = 0, double * rebins = 0){
      Plotter * p = new Plotter;      
      for(const auto& s : samps){
        TH1 * hd = 0;
        f->GetObject(TString::Format("%s_%s",s.first.Data(),den.Data()),hd);
        TH1 * hn = 0;
        f->GetObject(TString::Format("%s_%s",s.first.Data(),num.Data()),hn);
        if(hd == 0 || hn == 0) continue;
        hn = (TH1*)hn->Clone();
        hd = (TH1*)hd->Clone();

        if(rebin > 0){
          PlotTools::rebin(hn,rebin);
          PlotTools::rebin(hd,rebin);
        } else if(rebins){
          hn = PlotTools::rebin(hn,nR,rebins);
          hd = PlotTools::rebin(hd,nR,rebins);
        }
        PlotTools::toOverflow(hn);
        PlotTools::toOverflow(hd);
        hn->Divide(hn,hd,1,1,"b"); 
        p->addHist(hn,s.second);        
      }
      p->draw(false,TString::Format("%s.pdf",name.Data()));
    };
    
    // vector<TString> leps = {"e","mu"};
        // vector<TString> leps = {"e_lt50","e_50to200","e_gt200","mu_lt50","mu_50to200","mu_gt200"};
    vector<TString> leps = {"mu","mu_lt50","mu_50to100","mu_100to200","mu_200to500","mu_gt500"};
    // vector<TString> leps = {"e","e_lt50","e_50to100","e_100to200","e_200to500","e_gt500"};
    
    // vector<TString> dens = {"gen_dr0p4_incl","reco_dr0p4_incl","reco_dr0p4_ID"};
    // vector<TString> nums = {"gen_dr0p4_reco","reco_dr0p4_ID"  ,"reco_dr0p4_ISO"};
    vector<TString> dens = {"reco_dr0p4_ID"};
    vector<TString> nums = {"reco_dr0p4_ISO"};
    // vector<TString> dens = {"gen_dr0p4_incl"};
    // vector<TString> nums = {"gen_dr0p4_reco"};
    
    // vector<TString> vars = {"ptj_o_ptl","ptjpl_o_ptl"  ,"ptjml_o_ptl","ptjmlpl_o_ptl","jml_dr_l"};
    vector<TString> vars = {"pt0p5_jml_drN_l"};
    
    for(unsigned int iL = 0; iL < leps.size(); ++iL){
      for(unsigned int iT = 0; iT < nums.size(); ++iT){
              for(unsigned int iV = 0; iV < vars.size(); ++iV){
                TString num = TString::Format("%s_%s_%s",leps[iL].Data(),nums[iT].Data(),vars[iV].Data());
                TString den = TString::Format("%s_%s_%s",leps[iL].Data(),dens[iT].Data(),vars[iV].Data());
                distPlots(num,num,den,2);
        
              }
      }      
    }            
}

new TCanvas(); m1000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(m1000_mu_reco_dr0p4_ID_ptl_o_pt1q_v_dr);m1000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");
new TCanvas(); m2000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(m2000_mu_reco_dr0p4_ID_ptl_o_pt1q_v_dr);m2000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");

new TCanvas(); m3000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(m3000_mu_reco_dr0p4_ID_ptl_o_pt1q_v_dr);m3000_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");

new TCanvas(); m3000_mu_reco_dr0p4_ID_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); m3000_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Divide(m3000_mu_reco_dr0p4_ID_ptq_o_ptl_v_drN);m3000_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); m1000_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Divide(m1000_mu_reco_dr0p4_ID_ptq_o_ptl_v_drN);m1000_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); m1600_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Divide(m1600_mu_reco_dr0p4_ID_ptq_o_ptl_v_drN);m1600_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); ttbar_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Divide(ttbar_mu_reco_dr0p4_ID_ptq_o_ptl_v_drN);ttbar_mu_reco_dr0p4_ISO_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); ttbar_mu_reco_dr0p4_ID_ptq_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); m1000_mu_50to100_reco_dr0p4_ISO_ptq_o_ptl_v_dr->Divide(m1000_mu_50to100_reco_dr0p4_ID_ptq_o_ptl_v_dr);m1000_mu_50to100_reco_dr0p4_ISO_ptq_o_ptl_v_dr->Draw("COLZ");


new TCanvas(); ttbar_e_gen_dr0p4_reco_ptq_o_ptl_v_drN->Divide(ttbar_e_gen_dr0p4_incl_ptq_o_ptl_v_drN);ttbar_e_gen_dr0p4_reco_ptq_o_ptl_v_drN->Draw("COLZ");
new TCanvas(); ttbar_e_gen_dr0p4_reco_ptq_o_ptl_v_dr->Divide(ttbar_e_gen_dr0p4_incl_ptq_o_ptl_v_dr);ttbar_e_gen_dr0p4_reco_ptq_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); m3000_e_gen_dr0p4_reco_ptq_o_ptl_v_dr->Divide(m3000_e_gen_dr0p4_incl_ptq_o_ptl_v_dr);m3000_e_gen_dr0p4_reco_ptq_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); ttbar_e_gen_dr0p4_incl_ptq_o_ptl_v_dr->Draw("COLZ")
new TCanvas(); m3000_e_gen_dr0p4_incl_ptq_o_ptl_v_dr->Draw("COLZ")


  new TCanvas(); m3000_e_gen_dr0p4_reco_ptqt_o_ptl_v_dr->Divide(m3000_e_gen_dr0p4_incl_ptqt_o_ptl_v_dr);m3000_e_gen_dr0p4_reco_ptqt_o_ptl_v_dr->Draw("COLZ");


new TCanvas(); ttbar_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(ttbar_mu_reco_dr0p4_ID_ptl_o_pt1q_v_dr);ttbar_mu_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");

new TCanvas(); ttbar_mu_100to200_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(ttbar_mu_100to200_reco_dr0p4_ID_ptl_o_pt1q_v_dr);ttbar_mu_100to200_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");
  
  new TCanvas(); ttbar_mu_lt50_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(ttbar_mu_lt50_reco_dr0p4_ID_ptl_o_pt1q_v_dr);ttbar_mu_lt50_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");
  new TCanvas(); m600_mu_lt50_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Divide(m600_mu_lt50_reco_dr0p4_ID_ptl_o_pt1q_v_dr);m600_mu_lt50_reco_dr0p4_ISO_ptl_o_pt1q_v_dr->Draw("COLZ");
  

new TCanvas(); m3000_mu_reco_dr0p4_ID_jml_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); m3000_e_gen_dr0p4_incl_jml_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); ttbar_e_gen_dr0p4_incl_jml_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); m3000_e_gen_dr0p4_reco_jml_o_ptl_v_dr->Divide(m3000_e_gen_dr0p4_incl_jml_o_ptl_v_dr);m3000_e_gen_dr0p4_reco_jml_o_ptl_v_dr->Draw("COLZ");
new TCanvas(); ttbar_e_gen_dr0p4_reco_jml_o_ptl_v_dr->Divide(ttbar_e_gen_dr0p4_incl_jml_o_ptl_v_dr);ttbar_e_gen_dr0p4_reco_jml_o_ptl_v_dr->Draw("COLZ");




//T&P data/mc sfplot for AN
{
  TFile * fd  = new TFile("data_triggerTurnons.root");
  Plotter * pt = new Plotter();
  
  
  auto getEff = [&](TFile * f, TString pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
    TH1 * hd = 0;
    f->GetObject(TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()),hd);
    TH1 * hn = 0;
    f->GetObject(TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn);
    if(hn == 0){
      cout << TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
      return 0;
    }
    if(hd == 0){
      cout << TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
      return 0;
    }
    hn = (TH1*)hn->Clone();
    hd = (TH1*)hd->Clone();

    if(rebin > 0){
      PlotTools::rebin(hn,rebin);
      PlotTools::rebin(hd,rebin);
    } else if(rebins){
      hn = PlotTools::rebin(hn,nR,rebins);
      hd = PlotTools::rebin(hd,nR,rebins);
    }
    PlotTools::toOverflow(hn);
    PlotTools::toOverflow(hd);
    hn->Divide(hn,hd,1,1,"b"); return hn; 
    // return PlotTools::getBinomErrors(hn,hd);   
  };
  
  auto plotTurnons =[&](TString oHName, TString name, TString prefix,TString dataName,TString mcName, const TString& sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0 ){      
      Plotter * p = new Plotter();
        auto * mcEff = getEff(fd,mcName,prefix,sel,var,trig,rebin,nR,rebins);
        auto * dataEff = getEff(fd,dataName,prefix,sel,var,trig,rebin,nR,rebins);
        if(mcEff == 0 || dataEff == 0) return;
        p->addHist(mcEff,"MC",-1,1,4,20,1,true,true, false, "E X P");
        p->addHist(dataEff,"data",-1,1,4,20,1,true,true, false, "E X P");
        p->setUnderflow(false);
        p->setMinMax(0.68,1.09);
        p->setBotMinMax(0.925,1.075);
        p->setCMSLumi();
        p->setCMSLumiExtraText("Preliminary");
        p->setYTitle("trigger efficiency");
        p->setYTitleBot("data/MC");
        p->setLegendPos(0.7,0.4,0.95,0.65);
        // p->setXTitleBot("data/MC");
        p->drawSplitRatio(0,"stack",false,true,TString::Format("%s.pdf",name.Data()));

        of->cd();        
        TH1 * rat = (TH1*)dataEff->Clone(oHName);
        rat->SetDirectory(0);
        rat->Divide(mcEff);
        rat->Write();
  };
  
  std::vector<TString> muSels = {"mupt_25to30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
  std::vector<TString> elSels = {"elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
    std::vector<TString> htSels = {"ht_450","ht_475","ht_500","ht_600","ht_700","ht_800","ht_1000","ht_1200"};
    int nLepBins = 17;
    double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
    // int nLepBins = 8;
    // double lepBins[] = {5,10,15,20,25,30,35,50,75};
    int nHTBins = 28;
    double htBins[] = {100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,700,800,900,1000,1100,1200,1600,2000};
    

    plotTurnons("muonSF",TString::Format("turnOn_muon_ht"),"GL_passSE","singlee","ttbar","mupt_26","ht","passSMuoHtMuoBu"        ,0,nHTBins,htBins);
    
    
    delete of;
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

    TFile * f = new TFile("getHighPTSF.root","read");
  



  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres,bool doSig = false, float rebin = 0){
    for(auto v : vars) for(auto p : pres){
      Plotter * plots = new Plotter;      
      
      for(unsigned int iS = 0; bkgs[iS][0]; ++iS){
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s_%s",bkgs[iS].Data(),p.Data(),v.Data()),h);
        if(h == 0) continue;
        plots->addStackHist(h,bkgNamess[iS]);
      }
      TH1 * hd = 0;
      f->GetObject(TString::Format("data_%s_%s",p.Data(),v.Data()),hd);
      if(hd != 0) plots->addHist(hd,"data");
      
      if(doSig)
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s_%s",sigMasses[iM],p.Data(),v.Data()),h);
        if(h ==0)continue;
        h->Scale(sigNorm);        
        plots->addHistLine(h,TString::Format("#it{m}(X) %.1f TeV",float(sigMasses[iM])/1000.));
    }
        

    
    if(rebin > 0) plots->rebin(rebin);
    plots->setYTitle("events/bin width");
    plots->setBotMinMax(0.0001,1.999);
    plots->setYTitleBot("N/N(MC)");
    plots->setCMSLumi();
    plots->setCMSLumiPosition(0,1.1);
    plots->setCMSLumiExtraText("Preliminary");
    // plots->setLegendPos(0.6,0.6,0.85,0.95);
    // plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    plots->drawSplitRatio(-1,"stack",false,false,TString::Format("%s_%s_%s.pdf",name.Data(),p.Data(),v.Data()));
  }
};
  vector<TString> vars = { "pt0p5_jml_drN_l","pass_pt0p5_jml_drN_l","nGoodMuons"};
    vector<TString> pres = { "ltj","mtj"};
  distPlots("plots",vars,pres,false);

}