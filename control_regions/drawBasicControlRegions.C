/*
hadd data_e_basicCRPlots.root ../jobs_cr2/*data_e*.root &
hadd data_mu_basicCRPlots.root ../jobs_cr2/*data_mu*.root &
hadd data_jetht_basicCRPlots.root ../jobs_cr2/*data_jetht*.root &
hadd data_met_basicCRPlots.root ../jobs_cr2/*data_met*.root &
*/
{
	auto getH = [&](TString name, TString hist) ->TH1*{
		TFile *f = new TFile(name,"read");
		TH1 * h = 0;
		f->GetObject(hist,h);
		return h;
	};
	TString hist = "data_1el_leptonPT";
	Plotter * p = new Plotter();
	p->addHist  (getH("data_mu_basicCRPlots.root",hist),"singlemu");
	p->addHist  (getH("data_e_basicCRPlots.root",hist),"singlee");
	p->addHist  (getH("data_jetht_basicCRPlots.root",hist),"jetht");
	p->addHist  (getH("data_met_basicCRPlots.root",hist),"met");
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

    TFile * f = new TFile("basicCRPlots.root","read");
  



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
    plots->drawSplitRatio(-1,"stack",false,true,TString::Format("%s_%s_%s.pdf",name.Data(),p.Data(),v.Data()));
  }
};
  // vector<TString> vars = { "numvtx","woPUW_numvtx","ht"};
  // vector<TString> vars = { "ak4Wlep_nJets","ak4Wlep_nLooseBtags","ak4Wlep_nMedBtags","ak4Wlep_nTightBtags"};
  // vector<TString> vars = { "ak4Nolep_nJets","ak4Nolep_nLooseBtags","ak4Nolep_nMedBtags","ak4Nolep_nTightBtags"};
// vector<TString> vars = { "ak4Nolep_nJets","ak4Nolep_nMedBtags"};
  // vector<TString> vars = { "ak8_nJets","ak8_jetPT","ak8_jetETA","ak8_jetMass","ak8_jetSDMass","ak8_jetRSDMass","ak8_NSubjets"};
  // vector<TString> vars = { "ak8_leadingJetNSubjets","ak8_leadingJetMinSJCSV","ak8_leadingJetMaxSJCSV","ak8_leadingJetTau2oTau1","ak8_leadingJetSDMass"};
    // vector<TString> vars = {"ak8_leadingJetSDMass"};
  // vector<TString> vars = { "ak8_leadingJetNSubjets","ak8_leadingJetMinSJCSV","ak8_leadingJetMaxSJCSV","ak8_leadingJetTau2oTau1","ak8_leadingJetSDMass"};
  // vector<TString> vars = { "hhMass","hbbMass","hbbPT","hbbSDMass","hbbRSDMass"};
  vector<TString> vars = { "hhMass","hbbSDMass"};
  // vector<TString> vars = { "ak4Wlep_nJets","ak4Nolep_nJets"};
  // vector<TString> pres = { "1mu","1el","2lsf","2lof"};
    // vector<TString> pres = { "1mu_0b","1el_0b","1mu_2b","1el_2b"};
  // vector<TString> pres = { "1mu","1el","1mu_2b","1el_2b"};
    vector<TString> pres = { "1mu_HbbAntiBCR","1el_HbbAntiBCR"};
    // vector<TString> pres = { "1mu_wjjBtagTauCR","1el_wjjBtagTauCR"};
  distPlots("plots",vars,pres,false,5);

}