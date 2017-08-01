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
  vector<unsigned int> sigMasses = {
    600,
    // 800,
    1000,
    // 1200,
    // 1400,
    // 1600,
    // 1800,
    2000,
    // 2500,
    // 3000,
    // 3500,
    // 4000,
    4500

  };
  vector<TString> ecNs = {
    "inclusive",
    "gen lepton p_{T}> 20 GeV, |#eta| <2.4",
    "reco lepton p_{T}> 20 GeV, |#eta| <2.4",
    "+ medium ID",
    "+ medium ID and miniIso",
    "#it{H}_{T} < 1200 GeV: lep sel, #it{H}_{T} > 1200 GeV: no lep sel",
    "medium ID and miniIso for #it{H}_{T} < 1200 GeV, id for #it{H}_{T} > 1200 GeV",
    "iso for AK8 jet",
    "iso for AK8 jet, id for AK8 jet",
    "iso for AK8 sd jet",
    "iso for AK8 sd jet, id for AK8 sd jet",
    "medium ID and miniIso lepton > 20 GeV + #it{H}_{T} > 400 GeV",
    "(lepton > 30 GeV) or (#it{H}_{T} > 400 GeV, lepton > 20 GeV)"
  };

  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres){
    for(auto v : vars)     for(auto p : pres){
      Plotter * plots = new Plotter;
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TFile * f = new TFile(TString::Format("out_radion_hh_bbinc_m%u_0.root",sigMasses[iM]),"read");
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s",p.Data(),v.Data()),h);
        if(!h){
          cout <<   TString::Format("%s_%s",p.Data(),v.Data()) <<" "<< TString::Format("#it{m}(X) %u GeV",sigMasses[iM])<<endl;
        }
        plots->addHistLine(h,TString::Format("#it{m}(X) %u GeV",sigMasses[iM]));
    }
    plots->normalize();
    // plots->rebin(2);
    plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    plots->yAxis()->SetTitleOffset(1.5);
  }
};

  auto effPlots = [&](TString name, TString prefix, std::vector<unsigned int>& cuts){
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%u",name.Data(),iN),";#it{m}(X) [TeV]",40,0.550,4.550));
      }

      for(auto m : fullSigMasses){
        TFile * f = new TFile(TString::Format("out_radion_hh_bbinc_m%u_0.root",m),"read");
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_event_count",prefix.Data()),h);
        if(h==0){
          std::cout << TString::Format("%s_event_count",prefix.Data())<<endl;
          continue;
        }
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(float(m)/1000.),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHist(hists[iN],ecNs[cuts[iN]],-1,1,4,20,1,true,false);
      }
      p->setMinMax(0,1.0);
      p->drawRatio(false,"stack",false,false,name);

            // p->draw(false);
};

std::vector<TString> vars = {"genLepWDR","genlep_pt","ht"};
std::vector<TString> pres = {"all_incl"};

std::vector<unsigned int> cuts = {0,1,2,3,4};
// std::vector<unsigned int> cuts = {0,4,5};
// std::vector<unsigned int> cuts = {0,4,11,12};


// distPlots("plots",vars,pres);
effPlots("muRatio","mu",cuts);
effPlots("elRatio","el",cuts);

}