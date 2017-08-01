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
    1200,
    // 1400,
    1600
    // 1800,
    // 2000,
    // 2500,
    // 3000,
    // 3500,
    // 4000,
    // 4500

  };

  TFile * f = new TFile("getHbbJetCandSel_plots.root","read");
  

  auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres){
    for(auto v : vars)     for(auto p : pres){
      Plotter * plots = new Plotter;
      for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s%s",sigMasses[iM],p.Data(),v.Data()),h);
        // if(!h){
        //   cout <<   TString::Format("%s_%s",p.Data(),v.Data()) <<" "<< TString::Format("#it{m}(X) %u GeV",sigMasses[iM])<<endl;
        // }
        plots->addHistLine(h,TString::Format("#it{m}(X) %u GeV",sigMasses[iM]));
    }
    plots->normalize();
    // plots->rebin(2);
    plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
    plots->yAxis()->SetTitleOffset(1.5);
  }
};

  auto effPlots = [&](TString name, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%u",name.Data(),iN),";#it{m}(X) [GeV]",40,550,4550));
      }

      for(auto m : fullSigMasses){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_selection",m),h);
        if(h==0){continue;}
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(m),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHist(hists[iN],cutNames[cuts[iN]],-1,1,4,20,1,true,false);
      }
      p->setMinMax(0,1.0);
      p->drawRatio(false,"stack",false,false,name);

            // p->draw(false);
};

std::vector<TString> vars = {"nearestDR","ptRank","drlepRank","dphiWRank","dphilepRank","leadpt_drlepRank","leadpt_dphiWRank","leadpt_dphilepRank","pt","drlep","dphiW","dphilep"};
// std::vector<TString> pres = {""};
// std::vector<TString> vars = {"maxDRgenAK4Jet","drAK4Jets","min_dphiAK4Jetslep","min_dRAK4Jetslep","ak4jetptRank","ak4jetptPairRank","dphipairlep","drpairjj","bkg_dphipairlep","bkg_drpairjj","ak4jetptPairRank_byDR"};
std::vector<TString> pres = {""};
distPlots("plots",vars,pres);

std::vector<unsigned int> cuts = {1,3,4,6,7,8,10};
// std::vector<unsigned int> cuts = {10,11,12,13};
// std::vector<unsigned int> cuts = {1,4,9,10,13};
vector<TString> cutNames = {
  "inclusive",
  "good reco lepton",
  "gen h(bb) #it{p}_{T} > 50 GeV",
  "gen h(bb) #it{p}_{T} > 50 GeV (|#eta|<2.4)",
  "gen h(bb) #it{p}_{T},|#eta|, #DeltaR(bb)<0.8",
  "+ #geq1 #it{p}_{T} > 50 GeV |#eta|<2.4 fj",
  "+ matched fj",
  "+ not Wjj cand",
  "+ leading  fj",
  "+ greater |#Delta#phi(fj,l)|",
  "+ |#Delta#phi(fj,l)| > 2",
  "Good gen jets (wide angle)",
  "Good reco jets (wide angle)",
  "|#Delta#phi| > #pi/2",
  "leading pair"
    
};
effPlots("cutflow",cuts,cutNames);
// std::vector<unsigned int> cuts = {0,4,5};
// std::vector<unsigned int> cuts = {0,4,11,12};



// effPlots("muRatio","mu",cuts);
// effPlots("elRatio","el",cuts);

}