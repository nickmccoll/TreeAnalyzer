{
    bool print = false;
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
    // 600,
    // 800,
    // 1000,
    // 1200,
    // 1400,
    // 1600,
    // 1800,
    2000
    // 2500,
    // 3000,
    // 3500,
    // 4000,
    // 4500

  };
  
  vector<TString> bkgs = {
    // "ttbar",
    "qcd"
  };
  vector<TString> bkgNamess = {
    // "t#bar{t}",
    "QCD"
  };

  TFile * f = new TFile("lepsel_plots.root","read");
  
  double rebins[] = {0,5,10,15,20,25,30,40,50,75,100,150,200,300,500};
  int nBins = 14;
  
  double rebinsBKG[] = {0,20,30,50,75,100,150,200,300,500};
  int nBinsBKG = 9;
  
  auto sigDistPlots = [&](TString name, TString var, const std::vector<TString>& vars, const std::vector<TString>& pres){
for(auto p : pres){       for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
      Plotter * plots = new Plotter();
      int iF =0;
    for(unsigned int iV = 0; iV < vars.size();++iV) {
      const auto& v = vars[iV];
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s%s_%s",sigMasses[iM],p.Data(),v.Data(),var.Data()),h);
        if(h == 0) continue;
        ++iF;
        plots->addHistLine(h,TString::Format("#it{m}(X) %u GeV, %s",sigMasses[iM],v.Data()),StyleInfo::getLineColor(iV));
    }
    plots->rebin(nBins,rebins);
    plots->setMinMax(0,1);
    if(iF > 1){
    plots->drawRatio(0,"stack",true,print,TString::Format("sigDist_%s_%s_%s_%u.pdf",name.Data(),p.Data(),var.Data(),sigMasses[iM]));
    plots->yAxis()->SetTitleOffset(1.5);
  }
  }}
};

  auto effPlots = [&](TString name, const std::vector<TString>& pres, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
    for(auto pr : pres){  
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%s_%u",name.Data(),pr.Data(),iN),";#it{m}(X) [TeV]",40,0.550,4.550));
      }

      for(auto m : fullSigMasses){
        TH1 * h = 0;
        f->GetObject(TString::Format("m%u_%s%s_selection",m,pr.Data(),name.Data()),h);
        if(h==0){continue;}
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(float(m)/1000.),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHist(hists[iN],cutNames[cuts[iN]],-1,1,4,20,1,true,false);
      }
      p->setMinMax(0,1.0);
      p->drawRatio(0,"stack",false,print,TString::Format("sigEff_%s_%s.pdf",name.Data(),pr.Data()));
    }

            // p->draw(false);
};

  auto bkgDistPlots = [&](TString name, TString var, const std::vector<TString>& vars, const std::vector<TString>& pres){
for(auto p : pres){       for(unsigned int iM = 0; iM < bkgs.size(); ++iM){
      Plotter * plots = new Plotter();
      int iF =0;
      TH1 *hd = 0;
    for(unsigned int iV = 0; iV < vars.size(); ++iV) {
        TH1 * h = 0;
        f->GetObject(TString::Format("%s_%s%s_%s",bkgs[iM].Data(),p.Data(),vars[iV].Data(),var.Data()),h);
        if(h == 0) continue;
        if(iV == 0){
          hd = h;                      
        } else {
        ++iF;        
        plots->addHistLine(h,TString::Format("%s, %s",bkgNamess[iM].Data(),vars[iV].Data()),StyleInfo::getLineColor(iV));
      }
    }

    if(iF && hd != 0){    
      plots->rebin(nBinsBKG,rebinsBKG);
      plots->scale(1.0/hd->Integral(0,-1));
      // plots->setMinMax(0,1);
    plots->draw(print,TString::Format("bkgEff_%s_%s_%s.pdf",p.Data(),var.Data(),bkgs[iM].Data()));
    plots->yAxis()->SetTitleOffset(1.5);
  }
  }}
};
  
//
//   auto distPlots = [&](TString name, const std::vector<TString>& vars, const std::vector<TString>& pres){
//     for(auto v : vars)     for(auto p : pres){
//       Plotter * plots = new Plotter;
//       for(unsigned int iM = 0; iM < sigMasses.size(); ++iM){
//         TH1 * h = 0;
//         f->GetObject(TString::Format("m%u_%s%s",sigMasses[iM],p.Data(),v.Data()),h);
//         // if(!h){
//         //   cout <<   TString::Format("%s_%s",p.Data(),v.Data()) <<" "<< TString::Format("#it{m}(X) %u GeV",sigMasses[iM])<<endl;
//         // }
//         plots->addHistLine(h,TString::Format("#it{m}(X) %u GeV",sigMasses[iM]));
//     }
//     plots->normalize();
//     // plots->rebin(2);
//     plots->draw(false,TString::Format("%s_%s_%s",name.Data(),p.Data(),v.Data()));
//     plots->yAxis()->SetTitleOffset(1.5);
//   }
// };
//
//   auto effPlots = [&](TString name, std::vector<unsigned int>& cuts, std::vector<TString>& cutNames){
//       std::vector<TH1*> hists;
//       for(unsigned int iN = 0; iN < cuts.size(); ++iN){
//         hists.push_back(new TH1F(TString::Format("%s_%u",name.Data(),iN),";#it{m}(X) [GeV]",40,550,4550));
//       }
//
//       for(auto m : fullSigMasses){
//         TH1 * h = 0;
//         f->GetObject(TString::Format("m%u_selection",m),h);
//         if(h==0){continue;}
//         for(unsigned int iC = 0; iC < cuts.size(); ++iC){
//           hists[iC]->SetBinContent(hists[iC]->FindFixBin(m),h->GetBinContent(cuts[iC]+1));
//         }
//       }
//       Plotter * p = new Plotter;
//       for(unsigned int iN = 0; iN < cuts.size(); ++iN){
//         p->addHist(hists[iN],cutNames[cuts[iN]],-1,1,4,20,1,true,false);
//       }
//       p->setMinMax(0,1.0);
//       p->drawRatio(false,"stack",false,false,name);
//
//             // p->draw(false);
// };
std::vector<TString> sigElpres = {"","passTight_","matchedEl_","matchedEl_passTight_"};
std::vector<TString> sigMupres = {"","passTight_","matchedMu_","matchedMu_passTight_"};
std::vector<TString> bkgPres = {"","passTight_"};
//MU ID
// std::vector<TString> muIDvars = {"incl","soft","loose","med","tight","highPT"};
// std::vector<unsigned int> muIDcuts = {0,1,2,3,4,5};
// sigDistPlots("muID","muID",muIDvars,sigMupres);
// effPlots("muID",sigMupres,muIDcuts,muIDvars);
// bkgDistPlots("muID","muID",muIDvars,sigMupres);

// //Mu ISO
// std::vector<TString> muISOvars = {"incl","miniIso_0p1","miniIso_0p2","miniIso_0p3","pfIso_0p1","pfIso_0p2","pfIso_0p3"};
// std::vector<unsigned int> muISOcuts = {0,1,2,3,4,5,6};
// sigDistPlots("muISO","muISO",muISOvars,sigMupres);
// effPlots("muISO",sigMupres,muISOcuts,muISOvars);
// bkgDistPlots("muISO","muISO",muISOvars,sigMupres);

// //Mu D0
// std::vector<TString> muD0vars = {"incl","0p01","0p02","0p05","0p10","0p20"};
// std::vector<unsigned int> muD0cuts = {0,1,2,3,4,5};
// sigDistPlots("muD0","muD0",muD0vars,sigMupres);
// effPlots("muD0",sigMupres,muD0cuts,muD0vars);
// bkgDistPlots("muD0","muD0",muD0vars,sigMupres);

//Mu SIP3D
// std::vector<TString> muS3vars = {"incl","2","3","4","5","6"};
// std::vector<unsigned int> muS3cuts = {0,1,2,3,4,5};
// sigDistPlots("muS3","muS3",muS3vars,sigMupres);
// effPlots("muS3",sigMupres,muS3cuts,muS3vars);
// bkgDistPlots("muS3","muS3",muS3vars,sigMupres);

// //Mu DZ
// std::vector<TString> muDZvars = {"incl","0p1","0p2","0p3","0p4","0p5"};
// std::vector<unsigned int> muDZcuts = {0,1,2,3,4,5};
// sigDistPlots("muDZ","muDZ",muDZvars,sigMupres);
// effPlots("muDZ",sigMupres,muDZcuts,muDZvars);
// bkgDistPlots("muDZ","muDZ",muDZvars,sigMupres);

//El ID
std::vector<TString> elIDvars = {"incl","loose","med","tight","heep","mvaLoose","mva80","mva90"};

std::vector<unsigned int> elIDcuts = {0,1,2,3,4,5,6,7};
sigDistPlots("elID","elID",elIDvars,sigElpres);
effPlots("elID",sigElpres,elIDcuts,elIDvars);
bkgDistPlots("elID","elID",elIDvars,sigElpres);

// El ISO
// std::vector<TString> elISOvars = {"incl","miniIso_0p1","miniIso_0p2","miniIso_0p3","miniIso_0p4","miniIso_0p5",
//                                          "miniIsoFP_0p1","miniIsoFP_0p2","miniIsoFP_0p3","miniIsoFP_0p4","miniIsoFP_0p5",
//                                            "pfIso_0p1","pfIso_0p2","pfIso_0p3","pfIso_0p4","pfIso_0p5",
//                                            "trackerIso_0p1","trackerIso_0p2","trackerIso_0p3","tracekrIso_0p4","trackerIso_0p5"};
// std::vector<TString> elISOvarsD = {"incl","miniIso_0p1","miniIso_0p2","miniIso_0p3","pfIso_0p1","pfIso_0p2","pfIso_0p3"};
//
// std::vector<unsigned int> elISOcuts = {0,1,2,3,11,12,13};
// sigDistPlots("elISO","elISO",elISOvarsD,sigElpres);
// effPlots("elISO",sigElpres,elISOcuts,elISOvars);
// bkgDistPlots("elISO","elISO",elISOvarsD,sigElpres);

//El D0
// std::vector<TString> elD0vars = {"incl","0p01","0p02","0p05","0p10","0p20"};
// std::vector<unsigned int> elD0cuts = {0,1,2,3,4,5};
// sigDistPlots("elD0","elD0",elD0vars,sigElpres);
// effPlots("elD0",sigElpres,elD0cuts,elD0vars);
// bkgDistPlots("elD0","elD0",elD0vars,sigElpres);

// //El SIP3D
// std::vector<TString> elS3vars = {"incl","2","3","4","5","6","7","8"};
// std::vector<unsigned int> elS3cuts = {0,1,2,3,4,5,6,7};
// sigDistPlots("elS3","elS3",elS3vars,sigElpres);
// effPlots("elS3",sigElpres,elS3cuts,elS3vars);
// bkgDistPlots("elS3","elS3",elS3vars,sigElpres);
//

//El SIP3D
// std::vector<TString> elDZvars = {"incl","0p1","0p2","0p3","0p4","0p5"};
// std::vector<unsigned int> elDZcuts = {0,1,2,3,4,5};
// sigDistPlots("elDZ","elDZ",elDZvars,sigElpres);
// effPlots("elDZ",sigElpres,elDZcuts,elDZvars);
// bkgDistPlots("elDZ","elDZ",elDZvars,sigElpres);



}