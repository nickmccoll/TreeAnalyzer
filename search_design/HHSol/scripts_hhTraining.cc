rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2016_radion.root","2016_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2016_bulkgrav.root","2016_bulkgrav")' &
rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2017_radion.root","2017_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2017_bulkgrav.root","2017_bulkgrav")' &
rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2018_radion.root","2018_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(0,"hSolTrees/hhSolTrees_2018_bulkgrav.root","2018_bulkgrav")' &
hadd -f hSolTrees_HPTRes_2016.root hSolTrees_HPTRes_2016_*.root
hadd -f hSolTrees_HPTRes_2017.root hSolTrees_HPTRes_2017_*.root
hadd -f hSolTrees_HPTRes_2018.root hSolTrees_HPTRes_2018_*.root  

rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2016_radion.root","2016_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2016_bulkgrav.root","2016_bulkgrav")' &
rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2017_radion.root","2017_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2017_bulkgrav.root","2017_bulkgrav")' &
rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2018_radion.root","2018_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(3,"hSolTrees/hhSolTrees_2018_bulkgrav.root","2018_bulkgrav")' &
hadd -f hSolTrees_temp_2016.root hSolTrees_temp_2016_*.root
hadd -f hSolTrees_temp_2017.root hSolTrees_temp_2017_*.root
hadd -f hSolTrees_temp_2018.root hSolTrees_temp_2018_*.root
rr -b -q 'HHSol/getHHTraining.C+(4,"2016","")' &
rr -b -q 'HHSol/getHHTraining.C+(4,"2017","")' &
rr -b -q 'HHSol/getHHTraining.C+(4,"2018","")' &
  
  
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2016_bkg.root"     ,"2016_bkg")' &    
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2016_radion.root"  ,"2016_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2016_bulkgrav.root","2016_bulkgrav")' &
  
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2017_bkg.root"     ,"2017_bkg")' &    
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2017_radion.root"  ,"2017_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2017_bulkgrav.root","2017_bulkgrav")' &

rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2018_bkg.root"     ,"2018_bkg")' &    
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2018_radion.root"  ,"2018_radion")' &
rr -b -q 'HHSol/getHHTraining.C+(5,"hSolTrees/hhSolTrees_2018_bulkgrav.root","2018_bulkgrav")' &

hadd -f hSolTrees_test_2016.root hSolTrees_test_2016_*.root
hadd -f hSolTrees_test_2017.root hSolTrees_test_2017_*.root
hadd -f hSolTrees_test_2018.root hSolTrees_test_2018_*.root


rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2016_bkg.root"     ,"2016_bkg")' &    
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2016_radion.root"  ,"2016_radion")' &
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2016_bulkgrav.root","2016_bulkgrav")' &
  
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2017_bkg.root"     ,"2017_bkg")' &    
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2017_radion.root"  ,"2017_radion")' &
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2017_bulkgrav.root","2017_bulkgrav")' &

rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2018_bkg.root"     ,"2018_bkg")' &    
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2018_radion.root"  ,"2018_radion")' &
rr -b -q 'HHSol/getHHCuts.C+(0,"hSolTrees/hhSolTrees_2018_bulkgrav.root","2018_bulkgrav")' &
hadd -f hSolTrees_getCuts_2016.root hSolTrees_getCuts_2016_*.root
hadd -f hSolTrees_getCuts_2017.root hSolTrees_getCuts_2017_*.root
hadd -f hSolTrees_getCuts_2018.root hSolTrees_getCuts_2018_*.root
  
  
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2016_bkg.root"     ,"2016_bkg")' &    
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2016_radion.root"  ,"2016_radion")' &
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2016_bulkgrav.root","2016_bulkgrav")' &

rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2017_bkg.root"     ,"2017_bkg")' &    
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2017_radion.root"  ,"2017_radion")' &
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2017_bulkgrav.root","2017_bulkgrav")' &

rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2018_bkg.root"     ,"2018_bkg")' &    
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2018_radion.root"  ,"2018_radion")' &
rr -b -q 'HHSol/getHbbTaggingCutHistograms.C+("hSolTrees/hhSolTrees_2018_bulkgrav.root","2018_bulkgrav")' &
  hadd -f hSolTrees_hbbTaggingCutHistograms_2016.root hSolTrees_hbbTaggingCutHistograms_2016_*.root
hadd -f hSolTrees_hbbTaggingCutHistograms_2017.root hSolTrees_hbbTaggingCutHistograms_2017_*.root
hadd -f hSolTrees_hbbTaggingCutHistograms_2018.root hSolTrees_hbbTaggingCutHistograms_2018_*.root

//Fit HHRes
{
  TFile * f = new TFile("hSolTrees_HPTRes_2018.root");
  std::vector<float> ptBins {300,400,500,600,800,1000,1200,1400,1600,1800,2000};
  
  TGraphErrors * gF = new TGraphErrors();
  TGraphErrors * gC = new TGraphErrors();
  TGraphErrors * gCor = new TGraphErrors();
  int iGP = 0;
  for(unsigned int iP = 0; iP < ptBins.size()-1; ++iP){
    TString hName = TString::Format("signal_pt_%.0f_%.0f",ptBins[iP],ptBins[iP+1]);    
    TH1 * h = 0;
    TH1 * hc = 0;
    TH1 * hPT = 0;
    
    f->GetObject(hName+"_HPtRes",h);         
    f->GetObject(hName+"_true_Hpt",hPT);
    f->GetObject(hName+"_corrHPtRes",hc);
    std::cout <<hName<<" "<< h <<" "<< hPT << std::endl;
    if(h==0||hPT==0) continue;
    gF->SetPoint(iGP,hPT->GetMean(),h->GetMean() );
    gF->SetPointError(iGP,0,h->GetMeanError());
    h->GetXaxis()->SetRange(31,71);
    gC->SetPoint(iGP,hPT->GetMean(),h->GetMean() );
    gC->SetPointError(iGP,0,h->GetMeanError());
    
    if(hc){
      hc->GetXaxis()->SetRange(31,71);
      gCor->SetPoint(iGP,hPT->GetMean(),hc->GetMean() );
      gCor->SetPointError(iGP,0,hc->GetMeanError());              
    }

    iGP++;    
  }
  
  Plotter * p = new Plotter();
  p->addGraph(gF,"Full range");
  p->addGraph(gC,"0.6 #rightarrow 1.4");
  p->addGraph(gCor,"Corrected 0.6 #rightarrow 1.4");
  
  auto res = gC->Fit("pol1","S");
  cout << "RESULT: " << res->GetParams()[0] <<" "<< res->GetParams()[1] <<endl;
  
  
  
  p->draw();
}


///Look at dists
{
  TFile * f = new TFile("hhtrain_plots.root");
  std::vector<TString> sigs = {"radion","graviton"};
  std::vector<int> masses = {800,1000,3500};
  std::vector<TString> vars = {"hwwMag","true_hwwMag","hwwMagRes","extraMetPerp","extraMetParRelhwwMag","wqqPTRes","virtualSD_approxW","onshellSD_approxW","virtualSD_approxH","onshellSD_approxH","virtualSD_qqSDMassCoarse","onshellSD_qqSDMassCoarse"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    for(unsigned int iM = 0; iM < masses.size(); ++iM){
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = TString::Format("%s_m%i",sigs[iS].Data(),masses[iM]) +"_"+vars[iV];
      f->GetObject(hname,h);      
      p->addHistLine(h,TString::Format("%s: %.1f",sigs[iS].Data(),masses[iM]/1000.));
    }
  }
    p->normalize();
    p->drawSplitRatio(0,"stack",false,false,vars[iV]);
  }
  
}


///Figure out binning (fine)
{
  TFile * f = new TFile("hSolTrees_temp.root");
  std::vector<TString> sigs = {"signal_low","signal_test","signal_high"};
  std::vector<TString> shells = {"","_vqq","_osqq"};
  std::vector<double> avgs(3);
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    TH1 * h = 0;
    f->GetObject(sigs[iS]+"_true_hwwMag",h);
    if(h) avgs[iS] = h->GetMean();
  }      
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","approxW","approxH","approxJ","approxJ2","qqSDMassCoarse","hWW","Wlnu"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[1] +"_"+vars[iV];
      f->GetObject(hname,h);      
      if(h==0) continue;
      nH++;    
      p->addHistLine(h,sigs[iS]);
  }
  if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    p->draw(false,vars[iV]);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;    
      for(unsigned int iS = 0; iS < shells.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[0]+shells[iS] +"_"+vars[iV];
      f->GetObject(hname,h); 
      if(h==0) continue;
      nH++;         
      p->addHistLine(h,sigs[0]+shells[iS]);
  }
    if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]+"_shells");
    p->draw(false,vars[iV]+"_shells");
    // objs.push_back(c);
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}
///Figure out binning (final binning)
{
  TFile * f = new TFile("hSolTrees_temp.root");
  std::vector<TString> sigs = {"signal_low","signal_test","signal_high"};
  std::vector<TString> shells = {"_vqq","_osqq"};
    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","approxW","approxH","qqSDMassCoarse"};
  std::vector<int> vBinV      = {30  ,50  ,50  ,15,25,1};
  std::vector<double> vBinLEV = {-150,-0.6,-0.5,60,90,1};
  std::vector<double> vBinHEV = {150 ,0.4 ,0.5 ,90,190,1};
  
  std::vector<int> vBinO      = {30  ,50  ,50  ,15,20,1};
  std::vector<double> vBinLEO = {-150,-0.6,-0.5,0 ,110,1};
  std::vector<double> vBinHEO = {150 ,0.4 ,0.5 ,60,150,1};
  
  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
            Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[iSh] +"_"+vars[iV];                  
      f->GetObject(hname,h);      
      if(!h) continue;
      
      int nNB; double nL; double nH;
      if(iSh == 0){
        nNB = vBinV[iV];
        nL = vBinLEV[iV];
        nH = vBinHEV[iV];
      } else {
        nNB = vBinO[iV];
        nL = vBinLEO[iV];
        nH = vBinHEO[iV];
      }
      
      if(nL==nH) continue;
      TH1 * hR = new TH1F(TString(h->GetName()) +"_reb",h->GetTitle(),nNB,nL,nH);
      hR->SetBinContent(0,h->GetBinContent(0));
      hR->SetBinContent(nNB+1,h->GetBinContent(h->GetNbinsX()+1));      
      for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB ){
        int hRB = hR->FindFixBin(h->GetBinCenter(iB));
        
        hR->SetBinContent(hRB, hR->GetBinContent(hRB) + h->GetBinContent(iB) );
      }
      
      
      p->addHistLine(hR,sigs[iS]);
  
}
    p->normalize();
    p->setBotMinMax(0,2);
    TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]);
    objs.push_back(c);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}


//////////////////////////////////TEST RESULTS
  //Mass shape comparison
{
  TString year = "2016";
  TFile * f = new TFile("hSolTrees_test_"+year+".root");
  std::vector<TString> sigs = {"radion_m1000","radion_m2000","radion_m3000","ttbar","wjets","qcd"};
  // std::vector<TString> vars = {"simpleHH","chi2HH","likeliHH"};
    std::vector<TString> vars = {"simpleHH","likeliHH"};
  std::vector<TString> varNs = {"Old method","Likelihood method"};
  vector<TObject*> objs;

  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * p = new Plotter();
    for(unsigned int iV = 0; iV < vars.size(); ++iV){
      TH1* h = 0;
      f->GetObject(sigs[iS]+"_"+vars[iV],h);
      if(h==0) continue;
      p->addHistLine(h,varNs[iV]);
    }
    p->rebin(2);
    p->setBotMinMax(0,2);
    p->setXTitle(sigs[iS]+" #it{m}_{HH} [GeV]");
    p->setYTitle("Arbitary units");
    auto c = p->drawSplitRatio(0,"stack",false,false,sigs[iS]);
    objs.push_back(c);
    
  }
  Drawing::drawAll(objs,year);
  


}

//Efficiency curves of the likli/chi
{
  TFile * f = new TFile("hSolTrees_test_2018.root");
  // std::vector<TString> sigs = {"radion_m1000","radion_m2000","radion_m4000","ttbar_m1000","wjets_m1000","ttbar_m2000","wjets_m2000"};
    std::vector<TString> sigs = {"radion_m1000","radion_m2000","ttbar_m1000","wjets_m1000"};
    std::vector<TString> sigNs = {"1 TeV #it{m}_{X}","2 TeV #it{m}_{X}","t#bar{t}","W+jets"};
  // std::vector<TString> vars = {"chi2","likeli","OSQQ","OSQQPZ","VQQ","VQQPZ"};
  // std::vector<TString> vars = {"chi2","likeli","likeli_nAlt","alt"};
    // std::vector<TString> vars = {"chi2","Blikeli","Blikeli_nAlt","bAlt"};
        std::vector<TString> vars = {"likeli"};
        // std::vector<TString> vars = {"md"};
 
  
  // std::vector<TString> sigs = {"radion_m1000","radion_m2000","radion_m3000","ttbar_m1000","wjets_m1000"};
  // std::vector<TString> sigNs = {"#it{m}_{X} 1 TeV","#it{m}_{X} 2 TeV", "#it{m}_{X} 3 TeV","t#bar{t}","W+jets"};
  
  // std::vector<TString> sigs = {"radion_m1000","radion_m3000","ttbar_m1000"};
  // std::vector<TString> sigNs = {"#it{m}_{X} 1 TeV","#it{m}_{X} 3 TeV","t#bar{t}"};
  // std::vector<TString> vars = {"chi2","likeli","likeliNoSD","Blikeli"};


    for(unsigned int iV = 0; iV < vars.size(); ++iV){  
    Plotter * p = new Plotter();
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1* h = 0;
      f->GetObject(sigs[iS]+"_"+vars[iV],h);
      if(h==0) continue;
        // h = PlotTools::getIntegral(h,false,true);
      p->addHistLine(h,sigNs[iS]);
      // p->addHistLine(h,sigs[iS]);
    }
    p->setBotMinMax(0,2);
    // p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    p->normalize();
    p->rebin(2);
    // p->setXTitle("miniminized -2*log(L)");
      p->setXTitle("#it{D}");
      p->setMinMax(0.00001,1);
    auto c = p->draw(false,vars[iV]);
    c->SetLogy();
    c->Update();
  }
  
  
}


///Straight eff check
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets","qcd"};
  std::vector<TString> vars = {"simpleHH","chi2HH","likeliHH"};
  // std::vector<TString> cuts = {"chilt19","chilt11","mdlt125","lllt1p45","lllt1p35"};
    // std::vector<TString> cuts = {"mdlt125","chilt19","chilt11","lllt1p6","lllt1p45","lllt1p35"};
    std::vector<TString> cuts = {"mdlt125","chilt11","lllt1p45","lllt1p25","lllt1p2","lllt1p1"};


  //start With sigal
  Plotter * p = new Plotter();
  TH1* hd = 0;
  f->GetObject("signal_llltinf_sampParam",hd);
  if(hd!= 0){
    hd = (TH1*) hd->Clone();
    p->addHistLine(hd,"noCut");
    for(unsigned int iC = 0; iC < cuts.size(); ++iC){          
      TH1* h = 0;
      f->GetObject("signal_"+cuts[iC]+"_sampParam",h);      
      if(h==0){
        std::cout << "signal_"+cuts[iC]+"_sampParam" <<std::endl;
        continue;  
      } 
      p->addHistLine(h,cuts[iC]);
    }
    p->setMinMax(0.5,1);
    p->drawRatio(0,"stack",false,false,"signal_sampParam");
  } else {
    std::cout << "signal_sampParam" <<std::endl;
  }
  
  //Now BKG
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * p = new Plotter();
    TH1* hd = 0;
    f->GetObject(sigs[iS]+"_simpleHH",hd);      
    if(hd==0) continue;
    hd = (TH1*) hd->Clone();
    p->addHistLine(hd,"noC");    
      for(unsigned int iC = 0; iC < cuts.size(); ++iC){            
            TString v = vars[0];
            if(iC > 0) v = vars[1];
            if(iC >1) v = vars[2];
    TH1* h = 0;
    f->GetObject(sigs[iS]+"_"+cuts[iC]+"_"+v,h);      
      if(h==0) continue;      
      p->addHistLine(h,cuts[iC]);
    }

  p->setMinMax(0,1.5);
  p->rebin(10);
  p->drawRatio(0,"stack",false,false,sigs[iS]+"_cutComps");  
}
 
}


///BKG Shape
{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets"};
  std::vector<TString> vars = {"chi2HH","likeliHH"};
  std::vector<TString> dens = {"simpleHH","simpleHH"};

  std::vector<TString> cuts = {"chiltinf","chilt20","chilt9","chilt5","llltinf","lllt37","lllt32","lllt29p5"};


    for(unsigned int iV = 0; iV < vars.size(); ++iV){  
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      Plotter * p = new Plotter();
      TH1* hd = 0;
      f->GetObject(sigs[iS]+"_"+dens[iV],hd);
      
      if(hd==0) continue;
      hd = (TH1*) hd->Clone();
      hd->Scale(1./hd->Integral(hd->FindFixBin(1000),hd->FindFixBin(4000)));
              p->addHistLine(hd,"noC");
  for(unsigned int iC = 0; iC < cuts.size(); ++iC){
    
    
    
          
    TH1* h = 0;
    f->GetObject(sigs[iS]+"_"+cuts[iC]+"_"+vars[iV],h);      
      if(h==0) continue;      
            h->Scale(1./h->Integral(hd->FindFixBin(1000),hd->FindFixBin(4000)));
      p->addHistLine(h,cuts[iC]);
    }
    p->setBotMinMax(0,2);
    // p->normalize();
    p->drawSplitRatio(0,"stack",false,false,sigs[iS]+"_"+vars[iV]);

    // p->draw(false,vars[iV]);
  }
  
}
  
  
}


{
  TFile * f = new TFile("hSolTrees_test.root");
  std::vector<TString> sigs = {"ttbar","wjets"};
  std::vector<TString> vars = {"chi2HH","likeliHH"};
  std::vector<TString> dens = {"simpleHH","simpleHH"};

  std::vector<TString> cuts1 = {"chiltinf","chilt20","chilt9","chilt5"};
  std::vector<TString> cuts2 = {"llltinf","lllt37","lllt32","lllt29p5"};


    for(unsigned int iC = 0; iC < cuts1.size(); ++iC){  
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      Plotter * p = new Plotter();
      TH1* hd = 0;
      f->GetObject(sigs[iS]+"_"+dens[0],hd);
      
      if(hd==0){
        std::cout << sigs[iS]+"_"+dens[iC] << std::endl;
        continue;
      } 
      hd = (TH1*) hd->Clone();
              // p->addHistLine(hd,"noC");
              
              TH1* h1 = 0;
              f->GetObject(sigs[iS]+"_"+cuts1[iC]+"_"+vars[0],h1);      
                if(h1==0) {
                  std::cout << sigs[iS]+"_"+cuts1[iC]+"_"+vars[0]<<std::endl;
                  continue;      
                }
                p->addHistLine(h1,cuts1[iC]+"_"+vars[0]);        
                
                TH1* h2 = 0;
                  f->GetObject(sigs[iS]+"_"+cuts2[iC]+"_"+vars[1],h2);      
                  if(h2==0){
                                      std::cout << sigs[iS]+"_"+cuts2[iC]+"_"+vars[1]<<std::endl;
                                      continue;      
                  } 
                  p->addHistLine(h2,cuts2[iC]+"_"+vars[1]);     
                  
                  
                  p->setBotMinMax(0,2);
                  p->drawSplitRatio(0,"stack",false,false,sigs[iS]+"_"+cuts1[iC]);   
              
              

  }
  
}
  
  
}



///////. BACKGROUND LIKLI


///Figure out binning (fine)
{
  TFile * f = new TFile("hSolTrees_bkgTempStudy_bkg.root");
  std::vector<TString> sigs = {"ttbarPW","ttbar","wjets"};
  std::vector<TString> shells = {"","_low","_test","_high"};

  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","approxW","approxH2"};
  
  std::vector<TObject*> objs;
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;
      for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[0] +"_"+vars[iV];
      f->GetObject(hname,h);      
      if(h==0) continue;
      nH++;    
      p->addHistLine(h,sigs[iS]);
  }
  if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    p->draw(false,vars[iV]);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  
  for(unsigned int iV = 0; iV < vars.size();++iV){
    Plotter * p = new Plotter();
    int nH = 0;    
      for(unsigned int iS = 0; iS < shells.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[0]+shells[iS] +"_"+vars[iV];
      f->GetObject(hname,h); 
      if(h==0) continue;
      nH++;         
      p->addHistLine(h,sigs[0]+shells[iS]);
  }
    if(!nH) continue;
    p->normalize();
    p->setBotMinMax(0,2);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]+"_shells");
    p->draw(false,vars[iV]+"_shells");
    // objs.push_back(c);
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}
///Figure out binning (final binning)
{
  TFile * f = new TFile("hSolTrees_bkgTempStudy_bkg.root");
  std::vector<TString> sigs = {"ttbarPW"};
  std::vector<TString> shells = {"_low","_high"};
    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","approxW","approxH"};
  std::vector<int> vBinV      = {30  ,50  ,20,50};
  std::vector<double> vBinLEV = {-150,-0.6,60,100};
  std::vector<double> vBinHEV = {150 ,0.4 ,100,600};
  
  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
            Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = sigs[iS]+shells[iSh] +"_"+vars[iV];                  
      f->GetObject(hname,h);      
      if(!h) continue;
      
      int nNB; double nL; double nH;

        nNB = vBinV[iV];
        nL = vBinLEV[iV];
        nH = vBinHEV[iV];
  
  
      
      if(nL==nH) continue;
      TH1 * hR = new TH1F(TString(h->GetName()) +"_reb",h->GetTitle(),nNB,nL,nH);
      hR->SetBinContent(0,h->GetBinContent(0));
      hR->SetBinContent(nNB+1,h->GetBinContent(h->GetNbinsX()+1));      
      for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB ){
        int hRB = hR->FindFixBin(h->GetBinCenter(iB));
        
        hR->SetBinContent(hRB, hR->GetBinContent(hRB) + h->GetBinContent(iB) );
      }
      
      
      p->addHistLine(hR,sigs[iS]);
  
}
    p->normalize();
    p->setBotMinMax(0,2);
    TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]);
    objs.push_back(c);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  }
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}

///ROC CURVES
{
  TString year = "2018";
  TFile * f = new TFile("hSolTrees_test_"+year+".root");
  std::vector<TString> sigs = {"radion_m1000","radion_m2000","radion_m3000"};
  std::vector<TString> bkgCs = {"m1000","m2000","m2000"};
  
  
  
  
  // std::vector<TString> sigs = {"radion_m1000"};
  // std::vector<TString> bkgCs = {"m1000"};
  std::vector<TString> sigNs = {""};
  // std::vector<TString> varNss = {"signal L","bkg. L", "L ratio"};
  std::vector<TString> varNss = {"#it{D} (new)","#it{m}_{D} (old)"};
  
  // std::vector<TString> sigs = {"radion_m1000_lllt32"};
  // std::vector<TString> bkgCs = {"lllt32"};
  std::vector<TString> bkgs = {"ttbar","wjets","qcd"};
  
  // std::vector<TString> vars = {"chi2","likeli","Blikeli","SoBlikeli"};
  // std::vector<TString> vars = {"likeli","Blikeli","SoBlikeli"};
    // std::vector<TString> vars = {"likeliBHWW","likeliHWW"};
  
  
    // std::vector<TString> vars = {"Blikeli","Blikeli_nAlt"};
    // std::vector<TString> vars = {"SoBlikeli","SoBNlikeli"};
    std::vector<TString> vars = {"likeli","md"};
  
  

  int cutGT = 10;
  vector<TObject*> objs;
  for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
          for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * p = new Plotter();
    for(unsigned int iV = 0; iV < vars.size(); ++iV){

        TH1* hs = 0;
        f->GetObject(sigs[iS]+"_"+vars[iV],hs);     
        if(hs==0)continue;
        hs=(TH1*)hs->Clone();
        TH1* hb = 0;
        f->GetObject(bkgs[iB]+"_"+bkgCs[iS]+"_"+vars[iV],hb);     
        if(hb==0)continue;
        hb=(TH1*)hb->Clone();
        hb->SetName(bkgs[iB]+"_"+bkgCs[iS]+"_"+vars[iV] + "_"+iS );
        PlotTools::toUnderflow(hs);
        PlotTools::toOverflow(hs);
        PlotTools::toUnderflow(hb);
        PlotTools::toOverflow(hb);
        auto * roc = PlotTools::getRocCurve(hs,hb,iV==cutGT, "signal eff","bkg. eff");
        // p->addGraph(roc,sigs[iS]+"_"+vars[iV]);
        // p->addGraph(roc,sigs[iS]+": "+varNss[iV]);
                p->addGraph(roc,varNss[iV]);
                        // p->addGraph(roc,vars[iV]); 
      }                  
      // p->setXTitle("signal efficiency");
      p->setXTitle(sigs[iS]+" efficiency");
      // p->setYTitle("bkg. efficiency");
      p->setYTitle(bkgs[iB] +" efficiency");
      auto c = p->draw(false,sigs[iS]+"_"+bkgs[iB])      ;
      objs.push_back(c);
    }

  }
  auto c= Drawing::drawAll(objs,"rocCurves_"+year);
}


////QCD TESTING
{
  TFile * f = new TFile("hSolTrees_qcdTempStudy.root");
  std::vector<TString> sigs = {"qcd","radion_m1000","radion_m2000"};
  std::vector<TString> vars = {"extraMetParRelhwwMag","approxH","wPTRelhwwMag","extraMetParRelhwwMag2","approxH3"};
    // std::vector<TString> vars = {"deltaPhi","deltaPhi2","deltaTheta","deltaTheta2"};
    // std::vector<TString> vars = {"approxW","approxH","ptDR","approxW2","approxH2","ptDR2","DR","DR2"};

  std::vector<TString> cuts = {""};


    for(unsigned int iC = 0; iC < cuts.size(); ++iC){  
          for(unsigned int iV = 0; iV < vars.size(); ++iV){  
                  Plotter * p = new Plotter();
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){

              
              TH1* h1 = 0;
              f->GetObject(sigs[iS]+"_"+cuts[iC]+vars[iV],h1);      
                if(h1==0) {
                  std::cout << sigs[iS]+"_"+cuts[iC]+vars[iV]<<std::endl;
                  continue;      
                }
                p->addHistLine(h1,sigs[iS]);                                                    

  }
  p->rebin(4);
  p->normalize();
  p->draw(false,vars[iV]+"_"+cuts[iC]);   
  
}
  
  
}
}


//////////////. GET CUTS. //////////////. 
//First cut for high mass
  {
    TString year = "2017";
    TFile * f = new TFile("hSolTrees_getCuts_"+year+".root");
    std::vector<std::string> sigs = {"m800","m1000","m2000","m2500","m3000"};
  
          // std::vector<std::string> cuts = {"tau0p75_LMT_emu","tau0p75_MT_emu","tau0p45_T_emu","tau0p45_MT_emu","tau0p45to0p75_MT_emu","tau0p45to0p75_LMT_emu"};
                    // std::vector<std::string> cuts = {"tau0p75_LMT_emu","tau0p75_MT_emu","tau0p45_T_emu","tau0p45_MT_emu","tau0p35_T_emu","tau0p35_MT_emu"};
                          std::vector<std::string> cuts = {"tau0p45_T_emu","tau0p45_MT_emu","tau0p40_T_emu","tau0p40_MT_emu","tau0p55_T_emu","tau0p55_MT_emu"};
  
    // TString sigN ="radion";
    TString sigN =   "graviton";
    TString var = "likeli";
    vector<TObject*> objs_s;
    vector<TObject*> objs_b;
    vector<TObject*> objs_sob;
  
    for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      Plotter * pr = new Plotter();
      Plotter * ps = new Plotter();
      Plotter * pb = new Plotter();
    for(unsigned int iC = 0; iC < cuts.size(); ++iC){

    
      
        TH1* hs = 0;
        f->GetObject(sigN + "_"+sigs[iS]+"_"+cuts[iC]+"_mIncl_"+var,hs);
        TH1* hb = 0;
        f->GetObject("bkg_"+cuts[iC]+"_"+sigs[iS]+"_"+var,hb);      
        if(hs==0||hb==0) continue;      
        hs = PlotTools::getIntegral(hs,false,false);
        hb = PlotTools::getIntegral(hb,false,false);
        TH1* hr = (TH1*)hs->Clone();
        for(unsigned int iX = 1; iX <= hs->GetNbinsX(); ++iX){
          hr->SetBinContent(iX,hs->GetBinContent(iX)/(std::sqrt(3.+hb->GetBinContent(iX)) ));
        }
        ps->addHistLine(hs,cuts[iC]);
        pb->addHistLine(hb,cuts[iC]);
        pr->addHistLine(hr,cuts[iC]);        
      }
      pr->setXTitle(sigs[iS] +" likeli" );
      auto sob = pr->draw(false,"sosqrt3pb_" + sigs[iS]);
      objs_sob.push_back(sob);
      auto s = ps->draw(false,"sig_"+ sigs[iS]);
      objs_s.push_back(s);
      auto b = pb->draw(false,"bkg_"+ sigs[iS]);
      objs_b.push_back(b);
    
    }
  
    Drawing::drawAll(objs_sob,"sosqrt3pb_"+year);
    Drawing::drawAll(objs_s,"s_"+year);
    Drawing::drawAll(objs_b,"b_"+year);
    
  
  
  }

//////////////. Evalute impact on individual SRs. //////////////. 

{
  TFile * f = new TFile("hSolTrees_getCuts.root");
  
  std::vector<std::string> cuts = {"oHPMD_LMT_e","oHPMD_LMT_mu","oLPMD_LMT_e","oLPMD_LMT_mu","nHP_LMT_e","nHP_LMT_mu","nLP_LMT_e","nLP_LMT_mu"};
  // std::vector<std::string> cuts = {"nHP_L_e","nHP_M_e","nHP_T_e","nLP_L_e","nLP_M_e","nLP_T_e","nHP_L_mu","nHP_M_mu","nHP_T_mu","nLP_L_mu","nLP_M_mu","nLP_T_mu"};
std::vector<std::string> bkgs = {"ttbar","wjets","qcd"};
  TString var = "mhh";
  


  for(unsigned int iC = 0; iC < cuts.size(); ++iC){
    Plotter * p = new Plotter();
    
      for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
    
      
      TH1* hs = 0;
      f->GetObject(bkgs[iB] +"_"+cuts[iC]+"_"+var,hs);  
      if(hs==0) continue;      
      p->addHistLine(hs,bkgs[iB]);
      
    }
    p->setBotMinMax(0,1);
    p->drawSplitRatio(0,"stack",false,false,cuts[iC]);
    // ps->draw(false,"sig_"+ sigs[iS]);
    // pb->draw(false,"bkg_"+ sigs[iS]);
    
  }
  
}
{
  TFile * f = new TFile("hSolTrees_getCuts.root");
  

std::vector<std::string> bkgs = {"ttbar","wjets","qcd"};
std::vector<std::string> btags = {"L","M","T"};
std::vector<std::string> leps = {"e","mu"};
std::vector<std::string> regs = {"HP","LP"};

  TString var = "mhh";
  

for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
  for(unsigned int iT = 0; iT < btags.size(); ++iT){          
    Plotter * p = new Plotter();    
    
      for(unsigned int iL = 0; iL < leps.size(); ++iL){          
          for(unsigned int iR = 0; iR < regs.size(); ++iR){          
      
      TH1* ho = 0;
      f->GetObject(bkgs[iB] +"_o"+regs[iR]+"MD_"+btags[iT]+"_"+leps[iL]+"_"+var,ho);  
      TH1* hn = 0;
      f->GetObject(bkgs[iB] +"_n"+regs[iR]+"_"+btags[iT]+"_"+leps[iL]+"_"+var,hn);  
      if(ho==0) continue;      
            if(hn==0) continue;  
            
            hn = PlotTools::rebin(hn,4);
            ho = PlotTools::rebin(ho,4); 
            
      hn->Divide(ho);    
      p->addHistLine(hn,regs[iR]+"_"+btags[iT]+"_"+leps[iL]);
      
    }
  }
    p->setMinMax(0.01,10);
    auto *c = p->draw(false,bkgs[iB]+"_"+btags[iT]);
    c->SetLogy();
    c->Update();
    // ps->draw(false,"sig_"+ sigs[iS]);
    // pb->draw(false,"bkg_"+ sigs[iS]);
    
  }
}
  
}


////////////// PLOTS for the Meeting ///////////////////
{
  TFile * f = new TFile("hhSol_templates.root","read");
  std::vector<TString> hists = {
  "signal_low_vqq_hWW","signal_low_osqq_hWW",
  "signal_low_vqq_Wlnu","signal_low_osqq_Wlnu",
  "signal_low_vqq_extraMetPerp","signal_low_osqq_extraMetPerp",
  "signal_low_vqq_extraMetParRelhwwMag","signal_low_osqq_extraMetParRelhwwMag",
  "signal_low_vqq_wqqPTRes","signal_low_osqq_wqqPTRes",
  "signal_low_vqq_qqSDMassCoarse","signal_low_osqq_qqSDMassCoarse",

};
  
  std::vector<TString> xtitles = {
    "#it{m}_{#it{l}#nu qq'}","#it{m}_{#it{l}#nu qq'}",
    "#it{m}_{#it{l}#nu}","#it{m}_{#it{l}#nu}",
    "extra perpendicular #it{p}^{miss}_{T}","extra perpendicular #it{p}^{miss}_{T}",
    "relative extra parallel #it{p}^{miss}_{T}","relative extra parallel #it{p}^{miss}_{T}",
    "qq' jet #it{p}_{T} response","qq' jet #it{p}_{T} response",
    "#it{m}_{jet} [GeV]","#it{m}_{jet} [GeV]"
  };
    
    
  TString ytitle = "Prob. density";
  TString ytitle2 = "Probability";
  int yt2 = 10;
  
  for(unsigned int iH = 0; iH < hists.size(); ++iH){
    TH1* h = 0;
    f->GetObject(hists[iH],h);
    if(!h) continue;
    h->SetXTitle(xtitles[iH]);    
    if(iH>= yt2){
      h->GetYaxis()->SetRangeUser(0,1);
      h->SetYTitle(ytitle2);  
    } else {
      h->SetYTitle(ytitle);        
    }
    h->GetXaxis()->SetTitleOffset(1.0); 
    h->GetYaxis()->SetTitleOffset(1.45); 
    auto * c = new TCanvas();
    h->Draw("hist");
    c->Print("plots/"+hists[iH] +".pdf");
  }
}

{
  Plotter * p = new Plotter();
  p->addStackHist(radion_m1000_onshellSD_qqSDMass,"Virtual W#rightarrow#it{l}#nu");
  p->addStackHist(radion_m1000_virtualSD_qqSDMass,"Virtual W#rightarrowqq'");
  p->rebin(4);
  p->normalize();
  p->setXTitle("#it{m}_{jet} [GeV]");
  p->setYTitle("arbitrary units");
  p->draw();
  
}

{
  std::vector<double> masses = {800 ,900 ,1000,2000,2500,3000};
  // std::vector<double> oldV = {0.135254, 0.0766602, 0.0522461, 0.0112305, 0.00927734, 0.00732422};
  // std::vector<double> newV = { 0.0922852, 0.0532227, 0.0366211, 0.00927734, 0.0078125, 0.00732422}
  
  std::vector<double> oldV = {0.059082, 0.0356445, 0.0258789, 0.00634766, 0.00537109, 0.00439453};
  std::vector<double> newV = {0.0512695, 0.0307617, 0.0219727, 0.00634766, 0.00439453, 0.00439453};
  ;
  
  TGraph * go = new TGraph();
  TGraph * gn = new TGraph();
  TGraph * gr = new TGraph();
  for(unsigned int iP = 0; iP < masses.size(); ++iP){
    go->SetPoint(iP,masses[iP],oldV[iP]*1000);
    gn->SetPoint(iP,masses[iP],newV[iP]*1000);
    gr->SetPoint(iP,masses[iP],newV[iP]/oldV[iP]);
  }
  Plotter*p = new Plotter();
  p->addGraph(go,"Old");
  p->addGraph(gn,"New");
  p->setXTitle("#it{m}_{X} [GeV]");
  p->setYTitle("#sigma#it{B}(X#rightarrowHH) [fb]");
  p->setMinMax(1,300.);
  auto *c = p->draw(true,"plots/12l_diff.pdf");
  c->SetLogy();
  c->Update();
  
  
  Plotter*p2 = new Plotter();
  p2->addGraph(gr,"New/Old");
  p2->setXTitle("#it{m}_{X} [GeV]");
  p2->setYTitle("Relative limit (new/old)");
  // p->setMinMax(1,300.);
  p2->draw(true,"plots/12lrat.pdf");

  
  
  
}

// Compare different years for signal
{
  vector<TString> fileNs = {"hSolTrees_temp_2016_bulkgrav.root","hSolTrees_temp_2017_bulkgrav.root","hSolTrees_temp_2018_bulkgrav.root"};
  vector<TString> fileTs = {"2016","2017","2018"};
  std::vector<TString> sigs = {"signal_low","signal_test","signal_high"};
  std::vector<TString> shells = {"_vqq","_osqq"};    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","Wlnu","hWW","qqSDMassCoarse"};

  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
                for(unsigned int iS = 0; iS < sigs.size(); ++iS){
            Plotter * p = new Plotter();
            
            for(unsigned int iY = 0; iY< fileTs.size(); ++iY){
              auto * f = new TFile(fileNs[iY]);
              if(f == 0) continue;
              TH1 *h = 0;
              auto hname = sigs[iS]+shells[iSh] +"_"+vars[iV];
              f->GetObject(hname,h);
              if(!h) continue;
              h->SetDirectory(0);
              h->SetName(hname+"_"+fileTs[iY]);
              p->addHistLine(h,fileTs[iY]);
              // f->Close();
              // delete f;
            }       
            p->normalize();
            p->setBotMinMax(0,2);
            TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]+"_"+sigs[iS]);
            objs.push_back(c);
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    // objs.push_back(c);
  }
  }
}
  
  // auto c= Drawing::drawAll(objs,"templatesUnbinned");
  // c->Print("templatesUnbinned.pdf");
  
}

//Look at all templates
{
  TString year = "2018";
  
  TFile * f = new TFile("hhSol_templates_"+year+".root");
  std::vector<TString> sigs = {"low","high"};
  std::vector<TString> shells = {"vqq","osqq"};
    
  std::vector<TString> vars = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","Wlnu","hWW","qqSDMassCoarse"};
  std::vector<TObject*> objs;

  for(unsigned int iV = 0; iV < vars.size();++iV){    
      for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
            Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      TH1 *h = 0;
      auto hname = "signal_" + sigs[iS]+"_"+shells[iSh] +"_"+vars[iV];                  
      f->GetObject(hname,h);      
      if(!h) continue;
      for(unsigned int iB = 0; iB <= h->GetNbinsX() + 1; ++iB){
        if(h->GetBinContent(iB) <= 0) cout << sigs[iS]+"_"+shells[iSh] +"_"+vars[iV] <<" -> "<< iB <<" "<< 
          h->GetBinCenter(iB) <<endl;
      }
      p->addHistLine(h,shells[iSh]+": "+sigs[iS]);
  
}
    TCanvas * c = p->draw(false,vars[iV]+"_"+shells[iSh]);
    // c->SetLogy();
    // c->Update();
    // TCanvas * c = p->drawSplitRatio(0,"stack",false,false,vars[iV]);
    objs.push_back(c);
  }
  }
  
  auto c= Drawing::drawAll(objs,"templates_"+year);
  // c->Print("templatesUnbinned.pdf");
  
}

// Compare different years of likelihoods
{
  vector<TString> fileNs = {"hSolTrees_getCuts_2016.root","hSolTrees_getCuts_2017.root","hSolTrees_getCuts_2018.root"};
  vector<TString> fileTs = {"2016","2017","2018"};
  std::vector<TString> sigs = {"radion_m1000","radion_m2000","radion_m3000"};
  std::vector<std::string> cuts = {"tau0p75_LMT_emu","tau0p75_MT_emu","tau0p45_T_emu","tau0p45_MT_emu"};
  TString var = "likeli";
  std::vector<TObject*> objs;

      for(unsigned int iC = 0; iC < cuts.size(); ++iC){
                for(unsigned int iS = 0; iS < sigs.size(); ++iS){
            Plotter * p = new Plotter();
            
            for(unsigned int iY = 0; iY< fileTs.size(); ++iY){
              auto * f = new TFile(fileNs[iY]);
              if(f == 0) continue;
              TH1 *h = 0;
              auto hname = sigs[iS]+"_"+cuts[iC]+"_mIncl_"+var;
              f->GetObject(hname,h);
              if(!h) continue;
              h->SetDirectory(0);
              h->SetName(hname+"_"+fileTs[iY]);
              PlotTools::toUnderflow(h);
              PlotTools::toOverflow(h);
              h = PlotTools::getIntegral(h,false,true);
              p->addHistLine(h,fileTs[iY]);
              p->setYTitle(cuts[iC]);
              p->setXTitle(sigs[iS] +": likeli");
            }       
            TCanvas * c = p->draw(false,var+"_"+cuts[iC]+"_"+sigs[iS]);
            objs.push_back(c);
  }
  }
  Drawing::drawAll(objs,"efficiencies");
  
  
  
}


///////////. BTAGGGING SCRIPTS!!!!

{
  TString year = "2017";
  TFile * f = new TFile("hSolTrees_hbbTaggingCutHistograms_"+year+".root");
  std::vector<std::string> sigs = {"m800","m1000","m2000","m2500","m3000"};

                        std::vector<std::string> cuts = {"LHP_emu","HP_emu","LP_emu","HP_e","HP_mu"};

  // TString sigN ="radion";
  TString sigN =   "graviton";
  // TString var ="dCSVCat";
  TString var = "ak8Tag";
  vector<TObject*> objs_s;
  vector<TObject*> objs_b;
  vector<TObject*> objs_sob;

  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * pr = new Plotter();
    Plotter * ps = new Plotter();
    Plotter * pb = new Plotter();
  for(unsigned int iC = 0; iC < cuts.size(); ++iC){

  
    
      TH1* hs = 0;
      f->GetObject(sigN + "_"+sigs[iS]+"_"+cuts[iC]+"_mIncl_"+var,hs);
      TH1* hb = 0;
      f->GetObject("bkg_"+cuts[iC]+"_"+sigs[iS]+"_"+var,hb);      
      if(hs==0||hb==0) continue;      
      hs = PlotTools::getIntegral(hs,true,false);
      hb = PlotTools::getIntegral(hb,true,false);
      TH1* hr = (TH1*)hs->Clone();
      for(unsigned int iX = 1; iX <= hs->GetNbinsX(); ++iX){
        hr->SetBinContent(iX,hs->GetBinContent(iX)/(std::sqrt(3.+hb->GetBinContent(iX)) ));
      }
      ps->addHistLine(hs,cuts[iC]);
      pb->addHistLine(hb,cuts[iC]);
      pr->addHistLine(hr,cuts[iC]);        
    }
    pr->setXTitle(sigs[iS] +" likeli" );
    auto sob = pr->draw(false,"sosqrt3pb_" + sigs[iS]);
    objs_sob.push_back(sob);
    ps->normalize();
    auto s = ps->draw(false,"sig_"+ sigs[iS]);
    objs_s.push_back(s);
    pb->normalize();
    auto b = pb->draw(false,"bkg_"+ sigs[iS]);
    objs_b.push_back(b);
  
  }

  Drawing::drawAll(objs_sob,"sosqrt3pb_"+year);
  Drawing::drawAll(objs_s,"s_"+year);
  Drawing::drawAll(objs_b,"b_"+year);


}

//compare years

{

  std::vector<std::string> sigs = {"m800","m1000","m2000","m2500","m3000"};

  std::vector<std::string> cuts = {"LHP_emu"};
  std::vector<TString> years = {"2016","2017","2018"};

  // TString sigN ="radion";
  TString sigN =   "graviton";
  // TString var ="dCSVCat";
  TString var = "ak8Tag";
  vector<TObject*> objs_s;
  vector<TObject*> objs_b;

  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
    Plotter * pr = new Plotter();
    Plotter * ps = new Plotter();
    Plotter * pb = new Plotter();
  for(unsigned int iY = 0; iY < years.size(); ++iY){

        TFile * f = new TFile("hSolTrees_hbbTaggingCutHistograms_"+years[iY]+".root");
        
    
      TH1* hs = 0;
      f->GetObject(sigN + "_"+sigs[iS]+"_"+cuts[0]+"_mIncl_"+var,hs);
      TH1* hb = 0;
      f->GetObject("bkg_"+cuts[0]+"_"+sigs[iS]+"_"+var,hb);      
      if(hs==0||hb==0) continue;      
      
      
      hs = (TH1*)hs->Clone();
      hs->SetDirectory(0);
      hs->SetName(sigN + "_"+sigs[iS]+"_"+cuts[0]+"_mIncl_"+var +"_"+ years[iY]);
      
      hb = (TH1*)hb->Clone();
      hb->SetDirectory(0);
      hb->SetName("bkg_"+cuts[0]+"_"+sigs[iS]+"_"+var +"_"+ years[iY]);
      
      hs = PlotTools::getIntegral(hs,true,true);
      hb = PlotTools::getIntegral(hb,true,true);
      ps->addHistLine(hs,sigs[iS] +": "+years[iY]);
      pb->addHistLine(hb,sigs[iS] +": "+years[iY]);
    }
    pr->setXTitle(sigs[iS] +" likeli" );
    auto s = ps->draw(false,"sig_"+ sigs[iS]);
    objs_s.push_back(s);
    auto b = pb->draw(false,"bkg_"+ sigs[iS]);
    objs_b.push_back(b);
  
  }

  Drawing::drawAll(objs_s,"s");
  Drawing::drawAll(objs_b,"b");


  }

