Setup
https://twiki.cern.ch/twiki/bin/view/CMSPublic/ElectronTagAndProbe
Follow "Workflow Description in 80X"
  cmsrel CMSSW_8_0_27
  cd CMSSW_8_0_27/src
  cmsenv
  git cms-init
  git cms-merge-topic cms-egamma:EGM_gain_v1
  cd EgammaAnalysis/ElectronTools/data
  git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git
  cd $CMSSW_BASE/src
  git clone -b v2017.05.23_legacy80X_prelim https://github.com/nickmccoll/EgammaAnalysis-TnPTreeProducer EgammaAnalysis/TnPTreeProducer
  scram b -j8
    
    test:
    cmsRun TnPTreeProducer_cfg.py doEleID=True isMC=False maxEvents=5000 doPhoID=False doRECO=True

    submit:
    python ../CMSSW_8_0_27/src/EgammaAnalysis/TnPTreeProducer/scripts/crab/tnpCrabSubmitMiniAOD.py 



# Hadding trees

      hadd TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root `eos root://cmseos.fnal.gov find -f /eos/uscms/store/group/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/mc/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root `eos root://cmseos.fnal.gov find -f /eos/uscms/store/group/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/mc/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"`   &
      hadd TnPTree_SingleElectron_2016prompt_RunHv2.root `eos root://cmseos.fnal.gov find -f   /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016prompt_RunHv2/ | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunHv3.root `eos root://cmseos.fnal.gov find -f   /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016prompt_RunHv3/ | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunB.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunB/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunC.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunC/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunD.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunD/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunE.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunE/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunF.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunF/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &
      hadd TnPTree_SingleElectron_2016prompt_RunG.root   `eos root://cmseos.fnal.gov find -f  /eos/uscms/store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/data/SingleElectron/crab_2016rereco_RunG/    | egrep '.*root' | sed "s/.*\/eos\/uscms\(.*\)/root\:\/\/cmsxrootd.fnal.gov\/\1/"` &


      xrdcp TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root &
      xrdcp TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root &
      xrdcp TnPTree_SingleElectron_2016prompt_RunHv2.root root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunHv2.root &
      xrdcp TnPTree_SingleElectron_2016prompt_RunHv3.root root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunHv3.root &
      xrdcp TnPTree_SingleElectron_2016prompt_RunB.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunB.root   &
      xrdcp TnPTree_SingleElectron_2016prompt_RunC.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunC.root   &
      xrdcp TnPTree_SingleElectron_2016prompt_RunD.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunD.root   &
      xrdcp TnPTree_SingleElectron_2016prompt_RunE.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunE.root   &
      xrdcp TnPTree_SingleElectron_2016prompt_RunF.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunF.root   &
      xrdcp TnPTree_SingleElectron_2016prompt_RunG.root   root://cmseos.fnal.gov//store/user/lpchh/TagAndProbe/Electron_7_16_18/Moriond17_GainSwitch_newTnP_v5/combined/TnPTree_SingleElectron_2016prompt_RunG.root   &

# SKimming trees


      {
  
        std::vector<TString> inFiles = {
      "TnPTree_SingleElectron_2016prompt_RunG.root",
      "TnPTree_SingleElectron_2016prompt_RunB.root",
      "TnPTree_SingleElectron_2016prompt_RunE.root",
      "TnPTree_SingleElectron_2016prompt_RunHv2.root",
      "TnPTree_SingleElectron_2016prompt_RunD.root",
      "TnPTree_SingleElectron_2016prompt_RunF.root",
      "TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root",
      "TnPTree_SingleElectron_2016prompt_RunC.root",
      "TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
      "TnPTree_SingleElectron_2016prompt_RunHv3.root"
      };
      // TString outDir = "reco_skim";
      // TString treedir = "tnpEleReco";
      // TString sel = "tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && tag_Ele_nonTrigMVA >0.96 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 60";

      TString outDir = "id_skim";
      TString treedir = "tnpEleIDs";
      TString sel = "tag_Ele_pt > 27 && abs(tag_sc_eta) < 2.17 && el_q*tag_Ele_q < 0";

      for(const auto& iF : inFiles){
        TFile * fin = new TFile(iF,"read");
        TTree * tin = 0;
        fin->GetObject(treedir +"/fitter_tree",tin);
        if(tin == 0) return;
        TFile * fout = new TFile(outDir + "/"+iF ,"recreate");
        TDirectory * cdtof = fout->mkdir(treedir);
        cdtof->cd();
        TTree *tout = tin->CopyTree(sel);
        tout->Write();
        delete fin;
        delete fout;
      }
  
      }
    
# Now analyze the trees
#download this fork
https://github.com/michelif/nickmccoll/tree/egm_tnp_Moriond17_v2.0
---> need to link to libraries:
ln -s egm_tnp_analysis/libCpp .
  
  
  had to refit with special "doRefitFunc" on 8,21
  For --altSig --iBin 1, hadd to manually restrict "radFracFail[0.5,0.0,0.05]"]  
    
    
#examples
    https://indico.cern.ch/event/604907/contributions/2452909/attachments/1401631/2139468/Ele_IdSF_Moriond17_25Jan.pdf
    https://indico.cern.ch/event/604907/contributions/2452907/attachments/1401460/2139067/RecoSF_ApprovalMoriond17_25Jan2017.pdf
  
  //Eff for reco for signal in signaltriggerkine
//Compare efficiencies w/ new vs old method
  
  
{
TFile * nf = new TFile("results/reco_std_func/passingRECO/egammaEffi.txt_EGM2D.root","read");
TFile * of = new TFile("results/reco_std/passingRECO/egammaEffi.txt_EGM2D.root","read");
std::vector<TString> hnames = {"EGamma_SF2D","EGamma_EffData2D","EGamma_EffMC2D"};
std::vector<TString> ht = {"scale factor","data efficiency","MC efficiency"};

for(unsigned int iV = 0; iV < hnames.size();++iV){
  TH2*hn = 0;
  TH2*ho = 0;
  nf->GetObject(hnames[iV],hn);
  of->GetObject(hnames[iV],ho);
  
  auto hn1 = hn->ProjectionX(hnames[iV]+"_hnx");
  auto ho1 = ho->ProjectionX(hnames[iV]+"_hox");
  for(unsigned int iB = 1; iB <=hn1->GetNbinsX(); ++iB ) hn1->SetBinError(iB,0);
  
  Plotter * p = new Plotter();
  p->addHist(ho1,"method: template fit");
  p->addHist(hn1,"method: function fit");
  if(iV)p->setMinMax(0,1);
  else p->setMinMax(0.5,1.5);
  p->setBotMinMax(0.9,1.1);
  p->setYTitle(ht[iV]);
  p->setXTitle("supercluster #eta");
  p->setYTitleBot("template/function");  
  p->drawSplitRatio(1,"stack",false,false,hnames[iV]);
  
  
  
  
}

}

//Projection
{
TFile * nf = new TFile("results/recoISO/passingRECO/egammaEffi.txt_EGM2D.root","read");
std::vector<TString> hnames = {"EGamma_SF2D","EGamma_EffData2D","EGamma_EffMC2D"};
std::vector<TString> ht = {"scale factor","data efficiency","MC efficiency"};

for(unsigned int iV = 0; iV < hnames.size();++iV){
  TH2*hn = 0;
  nf->GetObject(hnames[iV],hn);
  Plotter * p = new Plotter();
  for(unsigned int iB = 6; iB <= hn->GetNbinsX(); ++iB){
    auto h1 = hn->ProjectionY(hnames[iV]+TString::Format("%u",iB),iB,iB);
    TString title = TString::Format("%.1f < #it{J}_{#it{A}} < %.1f",hn->GetXaxis()->GetBinLowEdge(iB),hn->GetXaxis()->GetBinLowEdge(iB)+hn->GetXaxis()->GetBinWidth(iB));
    p->addHist(h1,title);
    
  }  

  // if(iV)p->setMinMax(0,1);
  // else p->setMinMax(0.5,1.5);
  p->setYTitle(ht[iV]);
  p->setXTitle("#Delta#it{R}_{#it{A}}");
  p->draw(false,hnames[iV]);
  // p->drawSplitRatio(1,"stack",false,false,hnames[iV]);
  
  
  
  
}

}