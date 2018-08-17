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

# see electronEffPlots.C for scripts

# test on signal:

    nohup cmsRun run/searchRegionTreeMaker_cfg.py inputFiles="root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/Radion_hh_hVVhbb_inclusive_narrow_M2000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/B20383C9-295E-E711-A38B-008CFAFC0416.root" outputFile="testTree_2000.root" isCrab=False sample=signal  &
    nohup cmsRun run/searchRegionTreeMaker_cfg.py inputFiles="root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/Radion_hh_hVVhbb_inclusive_narrow_M4000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5C406E32-E45E-E711-804C-E0071B7A7830.root" outputFile="testTree_4000.root" isCrab=False sample=signal  &
    nohup cmsRun run/searchRegionTreeMaker_cfg.py inputFiles="root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/Radion_hh_hVVhbb_inclusive_narrow_M1000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/663F2191-DD5E-E711-9612-0CC47A4D767C.root" outputFile="testTree_1000.root" isCrab=False sample=signal &
    nohup cmsRun run/searchRegionTreeMaker_cfg.py inputFiles="root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/Radion_hh_hVVhbb_inclusive_narrow_M3000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8A905647-365F-E711-821E-7CD30ACE19BA.root" outputFile="testTree_3000.root" isCrab=False sample=signal  &


    rr -b -q 'understandElectronReco.C+("testTree_1000.root",1,"out_radion_hh_bbinc_m1000_0.root")' &
    rr -b -q 'understandElectronReco.C+("testTree_2000.root",1,"out_radion_hh_bbinc_m2000_0.root")' &
    rr -b -q 'understandElectronReco.C+("testTree_3000.root",1,"out_radion_hh_bbinc_m3000_0.root")' &
    rr -b -q 'understandElectronReco.C+("testTree_4000.root",1,"out_radion_hh_bbinc_m4000_0.root")' &



  
