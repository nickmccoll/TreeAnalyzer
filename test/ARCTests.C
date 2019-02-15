nohup combine -m 2500 -M AsymptoticLimits -v 2 --run expected  -n "Test_L_LP" datacard_std_e_L_LP_full_13TeV.txt &
nohup combine -m 1000 -M AsymptoticLimits -v 2 --run expected  -n "Test_L_LP" datacard_std_e_L_LP_full_13TeV.txt &
nohup combine -m 2500 -M AsymptoticLimits -v 2 --run expected  -n "Test_L_HP" datacard_std_e_L_HP_full_13TeV.txt &
nohup combine -m 1000 -M AsymptoticLimits -v 2 --run expected  -n "Test_L_HP" datacard_std_e_L_HP_full_13TeV.txt &


nohup combine -m 2500 -M AsymptoticLimits -v 2 --run expected  -n "Test_M_LP" datacard_std_e_M_LP_full_13TeV.txt &
nohup combine -m 1000 -M AsymptoticLimits -v 2 --run expected  -n "Test_M_LP" datacard_std_e_M_LP_full_13TeV.txt &
nohup combine -m 2500 -M AsymptoticLimits -v 2 --run expected  -n "Test_M_HP" datacard_std_e_M_HP_full_13TeV.txt &
nohup combine -m 1000 -M AsymptoticLimits -v 2 --run expected  -n "Test_M_HP" datacard_std_e_M_HP_full_13TeV.txt &


  rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/testATLASSel.C+("../7_13_18_eleSF_andMuTrig/eleRecoFollowUp/testTree_1000.root",2,1,"atlasTS_rad_m1000.root")' &
  rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/testATLASSel.C+("../7_13_18_eleSF_andMuTrig/eleRecoFollowUp/testTree_2000.root",2,1,"atlasTS_rad_m2000.root")' &
  rr -b -q '~/GitRepositories/TreeAnalyzer/search_design/testATLASSel.C+("../7_13_18_eleSF_andMuTrig/eleRecoFollowUp/testTree_4000.root",2,1,"atlasTS_rad_m4000.root")' &


{

    // std::vector<TString> files = {"higgsCombineTest_addTTHFS711.AsymptoticLimits.mH1000.root","higgsCombineTest_addTTHFS711.AsymptoticLimits.mH2500.root",
    //                               "higgsCombineTest_times2Syst.AsymptoticLimits.mH1000.root","higgsCombineTest_times2Syst.AsymptoticLimits.mH2500.root",
    //                               "higgsCombineTest_L_LP.AsymptoticLimits.mH1000.root","higgsCombineTest_L_LP.AsymptoticLimits.mH2500.root",
    //                               "higgsCombineTest_L_HP.AsymptoticLimits.mH1000.root","higgsCombineTest_L_HP.AsymptoticLimits.mH2500.root",
    //                               "higgsCombineTest_M_LP.AsymptoticLimits.mH1000.root","higgsCombineTest_M_LP.AsymptoticLimits.mH2500.root",
    //                               "higgsCombineTest_M_HP.AsymptoticLimits.mH1000.root","higgsCombineTest_M_HP.AsymptoticLimits.mH2500.root"};
    // std::vector<TString> names = {"std","std",
    //                               "x2Syst","x2Syst",
    //                               "L_LP","L_LP",
    //                               "L_HP","L_HP",
    //                               "M_LP","M_LP",
    //                               "M_HP","M_HP"};
  
  
  std::vector<TString> files = {"higgsCombineTest_mu_T_HP.Significance.mH1300.root","higgsCombineTest_mu_T_HP.Significance.mH2300.root",
"higgsCombineTest_mu_T_LP.Significance.mH1300.root","higgsCombineTest_mu_T_LP.Significance.mH2300.root",
"higgsCombineTest_mu_M_HP.Significance.mH1300.root","higgsCombineTest_mu_M_HP.Significance.mH2300.root",
"higgsCombineTest_mu_M_LP.Significance.mH1300.root","higgsCombineTest_mu_M_LP.Significance.mH2300.root",
"higgsCombineTest_mu_L_HP.Significance.mH1300.root","higgsCombineTest_mu_L_HP.Significance.mH2300.root",
"higgsCombineTest_mu_L_LP.Significance.mH1300.root","higgsCombineTest_mu_L_LP.Significance.mH2300.root",
"higgsCombineTest_e_T_HP.Significance.mH1300.root","higgsCombineTest_e_T_HP.Significance.mH2300.root",
"higgsCombineTest_e_T_LP.Significance.mH1300.root","higgsCombineTest_e_T_LP.Significance.mH2300.root",
"higgsCombineTest_e_M_HP.Significance.mH1300.root","higgsCombineTest_e_M_HP.Significance.mH2300.root",
"higgsCombineTest_e_M_LP.Significance.mH1300.root","higgsCombineTest_e_M_LP.Significance.mH2300.root",
"higgsCombineTest_e_L_HP.Significance.mH1300.root","higgsCombineTest_e_L_HP.Significance.mH2300.root",
"higgsCombineTest_e_L_LP.Significance.mH1300.root","higgsCombineTest_e_L_LP.Significance.mH2300.root"};
  std::vector<TString> names = {"mu_T_HP","mu_T_HP",
"mu_T_LP","mu_T_LP",
"mu_M_HP","mu_M_HP",
"mu_M_LP","mu_M_LP",
"mu_L_HP","mu_L_HP",
"mu_L_LP","mu_L_LP",
"e_T_HP","e_T_HP",
"e_T_LP","e_T_LP",
"e_M_HP","e_M_HP",
"e_M_LP","e_M_LP",
"e_L_HP","e_L_HP",
"e_L_LP","e_L_LP"};
  
  
  
  
  
  
    std::vector<float> vals;
    
    for(unsigned int iF = 0; iF < files.size();++iF){
    TFile * f = new TFile(files[iF],"read");
    TTree * t = 0;
    f->GetObject("limit",t);

    double limit;
    float quantileExpected;
    double mh;
    t->SetBranchAddress("limit",&limit);
    t->SetBranchAddress("quantileExpected",&quantileExpected);
    t->SetBranchAddress("mh",&mh);
    int eventNumber = 0;
    float val;
    while(t->GetEntry(eventNumber)){
        if(quantileExpected==-1) val = limit;
        ++eventNumber;
    }
    vals.push_back(val);
    f->Close();
    }
    
    std::cout <<std::endl<< vals.size()<<std::endl;
    for(unsigned int iF = 0; iF < files.size();++iF){
        if(iF%2)
        std::cout << vals[iF]<<"\n";
        else std::cout << names[iF]<<"\t"<<vals[iF]<<"\t";
    
    }
    std::cout << std::endl;


}


{


    Plotter * p = new Plotter();
    // p->addHistLine(m1000_qq_dr,"1 TeV #it{m}_{X}");
    // p->addHistLine(m2000_qq_dr,"2 TeV #it{m}_{X}");
    // p->addHistLine(m3000_qq_dr,"3 TeV #it{m}_{X}");
    
    auto add =[&](TH1* h, const char* name){
    p->addHistLine(h,name);
    int bin0 = h->FindFixBin(0.4);
    int bin1 = h->FindFixBin(0.8);
    float int0 = h->Integral(0,bin0 -1)/h->Integral(0,-1);
    float int1 = h->Integral(bin0,bin1-1)/h->Integral(0,-1);
    float int2 = h->Integral(bin1,-1)/h->Integral(0,-1);
    std::cout << name << "\t"<<int0<< "\t"<<int1<< "\t"<<int2<<"\t"<<int0+int1+int2<<"\n";
    } ;
    
    add(m1000_qq_dr,"1 TeV #it{m}_{X}");
    add(m2000_qq_dr,"2 TeV #it{m}_{X}");
    add(m3000_qq_dr,"3 TeV #it{m}_{X}");
    p->normalize();
    p->setXTitle("#DeltaR W#rightarrowqq quarks");
    p->setYTitle("a.u.");
    p->draw();
}

{
  ttbar_nevents->Add(singlet_nevents);
  ttbar_nevents->Add(ttX_nevents);
  ttbar_nevents->Add(wjets_njets);
  ttbar_nevents->Add(singlet_nevents);
  ttbar_nevents->Add(diboson_nevents);
  ttbar_nevents->Add(qcd_nevents);
  ttbar_nevents->Add(zjets_nevents);
  ttbar_nevents->Draw();
  
  
  ttbar_njets->Add(singlet_njets);
  ttbar_njets->Add(ttX_njets);
  ttbar_njets->Add(wjets_njets);
  ttbar_njets->Add(singlet_njets);
  ttbar_njets->Add(diboson_njets);
  ttbar_njets->Add(qcd_njets);
  ttbar_njets->Add(zjets_njets);
  ttbar_njets->Draw();
  
}





{
  Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>1500&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&hbbCSVCat>=4)"); h->Integral()
  Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>1500&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&hbbCSVCat>=4&&(met>50))"); h->Integral()
  Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>1500&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&hbbCSVCat>=4&&(hwwPT/hhMass>0.3))"); h->Integral()
    
    
    Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>700&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&wwDM<125.0&&wjjTau2o1<0.75&&process==3&&hbbCSVCat>=4)"); h->Integral()
    Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>700&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&wwDM<125.0&&wjjTau2o1<0.75&&process==3&&hbbCSVCat>=4&&met>50)"); h->Integral()
    Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>700&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&wwDM<125.0&&wjjTau2o1<0.75&&process==3&&hbbCSVCat>=4&&(hwwPT/hhMass>0.3))"); h->Integral()
      
      
    Events->Draw("1>>h","xsec*trig_N*pu_N*lep_N*btag_N*(passPre==1&&hhMass>1500&&hbbMass>30&&hbbMass<210&&nAK4Btags==0&&(hwwPT/hhMass>0.3)&&wjjTau2o1<0.75&&process==2&&hbbCSVCat>=4&&wwDM<125.0)"); h->Integral()      


}

std::vector<std::string> nms = {
"e_L_LP",
"e_L_HP",
"e_M_LP",
"e_M_HP",
"e_T_LP",
"e_T_HP",
"mu_L_LP",
"mu_L_HP",
"mu_M_LP",
"mu_M_HP",
"mu_T_LP",
"mu_T_HP"
};
  
  std::vector<float> counts;
  
  counts.push_back(prefit_e_L_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_e_L_HP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_e_M_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_e_M_HP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_e_T_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_e_T_HP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_L_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_L_HP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_M_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_M_HP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_T_LP_full__MJ_MR ->Integral(36,60,0,-1));
  counts.push_back(prefit_mu_T_HP_full__MJ_MR ->Integral(36,60,0,-1));
  
  for(unsigned int i = 0; i < counts.size(); ++i){
    std::cout << nms[i] <<"\t"<< counts[i] <<std::endl;
  }
  
e_L_LP	
e_L_HP	
e_M_LP	
e_M_HP	
e_T_LP	
e_T_HP	
mu_L_LP
mu_L_HP
mu_M_LP
mu_M_HP
mu_T_LP
mu_T_HP