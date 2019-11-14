{
  // TFile * f = new TFile("testHBBTagging_TT12LSig.root");
    TFile * f = new TFile("outPlots_HBBTag.root");
  // std::vector<TString> sigs = {"m800","m1000","m2000","m2500","m3000"};
    std::vector<TString> sigs = {"m1000"};
        std::vector<TString> bkgs = {"ttbar"};
        // std::vector<TString> cuts = {"emu_LHP_HbbIncl","emu_LHP_Hbb2SJ","emu_LHP_Hbb2SJTMass"};
                // std::vector<TString> cuts = {"emu_LHP_HbbTMass","emu_LHP_Hbb2SJ","emu_LHP_Hbb2SJTMass"};
                std::vector<TString> cuts = {"emu_LHP_HbbIncl"};
                
        // std::vector<TString> vars = {"bbt","mdZHbb","mdHbb","Hbb","sjMaxDeepCSV","sJBTagging"};
                std::vector<TString> vars = {"bbt","mdZHbb","Hbb","sJBTagging"};  
                std::vector<TString> varNs = {"H#rightarrow b#bar{b} tagger","H/Z#rightarrow b#bar{b} tagger"
                                              ,"H/Z#rightarrow b#bar{b} w/ mass info.","subjet b tagging"};
        TString sigName = "radion";

          // TString sigName = "blkgrav"
            
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
          for(unsigned int iC = 0; iC < cuts.size(); ++iC){
    Plotter * pr = new Plotter();
    Plotter * pe = new Plotter();
    Plotter * pb = new Plotter();    
    
          for(unsigned int iV = 0; iV < vars.size(); ++iV){
            TH1* hs = 0;
            f->GetObject(sigName + "_"+sigs[iS]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hs);

            TH1* hb = 0;
            f->GetObject(bkgs[iB]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hb);      
            
            if(hb==0){
              std::cout << bkgs[iB]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV]<< std::endl;
              continue;
            }
            if(hs==0){
              std::cout << sigName + "_"+sigs[iS]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV]<< std::endl;
              continue;
            }
            PlotTools::toUnderflow(hs);
            PlotTools::toOverflow(hs);
            PlotTools::toUnderflow(hb);
            PlotTools::toOverflow(hb);
            

            
            auto roc = PlotTools::getRocCurve(hs,hb,true,sigs[iS],bkgs[iB]);
            pr->addGraph(roc,varNs[iV]);           
            
          }
          pe->setMinMax(0.0001,1.0);
          auto c = pr->draw(false,bkgs[iB]+"_"+sigName +"_"+sigs[iS]+"_"+cuts[iC]);
          c->SetLogy(true);
          c->Update();
          
          
                    // ps->draw(false,sigs[iS]+cuts[iC]+"_"+"_"+hbbStat);
                              // pb->draw(false,bkgs[iB]+cuts[iC]+"_"+"_"+hbbStat);
        }
      }
    }
  
}


{
  TFile * f = new TFile("getHbbBtagging_plots.root");
  std::vector<TString> sigs = {"m800","m1000","m2000","m2500","m3000"};
        std::vector<TString> bkgs = {"bkg"};
        std::vector<TString> cuts = {"emu_LHP_HbbTMass"};        
        
                std::vector<TString> vars = {"bbt","mdZHbb","Hbb","sJBTagging"};  
        TString sigName = "blkgrav";


          for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          for(unsigned int iV = 0; iV < vars.size(); ++iV){            
            Plotter * pe = new Plotter();
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){

    TH1* hs = 0;
    f->GetObject(sigName + "_"+sigs[iS]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hs);
    PlotTools::toUnderflow(hs);
    PlotTools::toOverflow(hs);
    hs = PlotTools::getIntegral(hs,true,true);        
    if(hs==0) continue;   
    pe->addHistLine(hs,sigs[iS],StyleInfo::getLineColor(iS));  
      for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
            TH1* hb = 0;
            f->GetObject(bkgs[iB]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hb);      
            if(hb==0) continue;            
            PlotTools::toUnderflow(hb);
            PlotTools::toOverflow(hb);                                                
            hb = PlotTools::getIntegral(hb,true,true);                      
            // pe->addHistLine(hb,bkgs[iB]+": "+sigs[iS],StyleInfo::getLineColor(iS),9);            
          }          
        }
        pe->setYTitle("eff.");
        pe->setMinMax(0.001,100.0);
        auto c = pe->draw(false,cuts[iC]+"_"+vars[iV]);
        c->SetLogy(true);
        c->Update();
      }
    }
  
}


//// Wqq


{
  // TFile * f = new TFile("testHBBTagging_TT12LSig.root");
    TFile * f = new TFile("getWqqTagging_plots.root");
  // std::vector<TString> sigs = {"m800","m1000","m2000","m2500","m3000"};
    std::vector<TString> sigs = {"m1000","m2000","m3000"};
        std::vector<TString> bkgs = {"ttbar","wjets","qcd"};
        // std::vector<TString> cuts = {"emu_LHP_HbbIncl","emu_LHP_Hbb2SJ","emu_LHP_Hbb2SJTMass"};
                // std::vector<TString> cuts = {"emu_LHP_HbbTMass","emu_LHP_Hbb2SJ","emu_LHP_Hbb2SJTMass"};
                std::vector<TString> cuts = {"emu_LHP"};
                
                // std::vector<TString> vars = {"tau2otau1","tau21ddt","wqq",};
                // std::vector<TString> varNs = {"tau2otau1","tau21ddt","wqq",};


                std::vector<TString> vars = {"tau2otau1","oldTau21Bins","newTau21Bins"};
                std::vector<TString> varNs = {"tau2otau1","oldTau21Bins","newTau21Bins"};
        TString sigName = "radion";

          // TString sigName = "blkgrav"
            
  for(unsigned int iS = 0; iS < sigs.size(); ++iS){
      for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
          for(unsigned int iC = 0; iC < cuts.size(); ++iC){
    Plotter * pr = new Plotter();
    Plotter * pe = new Plotter();
    Plotter * pb = new Plotter();    
    
          for(unsigned int iV = 0; iV < vars.size(); ++iV){
            TH1* hs = 0;
            f->GetObject(sigName + "_"+sigs[iS]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hs);

            TH1* hb = 0;
            f->GetObject(bkgs[iB]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV],hb);      
            
            if(hb==0){
              std::cout << bkgs[iB]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV]<< std::endl;
              continue;
            }
            if(hs==0){
              std::cout << sigName + "_"+sigs[iS]+"_"+cuts[iC]+"_"+sigs[iS]+"_"+vars[iV]<< std::endl;
              continue;
            }
            PlotTools::toUnderflow(hs);
            PlotTools::toOverflow(hs);
            PlotTools::toUnderflow(hb);
            PlotTools::toOverflow(hb);
            

            
            auto roc = PlotTools::getRocCurve(hs,hb,iV>0,sigs[iS],bkgs[iB]);
            pr->addGraph(roc,varNs[iV]);           
            
          }
          pe->setMinMax(0.0001,1.0);
          auto c = pr->draw(false,bkgs[iB]+"_"+sigName +"_"+sigs[iS]+"_"+cuts[iC]);
          // c->SetLogy(true);
          // c->Update();
          
          
                    // ps->draw(false,sigs[iS]+cuts[iC]+"_"+"_"+hbbStat);
                              // pb->draw(false,bkgs[iB]+cuts[iC]+"_"+"_"+hbbStat);
        }
      }
    }
  
}


https://www.dropbox.com/sh/p2959irsplp1equ/AADFYNuUYBqD4wBlSJtbU-lLa?dl=0



