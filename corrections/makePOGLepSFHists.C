//Make electron scale factor files
{
  vector<TString> years = {"2016","2017","2018"};
  vector<TString> chans = {"1l","2l"};  
  
  auto makeCopy =[](TString fileName, TString histName, TString newName)->TH1* {
    TFile * f = new TFile(fileName,"read");    
    TH1 *h = 0;
    f->GetObject(histName,h);
    h->SetDirectory(0);
    h->SetName(newName);
    f->Close();
    delete f;
    return h;
  };
  
  auto fillHistogram = [](TH2 * inHist, TH2 * outHist){
    for(int iX = 1; iX <= inHist->GetNbinsX(); ++iX){
      int newX = outHist->GetXaxis()->FindFixBin(inHist->GetXaxis()->GetBinCenter(iX));
      for(int iY = 1; iY <= inHist->GetNbinsY(); ++iY){
          int newY = outHist->GetYaxis()->FindFixBin(inHist->GetYaxis()->GetBinCenter(iY));
          outHist->SetBinContent(newX,newY,inHist->GetBinContent(iX,iY));
          outHist->SetBinError(newX,newY,inHist->GetBinError(iX,iY));
      }
    }
  };
  auto fillHistogramSlice = [](TH2 * inHist, TH2 * outHist){
    for(int iX = 1; iX <= outHist->GetNbinsX(); ++iX){
      int oldX = inHist->GetXaxis()->FindFixBin(outHist->GetXaxis()->GetBinCenter(iX));
      outHist->SetBinContent(iX,1,inHist->GetBinContent(oldX,1));
      outHist->SetBinError(iX,1,inHist->GetBinError(oldX,1));
    }
  };
  
  for(int iY = 0; iY < years.size();++iY){
    for(int iC = 0; iC < chans.size();++iC){
      TFile * fo = new TFile(TString::Format("electron_%s_%s_SF.root",chans[iC].Data(),years[iY].Data()),"recreate");
            
      //reco
      TH2 * hr = 0;
      if(iY <  2){ //2016/7 is broken up into two files, let's merge them
        TH2 *hr1 = (TH2*)makeCopy(TString::Format("electron_%s_reco.root",years[iY].Data()),"EGamma_SF2D","reco_lowpt");
        TH2 *hr2 = (TH2*)makeCopy(TString::Format("electron_%s_reco_lowpt.root",years[iY].Data()),"EGamma_SF2D","reco_highpt");
        vector<double> xBins;
        vector<double> yBins;
        yBins.push_back(10);
        for(int iX = 1; iX <= hr1->GetNbinsX(); ++iX){
          xBins.push_back(hr1->GetXaxis()->GetBinLowEdge(iX));      
        }
        xBins.push_back(hr1->GetXaxis()->GetBinLowEdge(hr1->GetNbinsX())+hr1->GetXaxis()->GetBinWidth(hr1->GetNbinsX()) );              
        for(int iY = 1; iY <= hr1->GetNbinsY(); ++iY){
          yBins.push_back(hr1->GetYaxis()->GetBinLowEdge(iY));
        }
        yBins.push_back(hr1->GetYaxis()->GetBinLowEdge(hr1->GetNbinsY())+hr1->GetYaxis()->GetBinWidth(hr1->GetNbinsY()) );

        hr = new TH2F("reco","",xBins.size()-1, &xBins[0],yBins.size()-1, &yBins[0]);
        fillHistogram(hr1,hr);
        fillHistogramSlice(hr2,hr);
      } else {
        hr = (TH2*)makeCopy(TString::Format("electron_%s_reco.root",years[iY].Data()),"EGamma_SF2D","reco");
      }
      fo->cd();
      hr->Write();
      delete hr;            
      
      //ID
      TString idFileName;
      TString idHistName;
      if(chans[iC]=="1l"){
        idFileName = TString::Format("electron_%s_MVA90ID.root",years[iY].Data());
        idHistName = "EGamma_SF2D";
      } else {
        idFileName = TString::Format("electron_%s_SUSYSF.root",years[iY].Data());
        idHistName =TString::Format("Run%s_CutBasedMediumNoIso94XV2",years[iY].Data());            
      }
      TH2* hi =(TH2*)makeCopy(idFileName,idHistName,"id");
      fo->cd();
      hi->Write();
      delete hi;
      
      
      fo->Close();
    }    
  }
  
  
}




//MUON SCALE FACTORS
{
  vector<TString> years = {"2016","2017","2018"};
  vector<TString> chans = {"1l","2l"};  
  
  auto makeCopy =[](TString fileName, TString histName, TString newName)->TH1* {
    TFile * f = new TFile(fileName,"read");    
    TH1 *h = 0;
    f->GetObject(histName,h);
    h->SetDirectory(0);
    h->SetName(newName);
    f->Close();
    delete f;
    return h;
  };
  
  
  struct SpliceSettings{
    TH2* hist =0;
    bool isAbsEta = false;
    bool etaIsX = false;
    
    TAxis * ptAx() {
      if(etaIsX) return hist->GetYaxis();
      else return hist->GetXaxis();
    };

    
    TAxis * etaAx() {
      if(etaIsX) return hist->GetXaxis();
      else return hist->GetYaxis();
    };      
    
    int ptBins() {return ptAx()->GetNbins();}
    int etaBins() {return etaAx()->GetNbins();}
    
    float getETA(int iE){
      float eta = etaAx()->GetBinCenter(iE);
      if(isAbsEta) return fabs(eta);
      return eta;         
    }
    int getETABin(float eta){
      if(isAbsEta) eta = fabs(eta);
      return etaAx()->FindFixBin(eta);
    }
    
    int getXBin(int iP, int iE){
      return etaIsX ? iE : iP;
    }
    int getYBin(int iP, int iE){
      return etaIsX ? iP : iE;
    }
    
    float getContent(int iP, int iE){
      int iX = getXBin(iP,iE);
      int iY = getYBin(iP,iE);
      return hist->GetBinContent(iX,iY);
    }
    float getError(int iP, int iE){
      int iX = getXBin(iP,iE);
      int iY = getYBin(iP,iE);
      return hist->GetBinError(iX,iY);
    }
    void setBin(int iP, int iE, float content, float error){
      int iX = getXBin(iP,iE);
      int iY = getYBin(iP,iE);
      hist->SetBinContent(iX,iY,content);
      hist->SetBinError(iX,iY,error);
    }
  };
  
  

  auto spliceHistograms = [](TString newHistName ,SpliceSettings lowPT, SpliceSettings highPT)->TH2*{
    
    vector<double> ptBins;    
    for(int iP = 1; iP <= lowPT.ptBins(); ++iP){
      double lowEdge = lowPT.ptAx()->GetBinLowEdge(iP);
      if(lowEdge < 5 ) continue;
      if(lowEdge>=20) break;
      ptBins.push_back(lowEdge);
    }    
    int nlowPTBins = ptBins.size();
    
    for(int iP = 1; iP <= highPT.ptBins(); ++iP){
      double lowEdge = highPT.ptAx()->GetBinLowEdge(iP);
      if(lowEdge < 20 ) continue;
      ptBins.push_back(lowEdge);
    }
    ptBins.push_back(ptBins.back()+highPT.ptAx()->GetBinWidth(highPT.ptBins()) );
    
    for(auto y: ptBins) cout <<y <<", ";
    cout << nlowPTBins << endl;
    cout << endl;
    
    //output will always be x-axis is non-ABS eta
    //will do a simple +/-2.4 with .1 binning
    SpliceSettings spliced;
    spliced.isAbsEta = false;
    spliced.etaIsX = true;
    spliced.hist = new TH2F(newHistName,"",48, -2.4,2.4,ptBins.size()-1, &ptBins[0]);
    
    //add in the highPT            
    for(int newETA = 1; newETA <= spliced.etaBins(); ++newETA){
      int oldETA = highPT.getETABin(spliced.getETA(newETA));
      //so we don't overlap in the two areas
      for(int newPT = nlowPTBins+1; newPT <= spliced.ptBins(); ++newPT){
        int oldPT = highPT.ptAx()->FindFixBin(spliced.ptAx()->GetBinCenter(newPT)); 
        float content = highPT.getContent(oldPT,oldETA);      
        float error   = highPT.getError(oldPT,oldETA);
        spliced.setBin(newPT,newETA,content,error);
      }
    }
    //add in the lowPT            
    for(int newETA = 1; newETA <= spliced.etaBins(); ++newETA){
      int oldETA = lowPT.getETABin(spliced.getETA(newETA));
      //so we don't overlap in the two areas
      for(int newPT = 1; newPT <= nlowPTBins; ++newPT){
        int oldPT = lowPT.ptAx()->FindFixBin(spliced.ptAx()->GetBinCenter(newPT)); 
        float content = lowPT.getContent(oldPT,oldETA);      
        float error   = lowPT.getError(oldPT,oldETA);
        spliced.setBin(newPT,newETA,content,error);
      }
    }
    return spliced.hist;
  };
  
  
  
  
  for(int iY = 0; iY < years.size();++iY){
    for(int iC = 0; iC < chans.size();++iC){
      TFile * fo = new TFile(TString::Format("muon_%s_%s_SF.root",chans[iC].Data(),years[iY].Data()),"recreate");
            
      //ID
      TH2 * hi = 0;
      TString histNameLowPT;
      TString histNameHighPT;
      if(chans[iC]=="1l"){
        histNameLowPT = "NUM_MediumID_DEN_genTracks_pt_abseta";
        histNameHighPT ="NUM_MediumID_DEN_genTracks_pt_abseta";
        if(iY == 0) histNameHighPT = "NUM_MediumID_DEN_genTracks_eta_pt";
        if(iY == 2) histNameHighPT = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
      } else {
        histNameLowPT = "NUM_LooseID_DEN_genTracks_pt_abseta";
        histNameHighPT ="NUM_LooseID_DEN_genTracks_pt_abseta";
        if(iY == 0) histNameHighPT = "NUM_LooseID_DEN_genTracks_eta_pt";
        if(iY == 2) histNameHighPT = "NUM_LooseID_DEN_TrackerMuons_pt_abseta";
      }
      
      if(iY == 0){ //2016 has two seperate files for some reason        
        SpliceSettings lowPT;
        lowPT.isAbsEta = true;
        lowPT.etaIsX = false;
        lowPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_BCDEF_ID_lowpt.root",years[iY].Data()),histNameLowPT,"id_lowpt_1");
        
        SpliceSettings highPT;
        highPT.isAbsEta = false;
        highPT.etaIsX = true;
        highPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_BCDEF_ID.root",years[iY].Data()),histNameHighPT,"id_highpt_1");
        
        auto * h_1 = spliceHistograms("id_1",lowPT,highPT);
        
        lowPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_GH_ID_lowpt.root",years[iY].Data()),histNameLowPT,"id_lowpt_2");
        highPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_GH_ID.root",years[iY].Data()),histNameHighPT,"id_highpt_2");

        auto * h_2 = spliceHistograms("id_2",lowPT,highPT);

        if(h_1->GetNbinsX() != h_2->GetNbinsX() || h_1->GetNbinsY() != h_2->GetNbinsY()){
          cout <<"BAD BINNING!!!"<<endl;
          return;
        }        
        hi = (TH2*)h_1->Clone();
        for(int iX = 1; iX <= h_1->GetNbinsX(); ++iX ){
          if(h_1->GetXaxis()->GetBinCenter(iX) != h_2->GetXaxis()->GetBinCenter(iX)){
            cout <<"BAD BINNING in the XAXIS !!!"<<endl;
            return;
          }
          for(int iY = 1; iY <= h_1->GetNbinsY(); ++iY ){
            if(h_1->GetYaxis()->GetBinCenter(iY) != h_2->GetYaxis()->GetBinCenter(iY)){
              cout <<"BAD BINNING in the YAXIS !!!"<<endl;
              return;
            }
            // 35.9 total luminosity for 4b paper
            // G_H lumnosity 16.6
            float L_tot = 35.9;
            float L_2 = 16.6;
            float L_1 = L_tot-L_2;
            float newSF = L_1*h_1->GetBinContent(iX,iY)/L_tot + L_2*h_2->GetBinContent(iX,iY)/L_tot;
            float newErr = L_1*h_1->GetBinError(iX,iY)/L_tot + L_2*h_2->GetBinError(iX,iY)/L_tot;
            hi->SetBinContent(iX,iY,newSF);
            hi->SetBinError(iX,iY,newErr);
          }
        }
        delete h_1;
        delete h_2;
        
        
      } else {
        SpliceSettings lowPT;
        lowPT.isAbsEta = true;
        lowPT.etaIsX = false;
        lowPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_ID_lowpt.root",years[iY].Data()),histNameLowPT,"id_lowpt");
        
        SpliceSettings highPT;
        highPT.isAbsEta = true;
        highPT.etaIsX = false;
        highPT.hist = (TH2*)makeCopy(TString::Format("muon_%s_ID.root",years[iY].Data()),histNameHighPT,"id_highpt");
        
        hi = spliceHistograms("id",lowPT,highPT);
      }
      fo->cd();
      hi->Write();
      delete hi;                        
      
      fo->Close();
    }    
  }
  
  
}


