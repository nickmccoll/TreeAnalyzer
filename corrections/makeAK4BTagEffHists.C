{
  TString inFN = "ak4BTagSFEffHists.root";
  TString outFN = "ak4_csvEff.root";
  TFile * inF = new TFile(inFN,"read");
  TFile * outF = new TFile(outFN,"recreate");
  
  TString denN = "incl";
  std::vector<TString> numNs {"loose","med","tight"};
  std::vector<TString> flvs {"l","c","b"};
  
  double ptBins[] = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,700,800,900,1000};
  int nPTBins = 22;
  
  double ptBins2[] = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,800,1000};
  int nPTBins2 = 20;
  
  float maxX = 900;
  
  double etaBins[] = {0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,2.4};
  int nETABins = 9;
  
  for(unsigned int iF = 0; iF < flvs.size(); ++iF){
    TH2 * hd = 0; 
    inF->GetObject(TString::Format("%s_%s",flvs[iF].Data(),denN.Data()),hd);
    if(hd == 0) continue;
    
    for(unsigned int iN = 0; iN < numNs.size(); ++iN){
      TH2 * hn = 0; 
      TString name = TString::Format("%s_%s",flvs[iF].Data(),numNs[iN].Data());
      inF->GetObject(name,hn);
      if(hn == 0) continue;
      
      auto doRebin =[&](const TH2* oh, const TString& newName, int nptbins, double * ptbins) ->TH2* {
        TH2* nh = new TH2F(newName,";jet p_{T}[GeV];jet |#eta|",nptbins,ptbins,nETABins,etaBins);
        for(unsigned int iX = 1; iX <= oh->GetNbinsX()+1; ++iX){
          for(unsigned int iY = 1; iY <= oh->GetNbinsY(); ++iY){
            float xVal = oh->GetXaxis()->GetBinCenter( iX == oh->GetNbinsX()+1 ? oh->GetNbinsX() : iX );
            float yVal = oh->GetYaxis()->GetBinCenter( iY );
            float val = oh->GetBinContent(iX,iY);
            nh->Fill(xVal,yVal,val);
          }
        }
        return nh;
      };

      hn = doRebin(hn,name + "new_N", iF == 0 ? nPTBins: nPTBins2 , iF == 0 ? ptBins:ptBins2    );
      hn->Divide(doRebin(hd,name + "new_D", iF == 0 ? nPTBins: nPTBins2 , iF == 0 ? ptBins:  ptBins2   ));
            
      // hn->Divide(hd);
      outF->cd();
      hn->Write(name);
      
    }
    
    
  }
  
  outF->Close();
  
}

{
  TString inFN = "ak4_csvEff.root";
  TFile * inF = new TFile(inFN,"read");
  
  std::vector<TString> numNs {"loose","med","tight"};
  std::vector<TString> flvs {"l","c","b"};
  
  for(unsigned int iF = 0; iF < flvs.size(); ++iF){

    for(unsigned int iN = 0; iN < numNs.size(); ++iN){
      TString name = TString::Format("%s_%s",flvs[iF].Data(),numNs[iN].Data());
      TH2 * hn = 0; 
      inF->GetObject(name,hn);
      if(hn == 0) continue;
      new TCanvas(name,name); hn->Draw("COLZ"); 
    }
    
    
  }

  
}