{
  Plotter * p = new Plotter();
  p->addHist(signal_numTruePUInt,"signal");
  p->addHist(ttbar_numTruePUInt,"ttbar");
  p->normalize();
  p->draw();
}

{
  TFile * fm = new TFile("truePUDist.root","read");
  vector<TString> dataFiles = {"puDist_data_2016_69200.root","puDist_data_2016_66017.root","puDist_data_2016_72383.root"};
  vector<TString> dataPre = {"nominal","down","up"};
  TFile *fo = new TFile("puSF.root","recreate");
  TH1 * mcH = 0;
  fm->GetObject("ttbar_numTruePUInt",mcH);
  mcH->Scale(1.0/mcH->Integral(0,-1));
  for(unsigned int iF = 0; iF < dataFiles.size(); ++iF){
    TFile * fd = new TFile(dataFiles[iF],"read");
    TH1 *h = 0;
    fd->GetObject("pileup",h);
    h->Scale(1.0/h->Integral(0,-1));
    h->Divide(mcH);
    
    fo->cd();
    h->Write(TString::Format("puSF_%s",dataPre[iF].Data()));
      
    
  }
  
  
}

{
  Plotter * p = new Plotter();
  p->addHist(puSF_nominal,"puSF_nominal");
  p->addHist(puSF_down,"puSF_down");
  p->addHist(puSF_up,"puSF_up");
  p->draw();
}
