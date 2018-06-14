{
  bool makeWithHist = true;
  const float a = 0.0615;
  const float b = -0.0005;
  TH1 * h = new TH1F("weightConsts","; a / b / nf",3,-0.5,2.5);
  h->SetBinContent(1,a);
  h->SetBinContent(2,b);
  if(makeWithHist){
    TFile *fin = new TFile("topPTWeightInput.root","read");
    TH1 * hi = 0;
    fin->GetObject("counts",hi);
    h->SetBinContent(3,hi->GetBinContent(1)/hi->GetBinContent(2));
    fin->Close();
  } else {
    h->SetBinContent(3,1);
  }
    
    TFile * f = new TFile("topPTWeight.root","recreate");
    h->Write();
    f->Close();
    
}