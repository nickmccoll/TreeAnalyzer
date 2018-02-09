//CORRRECT VERSION...do not use other!!!!
{
  int n = 100000;
  float sigma = 0.1;
  float mu = 1.1;
  
  float gen = 100;
  
  float r = 1.1;
  float s = 1.3;
  
  auto ran = new TRandom3();
  
  TH1 * hn = new TH1F("orig","orig",300,0,gen*3);
  TH1 * h1 = new TH1F("v1","v1",300,0,gen*3);
  TH1 * h2 = new TH1F("v2","v2",300,0,gen*3);
  
  for(unsigned int i = 0; i < n; ++i){
    float reco = ran->Gaus(gen*mu,sigma*gen);
    hn->Fill(reco);
    // float cor1 = s*reco +(r-1)*(s*reco-gen);
    // float cor2 = s*(reco +(r-1)*(reco-gen));
    // float cor1 = gen +r*(s*reco-gen);
    float cor1 = gen +r*(reco-gen) + s*reco ;
    
    float cor2 = s*(gen +r*(reco-gen));
    h1->Fill(cor1);
        h2->Fill(cor2);
    
  }
  
  new TCanvas();
  hn->Draw();
  new TCanvas();
  h1->Draw();
  new TCanvas();
  h2->Draw();
  
  
  
}

{
  auto normAbove = [] (TH1 * h, float xv = 50.0) -> TH1 *{
    int b = h->FindFixBin(xv);
    float n = h->Integral(b,-1);
    h->Scale(1./n);
    return h;
  };
  
  Plotter * p = new Plotter();
  p->addHistLine(normAbove(res_emu_LMT_lWW_hbb_mass),"Resonant");
  p->addHistLine(normAbove(ttbar_emu_LMT_lWW_hbb_mass),"ttbar");
  p->addHistLine(normAbove(nonResWVV_emu_LMT_lWW_hbb_mass),"Non-resonant: MC w/o t quarks");
  p->addHistLine(normAbove(nonResT_emu_LMT_lWW_hbb_mass),"Non-resonant: MC w/ t quarks");

  p->setYTitle("arbitrary units");
    p->draw();
}


{
  auto normAbove = [] (TH1 * h, float xv = 50.0) -> TH1 *{
    int b = h->FindFixBin(xv);
    float n = h->Integral(b,-1);
    h->Scale(1./n);
    return h;
  };
  
  Plotter * p = new Plotter();
  p->addHistLine(normAbove(nonResWVV_emu_L_lWW_hbb_mass),"Loose Hbb tagging");
  p->addHistLine(normAbove(nonResWVV_emu_M_lWW_hbb_mass),"Medium Hbb tagging");
  p->addHistLine(normAbove(nonResWVV_emu_T_lWW_hbb_mass),"Tight Hbb tagging");
  p->setYTitle("arbitrary units");
  p->rebin(2);
  p->draw(false,"w");
  
  Plotter * p2 = new Plotter();
  p2->addHistLine(normAbove(nonResT_emu_L_lWW_hbb_mass),"Loose Hbb tagging");
  p2->addHistLine(normAbove(nonResT_emu_M_lWW_hbb_mass),"Medium Hbb tagging");
  p2->addHistLine(normAbove(nonResT_emu_T_lWW_hbb_mass),"Tight Hbb tagging");
  p2->setYTitle("arbitrary units");
    p2->rebin(2);
  p2->draw(false,"t");
  
  
  Plotter * p3 = new Plotter();
  p3->addHistLine(normAbove(nonRes_emu_L_lWW_hbb_mass),"Loose Hbb tagging");
  p3->addHistLine(normAbove(nonRes_emu_M_lWW_hbb_mass),"Medium Hbb tagging");
  p3->addHistLine(normAbove(nonRes_emu_T_lWW_hbb_mass),"Tight Hbb tagging");
  p3->setYTitle("arbitrary units");
    p3->rebin(2);
  p3->draw(false,"all");
}

{
  auto normAbove = [] (TH1 * h, float xv = 50.0) -> TH1 *{
    int b = h->FindFixBin(xv);
    float n = h->Integral(b,h->GetNbinsX());
    h->Scale(1./n);
    return h;
  };
  auto zeroError = [] (TH1 * h) -> TH1 *{
    for(unsigned int iB = 0; iB <= h->GetNbinsX()+1; ++iB)
      h->SetBinError(iB,0);
    return h;
  };
  
  
  Plotter * p3 = new Plotter();


  // p3->addHistLine(zeroError(normAbove(downRpdf__x)),"down resolution PDF",StyleInfo::getLineColor(3));
  // p3->addHistLine(zeroError(normAbove(upRpdf__x)),"up resolution PDF",StyleInfo::getLineColor(2));
  p3->addHistLine(zeroError(normAbove(downSpdf__x)),"down scale PDF",StyleInfo::getLineColor(3));
  p3->addHistLine(zeroError(normAbove(upSpdf__x)),"up scale PDF",StyleInfo::getLineColor(2));
    p3->addHistLine(zeroError(normAbove(nompdf__x)),"nominal PDF",StyleInfo::getLineColor(1));
    p3->addHist(normAbove(nomhist__x),"MC",1);
    
    
        p3->setUnderflow(false);
    p3->setOverflow(false);
  p3->setYTitle("arbitrary units");
  p3->setXTitle("Hbb soft-drop mass [GeV]");
  // p3->rebin(5);
  // p3->drawSplitRatio(2);
  p3->draw();
}


/// Explode a 2D
{
  std::vector<double> bins = {50,60,80,100,120,140,160,180,200,220,250};
  bool binInY = true;
    // std::vector<double> bins = {700,800,900,1000,1500,2000,3000,4000,5000};
      // bool binInY = false;
  TH2 * hm = mvvcond__nom_data;
  TH2 * hf = mvvcond__nom_pdfD;  
  hf->Scale(hm->Integral()/hf->Integral());
  
  const TAxis * ax = hm->GetXaxis();
  if(binInY) ax = hm->GetYaxis();
  
  for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
    int binL = ax->FindFixBin(bins[iB]);
    int binH = ax->FindFixBin(bins[iB+1]) -1;
    Plotter * p = new Plotter();
    
    TH1 * hm1D = 0;
    TH1 * hf1D = 0;
    if(binInY){
      hm1D = hm->ProjectionX(TString::Format("hm_proj_%u",iB),binL,binH);
      hf1D = hf->ProjectionX(TString::Format("hf_proj_%u",iB),binL,binH);
    } else {
      hm1D = hm->ProjectionY(TString::Format("hm_proj_%u",iB),binL,binH);
      hf1D = hf->ProjectionY(TString::Format("hf_proj_%u",iB),binL,binH);
    }
    
    p->addHist(hm1D,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    p->addHistLine(hf1D,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    p->setMinMax(.0001,hm1D->Integral());
    p->rebin(5);
    auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));    
    c->SetLogy();
    c->Update();
  }
  
}


  
  //Compare distributions;
  {
    TString type = "radHH";
    // std::vector<TString> sels = {"emu_LMT","emu_LMT_lWW","emu_LMT_lb","emu_LMT_ltb","emu_LMT_ltmb" };
    // std::vector<TString> sels = {"emu_L_full","emu_L_ltmb","emu_L_none","emu_M_ltmb","emu_M_full" };
    // std::vector<TString> sels = {"emu_LMT_ltmb","emu_LMT_full","emu_LMT_none" };
        // std::vector<TString> sels = {"emu_LMT_none","emu_L_none","emu_M_none","emu_T_none" };
        // std::vector<TString> sels = {"m1000_emu_L_none","m1000_emu_M_none","m1000_emu_T_none","m2500_emu_L_none","m2500_emu_M_none","m2500_emu_T_none" };
        std::vector<TString> sels = {"m4500_emu_T_none","m4500_emu_T_ltmb","m4500_emu_T_full" };

  // std::vector<double> bins = {30,210};
  // bool binInY = false;
    // std::vector<double> bins = {30,50,100,150,210};
    // bool binInY = false;
    // std::vector<double> bins = {800,1000,2000,3000,5000};
      // bool binInY = true;
  // std::vector<double> bins = {800,900,1000,1250,1500,2000,3000,4000,5000};
  //   bool binInY = true;
    // std::vector<double> bins = {800,900,1000,1200,1400,1600,2000};
      // bool binInY = true;
    std::vector<double> bins = {800,5000};
      bool binInY = true;
  std::vector<TH2*> hs;
  std::vector<TString> hNs;
  TFile *fY = new TFile(TString::Format("HHlnujj_%s_distributions.root",type.Data()),"read");
  for(unsigned int iS = 0; iS < sels.size(); ++iS){
    TH2* h = 0; 
    fY->GetObject(TString::Format("%s_%s_hbbMass_hhMass",type.Data(),sels[iS].Data()),h);hs.push_back(h);hNs.push_back(sels[iS]);    
  }    
  // for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(1.0/hs[iH]->Integral());
  
  const TAxis * ax = hs[0]->GetXaxis();
  if(binInY) ax = hs[0]->GetYaxis();

  for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
    int binL = ax->FindFixBin(bins[iB]);
    int binH = ax->FindFixBin(bins[iB+1]) -1;
    Plotter * p = new Plotter();
    
    for(unsigned int iH = 0; iH < hs.size(); ++iH){
      TH1 * h1D = 0;
      if(binInY){
        h1D = hs[iH]->ProjectionX(TString::Format("%s_proj_%u",hNs[iH].Data(),iB),binL,binH);
      } else {
        h1D = hs[iH]->ProjectionY(TString::Format("%s_proj_%u",hNs[iH].Data(),iB),binL,binH);
      }
      // if(iH == 0) p->addHist(h1D,TString::Format("%s %.0f-%.0f",hNs[iH].Data(),bins[iB],bins[iB+1]));
      // else  p->addHistLine(h1D,TString::Format("%s %.0f-%.0f",hNs[iH].Data(),bins[iB],bins[iB+1]));
      p->addHist(h1D,TString::Format("%s",hNs[iH].Data()));    
    }    
    // p->setMinMax(.00001,0.3);
    p->setBotMinMax(0,3);
    p->setUnderflow(false);
    p->setOverflow(false);
    // p->rebin( binInY ? 4 : 24);
        p->rebin( 2);
    // auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    // c->SetLogy();
    // c->Update();
    auto * c = p->drawSplitRatio(0,"stack",false,false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    // c->GetPad(1)->SetLogy();
    c->GetPad(1)->Update();
  }

}


//TEST CATEGORIES
{
  TFile * f = new TFile("test.root");
  //std::vector<TString> cats = {"w","nw","wqq","nwqq"};
    // std::vector<TString> cats = {"nw","wnwb","wb"};
        std::vector<TString> cats = {"nwqq","wqqnwqqb","wqqb"};
  
  Plotter * p = new Plotter();
  for(const auto& c : cats){
    TH1 * h = 0;
    f->GetObject(TString::Format("ttbar_%s_emu_LI_hbb_mass",c.Data()),h);
    if(h == 0) continue;
    p->addHistLine(h,c);
  }
  p->normalize();
  p->rebin(2);
  p->draw();
   
}


{
  // TFile * f = new TFile("test.root");
  TFile * f = new TFile("testResCats.root");
  // std::vector<TString> cats = {"w","nw","wqq","nwqq","wqqnwqqb","wqqb"};
  // std::vector<TString> cats = {"wq","nwq","wqnwqqb","wqqb"};
    // std::vector<TString> cats = {"nwq","wq","wqnwqqb","wqqb","w","nw","wqq","nwqq","wqqnwqqb","nwqqb"};
    // std::vector<TString> cats = {"incl","wq","nwq","wqnwqqb","wqqb","nw","w"};
        std::vector<TString> cats = {"nwq","incl","wqnwqqb","wqqb","nwqqb"};
        // std::vector<TString> cats = {"nwq","wq","wqnwqqb","wqqb","w","nw"};
    // std::vector<TString> cats = {"nw","wnwb","wb"};
        // std::vector<TString> cats = {"nwqq","wqqnwqqb","wqqb"};
        // std::vector<TString> cats = {"nw","wnwb","wb"};
        // std::vector<TString> recos = {"lW","lW_LI"};
        // std::vector<TString> recos = {"lW_L"};
                // std::vector<TString> recos = {"lW_L","L","lW_M","M","lW_T","T"};
                  // std::vector<TString> recos = {"lW_L","lW_M","lW_T"};
                                    // std::vector<TString> recos = {"lW_L_hh700to900","lW_M_hh700to900","lW_T_hh700to900","lW_LI_hh900to1100","lW_LI_hh1400to1800"};
                                    
                                                                        std::vector<TString> recos = {"lW_M_hh700to900","lW_M_hh900to1100","lW_M_hh1400to1800","lW_M_hh1400to1800"};
                                                                        // {"lW_L_hh700to900","lW_L_hh900to1100","lW_L_hh1400to1800","lW_L_hh1400to1800"};
        // std::vector<TString> recos = {"lW_L","lW_M","lW_T","lW_L_hh1400to1800","lW_M_hh1400to1800","lW_T_hh1400to1800"};
        
                // std::vector<TString> recos = {"lW","lW_L"};
        
        // std::vector<TString> recos = {"L","lW_L","L_hh900to1100","lW_L_hh900to1100","L_hh1400to1800","lW_L_hh1400to1800"};
  

  for(const auto& c : cats){
      Plotter * p = new Plotter();
      for(const auto& r : recos){
      
    TH1 * h = 0;
    f->GetObject(TString::Format("ttbar_%s_emu_%s_hbb_mass",c.Data(),r.Data()),h);
    if(h == 0) continue;
    p->addHist(h,TString::Format("%s %s",c.Data(), r.Data()));
  }
  p->normalize();
  p->rebin(2);
  p->setBotMinMax(0.,2.);
  p->drawSplitRatio(0,"stack",false,false,c);
  }
   
}

{
  // TFile * f = new TFile("test.root");
  TFile * f = new TFile("testResCats.root");
  // std::vector<TString> cats = {"nqq","qqnqqb","nwqqb"};
  std::vector<TString> cats = {"nw","nwq"};
  // std::vector<TString> bkgs = {"ttbar","other"};
    std::vector<TString> bkgs = {"bkgNoQCD","wjets","ttbar"};
  // std::vector<TString> cats = {"qq","nqq","qqnqqb","qqb"};
    // std::vector<TString> cats = {"nw","wnwb","wb"};
        // std::vector<TString> cats = {"nwqq","wqqnwqqb","wqqb"};
        // std::vector<TString> cats = {"nw","wnwb","wb"};
        std::vector<TString> recos = {"lW_M_hh900to1100","lW_L_hh900to1100","lW_L_hh700to900","lW_M_hh700to900"};
    // std::vector<TString> recos = {"L_hh700to900","M_hh700to900","T_hh700to900","L_hh900to1100","M_hh900to1100","T_hh900to1100","L_hh1400to1800","M_hh1400to1800"};
    //
    // std::vector<TString> recos = {"lW_L_hh700to900","lW_M_hh700to900","lW_T_hh700to900","lW_L_hh900to1100","lW_M_hh900to1100","lW_T_hh900to1100","lW_L_hh1400to1800","lW_M_hh1400to1800"};

  

  for(const auto& c : cats){
          for(const auto& r : recos){
      Plotter * p = new Plotter();

          for(const auto& b : bkgs){      
    TH1 * h = 0;
    f->GetObject(TString::Format("%s_%s_emu_%s_hbb_mass",b.Data(),c.Data(),r.Data()),h);
    if(h == 0) continue;
    p->addHist(h,b);
  }
  // p->normalize();
  p->rebin(2);
  p->setBotMinMax(0.,2.);
  p->drawSplitRatio(0,"stack",false,false,c+"_"+r);
  }
}
   
}


//Plot bkg categories
{
// std::vector<TString> bkgs = {"nonResH0","nonResHM","resW","resT"};
// std::vector<TString> bkgNs = {"lost t/W (0 q) bkg.","lost t/W (#geq1 q) bkg.","m_{W} bkg.","m_{t} bkg."};

std::vector<TString> bkgs = {"qg","losttw","mw","mt"};
std::vector<TString> bkgNs = {"q/g bkg.","lost t/W bkg.","m_{W} bkg.","m_{t} bkg."};
std::vector<TString> sels = {"emu_LMT_ltmb","emu_M_ltmb"};
TFile * f = new TFile("HHlnujj_distributions.root");
for(auto& s: sels){
Plotter * p = new Plotter;
for(unsigned int iB = 0; iB < bkgs.size(); ++iB){
TH1 * h= 0;
f->GetObject(bkgs[iB]+"_"+s+"_hbbMass",h);
if(h==0) continue;
p->addStackHist(h,bkgNs[iB]);

}
p->rebin(2);
p->draw(false,s);
}

}
  
