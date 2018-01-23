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

{
  TH2 * hm = mvvcond__coarse;
  TH2 * hf = mvvcond__smoothTail;  
  hf->Scale(hm->Integral()/hf->Integral());
  
  const TAxis * ax = hm->GetYaxis();
  
  for(unsigned int iB = 1; iB <= ax->GetNbins(); ++iB){
    int binL = iB;
    int binH = iB;
    Plotter * p = new Plotter();
    
    TH1 * hm1D = 0;
    TH1 * hf1D = 0;
      hm1D = hm->ProjectionX(TString::Format("hm_proj_%u",iB),binL,binH);
      hf1D = hf->ProjectionX(TString::Format("hf_proj_%u",iB),binL,binH);
    
    p->addHistLine(hm1D,TString::Format("%.0f-%.0f",ax->GetBinLowEdge(iB),ax->GetBinLowEdge(iB)+ax->GetBinWidth(iB)));
    p->addHistLine(hf1D,TString::Format("%.0f-%.0f",ax->GetBinLowEdge(iB),ax->GetBinLowEdge(iB)+ax->GetBinWidth(iB)));
    p->setMinMax(.0001,hm1D->Integral());
    // p->rebin(5);
    auto * c = p->draw(false,TString::Format("%.0f-%.0f",ax->GetBinLowEdge(iB),ax->GetBinLowEdge(iB)+ax->GetBinWidth(iB)));    
    c->SetLogy();
    c->Update();
  }
  
  

  /// Test 2D kernels
  {
    
    
    auto addComp =[](TString filename,TString name, std::vector<TH2*>& hs,std::vector<TString>& hNs,bool addMC = false){
          TFile *f = new TFile(filename,"read");
          TH2* h = 0; 
          if(addMC){
            f->GetObject("histo_fine_data",h);hs.push_back(h);
            hNs.push_back("mc");
          }
          f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back(name);                    
    };
    auto addComp2 =[](TString prefix, std::vector<TH2*>& hs,std::vector<TString>& hNs,bool addMC = false){
          TFile *f = new TFile(TString("hhLnuJJ_nonResW_")+prefix+"_COND2D_template.root","read");
          TH2* h = 0; 
          if(addMC){
            f->GetObject("histo_fine_data" ,h);hs.push_back(h);
            hNs.push_back("mc");
          }
          f->GetObject("histo_smooth",h);     hs.push_back(h);hNs.push_back(prefix);                    
    };
    
    
    // std::vector<double> bins = {30,40,50,60,80,100,120,140,160,180,200,220,250};
    // bool binInY = true;
      std::vector<double> bins = {700,800,900,1000,1500,2000,3000,4000,5000};
        bool binInY = false;
    std::vector<TH2*> hs;
    std::vector<TString> hNs;
    
    
    TFile *f = new TFile("hhLnuJJ_nonResT_COND2D_template.root","read");
    TH2* h = 0;
    f->GetObject("histo_fine_data",h);hs.push_back(h);hNs.push_back("mc");
    f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back("KDE");
    f->GetObject("histo_smooth",h);  hs.push_back(h);hNs.push_back("KDES");
    // addComp2("khxs_1_khxc_2_khys_1_khyc_2"     ,hs,hNs,true);
    // addComp2("khxs_1_khxc_1p5_khys_1_khyc_1p5"     ,hs,hNs);
        // addComp2("khxs_1_khxc_1_khys_1_khyc_1p5"     ,hs,hNs);
    // addComp2("khxs_0p75_khxc_2_khys_1_khyc_2"   ,hs,hNs);
    // addComp2("khxs_0p75_khxc_2_khys_0p75_khyc_2"     ,hs,hNs,true);
    // addComp2("khxs_0p75_khxc_1_khys_0p75_khyc_1"  ,hs,hNs);
    // addComp2("khxs_0p75_khxc_1_khys_0p75_khyc_2"   ,hs,hNs);
    // addComp2("khxs_0p75_khxc_1_khys_1_khyc_1"  ,hs,hNs);
    // addComp2("khxs_0p5_khxc_2_khys_1_khyc_2"   ,hs,hNs);
        // addComp2("khxs_0p5_khxc_2_khys_0p5_khyc_2"   ,hs,hNs);
        // addComp2("khxs_0p25_khxc_2_khys_0p25_khyc_2"   ,hs,hNs);

                // addComp2("khxs_0p25_khxc_2_khys_0p5_khyc_2"   ,hs,hNs);
        // addComp2("khxs_0p75_khxc_2_khys_0p75_khyc_2"   ,hs,hNs);
        // addComp2("khxs_0p75_khxc_1_khys_0p75_khyc_1"   ,hs,hNs);
        // addComp2("khxs_0p75_khxc_0p75_khys_0p75_khyc_0p75"   ,hs,hNs);
        // addComp2("khxs_0p75_khxc_0p5_khys_0p75_khyc_0p5"   ,hs,hNs);
        // addComp2("khxs_0p5_khxc_0p5_khys_0p5_khyc_0p5"   ,hs,hNs);
              // addComp2("khxs_0p5_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p25_khxc_0p5_khys_0p5_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p25_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p5_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p75_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                                            // addComp2("khxs_0p75_khxc_0p75_khys_0p25_khyc_0p5"   ,hs,hNs);
                                            // addComp2("khxs_0p75_khxc_1_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_1_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p75_khxc_0p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_1_khxc_1p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_0p75_khxc_1p5_khys_0p25_khyc_0p5"   ,hs,hNs);
                        // addComp2("khxs_1_khxc_2_khys_0p25_khyc_0p5"   ,hs,hNs);
    for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral()/hs[iH]->Integral());              
    for(unsigned int iH = 1; iH < hs.size(); ++iH) for(unsigned int iX = 1; iX <= hs[0]->GetNbinsX(); ++iX) hs[iH]->SetBinError(iX,0);
  
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
        if(iH == 0) p->addHist(h1D,TString::Format("%s",hNs[iH].Data()));
        // else  p->addHist(h1D,TString::Format("%s",hNs[iH].Data()),-1,1,4,20,1,false,true,false,"E");
        else  p->addHistLine(h1D,TString::Format("%s",hNs[iH].Data()));
      }    
      p->setMinMax(.0001,hs[0]->Integral()/4);
      p->setUnderflow(false);
      p->setOverflow(false);
      p->setBotMinMax(0,2);
      p->rebin(25);
      // auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
      // c->SetLogy();
      // c->Update();
      
      auto * c = p->drawSplitRatio(0,"stack",false,false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
      c->GetPad(1)->SetLogy();
      c->GetPad(1)->Update();
    }
  
  }
  
  
  
  /// Test 1D kernels
  {
    TFile *f = new TFile("HHlnujj_resT_incl_template.root","read");
    std::vector<TH1*> hs;
    std::vector<TString> hNs;
    TH1* h = 0; 
    f->GetObject("histo_data",h);hs.push_back(h);hNs.push_back("MC");
    // f->GetObject("histo_OPTUp",h);     hs.push_back(h);hNs.push_back("Up scale template");
    f->GetObject("histo",h);     hs.push_back(h);hNs.push_back("KDE w/ expo. tail smoothing");
    f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back("KDE");
    // f->GetObject("histo_PTUp",h);     hs.push_back(h);hNs.push_back("KDESU");
    // f->GetObject("histo_PTDown",h);     hs.push_back(h);hNs.push_back("KDESD");
    // f->GetObject("histo_OPTUp",h);     hs.push_back(h);hNs.push_back("KDERU");
    // f->GetObject("histo_OPTDown",h);     hs.push_back(h);hNs.push_back("KDERD");

    
    //For Hbb
    // for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral(16,125)/hs[iH]->Integral(16,125));
    //For HH
    for(unsigned int iH = 1; iH < hs.size(); ++iH){
      std::cout << hs[0]->Integral(33,201) <<" "<< hs[iH]->Integral(33,201) <<std::endl;
      hs[iH]->Scale(hs[0]->Integral(33,201)/hs[iH]->Integral(33,201));              
    } 
  
    Plotter * p = new Plotter();
    
    for(unsigned int iH = 0; iH < hs.size(); ++iH){
      if(iH == 0) p->addHist(hs[iH],TString::Format("%s",hNs[iH].Data()));
      else  p->addHistLine(hs[iH],TString::Format("%s",hNs[iH].Data()));
    }    
    p->setMinMax(.0001,hs[0]->Integral());
    p->setUnderflow(false);
    p->setOverflow(false);
    p->rebin(5);
    // auto * c = p->drawSplitRatio(1);    
        auto * c = p->draw();    
    c->SetLogy();
    c->Update();

  }
  
  //Test combined conditional kernel
  {
  // std::vector<double> bins = {30,40,50,60,80,100,120,140,160,180,200,220,250};
  // bool binInY = false;
    std::vector<double> bins = {700,800,900,1000,1500,2000,3000,4000,5000};
      bool binInY = true;
  std::vector<TH2*> hs;
  std::vector<TString> hNs;
  TFile *fY = new TFile("hhLNuJJ_nonResTCOND2D_template.root","read");
  TH2* h = 0; 
  fY->GetObject("histofine_data",h);hs.push_back(h);hNs.push_back("mc");
  fY->GetObject("histo_smooth",h);  hs.push_back(h);hNs.push_back("KDES");
  
  auto cutHistograms =[](const TH2* inH) -> TH2F*{
    TH2F * outH = new TH2F(TString(inH->GetName()) + "_cut","",90,30,210,168,800,5000);
    for(unsigned int iX =1; iX <= inH->GetNbinsX(); ++iX){
      const int outIY =outH->GetYaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(iX));
      if(outIY < 1 || outIY > outH->GetNbinsY() ) continue;
      for(unsigned int iY =1; iY <= inH->GetNbinsY(); ++iY){
        const int outIX = outH->GetXaxis()->FindFixBin(inH->GetYaxis()->GetBinCenter(iY));
        if(outIX < 1 || outIX > outH->GetNbinsX() ) continue;
        outH->SetBinContent(outIX,outIY,inH->GetBinContent(iX,iY));                
        outH->SetBinError(outIX,outIY,inH->GetBinError(iX,iY));
      }      
    }    
    return outH;    
  };
  
  for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH] = cutHistograms(hs[iH]);
  
  TFile * fY2 = new TFile("hhLNuJJ_nonResT_2D.root","read");
  fY2->GetObject("histo",h);hs.push_back(h);hNs.push_back("cond");
  
  for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral()/hs[iH]->Integral()); 
  
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
      if(iH == 0) p->addHist(h1D,TString::Format("%s",hNs[iH].Data()));
      else  p->addHistLine(h1D,TString::Format("%s",hNs[iH].Data()));
    }    
    p->setMinMax(.0001,hs[0]->Integral());
    p->setUnderflow(false);
    p->setOverflow(false);
    // p->rebin(5);
    auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));    
    c->SetLogy();
    c->Update();
  }

}


  //Test fits
  {
    TString type = "W";
    TString sel = "emu_LMT_ltmb";
  // std::vector<double> bins = {30,40,50,60,80,100,120,140,160,180,200,220,250};
  // bool binInY = false;
    std::vector<double> bins = {30,50,100,150,210};
    bool binInY = false;
    // std::vector<double> bins = {700,800,900,1000,1500,2000,3000,4000,5000};
      // bool binInY = true;
  std::vector<TH2*> hs;
  std::vector<TString> hNs;
  TFile *fY = new TFile(TString::Format("hhLNuJJ_nonRes%s_distributions.root",type.Data()),"read");
  TH2* h = 0; 
  fY->GetObject(TString::Format("nonRes%s_%s_hbb_mass_hh_mass",type.Data(),sel.Data()),h);hs.push_back(h);hNs.push_back("mc");
    
  TFile * fY1 = new TFile(TString::Format("hhLNuJJ_nonRes%s_2D.root",type.Data()),"read");
  fY1->GetObject("histo",h);hs.push_back(h);hNs.push_back("cond");
  
  TFile * fY2 = new TFile(TString::Format("nonRes%s_%s_explode.root",type.Data(),sel.Data()),"read");
  fY2->GetObject("histo__x_y",h);hs.push_back(h);hNs.push_back("fit");
  
  for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral()/hs[iH]->Integral()); 
  
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
      if(iH == 0) p->addHist(h1D,TString::Format("%s",hNs[iH].Data()));
      // else  p->addHistLine(h1D,TString::Format("%s",hNs[iH].Data()));
        else {
          for(unsigned int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
          p->addHist(h1D,TString::Format("%s",hNs[iH].Data()),-1,1,4,20,1,false,true,false,"E");
        } 
      
    }    
    p->setMinMax(.0001,hs[0]->Integral());
    p->setUnderflow(false);
    p->setOverflow(false);
    p->rebin(12);
    // auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    // c->SetLogy();
    // c->Update();
    
    p->setBotMinMax(0,2);    
    auto * c = p->drawSplitRatio(0,"stack",false,false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    c->GetPad(1)->SetLogy();
    c->GetPad(1)->Update();
  }

}

  //Test TW fits
  {
    std::vector<TString> types = {"resT"};
    std::vector<TString> typeNs = {"m_{t} bkg. model"};
    TString sel = "mu_M_full";
    // std::vector<double> bins = {30,50,100,150,210};
    // bool binInY = false;
    std::vector<double> bins = {800,900,1000,1500,2000,3000,4000,5000};
      bool binInY = true;
  std::vector<TH2*> hs;
  std::vector<TString> hNs;
  std::vector<TH2*> dhs;  
  
  
  TH2 * dh = 0;  
  for(unsigned int iT = 0; iT < types.size(); ++iT){
    TFile * fY = new TFile(TString::Format("HHlnujj_%s_distributions.root",types[iT].Data()),"read");
    TH2 * h = 0; 
    fY->GetObject(TString::Format("%s_%s_hbbMass_hhMass",types[iT].Data(),sel.Data()),h);
    if(h==0) continue;
    if(dh==0)dh= (TH2*)h->Clone();
    else dh->Add(h);
    dhs.push_back(h);
  }    
  hs.push_back(dh);
  hNs.push_back("MC");  
  
  for(unsigned int iT = 0; iT < types.size(); ++iT){
    TFile * fY = new TFile(TString::Format("HHlnujj_%s_debug_2DTemplate.root",types[iT].Data()),"read");
    TH2 * h = 0; 
    fY->GetObject(TString::Format("%s_%s",types[iT].Data(),sel.Data()),h);
    if(h==0) continue;
    h->Scale(dhs[iT]->Integral()/h->Integral());
    hs.push_back(h);
    hNs.push_back(typeNs[iT]);  
  }    
    
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
      if(iH == 0) p->addHist(h1D,TString::Format("%s",hNs[iH].Data()));
      // else  p->addHistLine(h1D,TString::Format("%s",hNs[iH].Data()));
        else {
          for(unsigned int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
          p->addStackHist(h1D,TString::Format("%s",hNs[iH].Data()));
        } 
      
    }    
    // p->setMinMax(.0001,hs[0]->Integral());
    p->setUnderflow(false);
    p->setOverflow(false);
    p->rebin(5);
    p->setXTitle("Hbb soft-drop mass [GeV]");
    p->setYTitle("N. events");
    auto * c = p->draw(false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    // c->SetLogy();
    // c->Update();
    
    // p->setBotMinMax(0,2);   
    // auto * c = p->drawSplitRatio(-1,"stack",false,false,TString::Format("%.0f-%.0f",bins[iB],bins[iB+1]));
    // c->GetPad(1)->SetLogy();
    // c->GetPad(1)->Update();
  }

}

  //Test 1D fits
  {
    TString type = "resT";

    // std::vector<TString> sels = {"emu_LMT_ltmb","emu_LMT_none","emu_I_full","emu_LMT_full","e_L_ltmb","e_M_ltmb","e_T_ltmb","e_L_full","e_M_full","e_T_full","mu_L_ltmb","mu_M_ltmb","mu_T_ltmb","mu_L_full","mu_M_full","mu_T_full"};
    
    
    std::vector<TString> sels = {"mu_L_full","mu_M_full","mu_T_full"};
    
        // std::vector<TString> sels = {"emu_LMT_lb","emu_L_full","emu_M_full","emu_T_full","e_L_full","e_M_full","e_T_full","mu_L_full","mu_M_full","mu_T_full"};
                // std::vector<TString> sels = {"emu_M_full","emu_M_lb","emu_M_ltmb","emu_T_ltmb","emu_T_full","e_M_full","mu_M_full"};
            // std::vector<TString> sels = {"emu_LMT_lb","e_L_lb","e_M_lb","e_T_lb","mu_L_lb","mu_M_lb","mu_T_lb"};
                // std::vector<TString> sels = {"emu_LMT_full","emu_LMT_lWW","emu_LMT_lb","emu_LMT_lm","emu_LMT_lt","emu_LMT_ltmb","emu_LMT_none"};
    

    TFile * fd = new TFile(TString::Format("HHlnujj_%s_distributions.root",type.Data()),"read");
    TFile * fo = new TFile(TString::Format("HHlnujj_%s_template.root",type.Data()),"read");
    TH1 * hof = 0;
    fo->GetObject("histo",hof);
    
    for(const auto& s : sels){      
      TH1 * hd = 0;
      fd->GetObject(TString::Format("%s_%s_hhMass",type.Data(),s.Data()),hd);
      if(hd==0) continue;
      std::vector<TH1*> hs;
      std::vector<TString> hNs;
      hs.push_back((TH1*)hof->Clone());
      hNs.push_back("Baseline KDE before MC fit");
      
      TFile * ff = new TFile(TString::Format("HHlnujj_%s_%s_fitTemplate.root",type.Data(),s.Data()),"read");
      TH1* h = 0; 
      ff->GetObject("histo",h);hs.push_back(h);hNs.push_back("Nominal SR PDF");
      
      for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(hd->Integral()/hs[iH]->Integral()); 
            
      Plotter * p = new Plotter();
      p->addHist(hd,"SR MC");
      for(unsigned int iH = 0; iH < hs.size(); ++iH){
        TH1 * h1D = hs[iH];      
        for(unsigned int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
        p->addHist(h1D,TString::Format("%s",hNs[iH].Data()),-1,1,4,20,1,false,true,false,"E");
      }
      p->setMinMax(.0001,hs[0]->Integral());
      p->setUnderflow(false);
      p->setOverflow(false);
      p->rebin(8);
      p->setXTitle("HH mass [GeV]");
      p->setYTitle("N. of events");
      p->setBotMinMax(0,2);    
      // auto * c = p->drawSplitRatio(0,"stack",false,false,s);
      // c->GetPad(1)->SetLogy();
      // c->GetPad(1)->Update();
      auto * c = p->draw(false,s);
      c->SetLogy();
      c->Update();  
    }
}

  //Compare distributions;
  {
    TString type = "resT";
    // std::vector<TString> sels = {"emu_LMT","emu_LMT_lWW","emu_LMT_lb","emu_LMT_ltb","emu_LMT_ltmb" };
    std::vector<TString> sels = {"emu_L","emu_L_ltmb","emu_M_ltmb" };

  // std::vector<double> bins = {30,40,50,60,80,100,120,140,160,180,200,220,250};
  // bool binInY = false;
    // std::vector<double> bins = {30,50,100,150,210};
    // bool binInY = false;
    // std::vector<double> bins = {800,1000,2000,3000,5000};
      // bool binInY = true;
  std::vector<double> bins = {800,900,1000,1250,1500,2000,2500,3000,3500,4000};
    bool binInY = true;
  std::vector<TH2*> hs;
  std::vector<TString> hNs;
  TFile *fY = new TFile(TString::Format("hhLNuJJ_%s_distributions.root",type.Data()),"read");
  for(unsigned int iS = 0; iS < sels.size(); ++iS){
    TH2* h = 0; 
    fY->GetObject(TString::Format("%s_%s_hbb_mass_hh_mass",type.Data(),sels[iS].Data()),h);hs.push_back(h);hNs.push_back(sels[iS]);    
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
    p->rebin( binInY ? 4 : 24);
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

std::vector<TString> bkgs = {"nonResAQ","nonResH","resW","resT"};
std::vector<TString> bkgNs = {"q/g bkg.","lost t/W bkg.","m_{W} bkg.","m_{t} bkg."};
std::vector<TString> sels = {"emu_L_ltmb","emu_M_ltmb"};
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
  
