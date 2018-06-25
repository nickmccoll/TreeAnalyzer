data_cr_hhMass->Rebin(4);
other_cr_hhMass->Rebin(4);
ttbar_cr_hhMass->Rebin(4);

data_cr_hhMass->Add(other_cr_hhMass,-1.0);
data_cr_hhMass->Divide(ttbar_cr_hhMass);
data_cr_hhMass->Draw();

hadd  -f HHlnujj_ttbarSF_all_inputPlots.root HHlnujj_ttbarSF_data_inputPlots.root HHlnujj_ttbarSF_mc_inputPlots.root


{

  TFile * f = new TFile("HHlnujj_ttbarSF_all_testPlots.root","read");
  TH1* hd;
  TH1*ho;
  TH1*ht;
  f->GetObject("data_cr_hhMass",hd);
  f->GetObject("other_cr_hhMass",ho);
  f->GetObject("ttbar_cr_hhMass",ht);

  std::vector<TString> plots ={"cr_hhMass","cr_hbbMass","noHbb_hhMass","noHbb_hbbMass","fullHbbHH_hhMass","fullHbbHH_hbbMass"};
  std::vector<TObject*> objects;
  for(unsigned int iP = 0; iP < plots.size(); ++iP){
    TH1* hd;
    TH1*ho;
    TH1*ht;
    TH1*hmv;
    f->GetObject( "data_"+plots[iP],hd);
    f->GetObject("other_"+plots[iP],ho);
    f->GetObject("ttbar_noSF_"+plots[iP],ht);

    Plotter * p = new Plotter();
    p->addHist(hd,"data");
    p->addStackHist(ht,"t#bar{bar} MC");
    p->addStackHist(ho,"non-t#bar{t} MC");
    p->rebin(5);
    p->setBotMinMax(0,2);
    p->setYTitle(plots[iP]);
    // p->normalize();
      // auto * c = p->drawSplitRatio(-1,"stack",false,false,plots[iP]);
      // objects.push_back(c);
    p->setYTitle("N. of events / 10 GeV");
    p->setCMSLumi();
    p->draw(true,"ttbarSFregion_"+plots[iP]+".pdf");
  }
  // Drawing::drawAll(objects,"all");


}


{
  TFile * f = new TFile("HHlnujj_ttbarSF_all_inputPlots.root","read");
    // TFile * f = new TFile("HHlnujj_ttbarSF_all_testPlots.root","read");
  // std::string selection = "cr_hhMass";
  TString selection = "cr_lmass_hhMass";
  TH1* hd;
  TH1*ho;
  TH1*ht;
  f->GetObject("data_"+selection,hd);
  f->GetObject("other_"+selection,ho);
  f->GetObject("ttbar_"+selection,ht);

  double bins[]={0,500,600,700,800,900,1000,1100,1200,1400,1600,2000,2500,3000,4000,5000};
  int nBins = 15;


  TH1 * hd_c = new TH1F("cut_data",";#it{m}_{HH} [GeV]",nBins,bins);
  TH1 * ho_c = new TH1F("cut_other",";#it{m}_{HH} [GeV]",nBins,bins);
  TH1 * ht_c = new TH1F("cut_ttbar",";#it{m}_{HH} [GeV]",nBins,bins);

  hd_c->Sumw2(true);
  ho_c->Sumw2(true);
  ht_c->Sumw2(true);
  TH1 * h_tot = new TH1F("tot",";#it{m}_{HH} [GeV]",nBins,bins);
  TH1 * h_avg = new TH1F("avg",";#it{m}_{HH} [GeV]",nBins,bins);

        auto fillHistogram =[](TH1 * inH, TH1* outH, TH1 * avgH = 0, TH1 * totH=0){
            const int nOutBins = outH->GetNbinsX();
            for(int iB = 1; iB <=inH->GetNbinsX(); ++iB ){
                const float val = inH->GetBinContent(iB);
                const float err = inH->GetBinError(iB);
                const float xval = inH->GetBinCenter(iB);
                const int oB = outH->FindFixBin(xval);
                if(oB < 1 || oB > nOutBins) continue;
                const float outV = outH->GetBinContent(oB);
                outH->SetBinContent(oB,outH->GetBinContent(oB) + val);
                (*outH->GetSumw2())[oB]+=err*err;

                if(avgH) {
                    avgH->SetBinContent(oB,avgH->GetBinContent(oB) + xval*val);
                    totH->SetBinContent(oB,totH->GetBinContent(oB) + val);
                }
            }
        };

        fillHistogram(hd,hd_c,h_avg,h_tot);
        fillHistogram(ht,ht_c);
        fillHistogram(ho,ho_c);


        TGraphErrors * rat = new TGraphErrors();
        TGraphErrors * ratM50 = new TGraphErrors();
        TGraphErrors * ratP50 = new TGraphErrors();
        int curPt = 0;
        int curPtM50 = 0;
        int curPtP50 = 0;



        for(int iB = 1; iB <= hd_c->GetNbinsX(); ++iB){
            if(ht_c->GetBinContent(iB) <= 0) continue;
            if(hd_c->GetBinContent(iB) <= 0) continue;

            double nD = hd_c->GetBinContent(iB);
            double nT = ht_c->GetBinContent(iB);
            double nO = ho_c->GetBinContent(iB);

            double nDE = hd_c->GetBinError(iB);
            double nTE = ht_c->GetBinError(iB);
            double nOE = ho_c->GetBinError(iB);

            auto addPoint =[&](TGraphErrors *g, double otherSF, int& cP ){
              double xV = h_avg->GetBinContent(iB)/h_tot->GetBinContent(iB);
              double yV = 0;
              double yE = 1;

              if(nO*otherSF<nD){
                yV = (nD-otherSF*nO)/nT;
                yE = std::sqrt(nDE*nDE+otherSF*otherSF*nOE*nOE+yV*yV*nTE*nTE)/nT;
              }

              g->SetPoint(curPt,xV,yV);
              g->SetPointError(curPt,0.0,yE);
              cP++;
            };

            addPoint(rat,1.0,curPt);
            addPoint(ratM50,0.5,curPtM50);
            addPoint(ratP50,1.5,curPtP50);
        }

        TCanvas * c = new TCanvas();
        rat->Draw("APE");
        ratP50->SetLineColor(857);
        ratP50->Draw("SAMELP");
        ratM50->SetLineColor(634);
        ratM50->Draw("SAMELP");

}

{

  TFile * f = new TFile("HHlnujj_ttbarSF_all_testPlots.root","read");
  TH1* hd;
  TH1*ho;
  TH1*ht;
  f->GetObject("data_cr_hhMass",hd);
  f->GetObject("other_cr_hhMass",ho);
  f->GetObject("ttbar_cr_hhMass",ht);

  std::vector<TString> plots ={"cr_hhMass","cr_hbbMass","noHbb_hhMass","noHbb_hbbMass","fullHbbHH_hhMass","fullHbbHH_hbbMass"};
  std::vector<TObject*> objects;
  for(unsigned int iP = 0; iP < plots.size(); ++iP){
    TH1* hd;
    TH1*ho;
    TH1*ht;
    TH1*hmv;
    f->GetObject( "data_"+plots[iP],hd);
    f->GetObject("other_"+plots[iP],ho);
    f->GetObject("ttbar_"+plots[iP],ht);
    f->GetObject("mc_noSF_"+plots[iP],hmv);

    Plotter * p = new Plotter();
    p->addHist(hd,"data");
    p->addHistLine(hmv,"vanMC");
    p->addStackHist(ht,"ttbar");
    p->addStackHist(ho,"other");
    p->rebin(5);
    p->setBotMinMax(0,2);
    p->setYTitle(plots[iP]);
    // p->normalize();
      auto * c = p->drawSplitRatio(-1,"stack",false,false,plots[iP]);
      objects.push_back(c);
  }
  Drawing::drawAll(objects,"all");


}
