// SKimming trees


      {
  
        std::vector<TString> inFiles = {
      "TnPTree_SingleElectron_2016prompt_RunG.root",
      "TnPTree_SingleElectron_2016prompt_RunB.root",
      "TnPTree_SingleElectron_2016prompt_RunE.root",
      "TnPTree_SingleElectron_2016prompt_RunHv2.root",
      "TnPTree_SingleElectron_2016prompt_RunD.root",
      "TnPTree_SingleElectron_2016prompt_RunF.root",
      "TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DYToLL_madgraph.root",
      "TnPTree_SingleElectron_2016prompt_RunC.root",
      "TnPTree_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
      "TnPTree_SingleElectron_2016prompt_RunHv3.root"
      };
      // TString outDir = "reco_skim";
      // TString treedir = "tnpEleReco";
      // TString sel = "tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && tag_Ele_nonTrigMVA >0.96 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 60";

      TString outDir = "id_skim";
      TString treedir = "tnpEleIDs";
      TString sel = "tag_Ele_pt > 27 && abs(tag_sc_eta) < 2.17 && el_q*tag_Ele_q < 0";

      for(const auto& iF : inFiles){
        TFile * fin = new TFile(iF,"read");
        TTree * tin = 0;
        fin->GetObject(treedir +"/fitter_tree",tin);
        if(tin == 0) return;
        TFile * fout = new TFile(outDir + "/"+iF ,"recreate");
        TDirectory * cdtof = fout->mkdir(treedir);
        cdtof->cd();
        TTree *tout = tin->CopyTree(sel);
        tout->Write();
        delete fin;
        delete fout;
      }
  
      }

  //Eff for reco for signal in signaltriggerkine
//Compare efficiencies w/ new vs old method
  
  
{
TFile * nf = new TFile("results/reco_std_func/passingRECO/egammaEffi.txt_EGM2D.root","read");
TFile * of = new TFile("results/reco_std/passingRECO/egammaEffi.txt_EGM2D.root","read");
std::vector<TString> hnames = {"EGamma_SF2D","EGamma_EffData2D","EGamma_EffMC2D"};
std::vector<TString> ht = {"scale factor","data efficiency","MC efficiency"};

for(unsigned int iV = 0; iV < hnames.size();++iV){
  TH2*hn = 0;
  TH2*ho = 0;
  nf->GetObject(hnames[iV],hn);
  of->GetObject(hnames[iV],ho);
  
  auto hn1 = hn->ProjectionX(hnames[iV]+"_hnx");
  auto ho1 = ho->ProjectionX(hnames[iV]+"_hox");
  for(unsigned int iB = 1; iB <=hn1->GetNbinsX(); ++iB ) hn1->SetBinError(iB,0);
  
  Plotter * p = new Plotter();
  p->addHist(ho1,"method: template fit");
  p->addHist(hn1,"method: function fit");
  if(iV)p->setMinMax(0,1);
  else p->setMinMax(0.5,1.5);
  p->setBotMinMax(0.9,1.1);
  p->setYTitle(ht[iV]);
  p->setXTitle("supercluster #eta");
  p->setYTitleBot("template/function");  
  p->drawSplitRatio(1,"stack",false,false,hnames[iV]);
  
  
  
  
}

}


//Projection
{
TFile * nf = new TFile("results/recoISO/passingRECO/egammaEffi.txt_EGM2D.root","read");
std::vector<TString> hnames = {"EGamma_SF2D","EGamma_EffData2D","EGamma_EffMC2D"};
std::vector<TString> ht = {"scale factor","data efficiency","MC efficiency"};

for(unsigned int iV = 0; iV < hnames.size();++iV){
  TH2*hn = 0;
  nf->GetObject(hnames[iV],hn);
  Plotter * p = new Plotter();
  for(unsigned int iB = 6; iB <= hn->GetNbinsX(); ++iB){
    auto h1 = hn->ProjectionY(hnames[iV]+TString::Format("%u",iB),iB,iB);
    TString title = TString::Format("%.1f < #it{J}_{#it{A}} < %.1f",hn->GetXaxis()->GetBinLowEdge(iB),hn->GetXaxis()->GetBinLowEdge(iB)+hn->GetXaxis()->GetBinWidth(iB));
    p->addHist(h1,title);
    
  }  

  // if(iV)p->setMinMax(0,1);
  // else p->setMinMax(0.5,1.5);
  p->setYTitle(ht[iV]);
  p->setXTitle("#Delta#it{R}_{#it{A}}");
  p->draw(false,hnames[iV]);
  // p->drawSplitRatio(1,"stack",false,false,hnames[iV]); 
}
}



//Understand reconstruction inefficiency
{
  vector<unsigned int> fullSigMasses = {
    1000,
    2000,
    3000,
    4000
  };
  vector<TString> ecNs = {
    "gen electron p_{T}>30 GeV, |#eta|<2.5",
    "reconstructed photon",
    "photon H/E <0.15",
    "reconstructed electron",
    "ECAL and tracker seed",
    "ECAL seed, not tracker seed",
    "tracker seed, not ECAL seed",
    "not ECAL or tracker seed",
    "H/E scale + 5%, #gamma H/E <0.15"    
  };

  auto effPlots = [&](TString name, TString prefix, std::vector<unsigned int>& cuts){
      std::vector<TH1*> hists;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        hists.push_back(new TH1F(TString::Format("%s_%u",name.Data(),iN),";#it{m}(X) [TeV]",40,0.550,4.550));
      }

      for(auto m : fullSigMasses){
        TFile * f = new TFile(TString::Format("out_radion_hh_bbinc_m%u_0.root",m),"read");
        TH1 * h = 0;
        f->GetObject(TString::Format("%sevent_count",prefix.Data()),h);
        if(h==0){
          std::cout << TString::Format("%sevent_count",prefix.Data())<<endl;
          continue;
        }
        for(unsigned int iC = 0; iC < cuts.size(); ++iC){
          hists[iC]->SetBinContent(hists[iC]->FindFixBin(float(m)/1000.),h->GetBinContent(cuts[iC]+1));
        }
      }
      Plotter * p = new Plotter;
      for(unsigned int iN = 0; iN < cuts.size(); ++iN){
        p->addHist(hists[iN],ecNs[cuts[iN]],-1,1,4,20,1,true,false);
      }
      p->setMinMax(0,1.0);
      p->drawRatio(false,"stack",false,false,name);

            // p->draw(false);
};

std::vector<TString> vars = {"genLepWDR","genlep_pt","ht","genW_pt"};
std::vector<TString> pres = {"all_incl"};

std::vector<unsigned int> cuts = {0,1,2,3};
effPlots("efficiency","",cuts);

std::vector<unsigned int> cuts2 = {2,8};
effPlots("systTest","",cuts2);

std::vector<unsigned int> cuts3 = {3,4,5,6,7};
effPlots("breakdown","",cuts3);

}

//Comparison of electron res.
{
  
    vector<unsigned int> fullSigMasses = {
      1000,
      2000,
      3000,
      4000
    };
    Plotter * p = new Plotter;
    for(unsigned int iM = 0; iM < fullSigMasses.size(); ++iM){
    
    TFile * f = new TFile(TString::Format("out_radion_hh_bbinc_m%u_0.root",fullSigMasses[iM]),"read");
    TH1 * h = 0;
    f->GetObject("energryRatio",h);
    if(h==0) continue;
    p->addHist(h,TString::Format("#it{m}_{X}=%0.1fTeV",float(fullSigMasses[iM])/1000.) );
  }
  p->normalize();
  p->draw();
    
}

// Run following as a seperate file:  rr fit.C
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

void go(){
  
  vector<unsigned int> fullSigMasses = {
    1000,
    2000,
    3000,
    4000
  };
  
  Plotter* p = new Plotter();
  std::vector<TObject*> plts;
  for(unsigned int iM = 0; iM < fullSigMasses.size(); ++iM){
    TFile * f = new TFile(TString::Format("out_radion_hh_bbinc_m%u_0.root",fullSigMasses[iM]),"read");
    TH1 * h = 0;
    f->GetObject("energryRatio",h);
    if(h==0) continue;
    auto w = new RooWorkspace("w",false);
    w->factory("res[1,.8,1.2]");    
    w->var("res")->setMin(0.8);
    w->var("res")->setMax(1.2);
    w->var("res")->setBins(100);
    w->factory("mean[1,0,2]");
    w->factory("sigma[.1,0,.5]");
    w->factory("alpha[1,0.001,100]");
    w->factory("alpha2[1,0.001,100]");
    w->factory("n[5,0,100]");
    w->factory("n2[5,0,100]");
    w->var("n")->setConstant(true);
    w->var("n2")->setConstant(true);

    RooDataHist dataHist("data","data",RooArgList(*w->var("res")),h);
    w->import(dataHist);
    RooDoubleCB model("model","model",*w->var("res")  ,
           *w->var("mean"),*w->var("sigma"),*w->var("alpha"),*w->var("n") ,*w->var("alpha2"),*w->var("n2"));
    w->import(model);
    w->pdf("model")->fitTo(*w->data("data"));
    auto frame = w->var("res")->frame();
    w->data("data")->plotOn(frame);    
    w->pdf("model")->plotOn(frame);
    // w->pdf("model")->paramOn(frame);
    
    auto c = new TCanvas();
    frame->Draw();
    frame->GetXaxis()->SetTitle("electron resolution");
    frame->GetYaxis()->SetTitle("a.u.");
    
    TPaveText *text = new TPaveText(.15,.8,.45,.9,"NDC");
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->AddText(TString::Format("#it{m}_{X}=%0.1fTeV",float(fullSigMasses[iM])/1000.) );
  text->AddText(TString::Format("mean = %.4f #pm %.4f",                
                w->var("mean")->getVal(),
                w->var("mean")->getError() ) );
  c->cd();
  text->Draw();        
    c->Update();
    plts.push_back(c);
    
  }
  
  Drawing::drawAll(plts,"all");
  
}

void fit(){
  go();
}

//