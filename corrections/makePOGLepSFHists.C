//Make POG param scalefactors...muons
{

TFile * fo = new TFile("Muon/muonSF_medID_mini0p2ISO.root","recreate");

TString recoFile = "Muon/Tracking_EfficienciesAndSF_BCDEFGH.root";
TString recoHist = "ratio_eff_eta3_dr030e030_corr";

TString idFile = "Muon/EffAndSF_Medium2016_IP_PtEta.root";
TString idHist = "Medium2016_IP_PtEta/abseta_pt_ratio";

TString isoFile = "Muon/EffAndSF_MiniIso_PtEta.root";
TString isoHist = "MiniIso_PtEta/abseta_pt_ratio";

TString isoActFile = "Muon/EffAndSF_MiniIso_JetAct.root";
TString isoActHist = "MiniIso_JetAct/dR_lepact_ratio";

TFile * fr = new TFile(recoFile,"read");
TGraphAsymmErrors * g = 0;
fr->GetObject(recoHist,g);
fo->cd();
g->Write("reco");

TFile * fid = new TFile(idFile,"read");
TH2*  hid = 0;
fid->GetObject(idHist,hid);
fo->cd();
hid->Write("id");


TFile * fiso = new TFile(isoFile,"read");
TH2*  hiso = 0;
fiso->GetObject(isoHist,hiso);
fo->cd();
hiso->Write("iso");

TFile * fisoact = new TFile(isoActFile,"read");
TH2*  hisoact = 0;
fisoact->GetObject(isoActHist,hisoact);
fo->cd();
hisoact->Write("iso_act");

fo->Close();

}
//Make POG param scalefactors...electrons
//Run in /Users/nmccoll/Dropbox/Work/Projects/HHbbWW/7_13_18_eleSF
{

TFile * fo = new TFile("electronSF_tightID_mini0p1ISO.root","recreate");

TString recoFile = "POGRecoSF.root";
TString recoHist = "EGamma_SF2D";

TString idFile = "results/id/passingTight80XNoIsoIP/egammaEffi.txt_EGM2D.root";
TString idHist = "EGamma_SF2D";

TString isoFile = "results/iso/passingMiniIsoTight/egammaEffi.txt_EGM2D.root";
TString isoHist = "EGamma_SF2D";

TFile * fr = new TFile(recoFile,"read");
TH2*  hr = 0;
fr->GetObject(recoHist,hr);
fo->cd();
hr->Write("reco");

TFile * fid = new TFile(idFile,"read");
TH2*  hid = 0;
fid->GetObject(idHist,hid);
fo->cd();
hid->Write("id");


TFile * fiso = new TFile(isoFile,"read");
TH2*  hiso = 0;
fiso->GetObject(isoHist,hiso);
fo->cd();
hiso->Write("iso");

fo->Close();

}

