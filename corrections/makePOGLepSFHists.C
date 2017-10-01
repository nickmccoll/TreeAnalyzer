//Make POG param scalefactors...muons
{

TFile * fo = new TFile("muonSF_medID_mini0p2ISO_POGParam.root","recreate");

TString recoFile = "Tracking_EfficienciesAndSF_BCDEFGH.root";
TString recoHist = "ratio_eff_eta3_dr030e030_corr";

TString idFile = "EfficienciesAndSF_GH.root";
TString idHist = "MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio";

TString isoFile = "TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root";
TString isoHist = "SF";

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

fo->Close();

}
//Make POG param scalefactors...electrons

{

TFile * fo = new TFile("electronSF_tightID_mini0p1ISO_POGParam.root","recreate");

TString recoFile = "reco_egammaEffi.txt_EGM2D.root";
TString recoHist = "EGamma_SF2D";

TString idFile = "id_egammaEffi.txt_EGM2D.root";
TString idHist = "EGamma_SF2D";

TString isoFile = "suse_scaleFactors.root";
TString isoHist = "MVAVLooseElectronToMini";

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

