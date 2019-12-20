#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

double getEffOR(const double eff1, const double eff2) {
	return (eff1 + eff2 - eff1*eff2);
}

TH2 * makeTriggerEffs2D(TString name, TH1 *h1, vector<double>& bins1, TH1 *h2, vector<double>& bins2) {

	static const int nbins1 = bins1.size() - 1;
	static const int nbins2 = bins2.size() - 1;

	TH2 *hh = new TH2D(name,";H_{T};p_{T}",nbins1,&bins1[0],nbins2,&bins2[0]);
	for(unsigned int i1=1; i1<=nbins1; i1++) for(unsigned int i2=1; i2<=nbins2; i2++) {
		double eff = getEffOR(h1->GetBinContent(i1),h2->GetBinContent(i2));
		hh->SetBinContent(i1,i2,eff);
	}

	return hh;

}

TH1 * getEff1D(TH1* num, TH1* den, vector<double> binedges) {

	int nbins = binedges.size() - 1;

	PlotTools::toOverflow(num);
	PlotTools::toUnderflow(num);
	PlotTools::toOverflow(den);
	PlotTools::toUnderflow(den);

	num = PlotTools::rebin(num,nbins,&binedges[0]);
	den = PlotTools::rebin(den,nbins,&binedges[0]);

	PlotTools::toOverflow(num);
	PlotTools::toUnderflow(num);
	PlotTools::toOverflow(den);
	PlotTools::toUnderflow(den);

	TH1 *rat = (TH1*)num->Clone("ratio");
	rat->Divide(num,den,1,1,"b");
	return rat;
}

TH1 * getTotHadHistMC(TFile *f, TString hS) {
	vector<TString> samps = {"ttbar1_","wjets_","singlet_"};
	TH1 *h = (TH1*)f->Get(samps[0]+hS);
	h = (TH1*)h->Clone("mc_"+hS);
	for (unsigned int i=1; i<samps.size(); ++i) {
		h->Add((TH1*)f->Get(samps[i]+hS),1);
	}
	return h;
}

TH2 *getCombPt2HistWithHT(TString name, double ht, TH2 *hist) {
	TH2 *hh = (TH2*)hist->Clone(name);
	for(unsigned int i1=1; i1<=hist->GetNbinsX(); ++i1) for(unsigned int i2=1; i2<=hist->GetNbinsY(); ++i2) {
		double cont = hist->GetBinContent(i1,i2);
		hh->SetBinContent(i1,i2,getEffOR(ht,cont));
	}
	return hh;
}

vector<TH1*> getEffsAndSFs_lnuqq(bool onlyForFile, int year, TString elptS, TString muptS) {
	vector<TH1*> hists;
	TString yS = TString::Format("%d",year);
	TString fPre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/";
	TString wtS = "";
	if(year == 2017) wtS = "lumiwt_";

	TFile *fm = new TFile(fPre+"triggerInfo_mc_"+yS+"_toppt.root");
	TFile *fd = new TFile(fPre+"triggerInfo_data_"+yS+"_toppt.root");

	vector<double> htbins = {100,200,250,300,350,400,450,500,550,600,650,
			700,800,900,1000,1100,1200,2000};

	TString elDenS_1 = "id1_"+yS+"_el_tt2_passSMu_";
	TString elNumS_1 = "id1_"+yS+"_el_tt2_passSMuAndFull_";
	TString muDenS_1 = "id1_"+yS+"_mu_tt2_passSEl_";
	TString muNumS_1 = "id1_"+yS+"_mu_tt2_passSElAndFull_";

	TH1 *hen_mc = (TH1*)fm->Get("ttbar2_"+elNumS_1+wtS+elptS+"_ht");
	TH1 *hed_mc = (TH1*)fm->Get("ttbar2_"+elDenS_1+elptS+"_ht");
	TH1 *hen_da = (TH1*)fd->Get("SingleMuon_"+elNumS_1+elptS+"_ht");
	TH1 *hed_da = (TH1*)fd->Get("SingleMuon_"+elDenS_1+elptS+"_ht");

	TH1 *hmn_mc = (TH1*)fm->Get("ttbar2_"+muNumS_1+wtS+muptS+"_ht");
	TH1 *hmd_mc = (TH1*)fm->Get("ttbar2_"+muDenS_1+muptS+"_ht");
	TH1 *hmn_da, *hmd_da;
	if(year == 2016 || year == 2017) {
		hmn_da = (TH1*)fd->Get("SingleElectron_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("SingleElectron_"+muDenS_1+muptS+"_ht");
	} else if (year == 2018) {
		hmn_da = (TH1*)fd->Get("EGamma_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("EGamma_"+muDenS_1+muptS+"_ht");
	}

	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins);
	elEffMC = (TH1*)elEffMC->Clone("effMC_e_"+elptS);
	elEffMC->GetXaxis()->SetTitle("H_{T}");

	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins);
	muEffMC = (TH1*)muEffMC->Clone("effMC_m_"+muptS);
	muEffMC->GetXaxis()->SetTitle("H_{T}");

	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins);
	elEffDA = (TH1*)elEffDA->Clone("effDA_e_"+elptS);
	elEffDA->GetXaxis()->SetTitle("H_{T}");

	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins);
	muEffDA = (TH1*)muEffDA->Clone("effDA_m_"+muptS);
	muEffDA->GetXaxis()->SetTitle("H_{T}");

	if (!onlyForFile) {
		hists.push_back(elEffMC);
		hists.push_back(muEffMC);
		hists.push_back(elEffDA);
		hists.push_back(muEffDA);
	}

	elEffDA = (TH1*)elEffDA->Clone("sf_e_"+elptS);
	elEffDA->Divide(elEffMC);
	muEffDA = (TH1*)muEffDA->Clone("sf_m_"+muptS);
	muEffDA->Divide(muEffMC);

	if(onlyForFile) {
		elEffDA = (TH1*)elEffDA->Clone("electronSF_lnuqq");
		muEffDA = (TH1*)muEffDA->Clone("muonSF_lnuqq");
		hists.push_back(elEffDA);
		hists.push_back(muEffDA);
		return hists;
	}

	hists.push_back(elEffDA);
	hists.push_back(muEffDA);

	return hists;
}

vector<TH2*> getEffsAndSFs_lnulnu(int year, TString elptS, TString muptS, bool includeCross) {
	vector<TH2*> hists2D;
	TString yS = TString::Format("%d",year);
	TString fPre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/";
	TString wtS = "";
	if(year == 2017) wtS = "lumiwt_";

	TFile *fm = new TFile(fPre+"triggerInfo_mc_"+yS+"_toppt.root");
	TFile *fd = new TFile(fPre+"triggerInfo_data_"+yS+"_toppt.root");

	vector<double> htbins = {100,200,250,300,350,400,450,500,550,600,650,
			700,800,900,1000,1100,1200,2000};
	vector<double> ptbins = {5,10,15,20,25,30,35,40,50,60,70,80,90,100,150,200,300,1000};

	//  get trigger efficiencies for first lepton and non-lepton triggers in di-lepton channel -----------------------------
	TString elDenS_1 = "id2_"+yS+"_el_tt2_passSMu_";
	TString elNumS_1 = "id2_"+yS+"_el_tt2_passSMuAndFull_";
	TString muDenS_1 = "id2_"+yS+"_mu_tt2_passSEl_";
	TString muNumS_1 = "id2_"+yS+"_mu_tt2_passSElAndFull_";

	TH1 *hen_mc = (TH1*)fm->Get("ttbar2_"+elNumS_1+wtS+elptS+"_ht");
	TH1 *hed_mc = (TH1*)fm->Get("ttbar2_"+elDenS_1+elptS+"_ht");
	TH1 *hen_da = (TH1*)fd->Get("SingleMuon_"+elNumS_1+elptS+"_ht");
	TH1 *hed_da = (TH1*)fd->Get("SingleMuon_"+elDenS_1+elptS+"_ht");

	TH1 *hmn_mc = (TH1*)fm->Get("ttbar2_"+muNumS_1+wtS+muptS+"_ht");
	TH1 *hmd_mc = (TH1*)fm->Get("ttbar2_"+muDenS_1+muptS+"_ht");
	TH1 *hmn_da, *hmd_da;
	if(year == 2016 || year == 2017) {
		hmn_da = (TH1*)fd->Get("SingleElectron_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("SingleElectron_"+muDenS_1+muptS+"_ht");
	} else if (year == 2018) {
		hmn_da = (TH1*)fd->Get("EGamma_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("EGamma_"+muDenS_1+muptS+"_ht");
	}

	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins);
	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins);
	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins);
	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins);

//  get single-lepton trigger efficiencies to be used for second lepton -----------------------------
	TString crossS = (includeCross ? "Cross" : "");
	TString cS = includeCross ? (year == 2017 ? "Cross_lumiwt" : "Cross") : "";

	TString elDenS  = "id2_"+yS+"_el_tt2_passSMu_";
	TString elNumS  = "id2_"+yS+"_el_tt2_passSMuAndSEl"+crossS+"_";
	TString muDenS  = "id2_"+yS+"_mu_tt2_passSEl_";
	TString muNumS  = "id2_"+yS+"_mu_tt2_passSElAndSMu"+crossS+"_";

	hen_mc = (TH1*)fm->Get("ttbar2_"+elNumS+(includeCross ? wtS.Data() : "")+"ht400_pt");
	hed_mc = (TH1*)fm->Get("ttbar2_"+elDenS+"ht400_pt");
	hmn_mc = (TH1*)fm->Get("ttbar2_"+muNumS+(includeCross ? wtS.Data() : "")+"ht400_pt");
	hmd_mc = (TH1*)fm->Get("ttbar2_"+muDenS+"ht400_pt");

	hen_da = (TH1*)fd->Get("SingleMuon_"+elNumS+"ht400_pt");
	hed_da = (TH1*)fd->Get("SingleMuon_"+elDenS+"ht400_pt");
	if(year == 2016 || year == 2017) {
		hmn_da = (TH1*)fd->Get("SingleElectron_"+muNumS+"ht400_pt");
		hmd_da = (TH1*)fd->Get("SingleElectron_"+muDenS+"ht400_pt");
	} else if (year == 2018) {
		hmn_da = (TH1*)fd->Get("EGamma_"+muNumS+"ht400_pt");
		hmd_da = (TH1*)fd->Get("EGamma_"+muDenS+"ht400_pt");
	}

	TH1 *el2EffMC = getEff1D(hen_mc,hed_mc,ptbins);
	TH1 *mu2EffMC = getEff1D(hmn_mc,hmd_mc,ptbins);
	TH1 *el2EffDA = getEff1D(hen_da,hed_da,ptbins);
	TH1 *mu2EffDA = getEff1D(hmn_da,hmd_da,ptbins);

	//  get 2D effs and scale factors for di-lepton channel
	TH2 *effMC_ee = makeTriggerEffs2D("effMC_ee_"+elptS+crossS,elEffMC,htbins,el2EffMC,ptbins); hists2D.push_back(effMC_ee);
	TH2 *effMC_me = makeTriggerEffs2D("effMC_me_"+muptS+crossS,muEffMC,htbins,el2EffMC,ptbins); hists2D.push_back(effMC_me);
	TH2 *effMC_em = makeTriggerEffs2D("effMC_em_"+elptS+crossS,elEffMC,htbins,mu2EffMC,ptbins); hists2D.push_back(effMC_em);
	TH2 *effMC_mm = makeTriggerEffs2D("effMC_mm_"+muptS+crossS,muEffMC,htbins,mu2EffMC,ptbins); hists2D.push_back(effMC_mm);

	TH2 *effDA_ee = makeTriggerEffs2D("effDA_ee_"+elptS+crossS,elEffDA,htbins,el2EffDA,ptbins); hists2D.push_back(effDA_ee);
	TH2 *effDA_me = makeTriggerEffs2D("effDA_me_"+muptS+crossS,muEffDA,htbins,el2EffDA,ptbins); hists2D.push_back(effDA_me);
	TH2 *effDA_em = makeTriggerEffs2D("effDA_em_"+elptS+crossS,elEffDA,htbins,mu2EffDA,ptbins); hists2D.push_back(effDA_em);
	TH2 *effDA_mm = makeTriggerEffs2D("effDA_mm_"+muptS+crossS,muEffDA,htbins,mu2EffDA,ptbins); hists2D.push_back(effDA_mm);

	TH2 *sfee = (TH2*)effDA_ee->Clone("sf_ee_"+elptS+crossS);
	TH2 *sfem = (TH2*)effDA_em->Clone("sf_em_"+elptS+crossS);
	TH2 *sfme = (TH2*)effDA_me->Clone("sf_me_"+muptS+crossS);
	TH2 *sfmm = (TH2*)effDA_mm->Clone("sf_mm_"+muptS+crossS);

	sfee->Divide(effMC_ee); hists2D.push_back(sfee);
	sfem->Divide(effMC_em); hists2D.push_back(sfem);
	sfme->Divide(effMC_me); hists2D.push_back(sfme);
	sfmm->Divide(effMC_mm); hists2D.push_back(sfmm);

	return hists2D;
}

vector<TH1*> getEffsForFile_2l(int year, TString elptS, TString muptS, bool includeCross) {
	vector<TH1*> hists;
	TString yS = TString::Format("%d",year);
	TString fPre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/";
	TString wtS = "";
	if(year == 2017) wtS = "lumiwt_";

	TFile *fm = new TFile(fPre+"triggerInfo_mc_"+yS+"_toppt.root");
	TFile *fd = new TFile(fPre+"triggerInfo_data_"+yS+"_toppt.root");

	vector<double> htbins = {100,200,250,300,350,400,450,500,550,600,650,
			700,800,900,1000,1100,1200,2000};
	vector<double> ptbins = {5,10,15,20,25,30,35,40,50,60,70,80,90,100,150,200,300,1000};

	//  get trigger efficiencies for first lepton and non-lepton triggers in di-lepton channel -----------------------------
	TString elDenS_1 = "id2_"+yS+"_el_tt2_passSMu_";
	TString elNumS_1 = "id2_"+yS+"_el_tt2_passSMuAndFull_";
	TString muDenS_1 = "id2_"+yS+"_mu_tt2_passSEl_";
	TString muNumS_1 = "id2_"+yS+"_mu_tt2_passSElAndFull_";

	TH1 *hen_mc = (TH1*)fm->Get("ttbar2_"+elNumS_1+wtS+elptS+"_ht");
	TH1 *hed_mc = (TH1*)fm->Get("ttbar2_"+elDenS_1+elptS+"_ht");
	TH1 *hen_da = (TH1*)fd->Get("SingleMuon_"+elNumS_1+elptS+"_ht");
	TH1 *hed_da = (TH1*)fd->Get("SingleMuon_"+elDenS_1+elptS+"_ht");

	TH1 *hmn_mc = (TH1*)fm->Get("ttbar2_"+muNumS_1+wtS+muptS+"_ht");
	TH1 *hmd_mc = (TH1*)fm->Get("ttbar2_"+muDenS_1+muptS+"_ht");
	TH1 *hmn_da, *hmd_da;
	if(year == 2016 || year == 2017) {
		hmn_da = (TH1*)fd->Get("SingleElectron_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("SingleElectron_"+muDenS_1+muptS+"_ht");
	} else if (year == 2018) {
		hmn_da = (TH1*)fd->Get("EGamma_"+muNumS_1+muptS+"_ht");
		hmd_da = (TH1*)fd->Get("EGamma_"+muDenS_1+muptS+"_ht");
	}

	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins);
	elEffMC = (TH1*)elEffMC->Clone("electronEffs_mcHT_lnulnu");
	elEffMC->GetXaxis()->SetTitle("H_{T}");

	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins);
	muEffMC = (TH1*)muEffMC->Clone("muonEffs_mcHT_lnulnu");
	muEffMC->GetXaxis()->SetTitle("H_{T}");

	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins);
	elEffDA = (TH1*)elEffDA->Clone("electronEffs_dataHT_lnulnu");
	elEffDA->GetXaxis()->SetTitle("H_{T}");

	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins);
	muEffDA = (TH1*)muEffDA->Clone("muonEffs_dataHT_lnulnu");
	muEffDA->GetXaxis()->SetTitle("H_{T}");


//  get single-lepton trigger efficiencies to be used for second lepton -----------------------------
	TString crossS = (includeCross ? "Cross" : "");

	TString elDenS  = "id2_"+yS+"_el_tt2_passSMu_";
	TString elNumS  = "id2_"+yS+"_el_tt2_passSMuAndSEl"+crossS+"_";
	TString muDenS  = "id2_"+yS+"_mu_tt2_passSEl_";
	TString muNumS  = "id2_"+yS+"_mu_tt2_passSElAndSMu"+crossS+"_";

	hen_mc = (TH1*)fm->Get("ttbar2_"+elNumS+(includeCross ? wtS.Data() : "")+"ht400_pt");
	hed_mc = (TH1*)fm->Get("ttbar2_"+elDenS+"ht400_pt");
	hmn_mc = (TH1*)fm->Get("ttbar2_"+muNumS+(includeCross ? wtS.Data() : "")+"ht400_pt");
	hmd_mc = (TH1*)fm->Get("ttbar2_"+muDenS+"ht400_pt");

	TString cS = includeCross ? (year == 2017 ? "Cross_lumiwt" : "Cross") : "";

	hen_da = (TH1*)fd->Get("SingleMuon_"+elNumS+"ht400_pt");
	hed_da = (TH1*)fd->Get("SingleMuon_"+elDenS+"ht400_pt");
	if(year == 2016 || year == 2017) {
		hmn_da = (TH1*)fd->Get("SingleElectron_"+muNumS+"ht400_pt");
		hmd_da = (TH1*)fd->Get("SingleElectron_"+muDenS+"ht400_pt");
	} else if (year == 2018) {
		hmn_da = (TH1*)fd->Get("EGamma_"+muNumS+"ht400_pt");
		hmd_da = (TH1*)fd->Get("EGamma_"+muDenS+"ht400_pt");
	}

	TH1 *el2EffMC = getEff1D(hen_mc,hed_mc,ptbins);
	el2EffMC = (TH1*)el2EffMC->Clone("electronEffs_mcPT_lnulnu");
	el2EffMC->GetXaxis()->SetTitle("p_{T}");

	TH1 *mu2EffMC = getEff1D(hmn_mc,hmd_mc,ptbins);
	mu2EffMC = (TH1*)mu2EffMC->Clone("muonEffs_mcPT_lnulnu");
	mu2EffMC->GetXaxis()->SetTitle("p_{T}");

	TH1 *el2EffDA = getEff1D(hen_da,hed_da,ptbins);
	el2EffDA = (TH1*)el2EffDA->Clone("electronEffs_dataPT_lnulnu");
	el2EffDA->GetXaxis()->SetTitle("p_{T}");

	TH1 *mu2EffDA = getEff1D(hmn_da,hmd_da,ptbins);
	mu2EffDA = (TH1*)mu2EffDA->Clone("muonEffs_dataPT_lnulnu");
	mu2EffDA->GetXaxis()->SetTitle("p_{T}");

	hists.push_back(elEffMC);
	hists.push_back(el2EffMC);
	hists.push_back(muEffMC);
	hists.push_back(mu2EffMC);
	hists.push_back(elEffDA);
	hists.push_back(el2EffDA);
	hists.push_back(muEffDA);
	hists.push_back(mu2EffDA);

	return hists;
}

void makeTriggerCorrections(int year, bool makePlotsForFileOnly = true) {

	if(makePlotsForFileOnly) {
		vector<TH1*> hists1 = getEffsAndSFs_lnuqq(true,year,"elpt30","mupt27");
		vector<TH1*> hists2 = getEffsForFile_2l(year,"elpt30","mupt27",false);

		TFile *fout = new TFile("triggerSF_"+TString::Format("%d",year)+"_toppt.root","RECREATE");
		fout->cd();
		for(const auto& h : hists1) h->Write();
		for(const auto& h : hists2) h->Write();
		fout->Close();
		delete fout;
		return;
	}

	vector<TH2*> hists_nom       = getEffsAndSFs_lnulnu(year,"elpt30","mupt27",false);
	vector<TH2*> hists_up        = getEffsAndSFs_lnulnu(year,"elpt32","mupt29",false);
	vector<TH2*> hists_down      = getEffsAndSFs_lnulnu(year,"elpt28","mupt25",false);
	vector<TH2*> hists_downdown  = getEffsAndSFs_lnulnu(year,"elpt25","mupt22",false);
	vector<TH2*> hists_upup      = getEffsAndSFs_lnulnu(year,"elpt35","mupt32",false);
	vector<TH2*> hists_upCross   = getEffsAndSFs_lnulnu(year,"elpt32","mupt29",true);
	vector<TH2*> hists_downCross = getEffsAndSFs_lnulnu(year,"elpt28","mupt25",true);
	vector<TH2*> hists_nomCross  = getEffsAndSFs_lnulnu(year,"elpt30","mupt27",true);
	vector<TH2*> hists_upupCross  = getEffsAndSFs_lnulnu(year,"elpt35","mupt32",true);
	vector<TH2*> hists_downdownCross  = getEffsAndSFs_lnulnu(year,"elpt25","mupt22",true);

	vector<TH1*> hist_nom = getEffsAndSFs_lnuqq(false,year,"elpt30","mupt27");
	vector<TH1*> hist_up1 = getEffsAndSFs_lnuqq(false,year,"elpt32","mupt29");
	vector<TH1*> hist_up2 = getEffsAndSFs_lnuqq(false,year,"elpt35","mupt32");
	vector<TH1*> hist_down1 = getEffsAndSFs_lnuqq(false,year,"elpt28","mupt25");
	vector<TH1*> hist_down2 = getEffsAndSFs_lnuqq(false,year,"elpt25","mupt22");

	TFile *fout = new TFile("effsAndSF_"+TString::Format("%d",year)+"_toppt.root","RECREATE");
	fout->cd();
	for(const auto& h : hists_nom)       h->Write();
	for(const auto& h : hists_up)        h->Write();
	for(const auto& h : hists_down)      h->Write();
	for(const auto& h : hists_downdown)      h->Write();
	for(const auto& h : hists_upup)      h->Write();
	for(const auto& h : hists_nomCross)  h->Write();
	for(const auto& h : hists_upCross)   h->Write();
	for(const auto& h : hists_downCross) h->Write();
	for(const auto& h : hists_upupCross) h->Write();
	for(const auto& h : hists_downdownCross) h->Write();

	for(const auto& h : hist_nom) h->Write();
	for(const auto& h : hist_up1) h->Write();
	for(const auto& h : hist_up2) h->Write();
	for(const auto& h : hist_down1) h->Write();
	for(const auto& h : hist_down2) h->Write();

	fout->Close();
	delete fout;

}
