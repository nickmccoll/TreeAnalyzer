#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/PlotHelp.h"
#include "Configuration/interface/FillerConstants.h"
#include "TFile.h"
#include "TH1.h"

using namespace std;

TString getTrigName(int k) {
	TString trigName;
	switch (k) {
		case 0: trigName = "PFHT500_PFMET100_PFMHT100_IDTight"                  ; break;
		case 1: trigName = "PFHT700_PFMET85_PFMHT85_IDTight"                    ; break;
		case 2: trigName = "PFHT800_PFMET75_PFMHT75_IDTight"                    ; break;
		case 3: trigName = "AK8PFHT850_TrimMass50"                              ; break;
		case 4: trigName = "AK8PFJet400_TrimMass30"                             ; break;
		case 5: trigName = "AK8PFJet500"                                        ; break;
		case 6: trigName = "PFHT1050"                                           ; break;
		case 7: trigName = "PFMET120_PFMHT120_IDTight"                          ; break;
		case 8: trigName = "PFMET120_PFMHT120_IDTight_PFHT60"                   ; break;
		case 9: trigName = "PFMET140_PFMHT140_IDTight"                          ; break;
		case 10: trigName = "PFMETNoMu120_PFMHTNoMu120_IDTight"                 ; break;
		case 11: trigName = "PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"          ; break;
		case 12: trigName = "PFMETNoMu140_PFMHTNoMu140_IDTight"                 ; break;
		case 13: trigName = "PFMETTypeOne120_PFMHT120_IDTight"                  ; break;
		case 14: trigName = "PFMETTypeOne120_PFMHT120_IDTight_PFHT60"           ; break;
		case 15: trigName = "PFMETTypeOne140_PFMHT140_IDTight"                  ; break;
		case 16: trigName = "Ele115_CaloIdVT_GsfTrkIdT"                         ; break;
		case 17: trigName = "Ele15_IsoVVVL_PFHT450"                             ; break;
		case 18: trigName = "Ele28_eta2p1_WPTight_Gsf_HT150"                    ; break;
		case 19: trigName = "Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"; break;
		case 20: trigName = "Ele32_WPTight_Gsf"                                 ; break;
		case 21: trigName = "Ele32_WPTight_Gsf_L1DoubleEG"                      ; break;
		case 22: trigName = "Ele35_WPTight_Gsf"                                 ; break;
		case 23: trigName = "Ele50_CaloIdVT_GsfTrkIdT_PFJet165"                 ; break;
		case 24: trigName = "IsoMu27"                                           ; break;
		case 25: trigName = "Mu15_IsoVVVL_PFHT450"                              ; break;
		case 26: trigName = "Mu50"                                              ; break;
		case 27: trigName = "Photon200"                                         ; break;
		case 28: trigName = "NTrig"                                             ; break;
	}
	return trigName;
}

void plotIndTriggers(TFile *f, TString prefix, const vector<int>& trigs, TString sel, TString var, TString canName, float rebin = 0) {
	Plotter *p = new Plotter();
	TH1F *hd = (TH1F*)f->Get(prefix+"_TrigIncl__"+sel+"_"+var);
    if (hd==0) {
    	cout<<"bad den: "<<prefix+"_TrigIncl__"+sel+"_"+var<<endl;
    	return;
    }
    hd = (TH1F*)hd->Clone();
    PlotTools::toOverflow(hd);
    if (rebin > 0) PlotTools::rebin(hd,rebin);
    p->addHist(hd,"incl",-1,1,4,20,1,true,false);

    for (const auto& trg : trigs) {
    	TH1F *hn = (TH1F*)f->Get(TString::Format("%s_passTrig_%i__%s_%s",prefix.Data(),trg,sel.Data(),var.Data()));
	    if (hn==0) {
	    	cout<<"bad num: "<<TString::Format("%s_passTrig_%i__%s_%s",prefix.Data(),trg,sel.Data(),var.Data())<<endl;
	    	continue;
	    }
	    hn = (TH1F*)hn->Clone();
	    PlotTools::toOverflow(hn);
	    if (rebin > 0) PlotTools::rebin(hn,rebin);
	    p->addHist(hn,getTrigName(trg),-1,1,4,20,1,true,false);
    }
    p->drawRatio(0,prefix+sel+var,false,false,canName);
}

void plotTurnons_Sel(TFile *f, TString prefix, const vector<TString>& sels, TString trig, TString var, TString canName, float rebin = 0, int nR = 0, double *rebins = 0) {
	Plotter *p = new Plotter();
	for (const auto& s : sels) {
		TH1 *hd = 0;
	    f->GetObject(TString::Format("%s__%s_%s",prefix.Data(),s.Data(),var.Data()),hd);
	    TH1 *hn = 0;
	    f->GetObject(TString::Format("%s_%s__%s_%s",prefix.Data(),trig.Data(),s.Data(),var.Data()),hn);

	    if (hn==0) {
	    	cout<<"bad num: "<<TString::Format("%s_%s__%s_%s",prefix.Data(),trig.Data(),s.Data(),var.Data())<<endl;
	    	continue;
	    }
	    if (hd==0) {
	    	cout<<"bad den: "<<TString::Format("%s__%s_%s",prefix.Data(),s.Data(),var.Data())<<endl;
	    	continue;
	    }
	    hn = (TH1*)hn->Clone();
	    hd = (TH1*)hd->Clone();
	    PlotTools::toOverflow(hn);
	    PlotTools::toOverflow(hd);
	    PlotTools::toUnderflow(hn);
	    PlotTools::toUnderflow(hd);
	    if(rebin > 0){
	      PlotTools::rebin(hn,rebin);
	      PlotTools::rebin(hd,rebin);
	    } else if(rebins){
	      hn = PlotTools::rebin(hn,nR,rebins);
	      hd = PlotTools::rebin(hd,nR,rebins);
	    }
//	    auto *ratio = PlotTools::getBinomErrors(hn,hd);
//	    p->addGraph(ratio,s);
	    hn->Divide(hn,hd,1,1,"b");
	    p->addHist(hn,s,-1,1,4,20,1,true,false);
	}
//	p->draw(false,canName);
	p->drawSplitRatio(0,"stack",false,false,canName);
}

void plotTurnons_Trig(TFile *f, TString prefix, const TString sel, const vector<TString>& trigs, TString var, TString canName, float rebin = 0, int nR = 0, double *rebins = 0) {
	Plotter *p = new Plotter();
	for (const auto& trg : trigs) {
		TH1 *hd = 0;
	    f->GetObject(TString::Format("%s__%s_%s",prefix.Data(),sel.Data(),var.Data()),hd);
	    TH1 *hn = 0;
	    f->GetObject(TString::Format("%s_%s__%s_%s",prefix.Data(),trg.Data(),sel.Data(),var.Data()),hn);

	    if (hn==0) {
	    	cout<<"bad num: "<<TString::Format("%s_%s__%s_%s",prefix.Data(),trg.Data(),sel.Data(),var.Data())<<endl;
	    	continue;
	    }
	    if (hd==0) {
	    	cout<<"bad den: "<<TString::Format("%s__%s_%s",prefix.Data(),sel.Data(),var.Data())<<endl;
	    	continue;
	    }
	    hn = (TH1*)hn->Clone();
	    hd = (TH1*)hd->Clone();
	    PlotTools::toOverflow(hn);
	    PlotTools::toOverflow(hd);
	    if(rebin > 0){
	      PlotTools::rebin(hn,rebin);
	      PlotTools::rebin(hd,rebin);
	    } else if(rebins){
	      hn = PlotTools::rebin(hn,nR,rebins);
	      hd = PlotTools::rebin(hd,nR,rebins);
	    }
//	    auto *ratio = PlotTools::getBinomErrors(hn,hd);
//	    p->addGraph(ratio,s);
	    hn->Divide(hn,hd,1,1,"b");
	    p->addHist(hn,trg,-1,1,4,20,1,true,false);
	}
	p->draw(false,canName);
}

void getSF(TFile *fd, TFile *fmc, TString dataPre, TString mcPre, TString effDen, TString trigSel, TString lepSel, TString htSel, TString canName) {

	int nLepBins = 33;
	double lepBins[] = {5,10,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,50,75,100,150,200,250,300,350,400,450,500};
	int nHTBins = 28;
	double htBins[] = {100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,700,800,900,1000,1100,1200,1600,2000};

	auto getEff = [&](TFile * f, TString pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
	    TH1 * hd = 0;
	    f->GetObject(TString::Format("%s_%s__%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()),hd);
	    TH1 * hn = 0;
	    f->GetObject(TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn);
	    if(hn == 0){
	      cout << "Bad num: " + TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
	      return 0;
	    }
	    if(hd == 0){
	      cout << "Bad den: " + TString::Format("%s_%s__%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
	      return 0;
	    }
	    hn = (TH1*)hn->Clone();
	    hd = (TH1*)hd->Clone();
	    PlotTools::toOverflow(hn);
	    PlotTools::toOverflow(hd);

	    if(rebin > 0){
	      PlotTools::rebin(hn,rebin);
	      PlotTools::rebin(hd,rebin);
	    } else if(rebins){
	      hn = PlotTools::rebin(hn,nR,rebins);
	      hd = PlotTools::rebin(hd,nR,rebins);
	    }
//	    PlotTools::toOverflow(hn);
//	    PlotTools::toOverflow(hd);
	    hn->Divide(hn,hd,1,1,"b"); return hn;
	    // return PlotTools::getBinomErrors(hn,hd);
	};

	auto plotSFTurnons =[&](TString name, TString prefix,TString dataName,TString mcName, TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0 ){
	    Plotter * p = new Plotter();
	    auto * mcEff = getEff(fmc,mcName,prefix,sel,var,trig,rebin,nR,rebins);
	    auto * dataEff = getEff(fd,dataName,prefix,sel,var,trig,rebin,nR,rebins);
	    if(mcEff == 0 || dataEff == 0) {
	    	cout << "mc or data hist is empty" << endl;
	    	return;
	    }

	    p->addHist(mcEff,mcPre,-1,1,4,20,1,true,true, false);
	    p->addHist(dataEff,dataPre,-1,1,4,20,1,true,true, false);
	      // p->draw(true,TString::Format("%s.pdf",name.Data()));
	    p->drawSplitRatio(0,"stack",true,false,TString::Format("%s.pdf",name.Data()));
//	    p->draw(false,TString::Format("%s.pdf",name.Data()));
	};

	TString lepVar = (lepSel.BeginsWith("m")) ? "mu_pt" : "el_pt";

	plotSFTurnons("SF_"+canName+"_HTpup",effDen,dataPre,mcPre,lepSel,"ht_pup",trigSel,0,nHTBins,htBins);
	plotSFTurnons("SF_"+canName+"_HT",effDen,dataPre,mcPre,lepSel,"ht",trigSel,0,nHTBins,htBins);
//	plotSFTurnons("SF_"+canName+"_"+lepVar,effDen,dataPre,mcPre,htSel,lepVar+"_pup",trigSel,0,nLepBins,lepBins);
}

TH1 *getSFHist(TFile *fd, TFile *fmc, TString dataPre, TString mcPre, TString effDen, TString trigSel, TString sel, TString var, float rebin = 0, int nR = 0, double * rebins = 0) {

	auto getEff = [&](TFile * f, TString pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
	    TH1 * hd = 0;
	    f->GetObject(TString::Format("%s_%s__%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()),hd);
	    TH1 * hn = 0;
	    f->GetObject(TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn);
	    if(hn == 0){
	      cout << "Bad num: " + TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
	      return 0;
	    }
	    if(hd == 0){
	      cout << "Bad den: " + TString::Format("%s_%s__%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
	      return 0;
	    }
	    hn = (TH1*)hn->Clone();
	    hd = (TH1*)hd->Clone();
	    PlotTools::toOverflow(hn);
	    PlotTools::toOverflow(hd);

	    if(rebin > 0){
	      PlotTools::rebin(hn,rebin);
	      PlotTools::rebin(hd,rebin);
	    } else if(rebins){
	      hn = PlotTools::rebin(hn,nR,rebins);
	      hd = PlotTools::rebin(hd,nR,rebins);
	    }
//	    PlotTools::toOverflow(hn);
//	    PlotTools::toOverflow(hd);
	    hn->Divide(hn,hd,1,1,"b"); return hn;
	    // return PlotTools::getBinomErrors(hn,hd);
	};

	TH1 *dataEff = getEff(fd,dataPre,effDen,sel,var,trigSel,rebin,nR,rebins);
	TH1 *mcEff   = getEff(fmc,mcPre,effDen,sel,var,trigSel,rebin,nR,rebins);
	dataEff->Divide(dataEff,mcEff,1,1,"b");
	return dataEff;
}

void drawTriggerTurnons() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/trigger17_lnuqq/";
	TFile *fd = new TFile(prePath+"trigger_data.root");
	TFile *fdb = new TFile(prePath+"trigger_dataB.root");
//	TFile *fdelse = new TFile(prePath+"trigger_dataCDEF.root");
	TFile *fmc = new TFile(prePath+"trigger_ttbar2L.root");
	TFile *fs = new TFile(prePath+"trigger_blk.root");

	int nLepBins = 33;
	double lepBins[] = {5,10,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,50,75,100,150,200,250,300,350,400,450,500};
	int nHTBins = 19;
	double htBins[] = {0,100,150,200,250,300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1600,2000};

    vector<TString> muSels = {"mupt_27","mupt_27to30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
    vector<TString> elSels = {"elpt_30","elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
    vector<TString> htSels = {"ht_400","ht_500","ht_500to600","ht_600to700","ht_700to800","ht_800to900","ht_900to1000","ht_1000"};

	getSF(fd,fmc,"SingleMuon","ttbar","GL_passSMu","passSEloHtEloHEoBu","elpt_30","ht_400","electron");
	getSF(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","muon");

	TH1 *dm_lep = (TH1*)fd->Get("SingleElectron_GL_passSE_passSMuoHtMuoHMoBu__ht_400_mu_pt_pup");
	TH1 *de_lep = (TH1*)fd->Get("SingleMuon_GL_passSMu_passSEloHtEloHEoBu__ht_400_el_pt_pup");
	TH1 *dm_ht = (TH1*)fd->Get("SingleElectron_GL_passSE_passSMuoHtMuoHMoBu__mupt_27_ht_pup");
	TH1 *de_ht = (TH1*)fd->Get("SingleMuon_GL_passSMu_passSEloHtEloHEoBu__elpt_30_ht_pup");
	vector<TH1*> dataHists = {dm_lep,de_lep,dm_ht,de_ht};
	vector<TString> names = {"Muon pt","Electron pt","Muon Ht","Electron Ht"};

//	for (unsigned int i=0; i<dataHists.size(); i++) {
//		Plotter *p = new Plotter();
//		TH1 *hi = PlotTools::getIntegral(dataHists[i],false,true);
//		p->addHist(hi,"");
//		p->draw(false,names[i]);
//	}

    Plotter *pe = new Plotter();
    Plotter *pm = new Plotter();
    Plotter *phm = new Plotter();
    Plotter *phe = new Plotter();

/*    for (const auto& sel : muSels) {
    	TH1 *h = getSFHist(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu",sel,"ht",0,nHTBins,htBins);
    	pm->addHist(h,sel);
    }
    pm->drawSplitRatio(0,"mu syst",false,false,"mu syst");
    for (const auto& sel : elSels) {
    	TH1 *h = getSFHist(fd,fmc,"SingleMuon","ttbar","GL_passSMu_maxEta1p479","passSEloHtEloHEoBu",sel,"ht",0,nHTBins,htBins);
    	pe->addHist(h,sel);
    }
    pe->drawSplitRatio(0,"el syst",false,false,"el syst");
*/
//    for (const auto& sel : htSels) {
//    	TH1 *h = getSFHist(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu",sel,"mu_pt",0,nLepBins,lepBins);
//    	phm->addHist(h,sel);
//    }
//    phm->drawSplitRatio(0,"ht syst - mu",false,false,"ht syst - mu");
//    for (const auto& sel : htSels) {
//    	TH1 *h = getSFHist(fd,fmc,"SingleMuon","ttbar","GL_passSMu_maxEta1p479","passSEloHtEloHEoBu",sel,"el_pt",0,nLepBins,lepBins);
//    	phe->addHist(h,sel);
//    }
//    phe->drawSplitRatio(0,"ht syst - e",false,false,"ht syst - e");

//    plotTurnons_Sel(fmc,"ttbar_GL_passSE",muSels,"passSMuoHtMuoHMoBu","ht","Muon pt thresh comp",0,nHTBins,htBins);
//    plotTurnons_Sel(fmc,"ttbar_GL_passSMu",elSels,"passSEloHtEloHEoBu","ht","Electron pt thresh comp",0,nHTBins,htBins);

    // comparing single triggers (consult function at top for trigger names)
    const vector<int> mu_trigs = {24,25,26,10,11,12};
    const vector<int> el_trigs = {18,19};
    const vector<int> metht_trigs = {3,4,5,6,7,8};
    const vector<int> oth_trigs = {0,1,2,6,7};

//    plotIndTriggers(fb,"ttbar",el_trigs,"elpt_25","ht",+"el trigs pt 25 ht",40);

    // turnons for different trigger combinations

    vector<TString> muTrigs = {"passNomMuoJet","passNomMuoJetoMHT","passNomMuoJetoMET","passSMuoHtMuoHMoBu"};
    vector<TString> elTrigs = {"passNomEloJet","passNomEloJetoMHT","passNomEloJetoMET","passNomEloJetoPh","passNomEloJetoMHToPh",
    		"passNomEloJetoMEToPh","passNomEloJetoMHToMET","passSEloHtEloHEoBuoSEl_ETA"};

//    plotTurnons_Trig(fmc,"ttbar_GL_passSE" ,"mupt_27",muTrigs,"ht","muon dataset comp",0,nHTBins,htBins);
//    plotTurnons_Trig(fmc,"ttbar_GL_passSMu","elpt_30",elTrigs,"ht","electron dataset comp",0,nHTBins,htBins);

    // Main electron triggers
    elTrigs = {"passEl32","passEl35","passEl32Dbl","passEl32SngoDbl","passEl32SngoDblo35"};
//    plotTurnons_Trig(fb,"ttbar_passSMu","elpt_15",elTrigs,"ht",40);

    elTrigs = {"passSEl","passSEloHE","passSEloHtEl","passSEloHtEloHE","passSEloHtEloHEoBu"};
    muTrigs = {"passSMu","passSMuoHM","passSMuoHtMu","passSMuoHtMuoHM","passSMuoHtMuoHMoBu"};

    // Study other triggers such as Jet, MHT, MET_NoMu
    elTrigs = {"passSEloHtEloHE","passNomEloJet","passNomEloJetoMHT1","passNomEloJetoMHT2","passNomEloJetoMHT"};
    muTrigs = {"passSMuoHtMuoHM","passNomMuoJet","passNomMuoJetoMetNoMu","passNomMuoJetoMHT","passNomMuoJetoMHToMetNoMu"};

    // Scale factors
    muSels = {"mupt_27","mupt_30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
    elSels = {"elpt_30","elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
    htSels = {"ht_475","ht_500"};
//    vector<TString> etas = {"","_maxEta1p5","_maxEta2p5"};

    // |||| SYSTEMATICS ||||||
    // How much different are trigger efficiencies between signal and ttbar MC?
    TString eltrig = "passSEloHtEloHEoBu";
    TString mutrig = "passSMuoHtMuoHMoBu";
//	getSF(fs,fmc,"m1000","ttbar","GL_passSMu","passSEloHtEloHEoBu","elpt_30","ht_400","electron - 1 TeV signal ttbar comp");
//	getSF(fs,fmc,"m3000","ttbar","GL_passSMu","passSEloHtEloHEoBu","elpt_30","ht_400","electron - 3 TeV signal ttbar comp");
//	getSF(fs,fmc,"m1000","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","muon - 1 TeV signal ttbar comp");
//	getSF(fs,fmc,"m3000","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","muon - 3 TeV signal ttbar comp");

	vector<TString> masses = {"m800","m900","m1000","m1200","m1400","m1600","m1800","m2000",
			"m2500","m3000","m3500","m4000","m4500"};
	TString intHistm = "_SYST__mupt_27_ht_pup";
	TString intHiste = "_SYST__elpt_27_ht_pup";
/*
	Plotter *pp = new Plotter();
	Plotter *ppe = new Plotter();

	double intgm, intge;
	TString mhst = "_SYST_passSMuoHtMuoHMoBu__mupt_27_ht_pup";
	TString ehst = "_SYST_passSEloHtEloHEoBu__elpt_27_ht_pup";
	TH1 *ttm = (TH1*)fmc->Get("ttbar"+mhst);
	TH1 *tte = (TH1*)fmc->Get("ttbar"+ehst);
	for (const auto& m : masses) {
		TH1 *hhm = (TH1*)fs->Get(m+intHistm);
		TH1 *hhe = (TH1*)fs->Get(m+intHiste);
		TH1 *him = PlotTools::getIntegral(hhm,false,true);
		TH1 *hie = PlotTools::getIntegral(hhe,false,true);
		pp->addHist(him,m);
		ppe->addHist(hie,m);

		intgm = him->GetBinContent(him->FindBin(1200));
		intge = hie->GetBinContent(hie->FindBin(1200));
//		cout<<"integ = "<<intg<<endl;

    	TH1 *sfm = getSFHist(fs,fmc,m,"ttbar","SYST","passSMuoHtMuoHMoBu","mupt_27","ht_pup",0,nHTBins,htBins);
    	TH1 *sfe = getSFHist(fs,fmc,m,"ttbar","SYST","passSEloHtEloHEoBu","elpt_30","ht_pup",0,nHTBins,htBins);

		double maxerr_m = 0, maxerr_e = 0;
		int its = 15;
		for (int i=0; i<its; i++) {
			double err_m = sfm->GetBinContent(sfm->FindBin(400)+i);
			double err_e = sfe->GetBinContent(sfe->FindBin(400)+i);

//			cout<<"slurm error "<<fabs(1-err)<<endl;
			if (fabs(1-err_m) > maxerr_m && err_m > 0.5) maxerr_m = fabs(1-err_m);
			if (fabs(1-err_e) > maxerr_e && err_e > 0.5) maxerr_e = fabs(1-err_e);

		}
		double syst_m = maxerr_m*intgm;
		double syst_e = maxerr_e*intge;

		printf("MU %s = %f\n",m.Data(),syst_m);
		printf("EL %s = %f\n",m.Data(),syst_e);

		getSF(fs,fmc,m,"ttbar","SYST","passSEloHtEloHEoBu","elpt_30","ht_400",m+" electron signal ttbar comp");
		getSF(fs,fmc,m,"ttbar","SYST","passSMuoHtMuoHMoBu","mupt_27","ht_400",m+" muon signal ttbar comp");

	}
	pp->draw(false,"slurm u");
	ppe->draw(false,"electrons");
*/

	// How do efficiencies look for different lepton pt regimes?
//	plotTurnons_Sel(fmc,"ttbar_GL_passSE",muSels,mutrig,"ht","Muons",40);
//	plotTurnons_Sel(fmc,"ttbar_GL_passSMu",elSels,eltrig,"ht","Electrons",40);


//    for (const auto& sel : elSels) {getSF(fd,fmc,"SingleMuon","ttbar","GL_passSMu","passSEloHtEloHEoBu",sel,"ht_500","ht500_"+sel);}
//    for (const auto& sel : muSels) {getSF(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu",sel,"ht_500","ht500_"+sel);}
//    for (const auto& sel : htSels) {
//    	getSF(fd,fmc,"SingleMuon","ttbar","GL_passSMu","passSEloHtEloHEoBu","elpt_30",sel,"elpt_30_"+sel);
//        getSF(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_30",sel,"mupt_30_"+sel);
//    }

//    for (const auto& eta:etas) {getSF(fd,fmc,"SingleElectron","ttbar","GL_passSE"+eta,"passSMuoHtMuoHMoBu","mupt_26","ht_500",eta);}
//    getSF(fdb,   fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_30","ht_400","RunB");
//    getSF(fdelse,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_30","ht_400","RunCDEF");

	return;
}
