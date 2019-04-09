#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
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

void plotTurnons_Sel(TFile *f, TString prefix, const vector<TString>& sels, TString trig, TString var, TString canName, float rebin = 0) {
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
	    if (rebin > 0) {
	    	PlotTools::rebin(hn,rebin);
	    	PlotTools::rebin(hd,rebin);
	    }
//	    auto *ratio = PlotTools::getBinomErrors(hn,hd);
//	    p->addGraph(ratio,s);
	    hn->Divide(hd);
	    p->addHist(hn,s,-1,1,4,20,1,true,false);
	}
//	p->draw(false,canName);
	p->drawSplitRatio(0,"stack",false,false,canName);
}

void plotTurnons_Trig(TFile *f, TString prefix, const TString sel, const vector<TString>& trigs, TString var, TString canName, float rebin = 0) {
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
	    if (rebin > 0) {
	    	PlotTools::rebin(hn,rebin);
	    	PlotTools::rebin(hd,rebin);
	    }
//	    auto *ratio = PlotTools::getBinomErrors(hn,hd);
//	    p->addGraph(ratio,s);
	    hn->Divide(hd);
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

	    p->addHist(mcEff,mcPre,-1,1,4,20,1,true,true, false, "E X P");
	    p->addHist(dataEff,dataPre,-1,1,4,20,1,true,true, false, "E X P");
	      // p->draw(true,TString::Format("%s.pdf",name.Data()));
	    p->drawSplitRatio(0,"stack",true,false,TString::Format("%s.pdf",name.Data()));
//	    p->draw(false,TString::Format("%s.pdf",name.Data()));
	};

	TString lepVar = (lepSel.BeginsWith("m")) ? "mu_pt" : "el_pt";

	plotSFTurnons("SF_"+canName+"_HT",effDen,dataPre,mcPre,lepSel,"ht",trigSel,0,nHTBins,htBins);
	plotSFTurnons("SF_"+canName+"_"+lepVar,effDen,dataPre,mcPre,htSel,lepVar,trigSel,0,nLepBins,lepBins);
}

void drawTriggerTurnons() {
	TString prePath = "/Users/brentstone/Dropbox/Physics/HHbbWW/plots/";
	TFile *fd = new TFile(prePath+"trigger_data.root");
	TFile *fdb = new TFile(prePath+"trigger_dataB.root");
	TFile *fdelse = new TFile(prePath+"trigger_dataCDEF.root");
	TFile *fb = new TFile(prePath+"out_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_0.root");
	TFile *fmc = new TFile(prePath+"trigger_ttbar2L.root");
	TFile *fs = new TFile(prePath+"trigger_blkGrav.root");

    vector<TString> muSels = {"mupt_10","mupt_15","mupt_20","mupt_25","mupt_26","mupt_30","mupt_35"};
    vector<TString> elSels = {"elpt_10","elpt_15","elpt_20","elpt_25","elpt_30","elpt_35","elpt_40"};
    vector<TString> htSels = {"ht_350","ht_375","ht_400","ht_450","ht_475","ht_500","ht_600"};
//    plotTurnons_Sel(fb,"ttbar_GL_passSE",muSels,"passSMuoHtMuoHMoBu","ht",40);
//    plotTurnons_Sel(fb,"ttbar_GL_passSE",htSels,"passSMuoHtMuoHMoBu","mu_pt",25);
//    plotTurnons_Sel(fb,"ttbar_GL_passSMu",elSels,"passSEloHtEloHEoBu","ht",40);
//    plotTurnons_Sel(fb,"ttbar_GL_passSMu",htSels,"passSEloHtEloHEoBu","el_pt",25);

    // comparing single triggers (consult function at top for trigger names)
    const vector<int> mu_trigs = {24,25,26,10,11,12};
    const vector<int> el_trigs = {18,19};
    const vector<int> metht_trigs = {3,4,5,6,7,8};
    const vector<int> oth_trigs = {0,1,2,6,7};

//    plotIndTriggers(fb,"ttbar",el_trigs,"elpt_25","ht",+"el trigs pt 25 ht",40);
//    plotIndTriggers(fmc,"ttbar",el_trigs,"ht_500","el_pt","jonson",10);
//    plotIndTriggers(fb,"ttbar",mu_trigs,"mupt_15","ht",40);
//    plotIndTriggers(fb,"ttbar",metht_trigs,"mupt_15","ht",40);
//    plotIndTriggers(fb,"ttbar",metht_trigs,"elpt_15","ht",40);
//    plotIndTriggers(fb,"ttbar",oth_trigs,"mupt_15","ht",40);
//    plotIndTriggers(fb,"ttbar",oth_trigs,"elpt_15","ht",40);

    // turnons for different trigger combinations
    const TString mu_WP = "mupt_26";
    const TString el_WP = "elpt_30";
    const vector<TString> muTrigCombs = {"passSMu","passSMuoHtMu","passSMuoHtMuoHMoBu","passMuDenNoCross"};
    const vector<TString> elTrigCombs = {"passSEl","passSEloHtEl","passSEloHtEloHEoBu","passElDenNoCross"};

//    plotTurnons_Trig(fb,"ttbar_GL_passSE",mu_WP,muTrigCombs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu",el_WP,elTrigCombs,"ht",40);

    // Triggers from MC using probes only
    vector<TString> muTrigs = {"passSMu","passSMuoHtMu","passSMuoHM","passSMuoHtMuoHM","passSMuoHtMuoHMoBu"};
    vector<TString> elTrigs = {"passSEl","passSEloHtEl","passSEloHE","passSEloHtEloHE","passSEloHtEloHEoBu"};

//    plotTurnons_Trig(fb,"ttbar_MC","mupt_incl",muTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_MC","elpt_incl",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_MC","mupt_26",muTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_MC","elpt_30",elTrigs,"ht",40);
//
//    plotTurnons_Trig(fb,"ttbar_MC","ht_incl",muTrigs,"mu_pt",25);
//    plotTurnons_Trig(fb,"ttbar_MC","ht_incl",elTrigs,"el_pt",25);
//    plotTurnons_Trig(fb,"ttbar_MC","ht_400",muTrigs,"mu_pt",25);
//    plotTurnons_Trig(fb,"ttbar_MC","ht_400",elTrigs,"el_pt",25);

    // Triggers using tag and probe approach with emu ttbar
//    plotTurnons_Sel(fb,"ttbar_GL_passSE",muSels,"passSMu","ht",40);
//    plotTurnons_Sel(fb,"ttbar_GL_passSMu",elSels,"passSEl","ht",40);
//    plotTurnons_Sel(fb,"ttbar_GL_passSE",htSels,"passSMu","mu_pt",25);
//    plotTurnons_Sel(fb,"ttbar_GL_passSMu",htSels,"passSEl","el_pt",25);
//
//    plotTurnons_Trig(fb,"ttbar_GL_passSE","mupt_15",muTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSE","ht_incl",muTrigs,"mu_pt",25);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","elpt_15",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","ht_incl",elTrigs,"el_pt",25);

    // Main electron triggers
    elTrigs = {"passEl32","passEl35","passEl32Dbl","passEl32SngoDbl","passEl32SngoDblo35"};
//    plotTurnons_Trig(fb,"ttbar_passSMu","elpt_15",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_passSMu","ht_350",elTrigs,"el_pt",25);

    // High energy electron triggers
    elTrigs = {"passSEl","passSEloPh200","passSEloEl115","passSEloHE"};
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","elpt_15",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","ht_350",elTrigs,"el_pt",25);

    // Electron cross trigger comparison
    elTrigs = {"passEl32Dbl","passElHT1","passElHT2","passElHT"};
//    plotTurnons_Trig(fb,"ttbar_passSMu","elpt_15",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_passSMu","ht_350",elTrigs,"el_pt",25);

    // High energy, HT, and cross triggers
    elTrigs = {"passSEl","passSEloHE","passSEloHtEl","passSEloHtEloHE","passSEloHtEloHEoBu"};
    muTrigs = {"passSMu","passSMuoHM","passSMuoHtMu","passSMuoHtMuoHM","passSMuoHtMuoHMoBu"};
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","elpt_15",elTrigs,"ht",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","ht_350",elTrigs,"el_pt",25);
//    plotTurnons_Trig(fmc,"ttbar_GL_passSE","mupt_15",muTrigs,"ht","muon trig comp (pt>15) - ht",40);
//    plotTurnons_Trig(fmc,"ttbar_GL_passSE","ht_350",muTrigs,"mu_pt","muon trig comp (ht>350)- pt",25);

    // Study other triggers such as Jet, MHT, MET_NoMu
    elTrigs = {"passSEloHtEloHE","passNomEloJet","passNomEloJetoMHT1","passNomEloJetoMHT2","passNomEloJetoMHT"};
    muTrigs = {"passSMuoHtMuoHM","passNomMuoJet","passNomMuoJetoMetNoMu","passNomMuoJetoMHT","passNomMuoJetoMHToMetNoMu"};

//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","elpt_30",elTrigs,"ht","ttbar_GL_passSMu_elpt30_MiscTrigs",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSMu","ht_350",elTrigs,"el_pt","ttbar_GL_passSMu_ht350_MiscTrigs",25);
//    plotTurnons_Trig(fb,"ttbar_GL_passSE","mupt_26",muTrigs,"ht","ttbar_GL_passSE_mupt26_MiscTrigs",40);
//    plotTurnons_Trig(fb,"ttbar_GL_passSE","ht_350",muTrigs,"mu_pt","ttbar_GL_passSE_ht350_MiscTrigs",25);

    // Scale factors
    muSels = {"mupt_27","mupt_30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
    elSels = {"elpt_30","elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
    htSels = {"ht_475","ht_500"};
//    vector<TString> etas = {"","_maxEta1p5","_maxEta2p5"};

    elTrigs = {/*"passSEloHtEloHEoBu",*/"passSEloHtEloHEoBuoSEl_ETA"};
    getSF(fd,fmc,"SingleMuon","ttbar","GL_passSMu_maxEta2p1","passSEloHtEloHEoBuoSEl_ETA","elpt_30","ht_400","Electron trigger SF");
    getSF(fd,fmc,"SingleElectron","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","Muon trigger SF");

//    for (const auto& eta : etas) {
//    	plotTurnons_Trig(fmc,"ttbar_GL_passSMu"+eta,"ht_400",elTrigs,"el_pt","PT compare eta restricted el trig "+eta);
//    	plotTurnons_Trig(fmc,"ttbar_GL_passSMu"+eta,"elpt_27",elTrigs,"ht","HT compare eta restricted el trig "+eta,20);
//    	for (const auto& sel : elSels) {
//    		getSF(fd,fmc,"SingleMuon","ttbar","GL_passSMu","passSEloHtEloHEoBuoSEl_ETA",sel,"ht_400","passSEloHtEloHEoBuoSEl_ETA "+sel);
//    	}
//    }

    // |||| SYSTEMATICS ||||||
    // How much different are trigger efficiencies between signal and ttbar MC?
    TString eltrig = "passSEloHtEloHEoBuoSEl_ETA";
    TString mutrig = "passSMuoHtMuoHMoBu";
//	getSF(fs,fmc,"m1200","ttbar","GL_passSMu","passSEloHtEloHEoBuoSEl_ETA","elpt_30","ht_400","electron - 1.2 TeV signal ttbar comp");
//	getSF(fs,fmc,"m2000","ttbar","GL_passSMu","passSEloHtEloHEoBuoSEl_ETA","elpt_30","ht_400","electron - 2 TeV signal ttbar comp");
//	getSF(fs,fmc,"m1200","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","muon - 1.2 TeV signal ttbar comp");
//	getSF(fs,fmc,"m2000","ttbar","GL_passSE","passSMuoHtMuoHMoBu","mupt_27","ht_400","muon - 2 TeV signal ttbar comp");

	// How do efficiencies look for different lepton pt regimes?
//	plotTurnons_Sel(fmc,"ttbar_GL_passSE",muSels,mutrig,"ht","Muons",40);
//	plotTurnons_Sel(fmc,"ttbar_GL_passSMu",elSels,eltrig,"ht","Electrons",40);

 /*
    TString elPre = "ttbar_GL_passSMu";
    TString muPre = "ttbar_GL_passSE";
    std::vector<TString> etas = {"_maxEta1p5_","_maxEta2p5_"};

    TString eltrig = "passSEloHtEloHEoBuoSEl_ETA";
    TString mutrig = "passSMuoHtMuoHMoBu";

    Plotter *pe1 = new Plotter();
    Plotter *pe2 = new Plotter();
    Plotter *pm1 = new Plotter();
    Plotter *pm2 = new Plotter();

	int nLepBins = 33;
	double lepBins[] = {5,10,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,50,75,100,150,200,250,300,350,400,450,500};
	int nHTBins = 28;
	double htBins[] = {100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,700,800,900,1000,1100,1200,1600,2000};



    for (const auto& eta : etas) {
		std::cout<<elPre+eta+eltrig+"__ht_400_el_pt"<<std::endl;
		TH1 *pte = (TH1*)fmc->Get(elPre+eta+eltrig+"__ht_400_el_pt");
		TH1 *hte = (TH1*)fmc->Get(elPre+eta+eltrig+"__elpt_30_ht");
		TH1 *ptm = (TH1*)fmc->Get(muPre+eta+mutrig+"__ht_400_mu_pt");
		TH1 *htm = (TH1*)fmc->Get(muPre+eta+mutrig+"__mupt_27_ht");

		TH1 *den_pte = (TH1*)fmc->Get(elPre+eta+"_ht_400_el_pt");
		TH1 *den_hte = (TH1*)fmc->Get(elPre+eta+"_elpt_30_ht");
		TH1 *den_ptm = (TH1*)fmc->Get(muPre+eta+"_ht_400_mu_pt");
		TH1 *den_htm = (TH1*)fmc->Get(muPre+eta+"_mupt_27_ht");

	    PlotTools::toOverflow(den_pte);
	    PlotTools::toOverflow(den_hte);
	    PlotTools::toOverflow(den_ptm);
	    PlotTools::toOverflow(den_htm);

	    den_pte = PlotTools::rebin(den_pte,nLepBins,lepBins);
	    den_hte = PlotTools::rebin(den_hte,nHTBins,htBins);
	    den_ptm = PlotTools::rebin(den_ptm,nLepBins,lepBins);
	    den_htm = PlotTools::rebin(den_htm,nHTBins,htBins);

//		pe1->addHist(den_pte,"den",-1,1,4,20,1,true,false);
//		pe2->addHist(den_hte,"den",-1,1,4,20,1,true,false);
//		pm1->addHist(den_ptm,"den",-1,1,4,20,1,true,false);
//		pm2->addHist(den_htm,"den",-1,1,4,20,1,true,false);

	    pte = (TH1*)pte->Clone();
	    hte = (TH1*)hte->Clone();
	    ptm = (TH1*)ptm->Clone();
	    htm = (TH1*)htm->Clone();

	    PlotTools::toOverflow(pte);
	    PlotTools::toOverflow(hte);
	    PlotTools::toOverflow(ptm);
	    PlotTools::toOverflow(htm);

	    pte = PlotTools::rebin(pte,nLepBins,lepBins);
	    hte = PlotTools::rebin(hte,nHTBins,htBins);
	    ptm = PlotTools::rebin(ptm,nLepBins,lepBins);
	    htm = PlotTools::rebin(htm,nHTBins,htBins);

	    pte->Divide(den_pte);
	    hte->Divide(den_hte);
	    ptm->Divide(den_ptm);
	    htm->Divide(den_htm);

		pe1->addHist(pte,eta,-1,1,4,20,1,true,false);
		pe2->addHist(hte,eta,-1,1,4,20,1,true,false);
		pm1->addHist(ptm,eta,-1,1,4,20,1,true,false);
		pm2->addHist(htm,eta,-1,1,4,20,1,true,false);
    }
    pe1->drawSplitRatio(0,"stack",false,false,"electron pt");
    pe2->drawSplitRatio(0,"stack",false,false,"ht electrons");
    pm1->drawSplitRatio(0,"stack",false,false,"muon pt");
    pm2->drawSplitRatio(0,"stack",false,false,"ht muons");
*/
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

//
//
////T&P data/mc sf
//{
//  TFile * fd  = new TFile("data_h_triggerTurnons.root");
//  TFile * fMC = new TFile("all_triggerTurnons.root");
//  Plotter * pt = new Plotter();
//
//  // auto getEff = [&](TFile * f, TString pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0)->TGraphAsymmErrors*{
//  //   TH1 * hd = 0;
//  //   f->GetObject(TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()),hd);
//  //   TH1 * hn = 0;
//  //   f->GetObject(TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn);
//  //   if(hn == 0){
//  //     cout << TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
//  //     return 0;
//  //   }
//  //   if(hd == 0){
//  //     cout << TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
//  //     return 0;
//  //   }
//  //   hn = (TH1*)hn->Clone();
//  //   hd = (TH1*)hd->Clone();
//  //   PlotTools::toOverflow(hn);
//  //   PlotTools::toOverflow(hd);
//  //   if(rebin > 0){
//  //     int nLepBins = 10;
//  //     double lepBins[] = {5,10,15,20,25,30,35,50,75,100,500};
//  //     hn = PlotTools::rebin(hn,nLepBins,lepBins);
//  //     hd = PlotTools::rebin(hd,nLepBins,lepBins);
//  //
//  //     // PlotTools::rebin(hn,rebin);
//  //     // PlotTools::rebin(hd,rebin);
//  //   }
//  //   return PlotTools::getBinomErrors(hn,hd);
//  // };
//  auto getEff = [&](TFile * f, TString pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
//    TH1 * hd = 0;
//    f->GetObject(TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()),hd);
//    TH1 * hn = 0;
//    f->GetObject(TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn);
//    if(hn == 0){
//      cout << TString::Format("%s_%s_%s__%s_%s",pn.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
//      return 0;
//    }
//    if(hd == 0){
//      cout << TString::Format("%s_%s_%s_%s",pn.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
//      return 0;
//    }
//    hn = (TH1*)hn->Clone();
//    hd = (TH1*)hd->Clone();
//
//    if(rebin > 0){
//      PlotTools::rebin(hn,rebin);
//      PlotTools::rebin(hd,rebin);
//    } else if(rebins){
//      hn = PlotTools::rebin(hn,nR,rebins);
//      hd = PlotTools::rebin(hd,nR,rebins);
//    }
//    PlotTools::toOverflow(hn);
//    PlotTools::toOverflow(hd);
//          hn->Divide(hn,hd,1,1,"b"); return hn;
//    // return PlotTools::getBinomErrors(hn,hd);
//  };
//
//  auto plotTurnons =[&](TString name, TString prefix,TString dataName,TString mcName, TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0 ){
//      Plotter * p = new Plotter();
//      auto * mcEff = getEff(fMC,mcName,prefix,sel,var,trig,rebin,nR,rebins);
//      auto * dataEff = getEff(fd,dataName,prefix,sel,var,trig,rebin,nR,rebins);
//      if(mcEff == 0 || dataEff == 0) return;
//      // p->addGraph(mcEff,"MC");
//      // p->addGraph(dataEff,"data");
//      p->addHist(mcEff,"MC",-1,1,4,20,1,true,true, false, "E X P");
//      p->addHist(dataEff,"data",-1,1,4,20,1,true,true, false, "E X P");
//      // p->draw(true,TString::Format("%s.pdf",name.Data()));
//      p->drawRatio(0,"stack",false,true,TString::Format("%s.pdf",name.Data()));
//  };
//    std::vector<TString> muSels = {"mupt_20","mupt_25","mupt_30","mupt_50","mupt_100"};
//    std::vector<TString> elSels = {"elpt_20","elpt_25","elpt_30","elpt_50","elpt_100"};
//    std::vector<TString> htSels = {"ht_450","ht_475","ht_500","ht_525"};
//    // int nLepBins = 12;
//    // double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250};
//    int nLepBins = 8;
//    double lepBins[] = {5,10,15,20,25,30,35,50,75};
//    int nHTBins = 16;
//    double htBins[] = {0,50,100,150,200,250,300,350,400,450,500,550,600,800,1200,1600,2000};
//    //ht cuts
//    // for(const auto& s : htSels){
//    //   plotTurnons(TString::Format("turnOn_%s_muon",s.Data()),"GL_passSE","singlee","ttbar",s,"mu_pt","passSMuoHtMuoBu"        ,0,nLepBins,lepBins);
//    // }
//    // for(const auto& s : htSels){
//    //   plotTurnons(TString::Format("turnOn_%s_electron",s.Data()),"GL_passSMu","singlemu","ttbar",s,"el_pt","passSEloHtEloBu"        ,0,nLepBins,lepBins);
//    // }
//    for(const auto& s : muSels){
//      plotTurnons(TString::Format("turnOn_%s_muon_ht",s.Data()),"GL_passSE","singlee","ttbar",s,"ht","passSMuoHtMuoBu"        ,0,nHTBins,htBins);
//    }
//    for(const auto& s : elSels){
//      plotTurnons(TString::Format("turnOn_%s_electron_ht",s.Data()),"GL_passSMu","singlemu","ttbar",s,"ht","passSEloHtEloBu"        ,0,nHTBins,htBins);
//    }
//
//    //El denom
//    // plotTurnons("elDenom_selHT_muLeg"      ,"singlee_GL_passSE",htSels,"mu_pt","passSMu"        ,5);
//    // plotTurnons("elDenom_selHT_muLeg_di"   ,"singlee_GL_passSE",htSels,"mu_pt","passSMuoHtMu"   ,5);
//    // plotTurnons("elDenom_selHT_muLeg_diBu" ,"singlee_GL_passSE",htSels,"mu_pt","passSMuoHtMuoBu",5);
//    // plotTurnons("elDenom_selMu_htLeg"      ,"singlee_GL_passSE",muSels,"ht","passSMu"        ,25);
//    // plotTurnons("elDenom_selMu_htLeg_di"   ,"singlee_GL_passSE",muSels,"ht","passSMuoHtMu"   ,25);
//    // plotTurnons("elDenom_selMu_htLeg_diBu" ,"singlee_GL_passSE",muSels,"ht","passSMuoHtMuoBu",25);
//    // plotTurnons("muDenom_selHT_elLeg"      ,"singlemu_GL_passSMu",htSels,"el_pt","passSEl"        ,5);
//    // plotTurnons("muDenom_selHT_elLeg_di"   ,"singlemu_GL_passSMu",htSels,"el_pt","passSEloHtEl"   ,5);
//    // plotTurnons("muDenom_selHT_elLeg_diBu" ,"singlemu_GL_passSMu",htSels,"el_pt","passSEloHtEloBu",5);
//    // plotTurnons("muDenom_selMu_htLeg"      ,"singlemu_GL_passSMu",elSels,"ht","passSEl"        ,25);
//    // plotTurnons("muDenom_selMu_htLeg_di"   ,"singlemu_GL_passSMu",elSels,"ht","passSEloHtEl"   ,25);
//    // plotTurnons("muDenom_selMu_htLeg_diBu" ,"singlemu_GL_passSMu",elSels,"ht","passSEloHtEloBu",25);
//
//    // plotTurnons("selHT_muLeg"      ,"m1000_MC",htSels,"mu_pt","passSMu"        ,5);
//    // plotTurnons("selHT_muLeg_di"   ,"m1000_MC",htSels,"mu_pt","passSMuoHtMu"   ,5);
//    // plotTurnons("selHT_muLeg_diBu" ,"m1000_MC",htSels,"mu_pt","passSMuoHtMuoBu",5);
//    // plotTurnons("selMu_htLeg"      ,"m1000_MC",muSels,"ht","passSMu"        ,25);
//    // plotTurnons("selMu_htLeg_di"   ,"m1000_MC",muSels,"ht","passSMuoHtMu"   ,25);
//    // plotTurnons("selMu_htLeg_diBu" ,"m1000_MC",muSels,"ht","passSMuoHtMuoBu",25);
//    // plotTurnons("selHT_elLeg"      ,"m1000_MC",htSels,"el_pt","passSEl"        ,5);
//    // plotTurnons("selHT_elLeg_di"   ,"m1000_MC",htSels,"el_pt","passSEloHtEl"   ,5);
//    // plotTurnons("selHT_elLeg_diBu" ,"m1000_MC",htSels,"el_pt","passSEloHtEloBu",5);
//    // plotTurnons("selEl_htLeg"      ,"m1000_MC",elSels,"ht","passSEl"        ,25);
//    // plotTurnons("selEl_htLeg_di"   ,"m1000_MC",elSels,"ht","passSEloHtEl"   ,25);
//    // plotTurnons("selEl_htLeg_diBu" ,"m1000_MC",elSels,"ht","passSEloHtEloBu",25);
//
//
//
//}
//
//
////T&P data/mc sf2
//{
//  TFile * fd  = new TFile("data_triggerTurnons.root");
//  TFile * fMC = new TFile("all_triggerTurnons.root");
//  Plotter * pt = new Plotter();
//
//  auto getEff = [&](TFile * f, std::vector<TString> pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
//    TH1 * hd = 0;
//    TH1 * hn = 0;
//
//    for(const auto&  n : pn){
//      TH1 * hd1 = 0;
//      f->GetObject(TString::Format("%s_%s_%s_%s",n.Data(),prefix.Data(),sel.Data(),var.Data()),hd1);
//      TH1 * hn1 = 0;
//      f->GetObject(TString::Format("%s_%s_%s__%s_%s",n.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn1);
//      if(hn1 == 0){
//        cout << TString::Format("%s_%s_%s__%s_%s",n.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
//        return 0;
//      }
//      if(hd1 == 0){
//        cout << TString::Format("%s_%s_%s_%s",n.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
//        return 0;
//      }
//      if(hd == 0) hd = (TH1*)hd1->Clone();
//      else hd->Add(hd1);
//      if(hn == 0) hn = (TH1*)hn1->Clone();
//      else hn->Add(hn1);
//    }
//
//
//    if(rebin > 0){
//      PlotTools::rebin(hn,rebin);
//      PlotTools::rebin(hd,rebin);
//    } else if(rebins){
//      hn = PlotTools::rebin(hn,nR,rebins);
//      hd = PlotTools::rebin(hd,nR,rebins);
//    }
//    PlotTools::toOverflow(hn);
//    PlotTools::toOverflow(hd);
//          hn->Divide(hn,hd,1,1,"b"); return hn;
//    // return PlotTools::getBinomErrors(hn,hd);
//  };
//
//
//  auto plotTurnons =[&](TString name, TString prefix,std::vector<TString> dataName,std::vector<TString> mcName, const vector<TString>& sels,const vector<TString>& selNames, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0 ){
//      Plotter * p = new Plotter();
//      for(unsigned int iS = 0; iS <sels.size(); ++iS){
//        TString sel = sels[iS];
//        TString selN = selNames[iS];
//        auto * mcEff = getEff(fMC,mcName,prefix,sel,var,trig,rebin,nR,rebins);
//        auto * dataEff = getEff(fd,dataName,prefix,sel,var,trig,rebin,nR,rebins);
//        if(mcEff == 0 || dataEff == 0) return;
//        dataEff->Divide(mcEff);
//        p->addHist(dataEff,selN,-1,1,4,20,1,true,true, false, "E X P");
//      }
//      p->setUnderflow(false);
//      p->setMinMax(0.68,1.09);
//      p->setBotMinMax(0.88,1.12);
//      p->setCMSLumi();
//      p->setCMSLumiExtraText("Preliminary");
//      p->setYTitle("data/MC trigger SF");
//      p->setYTitleBot("SF(excl)/SF(incl)");
//      p->setLegendPos(0.45,0.2,0.9,0.65);
//      // p->draw(true,TString::Format("%s.pdf",name.Data()));
//      p->drawSplitRatio(0,"stack",false,true,TString::Format("%s.pdf",name.Data()));
//  };
//    // std::vector<TString> muSels = {"mupt_20to25","mupt_25to30","mupt_30to50","mupt_50to100","mupt_100"};
//    // std::vector<TString> elSels = {"elpt_20to25","elpt_25to30","elpt_30to50","elpt_50to100","elpt_100"};
//
//  std::vector<TString> muSels = {"mupt_26","mupt_26to30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
//  std::vector<TString> elSels = {"elpt_30","elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
//
//  std::vector<TString> muSelNs = {"muon #it{p}_{T}>26 GeV","muon #it{p}_{T} 26-30 GeV","muon #it{p}_{T} 20-35 GeV","muon #it{p}_{T} 25-40 GeV","muon #it{p}_{T} 40-50 GeV","muon #it{p}_{T} 50-100 GeV","muon #it{p}_{T}>100 GeV"};
//  std::vector<TString> elSelNs = {"electron #it{p}_{T}>30 GeV","electron #it{p}_{T} 30-35 GeV","electron #it{p}_{T} 35-40 GeV","electron #it{p}_{T} 40-50 GeV","electron #it{p}_{T} 50-100 GeV","electron #it{p}_{T}>100 GeV"};
//
//    std::vector<TString> htSels = {"ht_500","ht_500to600","ht_600to700","ht_700to800","ht_800to900","ht_900to1000","ht_1000to1200","ht_1200"};
//    int nLepBins = 13;
//    double lepBins[] = {25,30,35,50,75,100,150,200,250,300,350,400,450,500};
//    // int nLepBins = 8;
//    // double lepBins[] = {5,10,15,20,25,30,35,50,75};
//    int nHTBins = 14;
//    double htBins[] = {100,150,200,250,300,350,400,450,500,550,600,800,1200,1600,2000};
//
//    std::vector<TString> bkgNames = {"ttbar","diboson","wjets","ttX","singlet","zjets"};
//
//
//    plotTurnons(TString::Format("turnOn_vary_muon_ht"),"GL_passSE",{"singlee"},bkgNames,muSels,muSelNs,"ht","passSMuoHtMuoBu"        ,0,nHTBins,htBins);
//    plotTurnons(TString::Format("turnOn_vary_electron_ht"),"GL_passSMu",{"singlemu"},bkgNames,elSels,elSelNs,"ht","passSEloHtEloBu"        ,0,nHTBins,htBins);
//
//    plotTurnons(TString::Format("turnOn_vary_ht_muon"),"GL_passSE",{"singlee"},bkgNames,htSels,htSels,"mu_pt","passSMuoHtMuoBu"        ,0,nLepBins,lepBins);
//    plotTurnons(TString::Format("turnOn_vary_ht_electron"),"GL_passSMu",{"singlemu"},bkgNames,htSels,htSels,"el_pt","passSEloHtEloBu"        ,0,nLepBins,lepBins);
//
//
//
//}
//
//
////T&P data/mc sfplot for AN
//{
//  TFile * fd  = new TFile("data_triggerTurnons.root");
//  TFile * fMC = new TFile("all_triggerTurnons.root");
//  Plotter * pt = new Plotter();
//
//  TFile * of = new TFile("triggerSF.root","recreate");
//
//  auto getEff = [&](TFile * f, std::vector<TString> pn, TString prefix,TString sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0)->TH1*{
//    TH1 * hd = 0;
//    TH1 * hn = 0;
//
//    for(const auto&  n : pn){
//      TH1 * hd1 = 0;
//      f->GetObject(TString::Format("%s_%s_%s_%s",n.Data(),prefix.Data(),sel.Data(),var.Data()),hd1);
//      TH1 * hn1 = 0;
//      f->GetObject(TString::Format("%s_%s_%s__%s_%s",n.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()),hn1);
//      if(hn1 == 0){
//        cout << TString::Format("%s_%s_%s__%s_%s",n.Data(),prefix.Data(),trig.Data(),sel.Data(),var.Data()) << endl;
//        return 0;
//      }
//      if(hd1 == 0){
//        cout << TString::Format("%s_%s_%s_%s",n.Data(),prefix.Data(),sel.Data(),var.Data()) << endl;
//        return 0;
//      }
//      if(hd == 0) hd = (TH1*)hd1->Clone();
//      else hd->Add(hd1);
//      if(hn == 0) hn = (TH1*)hn1->Clone();
//      else hn->Add(hn1);
//    }
//
//
//    if(rebin > 0){
//      PlotTools::rebin(hn,rebin);
//      PlotTools::rebin(hd,rebin);
//    } else if(rebins){
//      hn = PlotTools::rebin(hn,nR,rebins);
//      hd = PlotTools::rebin(hd,nR,rebins);
//    }
//    PlotTools::toOverflow(hn);
//    PlotTools::toOverflow(hd);
//          hn->Divide(hn,hd,1,1,"b"); return hn;
//    // return PlotTools::getBinomErrors(hn,hd);
//  };
//
//
//  auto plotTurnons =[&](TString oHName, TString name, TString prefix,std::vector<TString> dataName,std::vector<TString> mcName, const TString& sel, TString var, TString trig, float rebin = 0, int nR = 0, double * rebins = 0 ){
//      Plotter * p = new Plotter();
//        auto * mcEff = getEff(fMC,mcName,prefix,sel,var,trig,rebin,nR,rebins);
//        auto * dataEff = getEff(fd,dataName,prefix,sel,var,trig,rebin,nR,rebins);
//        if(mcEff == 0 || dataEff == 0) return;
//        p->addHist(mcEff,"MC",-1,1,4,20,1,true,true, false, "E X P");
//        p->addHist(dataEff,"data",-1,1,4,20,1,true,true, false, "E X P");
//        p->setUnderflow(false);
//        p->setMinMax(0.68,1.09);
//        p->setBotMinMax(0.925,1.075);
//        p->setCMSLumi();
//        p->setCMSLumiExtraText("Preliminary");
//        p->setYTitle("trigger efficiency");
//        p->setYTitleBot("data/MC");
//        p->setLegendPos(0.7,0.4,0.95,0.65);
//        // p->setXTitleBot("data/MC");
//        p->drawSplitRatio(0,"stack",false,true,TString::Format("%s.pdf",name.Data()));
//
//        of->cd();
//        TH1 * rat = (TH1*)dataEff->Clone(oHName);
//        rat->SetDirectory(0);
//        rat->Divide(mcEff);
//        rat->Write();
//  };
//
//  std::vector<TString> muSels = {"mupt_25to30","mupt_30to35","mupt_35to40","mupt_40to50","mupt_50to100","mupt_100"};
//  std::vector<TString> elSels = {"elpt_30to35","elpt_35to40","elpt_40to50","elpt_50to100","elpt_100"};
//    std::vector<TString> htSels = {"ht_450","ht_475","ht_500","ht_600","ht_700","ht_800","ht_1000","ht_1200"};
//    int nLepBins = 17;
//    double lepBins[] = {5,10,15,20,25,30,35,50,75,100,150,200,250,300,350,400,450,500};
//    // int nLepBins = 8;
//    // double lepBins[] = {5,10,15,20,25,30,35,50,75};
//    int nHTBins = 28;
//    double htBins[] = {100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,700,800,900,1000,1100,1200,1600,2000};
//
//    std::vector<TString> bkgNames = {"ttbar","diboson","wjets","ttX","singlet","zjets"};
//
//    plotTurnons("muonSF",TString::Format("turnOn_muon_ht"),"GL_passSE",{"singlee"},bkgNames,"mupt_26","ht","passSMuoHtMuoBu"        ,0,nHTBins,htBins);
//    plotTurnons("electronSF",TString::Format("turnOn_electron_ht"),"GL_passSMu",{"singlemu"},bkgNames,"elpt_30","ht","passSEloHtEloBu"        ,0,nHTBins,htBins);
//
//    // plotTurnons(TString::Format("turnOn_ht_muon"),"GL_passSE","singlee","ttbar","ht_500","mu_pt","passSMuoHtMuoBu"        ,0,nLepBins,lepBins);
//    // plotTurnons(TString::Format("turnOn_ht_electron"),"GL_passSMu","singlemu","ttbar","ht_500","el_pt","passSEloHtEloBu"        ,0,nLepBins,lepBins);
//
//    delete of;
//}
