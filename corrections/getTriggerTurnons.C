#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"


#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/SignalHelper.h"

using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){

        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    void loadVariables() override {
        reader_event       =loadReader<EventReader>   ("event",isRealData());
        reader_jet         =loadReader<JetReader>     ("ak4Jet",isRealData());
        reader_electron    =loadReader<ElectronReader>("electron");
        reader_muon        =loadReader<MuonReader>    ("muon",isRealData());
        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData(),true,true);

        if(!isRealData()){
            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
        }

        checkConfig();
    }
    bool passTrig(Triggers_2017 trig) {return doesPass(triggerAccepts,trig);} // !! Need to use correct trigger enum for whatever era it is !!

    bool passTrig16(Triggers_2016 trig) {return doesPass(triggerAccepts,trig);} // !! Need to use correct trigger enum for whatever era it is !!
    bool passTrig17(Triggers_2017 trig) {return doesPass(triggerAccepts,trig);} // !! Need to use correct trigger enum for whatever era it is !!
    bool passTrig18(Triggers_2018 trig) {return doesPass(triggerAccepts,trig);} // !! Need to use correct trigger enum for whatever era it is !!

    void makeHTPlots(const TString& prefix, const TString& varname, float sel, float pt){

    	// !! input variable sel no longer used, reinstate it or get rid of it!!  ///////////////////

    	auto go = [&](float htsel, TString htid) {
            plotter.getOrMake1DPre(prefix,TString::Format("ht_incl_%s",(varname+htid).Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
            auto mkht = [&](float htCut){
                if(htsel >= htCut) plotter.getOrMake1DPre(prefix,TString::Format("ht_%.0f_%s",htCut,(varname+htid).Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
            };
            auto mkhtr = [&](float htCutmin, float htcutMax){
                if(htsel >= htCutmin && sel < htcutMax ) plotter.getOrMake1DPre(prefix,TString::Format("ht_%.0fto%.0f_%s",htCutmin,htcutMax,(varname+htid).Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
            };
            mkht(350);
            mkht(375);
            mkht(400);
            mkht(425);
            mkht(450);
            mkht(475);
            mkht(500);
            mkht(525);
            mkht(550);
            mkht(575);
            mkht(600);
            mkht(650);
            mkht(700);
            mkht(800);
            mkht(1000);
            mkht(1200);

            mkhtr(500,600);
            mkhtr(600,700);
            mkhtr(700,800);
            mkhtr(800,900);
            mkhtr(900,1000);
            mkhtr(1000,1200);
    	};

//    	go(ht_chs,"_chs");
    	go(ht,"");
    }

    void makeLepPlots(const TString& prefix, const TString& varname, TString selname, float sel){

    	auto go = [&](float HT, TString htid) {
        	plotter.getOrMake1DPre(prefix,TString::Format("%s_incl_%s",selname.Data(),(varname+htid).Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(HT,weight);

            auto mklp = [&](float lepcut){
                if(sel >= lepcut) plotter.getOrMake1DPre(prefix,TString::Format("%s_%.0f_%s",selname.Data(),lepcut,(varname+htid).Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(HT,weight);
            };
            auto mklpr = [&](float lepcutmin, float lepcutmax){
                if(sel >= lepcutmin && sel < lepcutmax) plotter.getOrMake1DPre(prefix,TString::Format("%s_%.0fto%.0f_%s",selname.Data(), lepcutmin,lepcutmax,(varname+htid).Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(HT,weight);
            };

            mklp(10);
            mklp(15);
            mklp(20);
            mklp(25);
            mklp(26);
            mklp(27);
            mklp(30);
            mklp(31);
            mklp(35);
            mklp(40);
            mklp(50);
            mklp(75);
            mklp(100);

            mklpr(15,20);
            mklpr(20,25);
            mklpr(25,30);
            mklpr(26,30);
            mklpr(27,30);
            mklpr(30,35);
            mklpr(35,40);
            mklpr(40,50);
            mklpr(50,75);
            mklpr(75,100);

            mklpr(30,40);
            mklpr(30,50);

            mklpr(40,75);
            mklpr(50,100);
    	};

//    	go(ht_chs,"_chs");
    	go(ht,"");

    }

    void testEachTriggerIndividually(const TString& prefix, bool doMuon, std::vector<const Muon*> tagMuons, std::vector<const Muon*> probeMuons, std::vector<const Electron*> tagElectrons, std::vector<const Electron*> probeElectrons) {
    	float maxLepPt, eta;
    	if (doMuon) {
    		if (!tagElectrons.size()) return;
    		maxLepPt = probeMuons.size() ? probeMuons.front()->pt() : 0;
    		eta = probeMuons.size() ? probeMuons.front()->absEta() : 99;
    	} else {
    		if (!tagMuons.size()) return;
    		maxLepPt = probeElectrons.size() ? probeElectrons.front()->pt() : 0;
    		eta = probeElectrons.size() ? probeElectrons.front()->absEta() : 99;
    	}

		makeHTPlots(prefix+"_TrigIncl_", doMuon ? "mu_pt" : "el_pt" ,ht,maxLepPt);
		makeLepPlots(prefix+"_TrigIncl_","ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
		if (eta < 2.1) {
			makeHTPlots(prefix+"_TrigIncl_eta2p1_", doMuon ? "mu_pt" : "el_pt" ,ht,maxLepPt);
			makeLepPlots(prefix+"_TrigIncl_eta2p1_","ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
		}
		if (eta < 1.479) {
			makeHTPlots(prefix+"_TrigIncl_eta1p479_", doMuon ? "mu_pt" : "el_pt" ,ht,maxLepPt);
			makeLepPlots(prefix+"_TrigIncl_eta1p479_","ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
		}

    	for (Triggers_2017 trg=(Triggers_2017)0; trg != HLT17_NTrig; trg=(Triggers_2017)(trg+1)) {
    		if (passTrig(trg)) {
        		TString preName = prefix + "_passTrig_"+TString::Format("%i_",int(trg));
        		TString varname = doMuon ? "mu_pt" : "el_pt";
    			makeHTPlots(preName,varname,ht,maxLepPt);
    			makeLepPlots(preName,"ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
    			if (eta < 2.1) {
        			makeHTPlots(preName+"eta2p1_",varname,ht,maxLepPt);
        			makeLepPlots(preName+"eta2p1_","ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
    			}
    			if (eta < 1.5) {
        			makeHTPlots(preName+"eta1p5_",varname,ht,maxLepPt);
        			makeLepPlots(preName+"eta1p5_","ht", doMuon ? "mupt" : "elpt" ,maxLepPt);
    			}
    		}
    	}
    }

    void doElectronLeg(const TString& prefix, std::vector<const Muon*> tagMuons, std::vector<const Electron*> probeElectrons){
    	// trigger to take tag muon
    	if( !passTrig(HLT17_IsoMu27) ) return;
        if(!tagMuons.size()) return;

        float maxLepPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;
        TString preName = prefix + "_passSMu_";
        makeHTPlots(preName,"el_pt",ht,maxLepPT);

        bool passEle35 = passTrig(FillerConstants::HLT17_Ele35_WPTight_Gsf);
        bool passEle32 = passTrig(FillerConstants::HLT17_Ele32_WPTight_Gsf);
        bool passEle32Double = passTrig(FillerConstants::HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross1 = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        bool passECross2 = passTrig(FillerConstants::HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        bool passECross = passECross1 || passECross2;

        if(passEle35) makeHTPlots(preName + "passEl35_","el_pt",ht,maxLepPT);
        if(passEle32) makeHTPlots(preName + "passEl32_","el_pt",ht,maxLepPT);
        if(passEle32Double) makeHTPlots(preName + "passEl32Dbl_","el_pt",ht,maxLepPT);
        if(passEle32Double || passEle32) makeHTPlots(preName + "passEl32SngoDbl_","el_pt",ht,maxLepPT);
        if(passEle32Double || passEle32 || passEle35) makeHTPlots(preName + "passEl32SngoDblo35_","el_pt",ht,maxLepPT);
        if(passECross1) makeHTPlots(preName + "passElHT1_","el_pt",ht,maxLepPT);
        if(passECross2) makeHTPlots(preName + "passElHT2_","el_pt",ht,maxLepPT);
        if(passECross) makeHTPlots(preName + "passElHT_","el_pt",ht,maxLepPT);
    }

    void doMuonLeg(const TString& prefix, std::vector<const Electron*> tagElectrons, std::vector<const Muon*> probeMuons){
    	// triggers to take a tag electron
        if( !passTrig(HLT17_Ele35_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG) ) return;
        if(!tagElectrons.size()) return;

        float maxLepPT = probeMuons.size() ? probeMuons.front()->pt() : 0;
        TString preName = prefix + "_passSE_";
        makeHTPlots(preName,"mu_pt",ht,maxLepPT);

        bool passMu27 = passTrig(FillerConstants::HLT17_IsoMu27);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);

        if(passMu27) makeHTPlots(preName + "passMu27_","mu_pt",ht,maxLepPT);
        if(passMCross) makeHTPlots(preName + "passMuHT_","mu_pt",ht,maxLepPT);
    }

    void doHTLegWithMuonDenom(const TString& prefix,std::vector<const Muon*> tagMuons, std::vector<const Electron*> probeElectrons){
    	if( !passTrig(HLT17_IsoMu27) ) return;
        if(!tagMuons.size()) return;

        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_passSMu_";
        makeLepPlots(preName,"ht","elpt",maxOtherPT);
        makeLepPlots(preName,"ht","mupt",maxSamePT);

        bool passEle35 = passTrig(FillerConstants::HLT17_Ele35_WPTight_Gsf);
        bool passEle32 = passTrig(FillerConstants::HLT17_Ele32_WPTight_Gsf);
        bool passEle32Double = passTrig(FillerConstants::HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross1 = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        bool passECross2 = passTrig(FillerConstants::HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        bool passECross = passECross1 || passECross2;

        bool passHighE = passTrig(FillerConstants::HLT17_Photon200) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
        bool passMu50 = passTrig(FillerConstants::HLT17_Mu50);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passJet = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);

        if(passEle35) makeLepPlots(preName + "passEl35_","ht","elpt",maxOtherPT);
        if(passEle32) makeLepPlots(preName + "passEl32_","ht","elpt",maxOtherPT);
        if(passEle32Double) makeLepPlots(preName + "passEl32Dbl_","ht","elpt",maxOtherPT);
        if(passEle32Double || passEle32) makeLepPlots(preName + "passEl32SngoDbl_","ht","elpt",maxOtherPT);
        if(passEle32Double || passEle32 || passEle35) makeLepPlots(preName + "passEl32SngoDblo35_","ht","elpt",maxOtherPT);

        if(passECross1) makeLepPlots(preName + "passElHT1_","ht","elpt",maxOtherPT);
        if(passECross2) makeLepPlots(preName + "passElHT2_","ht","elpt",maxOtherPT);
        if(passECross) makeLepPlots(preName + "passElHT_","ht","elpt",maxOtherPT);
        if(passMCross) makeLepPlots(preName + "passMuHT_","ht","mupt",maxSamePT);

        if(passHighE) makeLepPlots(preName + "passHighE_","ht","mupt",maxSamePT);
        if(passMu50) makeLepPlots(preName + "passMu50_","ht","mupt",maxSamePT);

        if(passECross|| passJet) makeLepPlots(preName + "passElHToHad_","ht","elpt",maxOtherPT);
        if(passMCross|| passJet) makeLepPlots(preName + "passMuHToHad_","ht","mupt",maxSamePT);

        if(passECross || passHighE || passJet) makeLepPlots(preName + "passElHToHEoHad_","ht","elpt",maxOtherPT);
        if(passMCross || passMu50 || passJet) makeLepPlots(preName + "passMuHToHMoHad_","ht","mupt",maxSamePT);
    }

    void doHTLegWithElDenom(const TString& prefix,std::vector<const Electron*> tagElectrons, std::vector<const Muon*> probeMuons){
        if( !passTrig(HLT17_Ele35_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG) ) return;
        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_passSE_";
        makeLepPlots(preName,"ht","mupt",maxOtherPT);
        makeLepPlots(preName,"ht","elpt",maxSamePT);

        bool passHighE = passTrig(FillerConstants::HLT17_Photon200) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
        bool passMu50 = passTrig(FillerConstants::HLT17_Mu50);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        bool passJet = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);

        if(passECross) makeLepPlots(preName + "passElHT_","ht","elpt",maxOtherPT);
        if(passMCross) makeLepPlots(preName + "passMuHT_","ht","mupt",maxSamePT);

        if(passHighE) makeLepPlots(preName + "passHighE_","ht","mupt",maxSamePT);
        if(passMu50) makeLepPlots(preName + "passMu50_","ht","mupt",maxSamePT);

        if(passECross|| passJet) makeLepPlots(preName + "passElHToHad_","ht","elpt",maxOtherPT);
        if(passMCross|| passJet) makeLepPlots(preName + "passMuHToHad_","ht","mupt",maxSamePT);

        if(passECross || passHighE || passJet) makeLepPlots(preName + "passElHToHEoHad_","ht","elpt",maxOtherPT);
        if(passMCross || passMu50 || passJet) makeLepPlots(preName + "passMuHToHMoHad_","ht","mupt",maxSamePT);
    }

    void doGrandLeptonWElDenom(const TString& prefix,std::vector<const Electron*> tagElectrons, std::vector<const Muon*> probeMuons){
        bool passEl = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
        		|| passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);

        if( !passEl ) return;
        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSE_";

        bool passSMu = passTrig(HLT17_IsoMu27);
        bool passHighMu = passTrig(HLT17_Mu50);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);

        bool passJetHT = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);
        bool passMet_1 = passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight) || passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
        bool passMet_2 = passTrig(HLT17_PFMET120_PFMHT120_IDTight) || passTrig(HLT17_PFMET120_PFMHT120_IDTight_PFHT60);
        bool passMetNoMu = passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

        bool passMET = passMet_1 || passMet_2 || passMetNoMu;
        bool passBu = passJetHT || passMET;

        bool passSMuoHM = passSMu || passHighMu;
        bool passSMuoHtMu = passSMu || passMCross;
        bool passSMuoHtMuoHM = passSMuoHtMu || passHighMu;
        bool passSMuoHtMuoHMoBu = passSMuoHtMuoHM || passBu;
        bool passMuDenNoCross = passBu || passSMu;

        bool passNomMuoJet = passSMuoHtMuoHM || passJetHT;
        bool passNomMuoJetoMetNoMu = passNomMuoJet || passMetNoMu;
        bool passNomMuoJetoMET = passNomMuoJet || passMET;
        bool passNomMuoJetoMet1 = passNomMuoJet || passMet_1;
        bool passNomMuoJetoMet2 = passNomMuoJet || passMet_2;

        auto makePlots = [&](TString pre, float pt) {
            //HT side
            makeLepPlots(pre,"ht","mupt",pt);
            if(passSMu )           makeLepPlots(pre +"passSMu_","ht","mupt",pt);
            if(passSMuoHtMu)       makeLepPlots(pre +"passSMuoHtMu_","ht","mupt",pt);
            if(passSMuoHM)         makeLepPlots(pre +"passSMuoHM_","ht","mupt",pt);
            if(passSMuoHtMuoHM)    makeLepPlots(pre +"passSMuoHtMuoHM_","ht","mupt",pt);
            if(passSMuoHtMuoHMoBu) makeLepPlots(pre +"passSMuoHtMuoHMoBu_","ht","mupt",pt);
            if(passMuDenNoCross)   makeLepPlots(pre +"passMuDenNoCross_","ht","mupt",pt);

            if(passNomMuoJet)             makeLepPlots(pre +"passNomMuoJet_","ht","mupt",pt);
            if(passNomMuoJetoMetNoMu)     makeLepPlots(pre +"passNomMuoJetoMetNoMu_","ht","mupt",pt);
            if(passNomMuoJetoMet1)        makeLepPlots(pre +"passNomMuoJetoMet1_","ht","mupt",pt);
            if(passNomMuoJetoMet2)        makeLepPlots(pre +"passNomMuoJetoMet2_","ht","mupt",pt);
            if(passNomMuoJetoMET)         makeLepPlots(pre +"passNomMuoJetoMET_","ht","mupt",pt);


            //Mu side
            makeHTPlots(pre,"mu_pt",ht,pt);
            if(passSMu )           makeHTPlots(pre+"passSMu_"         ,"mu_pt",ht,pt);
            if(passSMuoHtMu)       makeHTPlots(pre+"passSMuoHtMu_"    ,"mu_pt",ht,pt);
            if(passSMuoHM)         makeHTPlots(pre+"passSMuoHM_"    ,"mu_pt",ht,pt);
            if(passSMuoHtMuoHM)    makeHTPlots(pre+"passSMuoHtMuoHM_"    ,"mu_pt",ht,pt);
            if(passSMuoHtMuoHMoBu) makeHTPlots(pre+"passSMuoHtMuoHMoBu_","mu_pt",ht,pt);
            if(passMuDenNoCross)   makeHTPlots(pre+"passMuDenNoCross_","mu_pt",ht,pt);

            if(passNomMuoJet)             makeHTPlots(pre+"passNomMuoJet_","mu_pt",ht,pt);
            if(passNomMuoJetoMetNoMu)     makeHTPlots(pre+"passNomMuoJetoMetNoMu_","mu_pt",ht,pt);
            if(passNomMuoJetoMet1)        makeHTPlots(pre+"passNomMuoJetoMet1_","mu_pt",ht,pt);
            if(passNomMuoJetoMet2)        makeHTPlots(pre+"passNomMuoJetoMet2_","mu_pt",ht,pt);
            if(passNomMuoJetoMET)         makeHTPlots(pre+"passNomMuoJetoMET_","mu_pt",ht,pt);


            //2D
            if(pt > 0){
                plotter.getOrMake2DPre(pre,"mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSMu )            plotter.getOrMake2DPre(pre+"passSMu_"         ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSMuoHtMu)        plotter.getOrMake2DPre(pre+"passSMuoHtMu_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSMuoHM)          plotter.getOrMake2DPre(pre+"passSMuoHM_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSMuoHtMuoHM)     plotter.getOrMake2DPre(pre+"passSMuoHtMuoHM_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSMuoHtMuoHMoBu)  plotter.getOrMake2DPre(pre+"passSMuoHtMuoHMoBu_","mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
            }
        };

        makePlots(preName,maxOtherPT);
        if (probeMuons.size()) {
            if (probeMuons.front()->absEta() < 1.479) makePlots(preName+"maxEta1p479_",maxOtherPT);
            if (probeMuons.front()->absEta() < 2.1) makePlots(preName+"maxEta2p1_",maxOtherPT);
            if (probeMuons.front()->absEta() < 2.5) makePlots(preName+"maxEta2p5_",maxOtherPT);
        }
    }

    void doGrandLeptonWMuDenom(const TString& prefix,std::vector<const Muon*> tagMuons, std::vector<const Electron*> probeElectrons){
        if( !passTrig(HLT17_IsoMu27) ) return;
        if(!tagMuons.size()) return;

        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSMu_";

        bool passSEl = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
                    || passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450) || passTrig(FillerConstants::HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        bool passPh200 = passTrig(FillerConstants::HLT17_Photon200);
        bool passEle115 = passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
        bool passHighE = passPh200 || passEle115;

        bool passJetHT = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);
        bool passMet_1 = passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight) || passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
        bool passMet_2 = passTrig(HLT17_PFMET120_PFMHT120_IDTight) || passTrig(HLT17_PFMET120_PFMHT120_IDTight_PFHT60);
        bool passMET = passMet_1 || passMet_2;
        bool passSElSuite = passSEl || passECross || passEle115;

        bool passBu = passJetHT || passMET;

        bool passSEloHtEl = passSEl  || passECross;
        bool passSEloHE = passSEl  || passHighE;
        bool passSEloHtEloHE = passSEloHtEl || passSEloHE;
        bool passSEloHtEloHEoBu = passSEloHtEloHE || passBu;
        bool passElDenNoCross = passBu || passSEl;

        bool passNomEloJet = passSElSuite || passJetHT;
        bool passNomEloJetoMet1 = passNomEloJet || passMet_1;
        bool passNomEloJetoMet2 = passNomEloJet || passMet_2;
        bool passNomEloJetoMET = passNomEloJet || passMET;
        bool passNomEloJetoPh = passNomEloJet || passPh200;
        bool passNomEloJetoMEToPh = passNomEloJetoMET || passPh200;

        auto makePlots = [&](TString pre, float pt) {
            //HT side
            makeLepPlots(pre,"ht","elpt",pt);
            if(passSEl )           makeLepPlots(pre +"passSEl_","ht","elpt",pt);
            if(passSEloHtEl)       makeLepPlots(pre +"passSEloHtEl_","ht","elpt",pt);
            if(passSEloHE)         makeLepPlots(pre +"passSEloHE_","ht","elpt",pt);
            if(passSEloHtEloHE)    makeLepPlots(pre +"passSEloHtEloHE_","ht","elpt",pt);
            if(passSEloHtEloHEoBu) makeLepPlots(pre +"passSEloHtEloHEoBu_","ht","elpt",pt);

            if(passElDenNoCross)   makeLepPlots(pre +"passElDenNoCross_","ht","elpt",pt);
            if(passSEl || passPh200) makeLepPlots(pre +"passSEloPh200_","ht","elpt",pt);
            if(passSEl || passEle115) makeLepPlots(pre +"passSEloEl115_","ht","elpt",pt);

            if(passNomEloJet)         makeLepPlots(pre +"passNomEloJet_","ht","elpt",pt);
            if(passNomEloJetoMET)     makeLepPlots(pre +"passNomEloJetoMET_","ht","elpt",pt);
            if(passNomEloJetoPh)      makeLepPlots(pre +"passNomEloJetoPh_","ht","elpt",pt);
            if(passNomEloJetoMEToPh)  makeLepPlots(pre +"passNomEloJetoMEToPh_","ht","elpt",pt);


            //El side
            makeHTPlots(pre,"el_pt",ht,pt);
            if(passSEl )           makeHTPlots(pre+"passSEl_"        ,"el_pt",ht,pt);
            if(passSEloHtEl)       makeHTPlots(pre+"passSEloHtEl_"   ,"el_pt",ht,pt);
            if(passSEloHE)         makeHTPlots(pre+"passSEloHE_"   ,"el_pt",ht,pt);
            if(passSEloHtEloHE)    makeHTPlots(pre+"passSEloHtEloHE_"   ,"el_pt",ht,pt);
            if(passSEloHtEloHEoBu) makeHTPlots(pre+"passSEloHtEloHEoBu_","el_pt",ht,pt);

            if(passElDenNoCross)   makeHTPlots(pre+"passElDenNoCross_","el_pt",ht,pt);
            if(passSEl || passPh200) makeHTPlots(pre +"passSEloPh200_","el_pt",ht,pt);
            if(passSEl || passEle115) makeHTPlots(pre +"passSEloEl115_","el_pt",ht,pt);

            if(passNomEloJet)          makeHTPlots(pre+"passNomEloJet_","el_pt",ht,pt);
            if(passNomEloJetoMET)      makeHTPlots(pre+"passNomEloJetoMET_","el_pt",ht,pt);
            if(passNomEloJetoPh)       makeHTPlots(pre+"passNomEloJetoPh_","el_pt",ht,pt);
            if(passNomEloJetoMEToPh)   makeHTPlots(pre+"passNomEloJetoMEToPh_","el_pt",ht,pt);

            //2D
            if(pt > 0){
                plotter.getOrMake2DPre(pre,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSEl )           plotter.getOrMake2DPre(pre+"passSEl_"        ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSEloHtEl)       plotter.getOrMake2DPre(pre+"passSEloHtEl_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSEloHE)         plotter.getOrMake2DPre(pre+"passSEloHE_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSEloHtEloHE)    plotter.getOrMake2DPre(pre+"passSEloHtEloHE_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);
                if(passSEloHtEloHEoBu) plotter.getOrMake2DPre(pre+"passSEloHtEloHEoBu_","el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(pt,ht,weight);

            }
        };

        makePlots(preName,maxOtherPT);
        if (probeElectrons.size()) {
            if (probeElectrons.front()->absEta() < 1.479) makePlots(preName+"maxEta1p479_",maxOtherPT);
            if (probeElectrons.front()->absEta() < 2.1) makePlots(preName+"maxEta2p1_",maxOtherPT);
            if (probeElectrons.front()->absEta() < 2.5) makePlots(preName+"maxEta2p5_",maxOtherPT);
            if (probeElectrons.front()->absEta() > 2.1 && probeElectrons.front()->absEta() < 2.5) makePlots(preName+"eta2p1to2p5_",maxOtherPT);
        }
    }

    void doMCLepton(const TString& prefix, std::vector<const Muon*> probeMuons, std::vector<const Electron*> probeElectrons){
        float maxMu = probeMuons.size() ? probeMuons.front()->pt() : 0;
        float maxEl = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_MC_";

        bool passJetHT = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);
        bool passMet_1 = passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight) || passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
        bool passMet_2 = passTrig(HLT17_PFMET120_PFMHT120_IDTight) || passTrig(HLT17_PFMET120_PFMHT120_IDTight_PFHT60);
        bool passMetNoMu = passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

        bool passBu = passJetHT || passMet_1 || passMet_2 || passMetNoMu;

        bool passSMu = passTrig(HLT17_IsoMu27);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passHighMu = passTrig(FillerConstants::HLT17_Mu50);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passSMuoHM = passSMu || passHighMu;
        bool passSMuoHtMuoHM = passSMuoHtMu || passHighMu;
        bool passSMuoHtMuoHMoBu = passSMuoHtMuoHM || passBu;

        bool passSEl_ETA = passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
        bool passSEl = passSEl_ETA || passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450) || passTrig(FillerConstants::HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        bool passHighE = passTrig(FillerConstants::HLT17_Photon200) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
        bool passSEloHtEl = passSEl || passECross;
        bool passSEloHE = passSEl  || passHighE;
        bool passSEloHtEloHE = passSEloHtEl || passSEloHE;
        bool passSEloHtEloHEoBu = passSEloHtEloHE || passBu;

        bool passMuDenNoCross = passBu || passSMu;
        bool passElDenNoCross = passBu || passSEl;


        // el eta trig comparison
        bool passEl28 = passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150);
        bool passEl30 = passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
        bool passElEta = passEl28 || passEl30;

        //HT side
        makeLepPlots(preName,"ht","mupt",maxMu);

//        if(passSMu ) makeLepPlots(preName +"passSMu_","ht","mupt",maxMu);
//        if(passSMuoHtMu) makeLepPlots(preName +"passSMuoHtMu_","ht","mupt",maxMu);
//        if(passSMuoHM) makeLepPlots(preName +"passSMuoHM_","ht","mupt",maxMu);
//        if(passSMuoHtMuoHM) makeLepPlots(preName +"passSMuoHtMuoHM_","ht","mupt",maxMu);
//        if(passSMuoHtMuoHMoBu) makeLepPlots(preName +"passSMuoHtMuoHMoBu_","ht","mupt",maxMu);
//        if(passMuDenNoCross) makeLepPlots(preName +"passMuDenNoCross_","ht","mupt",maxMu);

        makeLepPlots(preName,"ht","elpt",maxEl);
        if(passEl28) makeLepPlots(preName +"passEl28_","ht","elpt",maxEl);
        if(passEl30) makeLepPlots(preName +"passEl30_","ht","elpt",maxEl);
        if(passElEta) makeLepPlots(preName +"passElEta_","ht","elpt",maxEl);

//        if(passSEl ) makeLepPlots(preName +"passSEl_","ht","elpt",maxEl);
//        if(passSEloHtEl) makeLepPlots(preName +"passSEloHtEl_","ht","elpt",maxEl);
//        if(passSEloHE) makeLepPlots(preName +"passSEloHE_","ht","elpt",maxEl);
//        if(passSEloHtEloHE) makeLepPlots(preName +"passSEloHtEloHE_","ht","elpt",maxEl);
//        if(passSEloHtEloHEoBu) makeLepPlots(preName +"passSEloHtEloHEoBu_","ht","elpt",maxEl);
//        if(passElDenNoCross) makeLepPlots(preName +"passElDenNoCross_","ht","elpt",maxEl);

        //Lep side
        makeHTPlots(preName,"mu_pt",ht,maxMu);

//        if(passSMu )           makeHTPlots(preName+"passSMu_"         ,"mu_pt",ht_puppi,maxMu);
//        if(passSMuoHtMu)       makeHTPlots(preName+"passSMuoHtMu_"    ,"mu_pt",ht_puppi,maxMu);
//        if(passSMuoHM)         makeHTPlots(preName+"passSMuoHM_"    ,"mu_pt",ht_puppi,maxMu);
//        if(passSMuoHtMuoHM)    makeHTPlots(preName+"passSMuoHtMuoHM_"    ,"mu_pt",ht_puppi,maxMu);
//        if(passSMuoHtMuoHMoBu) makeHTPlots(preName+"passSMuoHtMuoHMoBu_","mu_pt",ht_puppi,maxMu);
//        if(passMuDenNoCross)   makeHTPlots(preName+"passMuDenNoCross_","mu_pt",ht_puppi,maxMu);

        makeHTPlots(preName,"el_pt",ht,maxEl);
//        if(passSEl )           makeHTPlots(preName+"passSEl_"        ,"el_pt",ht_puppi,maxEl);
//        if(passSEloHtEl)       makeHTPlots(preName+"passSEloHtEl_"   ,"el_pt",ht_puppi,maxEl);
//        if(passSEloHE)         makeHTPlots(preName+"passSEloHE_"   ,"el_pt",ht_puppi,maxEl);
//        if(passSEloHtEloHE)    makeHTPlots(preName+"passSEloHtEloHE_"   ,"el_pt",ht_puppi,maxEl);
//        if(passSEloHtEloHEoBu) makeHTPlots(preName+"passSEloHtEloHEoBu_","el_pt",ht_puppi,maxEl);
//        if(passElDenNoCross)   makeHTPlots(preName+"passElDenNoCross_","el_pt",ht_puppi,maxEl);
/*
        if(maxEl > 0){
            plotter.getOrMake2DPre(preName,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
            if(passSMu )            plotter.getOrMake2DPre(preName+"passSMu_"         ,"mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
            if(passSMuoHtMu)        plotter.getOrMake2DPre(preName+"passSMuoHtMu_"    ,"mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
            if(passSMuoHM)          plotter.getOrMake2DPre(preName+"passSMuoHM_"    ,"mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
            if(passSMuoHtMuoHM)     plotter.getOrMake2DPre(preName+"passSMuoHtMuoHM_"    ,"mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
            if(passSMuoHtMuoHMoBu)  plotter.getOrMake2DPre(preName+"passSMuoHtMuoHMoBu_","mu_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_puppi,weight);
        }

        if(maxMu > 0){
            plotter.getOrMake2DPre(preName,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
            if(passSEl )           plotter.getOrMake2DPre(preName+"passSEl_"        ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
            if(passSEloHtEl)       plotter.getOrMake2DPre(preName+"passSEloHtEl_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
            if(passSEloHE)         plotter.getOrMake2DPre(preName+"passSEloHE_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
            if(passSEloHtEloHE)    plotter.getOrMake2DPre(preName+"passSEloHtEloHE_"   ,"el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
            if(passSEloHtEloHEoBu) plotter.getOrMake2DPre(preName+"passSEloHtEloHEoBu_","el_pt_v_ht",";lepton p_{T} [GeV]; #it{H}_{T} [GeV]",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_puppi,weight);
        }
*/
    }

    void doSystLepton(const TString& prefix){

    	if (!selectedLepton) return;
        TString preName = prefix + "_SYST_";

        bool passJetHT = passTrig(HLT17_PFHT1050)|| passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30);
        bool passMet_1 = passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight) || passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
        bool passMet_2 = passTrig(HLT17_PFMET120_PFMHT120_IDTight) || passTrig(HLT17_PFMET120_PFMHT120_IDTight_PFHT60);
        bool passMetNoMu = passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

        bool passBu = passJetHT || passMet_1 || passMet_2 || passMetNoMu;

        bool passSMu = passTrig(HLT17_IsoMu27);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passHighMu = passTrig(FillerConstants::HLT17_Mu50);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passSMuoHM = passSMu || passHighMu;
        bool passSMuoHtMuoHM = passSMuoHtMu || passHighMu;
        bool passSMuoHtMuoHMoBu = passSMuoHtMuoHM || passBu;

        bool passSEl_ETA = passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
        bool passSEl = passSEl_ETA || passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450) || passTrig(FillerConstants::HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        bool passHighE = passTrig(FillerConstants::HLT17_Photon200) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
        bool passSEloHtEl = passSEl || passECross;
        bool passSEloHE = passSEl  || passHighE;
        bool passSEloHtEloHE = passSEloHtEl || passSEloHE;
        bool passSEloHtEloHEoBu = passSEloHtEloHE || passBu;

        if (selectedLepton->isMuon()) {
        	if (selectedLepton->pt() < 27 || selectedLepton->absEta() > 2.4) return;
            makeLepPlots(preName,"ht","mupt",selectedLepton->pt());
            makeHTPlots(preName,"mu_pt",ht,selectedLepton->pt());

            if (selectedLeptons.size() == 1) {
                makeLepPlots(preName+"1lep_","ht","mupt",selectedLepton->pt());
                makeHTPlots(preName+"1lep_","mu_pt",ht,selectedLepton->pt());
            }
            if (passSMuoHtMuoHMoBu) {
            	makeLepPlots(preName +"passSMuoHtMuoHMoBu_","ht","mupt",selectedLepton->pt());
            	makeHTPlots(preName+"passSMuoHtMuoHMoBu_","mu_pt",ht,selectedLepton->pt());
            }
        } else {
        	if (selectedLepton->pt() < 30 || selectedLepton->absEta() > 1.479) return;
            makeLepPlots(preName,"ht","elpt",selectedLepton->pt());
            makeHTPlots(preName,"el_pt",ht,selectedLepton->pt());

            if (selectedLeptons.size() == 1) {
                makeLepPlots(preName+"1lep_","ht","elpt",selectedLepton->pt());
                makeHTPlots(preName+"1lep_","el_pt",ht,selectedLepton->pt());
            }
            if (passSEloHtEloHEoBu) {
            	makeLepPlots(preName +"passSEloHtEloHEoBu_","ht","elpt",selectedLepton->pt());
            	makeHTPlots(preName+"passSEloHtEloHEoBu_","el_pt",ht,selectedLepton->pt());
            }
        }
    }

    void doMuDenomPlotsForSF(const TString& prefix,std::vector<const Muon*> tagMuons, std::vector<const Electron*> probeElectrons){
    	if (!tagMuons.size()) return;
    	bool passDenomTrig = false;

    	if (FillerConstants::DataEra(*reader_event->dataEra) == ERA_2016) {
    		passDenomTrig = passTrig16(HLT16_IsoMu24);
    	} else if(FillerConstants::DataEra(*reader_event->dataEra) == ERA_2017) {
    		passDenomTrig = passTrig17(HLT17_IsoMu27);
    	} else if(FillerConstants::DataEra(*reader_event->dataEra) == ERA_2016) {
    		passDenomTrig = passTrig18(HLT18_IsoMu24);
    	}

    	if (!passDenomTrig) return;

        float maxTagPT  = tagMuons.front()->pt();
        float maxProbePT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSMu_";

        bool passProbePath = false;
        bool passJet = false;
        bool passMET = false;

        switch(FillerConstants::DataEra(*reader_event->dataEra)) {
        case ERA_2016:
        	passJet = passTrig16(HLT16_AK8PFJet450) || passTrig16(HLT16_AK8PFJet360_TrimMass30) ||
				passTrig16(HLT16_PFHT800) || passTrig16(HLT16_PFHT900);
        	passMET = true;
        	break;
        case ERA_2017:
        	passJet = passTrig17(HLT17_AK8PFHT850_TrimMass50) || passTrig17(HLT17_PFHT1050) ||
				passTrig17(HLT17_AK8PFJet400_TrimMass30) || passTrig17(HLT17_AK8PFJet500);

        	break;
        case ERA_2018:
        	passJet = passTrig18(HLT18_AK8PFJet400_TrimMass30) || passTrig18(HLT18_PFHT1050) ||
				passTrig18(HLT18_AK8PFHT800_TrimMass50) || passTrig18(HLT18_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) ||
				passTrig18(HLT18_AK8PFJet330_TrimMass30_PFAK8BTagCSV_p1);
        	break;
        default:
        	throw std::invalid_argument("data era needs to be set properly");
        }


    }

    void doElDenomPlotsForSF(const TString& prefix,std::vector<const Electron*> tagElectrons, std::vector<const Muon*> probeMuons){
    	if (!tagElectrons.size()) return;
    	bool passDenomTrig = false;

    	if (FillerConstants::DataEra(*reader_event->dataEra) == ERA_2016) {
    		passDenomTrig = passTrig16(HLT16_Ele27_WPTight_Gsf) ||
    				passTrig16(HLT16_Ele25_eta2p1_WPTight_Gsf);
    	} else if(FillerConstants::DataEra(*reader_event->dataEra) == ERA_2017) {
    		passDenomTrig = passTrig17(HLT17_Ele27_WPTight_Gsf) ||
    				passTrig17(HLT17_Ele32_WPTight_Gsf) ||
					passTrig17(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) ||
					passTrig17(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
    	} else if(FillerConstants::DataEra(*reader_event->dataEra) == ERA_2016) {
    		passDenomTrig = passTrig18(HLT18_Ele28_eta2p1_WPTight_Gsf_HT150) ||
    				passTrig18(HLT18_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned) ||
					passTrig18(HLT18_Ele32_WPTight_Gsf);
    	}

    	if (!passDenomTrig) return;

        float maxTagPT  = tagElectrons.front()->pt();
        float maxProbePT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_SF_passSEl_";

        bool passProbePath = false;
        bool passJet = false;
        bool passMET = false;
    }

    void studyDileptonTriggerEffs(TString prefix) {
    	if (lepChan != DILEP) return;
    	if(ht < 400) return;
    	TString llS = DileptonProcessor::getDilepStr(dilep1,dilep2);

    	int nBins = 12;
    	double bins[] = {0,10,20,30,40,50,60,70,80,90,100,300,1000};

    	bool passJet = passTrig17(HLT17_AK8PFHT850_TrimMass50) || passTrig17(HLT17_PFHT1050) ||
				passTrig17(HLT17_AK8PFJet400_TrimMass30) || passTrig17(HLT17_AK8PFJet500);
    	if (!passJet) return;

    	auto plt = [&](TString idS) {
        	plotter.getOrMake2DPre(prefix+"_"+idS,"pt1pt2",";pt1;pt2",nBins,bins,nBins,bins)->Fill(dilep1->pt(),dilep2->pt(),weight);
        	plotter.getOrMake2DPre(prefix+llS+idS,"pt1pt2",";pt1;pt2",nBins,bins,nBins,bins)->Fill(dilep1->pt(),dilep2->pt(),weight);
        	plotter.getOrMake2DPre(prefix+"_"+idS,"pt1ht",";pt1;ht",nBins,bins,nHTBins,htBins)->Fill(dilep1->pt(),ht,weight);
        	plotter.getOrMake2DPre(prefix+llS+idS,"pt1ht",";pt1;ht",nBins,bins,nHTBins,htBins)->Fill(dilep1->pt(),ht,weight);
        	plotter.getOrMake2DPre(prefix+"_"+idS,"pt2ht",";pt2;ht",nBins,bins,nHTBins,htBins)->Fill(dilep2->pt(),ht,weight);
        	plotter.getOrMake2DPre(prefix+llS+idS,"pt2ht",";pt2;ht",nBins,bins,nHTBins,htBins)->Fill(dilep2->pt(),ht,weight);
    	};

    	plt("passHT");

    	bool passSMu = passTrig17(HLT17_IsoMu27) || passTrig17(HLT17_Mu50);
    	bool passSEl = passTrig17(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig17(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned)
    			|| passTrig17(HLT17_Ele27_WPTight_Gsf) || passTrig17(HLT17_Ele32_WPTight_Gsf) || passTrig17(HLT17_Ele35_WPTight_Gsf);
    	bool passMET = passTrig17(HLT17_PFMET120_PFMHT120_IDTight_PFHT60) || passTrig17(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)
    			|| passTrig17(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);

    	bool passElCross = passTrig17(HLT17_Ele15_IsoVVVL_PFHT450) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT)
						|| passTrig17(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
    	bool passMuCross = passTrig17(HLT17_Mu15_IsoVVVL_PFHT450);

    	bool passPath = passSMu || passSEl || passMET || passTrig17(HLT17_Photon200);
    	bool passPath_wCross = passPath || passElCross || passMuCross;

    	if(passPath) plt("passLep");
    	if(passPath_wCross) plt("passLep_wCross");

/*    	if (!passPath) {
    		std::cout<<"does pass with cross trig: "<<passPath_wCross<<std::endl;
    		std::cout<<"ht = "<<ht<<std::endl;
    		printf("recolep1: %s, pt = %f, eta = %f, phi = %f\n",dilep1->isMuon() ? "mu":"el",
    				dilep1->pt(),dilep1->eta(),dilep1->phi());
    		printf("recolep2: %s, pt = %f, eta = %f, phi = %f\n",dilep2->isMuon() ? "mu":"el",
    				dilep2->pt(),dilep2->eta(),dilep2->phi());

    		if (diHiggsEvt.type == DiHiggsEvent::DILEP) {
        		printf("genlep1: %s, pt = %f, eta = %f, phi = %f\n",diHiggsEvt.w1_d1->absPdgId()==13 ? "mu":"el",
        				diHiggsEvt.w1_d1->pt(),diHiggsEvt.w1_d1->eta(),diHiggsEvt.w1_d1->phi());
        		printf("genlep2: %s, pt = %f, eta = %f, phi = %f\n",diHiggsEvt.w2_d1->absPdgId()==13 ? "mu":"el",
        				diHiggsEvt.w2_d1->pt(),diHiggsEvt.w2_d1->eta(),diHiggsEvt.w2_d1->phi());
    		} else {
    			std::cout<<"NOT DILEP EVT AT GEN LEVEL"<<std::endl;
            	plotter.getOrMake2DPre(prefix+"_notGenDilep_failLep","dilepPt",";pt1;pt2",100,0,1000,100,0,1000)->Fill(dilep1->pt(),dilep2->pt(),weight);
            	plotter.getOrMake2DPre(prefix+llS+"notGenDilep_failLep","dilepPt",";pt1;pt2",100,0,1000,100,0,1000)->Fill(dilep1->pt(),dilep2->pt(),weight);
    		}

    		std::cout<<std::endl;
    	} */
    }

    void testThings(TString prefix) {
    	if (diHiggsEvt.type < DiHiggsEvent::DILEP || diHiggsEvt.type == DiHiggsEvent::TAU_HAD) return;
    	if (ht < 400) return;

    	SignalHelper sigInfo(diHiggsEvt,reader_muon,reader_electron);
    	sigInfo.minElRecoPt = 10;
    	sigInfo.minMuRecoPt = 10;
    	sigInfo.maxMuRecoEta = 2.4;

    	TString chS;
    	if (diHiggsEvt.type == DiHiggsEvent::DILEP) {
    		chS = "2l";
    		sigInfo.maxElRecoEta = 2.5;
    		sigInfo.setRecoLeptons(0.1);
    		if (!sigInfo.hasMatchedDileps()) return;
    		if (sigInfo.genlep1->pt() < (sigInfo.genlep1->absPdgId()==13 ? 27 : 30) &&
    			sigInfo.genlep2->pt() < (sigInfo.genlep2->absPdgId()==13 ? 27 : 30) ) return;
    		if (sigInfo.recolep1->pt() < (sigInfo.recolep1->isMuon() ? 27 : 30) &&
    			sigInfo.recolep2->pt() < (sigInfo.recolep2->isMuon() ? 27 : 30) ) return;
    	} else {
    		chS = "1l";
    		sigInfo.maxElRecoEta = 1.479;
    		sigInfo.setRecoLeptons(0.1);
    		if (!sigInfo.hasMatchedSingleLep()) return;
    		if (sigInfo.genlep1->pt() < (sigInfo.genlep1->absPdgId()==13 ? 27 : 30)) return;
    		if (sigInfo.recolep1->pt() < (sigInfo.recolep1->isMuon() ? 27 : 30)) return;
    	}

    	sigInfo.setRecoHbb(reader_fatjet->jets,0.2);
    	if (!sigInfo.recoHbb) return;

    	bool passCross = passTrig17(HLT17_Ele15_IsoVVVL_PFHT450)
						|| passTrig17(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165)
						|| passTrig17(HLT17_Mu15_IsoVVVL_PFHT450);

    	bool passJet = passTrig17(HLT17_AK8PFHT850_TrimMass50) || passTrig17(HLT17_PFHT1050) ||
				passTrig17(HLT17_AK8PFJet400_TrimMass30) || passTrig17(HLT17_AK8PFJet500);

    	bool passSMu = passTrig17(HLT17_IsoMu27) || passTrig17(HLT17_Mu50);
    	bool passSEl = passTrig17(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150) || passTrig17(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned)
    			|| passTrig17(HLT17_Ele27_WPTight_Gsf) || passTrig17(HLT17_Ele32_WPTight_Gsf) || passTrig17(HLT17_Ele35_WPTight_Gsf)
				|| passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);
    	bool passMET = passTrig17(HLT17_PFMET120_PFMHT120_IDTight_PFHT60) || passTrig17(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)
    			|| passTrig17(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);

    	bool passSingleLep = passSMu || passSEl;

//    	std::cout<<"size of fatjets = "<<reader_fatjet->jets.size()<<std::endl;
//    	std::cout<<"reco Hbb is set"<<std::endl;

    	plotter.getOrMake1DPre(prefix,chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);

    	if(passCross) plotter.getOrMake1DPre(prefix+"_Cross",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passJet)   plotter.getOrMake1DPre(prefix+"_Jet",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passSingleLep) plotter.getOrMake1DPre(prefix+"_SLep",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passMET) plotter.getOrMake1DPre(prefix+"_Met",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passSingleLep || passCross) plotter.getOrMake1DPre(prefix+"_SLepoCross",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passSingleLep || passJet)   plotter.getOrMake1DPre(prefix+"_SLepoJet",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passSingleLep || passMET) plotter.getOrMake1DPre(prefix+"_SLepoMet",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passCross || passMET) plotter.getOrMake1DPre(prefix+"_CrossoMet",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if(passSingleLep || passCross || passMET || passJet) {
    		plotter.getOrMake1DPre(prefix+"_passAny",chS+"_ht",";ht",1000,0,4000)->Fill(ht,weight);
    	}
    }

    void testInSemileptonicTTBar(TString prefix) {

    	// dilepton ttbar handling in this region
    	if (!isRealData() && smDecayEvt.nLepsTT == 2 && smDecayEvt.topDecays.size() == 2) {
        	if (smDecayEvt.topDecays[0].type < TopDecay::MU) return;
        	if (smDecayEvt.topDecays[1].type < TopDecay::MU) return;
    	}
    	// can use this region to test non-lepton triggers by tagging the event with single high
    	// quality lepton (veto existence of other leptons) and testing the HT, cross (dif flavor
    	// than tag lep), and MET triggers. Also want to require AK8 b-jet to increase purity.

        LeptonParameters tagLeptonParam = parameters.leptons;
    	tagLeptonParam.mu_minPT = 27;
        tagLeptonParam.mu_getID = &Muon::passTightID;
        tagLeptonParam.mu_getISO = &Muon::pfIso;
        tagLeptonParam.mu_maxISO = 0.15;

        tagLeptonParam.el_minPT = 30;
        tagLeptonParam.el_getISO = &Electron::pfIso;
        tagLeptonParam.el_maxISO = 0.15;

        parameters.leptons.el_minPT = 5;
        parameters.leptons.mu_minPT = 5;

        auto tagElectrons = LeptonProcessor::getElectrons(tagLeptonParam,*reader_electron);
        auto tagMuons     = LeptonProcessor::getMuons(tagLeptonParam,*reader_muon);

        auto vetoElectrons = LeptonProcessor::getElectrons(parameters.leptons,*reader_electron);
        auto vetoMuons     = LeptonProcessor::getMuons(parameters.leptons,*reader_muon);

        if (tagElectrons.size() + tagMuons.size() != 1) return;
        if (vetoElectrons.size() + vetoMuons.size() > 1) return;

        const Lepton* tagLepton = tagMuons.size() ? (Lepton*)tagMuons.front() : (Lepton*)tagElectrons.front();

        if (vetoElectrons.size() + vetoMuons.size() == 1) {
        	const Lepton* vetoLep = vetoMuons.size() ? (Lepton*)vetoMuons[0] : (Lepton*)vetoElectrons[0];
        	if (vetoLep->isMuon() != tagLepton->isMuon()) return;
        	if (vetoLep->index() != tagLepton->index()) return;
        }

        // Now should have just one isolated lepton in the event
        plotter.getOrMake1DPre(prefix,"oneLep","ht",4000,0,4000)->Fill(ht,weight);

        // select fatjet with b-tagging
        if (reader_fatjet->jets.size() != 1) return;
        const FatJet* fatjet = &reader_fatjet->jets[0];
        if (BTagging::getCSVSJCat(parameters.jets,fatjet->subJets()) < BTagging::CSVSJ_MF) return;
        if (PhysicsUtilities::deltaR2(*fatjet,*tagLepton) < 1.6*1.6) return;
        if (fatjet->sdMom().mass() < 30) return;

        plotter.getOrMake1DPre(prefix,"oneLep_FatBjet","ht",4000,0,4000)->Fill(ht,weight);

        // test the triggers
        bool passSingleLep = false;
        bool passCross = false;
        bool isMuon = false;
        if (tagLepton->isMuon()) {
        	passSingleLep = passTrig(HLT17_IsoMu27) || passTrig(HLT17_Mu50);
//        	passCross = passTrig17(HLT17_Ele15_IsoVVVL_PFHT450) || passTrig17(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
        	passCross = passTrig17(HLT17_Mu15_IsoVVVL_PFHT450);
        	isMuon = true;
        } else {
        	passSingleLep = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
        			|| passTrig17(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150);
//        	passCross = passTrig17(HLT17_Mu15_IsoVVVL_PFHT450);
        	passCross = passTrig17(HLT17_Ele15_IsoVVVL_PFHT450) || passTrig17(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);        }

    	bool passJet = passTrig17(HLT17_AK8PFHT850_TrimMass50) || passTrig17(HLT17_PFHT1050) ||
				passTrig17(HLT17_AK8PFJet400_TrimMass30) || passTrig17(HLT17_AK8PFJet500);
    	bool passMET = passTrig17(HLT17_PFMET120_PFMHT120_IDTight_PFHT60) || passTrig17(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)
    			|| passTrig17(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);

    	bool passNonLepPath = passJet || passMET;
    	bool passNonLepAndCrossPath = passNonLepPath || passCross;

    	auto plt = [&](TString idS) {
    		TString name = (isMuon ? "_mu_" : "_el_") + idS;
    		plotter.getOrMake1DPre(prefix+name,"pt",";pt",2000,0,2000)->Fill(tagLepton->pt(),weight);
    		plotter.getOrMake1DPre(prefix+name,"ht",";ht",4000,0,4000)->Fill(ht,weight);

    		if(ht > 400) {
        		plotter.getOrMake1DPre(prefix+name+"_ht400","pt",";pt",2000,0,2000)->Fill(tagLepton->pt(),weight);
    		}
    	};

    	if(passSingleLep) plt("tt1_SingleLep");
    	if(passCross)     plt("tt1_Cross");
    	if(passSingleLep || passCross) plt("tt1_SingleLepoCross");

    	if(!passSingleLep) return;
    	plt("tt1_passSLep");
    	if(passNonLepPath) plt("tt1_passSLepAndMetHT");
    	if(passNonLepAndCrossPath) plt("tt1_passSLepAndMetHTCross");

    }

    void testInDrellYan(TString prefix) {
    	// can use this region to test non-lepton triggers by tagging the event with single high
    	// quality lepton (veto existence of other leptons) and testing the HT, cross (dif flavor
    	// than tag lep), and MET triggers. Also want to require AK8 b-jet to increase purity.

        LeptonParameters tagLeptonParam = parameters.leptons;
    	tagLeptonParam.mu_minPT = 10;
        tagLeptonParam.mu_getID = &Muon::passTightID;
        tagLeptonParam.mu_getISO = &Muon::pfIso;
        tagLeptonParam.mu_maxISO = 0.15;

        tagLeptonParam.el_minPT = 10;
        tagLeptonParam.el_getID = &Electron::passTightID_noIso;
        tagLeptonParam.el_getISO = &Electron::pfIso;
        tagLeptonParam.el_maxISO = 0.15;

        auto tagElectrons = LeptonProcessor::getElectrons(tagLeptonParam,*reader_electron);
        auto tagMuons     = LeptonProcessor::getMuons(tagLeptonParam,*reader_muon);

        if (tagElectrons.size() + tagMuons.size() != 2) return;
        if (tagElectrons.size() - tagMuons.size() == 0) return;

        const Lepton* lep1=0;
        const Lepton* lep2=0;
        bool hasHighPt = false;
        bool isDimuon = false;

        if (tagMuons.size()) {
        	lep1 = tagMuons[0];
        	lep2 = tagMuons[1];
        	hasHighPt = (lep1->pt() > 27);
        	isDimuon = true;
        } else {
        	lep1 = tagElectrons[0];
        	lep2 = tagElectrons[1];
        	hasHighPt = (lep1->pt() > 30);
        }

    	double dr2 = PhysicsUtilities::deltaR2(*lep1,*lep2);
    	MomentumF dilepMom = lep1->p4() + lep2->p4();

        if(!hasHighPt) return;
        if(dr2 > 1.2*1.2) return;
        if(dilepMom.mass() < 70 || dilepMom.mass() > 110) return;

        // Now should have just two nearby isolated leptons in the event
        plotter.getOrMake1DPre(prefix,"twoLep","ht",1000,0,4000)->Fill(ht,weight);

        // select fatjet away from leptons
        if (reader_fatjet->jets.size() != 1) return;
        const FatJet* fatjet = &reader_fatjet->jets[0];
        if (fatjet->pt() < 30) return;
        if (PhysicsUtilities::deltaR2(*fatjet,dilepMom) < 1.6*1.6) return;
        if (PhysicsUtilities::deltaR2(*fatjet,*lep1) < 0.8*0.8) return;
        if (PhysicsUtilities::deltaR2(*fatjet,*lep2) < 0.8*0.8) return;

        plotter.getOrMake1DPre(prefix,"twoLep_FatJet","ht",4000,0,4000)->Fill(ht,weight);

    	auto plt = [&](TString idS) {
    		TString name = (isDimuon ? "_mu_" : "_el_") + idS;
    		plotter.getOrMake1DPre(prefix+name,"pt1",";pt",2000,0,2000)->Fill(lep1->pt(),weight);
    		plotter.getOrMake1DPre(prefix+name,"pt2",";pt",2000,0,2000)->Fill(lep2->pt(),weight);
    		plotter.getOrMake1DPre(prefix+name,"ht",";ht",4000,0,4000)->Fill(ht,weight);

    		if(ht > 400) {
        		plotter.getOrMake1DPre(prefix+name+"_ht400","pt1",";pt",2000,0,2000)->Fill(lep1->pt(),weight);
        		plotter.getOrMake1DPre(prefix+name+"_ht400","pt2",";pt",2000,0,2000)->Fill(lep2->pt(),weight);
    		}
    	};

    	bool passSLep = false;
    	bool passCross = false;
    	if (isDimuon) {
        	passSLep = passTrig(HLT17_IsoMu27) || passTrig(HLT17_Mu50);
        	passCross = passTrig17(HLT17_Mu15_IsoVVVL_PFHT450);
    	} else {
        	passSLep = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
        			|| passTrig17(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150);
        	passCross = passTrig17(HLT17_Ele15_IsoVVVL_PFHT450) || passTrig17(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
    	}

    	bool passJet = passTrig17(HLT17_AK8PFHT850_TrimMass50) || passTrig17(HLT17_PFHT1050) ||
				passTrig17(HLT17_AK8PFJet400_TrimMass30) || passTrig17(HLT17_AK8PFJet500);
    	bool passMET = passTrig17(HLT17_PFMET120_PFMHT120_IDTight_PFHT60) || passTrig17(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)
    			|| passTrig17(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);

    	bool passNonLepPath = passJet || passMET;
    	bool passNonLepAndCrossPath = passNonLepPath || passCross;

    	if(passSLep) plt("dy_SingleLep");
    	if(passCross)     plt("dy_Cross");
    	if(passSLep || passCross) plt("dy_SingleLepoCross");

    	if(!passSLep) return;
    	plt("dy_passSLep");
    	if(passNonLepPath) plt("dy_passSLepAndMetHT");
    	if(passNonLepAndCrossPath) plt("dy_passSLepAndMetHTCross");
    }

    void testMuonInDileptonTTBar(TString sn, std::vector<const Electron*> tagElectrons, std::vector<const Muon*> probeMuons) {
    	// filter out dilepton ttbar events with taus
//    	if (smDecayEvt.nLepsTT != 2) return;
//    	if (smDecayEvt.topDecays.size() != 2) return;
//    	if (smDecayEvt.topDecays[0].type < TopDecay::MU) return;
//    	if (smDecayEvt.topDecays[1].type < TopDecay::MU) return;

    	if(!tagElectrons.size()) return;
    	if(!probeMuons.size()) return;

    	bool passSingleEl = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
    			|| passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150);
    	if(!passSingleEl) return;

    	bool passMuPath = passTrig(HLT17_IsoMu27) || passTrig(HLT17_Mu50);
    	bool passCross = passTrig(HLT17_Mu15_IsoVVVL_PFHT450);
        bool passMetNoMu = passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight);

        passMuPath = passMuPath || passMetNoMu;
        bool passSingleCross = passMuPath || passCross;

    	auto plt = [&](TString idS) {
    		TString name = "_mu_" + idS;
    		plotter.getOrMake1DPre(sn+name,"ht",";ht",4000,0,4000)->Fill(ht,weight);
    		plotter.getOrMake1DPre(sn+name,"pt",";pt",2000,0,2000)->Fill(probeMuons.front()->pt(),weight);

    		if(ht > 400) {
        		plotter.getOrMake1DPre(sn+name+"_ht400","pt",";pt",2000,0,2000)->Fill(probeMuons.front()->pt(),weight);
    		}
    	};

    	plt("tt2_passSEl");
    	if(passMuPath) plt("tt2_passSElAndSMu");
    	if(passSingleCross) plt("tt2_passSElAndSMuCross");
    }

    void testElectronInDileptonTTBar(TString sn, std::vector<const Muon*> tagMuons, std::vector<const Electron*> probeElectrons) {
    	if(!tagMuons.size()) return;
    	if(!probeElectrons.size()) return;

    	bool passSingleMu = passTrig(HLT17_IsoMu27);
    	if(!passSingleMu) return;

    	bool passElPath = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
    			|| passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150)
				|| passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);

    	bool passCross = passTrig(HLT17_Ele15_IsoVVVL_PFHT450) || passTrig(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);

    	auto plt = [&](TString idS) {
    		TString name = "_el_" + idS;
    		plotter.getOrMake1DPre(sn+name,"ht",";ht",4000,0,4000)->Fill(ht,weight);
    		plotter.getOrMake1DPre(sn+name,"pt",";pt",2000,0,2000)->Fill(probeElectrons.front()->pt(),weight);

    		if(ht > 400) {
        		plotter.getOrMake1DPre(sn+name+"_ht400","pt",";pt",2000,0,2000)->Fill(probeElectrons.front()->pt(),weight);
    		}
    	};

    	plt("tt2_passSMu");
    	if(passElPath) plt("tt2_passSMuAndSEl");
    	if(passElPath || passCross) plt("tt2_passSMuAndSElCross");
    }

    bool runEvent() override {

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(isRealData()) smpName = FillerConstants::DatasetNames[reader_event->dataset.val()];
/*
        plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(0.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_goodVertices) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(1.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_globalTightHalo2016Filter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(2.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_HBHENoiseFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(3.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_HBHENoiseIsoFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(4.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(5.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_eeBadScFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(6.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_BadPFMuonFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(7.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_muonBadTrackFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(8.0,weight);
        if(FillerConstants::doesPass(reader_event->metFilters.val(),FillerConstants::Flag_BadChargedCandidateFilter) ) plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(9.0,weight);
        if(reader_event->goodVtx.val() != 0)  plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(10.0,weight);
        if(passEventFilters)  plotter.getOrMake1DPre(smpName,"checkFilters",";checkFilters",15,-0.5,14.5)->Fill(11.0,weight);
*/
        if(!passEventFilters) return false;

        triggerAccepts = reader_event->triggerAccepts.val();
//        studyDileptonTriggerEffs(smpName);
        if(isSignal()) testThings(smpName);

        TString sn = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	sn += TString::Format("%d",smDecayEvt.nLepsTT);
        }

        testInSemileptonicTTBar(sn);
        testInDrellYan(sn);

//        std::cout<<"event num muons (electrons) = "<<reader_muon->muons.size()<<" ("<<reader_electron->electrons.size()<<")"<<std::endl;
//
//        if (reader_muon->muons.size()) {
//        	for (const auto& mu : reader_muon->muons) {
//        		printf("%i: pt = %.2f, eta = %.2f, mIso = %.2f, rIso = %.2f, d0 = %.2f, dz = %.2f, sip = %.2f, ",mu.index(),mu.pt(),mu.eta(),mu.miniIso(),mu.pfIso(),mu.d0(),mu.dz(),mu.sip3D());
//        		std::cout<<"passTight = "<<mu.passTightID()<<std::endl;
//        	}
//        }
//        if (reader_electron->electrons.size()) {
//        	for (const auto& el : reader_electron->electrons) {
//        		printf("%i: pt = %.2f, eta = %.2f, mIso = %.2f, rIso = %.2f, d0 = %.2f, dz = %.2f, sip = %.2f, ",el.index(),el.pt(),el.eta(),el.miniIso(),el.pfIso(),el.d0(),el.dz(),el.sip3D());
//        		std::cout<<"passTight = "<<el.passTightID_noIso()<<std::endl;
//        	}
//        }
        LeptonParameters tagLeptonParam = parameters.leptons;
    	tagLeptonParam.mu_minPT = 27;
        tagLeptonParam.mu_getID = &Muon::passTightID;
        tagLeptonParam.mu_getISO = &Muon::pfIso;
        tagLeptonParam.mu_maxISO = 0.15;

        tagLeptonParam.el_minPT = 30;
        tagLeptonParam.el_getISO = &Electron::pfIso;
        tagLeptonParam.el_maxISO = 0.15;

        parameters.leptons.el_minPT = 5;
        parameters.leptons.mu_minPT = 5;

        auto tagElectrons = LeptonProcessor::getElectrons(tagLeptonParam,*reader_electron);
        auto tagMuons     = LeptonProcessor::getMuons(tagLeptonParam,*reader_muon);
//        std::cout<<"num tags for muons (electrons) = "<<tagMuons.size()<<" ("<<tagElectrons.size()<<")"<<std::endl;

        auto probeElectrons = LeptonProcessor::getElectrons(parameters.leptons,*reader_electron);
        auto probeMuons     = LeptonProcessor::getMuons(parameters.leptons,*reader_muon);
//        std::cout<<"num probes for muons (electrons) = "<<probeMuons.size()<<" ("<<probeElectrons.size()<<")"<<std::endl;

        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleElectron) testMuonInDileptonTTBar(sn,tagElectrons,probeMuons);
        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleMuon)     testElectronInDileptonTTBar(sn,tagMuons,probeElectrons);


        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleElectron){
//        	testEachTriggerIndividually(smpName,false,tagMuons,probeMuons,tagElectrons,probeElectrons);
//            doMuonLeg(smpName,tagElectrons,probeMuons);
//            doHTLegWithElDenom(smpName,tagElectrons,probeMuons);
//            doGrandLeptonWElDenom(smpName,tagElectrons,probeMuons);
        }
        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleMuon){
//        	testEachTriggerIndividually(smpName,true,tagMuons,probeMuons,tagElectrons,probeElectrons);
//            doElectronLeg(smpName,tagMuons,probeElectrons);
//            doHTLegWithMuonDenom(smpName,tagMuons,probeElectrons);
//            doGrandLeptonWMuDenom(smpName,tagMuons,probeElectrons);
        }
        if(!isRealData()) {
        	if (isSignal() && diHiggsEvt.type < DiHiggsEvent::MU) return false;
//        	doMCLepton(smpName,probeMuons,probeElectrons);
//        	doSystLepton(smpName);
        }

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    size64 triggerAccepts=0;

    static const int nHTBins = 16;
    const double htBins[nHTBins+1] = {0,50,100,150,200,250,300,350,400,450,500,550,600,800,1200,1600,2000};
    static const int nLepBins = 10;
    const double lepBins[nLepBins+1] = {5,10,15,20,25,30,35,50,75,100,500};

};

#endif

void getTriggerTurnons(std::string fileName, int treeInt,  int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
