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

#include "Processors/Corrections/interface/EventWeights.h"

using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
        tagLeptonParam = parameters.leptons;
    	tagLeptonParam.mu_minPT = 26;
        tagLeptonParam.mu_getID = &Muon::passTightID;
        tagLeptonParam.mu_getISO = &Muon::pfIso;
        tagLeptonParam.mu_maxISO = 0.15;

        tagLeptonParam.el_minPT = 30;
        tagLeptonParam.el_getISO = &Electron::pfIso;
        tagLeptonParam.el_maxISO = 0.15;

        parameters.leptons.el_minPT = 5;
        parameters.leptons.mu_minPT = 5;

        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    bool passTrig(Triggers_2017 trig) {return doesPass(triggerAccepts,trig);} // !! Need to use correct trigger enum for whatever era it is !!

    void makeHTPlots(const TString& prefix, const TString& varname, float sel, float pt){
        plotter.getOrMake1DPre(prefix,TString::Format("ht_incl_%s",varname.Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
        auto mkht = [&](float htCut){
            if(sel >= htCut) plotter.getOrMake1DPre(prefix,TString::Format("ht_%.0f_%s",htCut,varname.Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
        };
        auto mkhtr = [&](float htCutmin, float htcutMax){
            if(sel >= htCutmin && sel < htcutMax ) plotter.getOrMake1DPre(prefix,TString::Format("ht_%.0fto%.0f_%s",htCutmin,htcutMax,varname.Data()),";lepton p_{T} [GeV]; arbitrary units",500,0,500 )->Fill(pt,weight);
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
    }

    void makeLepPlots(const TString& prefix, const TString& varname, TString selname, float sel, float ht){
        plotter.getOrMake1DPre(prefix,TString::Format("%s_incl_%s",selname.Data(),varname.Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(ht,weight);
        auto mklp = [&](float lepcut){
            if(sel >= lepcut) plotter.getOrMake1DPre(prefix,TString::Format("%s_%.0f_%s",selname.Data(), lepcut,varname.Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(ht,weight);
        };
        auto mklpr = [&](float lepcutmin, float lepcutmax){
            if(sel >= lepcutmin && sel < lepcutmax) plotter.getOrMake1DPre(prefix,TString::Format("%s_%.0fto%.0f_%s",selname.Data(), lepcutmin,lepcutmax,varname.Data()),";#it{H}_{T} [GeV]; arbitrary units",2000,0,2000 )->Fill(ht,weight);
        };
        mklp(15);
        mklp(20);
        mklp(25);
        mklp(26);
        mklp(30);
        mklp(35);
        mklp(40);
        mklp(50);
        mklp(75);
        mklp(100);

        mklpr(15,20);
        mklpr(20,25);
        mklpr(25,30);
        mklpr(26,30);
        mklpr(30,35);
        mklpr(35,40);
        mklpr(40,50);
        mklpr(50,75);
        mklpr(75,100);

        mklpr(30,40);
        mklpr(30,50);

        mklpr(40,75);
        mklpr(50,100);
    }

    void testEachTriggerIndividually(const TString& prefix, bool doMuon) {
    	float maxLepPt;
    	if (doMuon) {
    		if (!tagElectrons.size()) return;
    		maxLepPt = probeMuons.size() ? probeMuons.front()->pt() : 0;
    	} else {
    		if (!tagMuons.size()) return;
    		maxLepPt = probeElectrons.size() ? probeElectrons.front()->pt() : 0;
    	}

    	for (Triggers_2017 trg=(Triggers_2017)0; trg != HLT17_NTrig; trg=(Triggers_2017)(trg+1)) {
    		std::cout<<trg<<std::endl;
    		TString preName = prefix + "_passTrig_"+TString::Format("%i",trg);
    		TString varname = doMuon ? "mu_pt" : "el_pt";
    		if (passTrig(trg)) makeHTPlots(preName,varname,ht_chs,maxLepPt);
    	}
    }

    void doElectronLeg(const TString& prefix){
    	// trigger to take tag muon
    	if( !passTrig(HLT17_IsoMu27) ) return;
        if(!tagMuons.size()) return;

        float maxLepPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;
        TString preName = prefix + "_passSMu";
        makeHTPlots(preName,"el_pt",ht_chs,maxLepPT);

        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        if(passECross)
            makeHTPlots(preName + "_passElHT_","el_pt",ht_chs,maxLepPT);
    }

    void doMuonLeg(const TString& prefix){
    	// triggers to take a tag electron
        if( !passTrig(HLT17_Ele35_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG) ) return;
        if(!tagElectrons.size()) return;

        bool passOfflineE = false;
        float maxLepPT = probeMuons.size() ? probeMuons.front()->pt() : 0;
        TString preName = prefix + "_passSE";
        makeHTPlots(preName,"mu_pt",ht_chs,maxLepPT);

        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        if(passMCross)
            makeHTPlots(preName + "_passMuHT_","mu_pt",ht_chs,maxLepPT);
    }

    void doHTLegWithMuonDenom(const TString& prefix){
    	if( !passTrig(HLT17_IsoMu27) ) return;

        if(!tagMuons.size()) return;
        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_passSMu";
        makeLepPlots(preName,"ht","elpt",maxOtherPT,ht_chs);
        makeLepPlots(preName,"ht","mupt",maxSamePT,ht_chs);

        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);

        if(passECross) makeLepPlots(preName + "_passElHT_","ht","elpt",maxOtherPT,ht_chs);
        if(passMCross) makeLepPlots(preName + "_passMuHT_","ht","mupt",maxSamePT,ht_chs);

        if(passECross|| passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30))
            makeLepPlots(preName + "_passElHToHad_","ht","elpt",maxOtherPT,ht_chs);
        if(passMCross|| passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30))
            makeLepPlots(preName + "_passMuHToHad_","ht","mupt",maxSamePT,ht_chs);
    }
    void doHTLegWithElDenom(const TString& prefix){
        if( !passTrig(HLT17_Ele35_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG) ) return;

        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_passSE";
        makeLepPlots(preName,"ht","mupt",maxOtherPT,ht_chs);
        makeLepPlots(preName,"ht","elpt",maxSamePT,ht_chs);

        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);

        if(passECross)
            makeLepPlots(preName + "_passElHT_","ht","elpt",maxSamePT,ht_chs);
        if(passMCross)
            makeLepPlots(preName + "_passMuHT_","ht","mupt",maxOtherPT,ht_chs);

        if(passECross|| passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30))
            makeLepPlots(preName + "_passElHToHad_","ht","elpt",maxSamePT,ht_chs);

        if(passMCross|| passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30))
            makeLepPlots(preName + "_passMuHToHad_","ht","mupt",maxOtherPT,ht_chs);
    }

    void doGrandLeptonWElDenom(const TString& prefix){
        if( !passTrig(HLT17_Ele35_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf) && !passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG) ) return;
        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSE";

        bool passSMu = passTrig(HLT17_IsoMu27);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passBu = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30)
        		|| passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) || passTrig(HLT17_PFMETNoMu140_PFMHTNoMu140_IDTight);
        bool passSMuoHtMuoBu = passSMuoHtMu || passBu || passTrig(HLT17_Mu50);

        bool passMuDenNoCross = passBu || passSMu;

        //HT side
        makeLepPlots(preName,"ht","mupt",maxOtherPT,ht_chs);
        if(passSMu ) makeLepPlots(preName +"_passSMu_","ht","mupt",maxOtherPT,ht_chs);
        if(passSMuoHtMu) makeLepPlots(preName +"_passSMuoHtMu_","ht","mupt",maxOtherPT,ht_chs);
        if(passSMuoHtMuoBu) makeLepPlots(preName +"_passSMuoHtMuoBu_","ht","mupt",maxOtherPT,ht_chs);
        if(passMuDenNoCross) makeLepPlots(preName +"_passMuDenNoCross_","ht","mupt",maxOtherPT,ht_chs);

        //Mu side
        makeHTPlots(preName,"mu_pt",ht_chs,maxOtherPT);
        if(passSMu )        makeHTPlots(preName+"_passSMu_"         ,"mu_pt",ht_chs,maxOtherPT);
        if(passSMuoHtMu)    makeHTPlots(preName+"_passSMuoHtMu_"    ,"mu_pt",ht_chs,maxOtherPT);
        if(passSMuoHtMuoBu) makeHTPlots(preName+"_passSMuoHtMuoBu_","mu_pt",ht_chs,maxOtherPT);

        //2D
        if(maxOtherPT > 0){
            plotter.getOrMake2DPre(preName,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSMu )         plotter.getOrMake2DPre(preName+"_passSMu_"         ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSMuoHtMu)     plotter.getOrMake2DPre(preName+"_passSMuoHtMu_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSMuoHtMuoBu)  plotter.getOrMake2DPre(preName+"_passSMuoHtMuoBu_","mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
        }
    }

    void doGrandLeptonWMuDenom(const TString& prefix){
        if( !passTrig(HLT17_IsoMu27) ) return;
        if(!tagMuons.size()) return;

        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSMu";

        bool passSEl = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        bool passSEloHtEl = passSEl  || passECross;
        bool passBu = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30)
        		|| passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) || passTrig(HLT17_PFMETNoMu140_PFMHTNoMu140_IDTight);
        bool passSEloHtEloBu = passSEloHtEl || passBu ||  passTrig(HLT17_Photon200)|| passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);

        bool passElDenNoCross = passBu || passSEl;

        //HT side
        makeLepPlots(preName,"ht","elpt",maxOtherPT,ht_chs);
        if(passSEl ) makeLepPlots(preName +"_passSEl_","ht","elpt",maxOtherPT,ht_chs);
        if(passSEloHtEl) makeLepPlots(preName +"_passSEloHtEl_","ht","elpt",maxOtherPT,ht_chs);
        if(passSEloHtEloBu) makeLepPlots(preName +"_passSEloHtEloBu_","ht","elpt",maxOtherPT,ht_chs);
        if(passElDenNoCross) makeLepPlots(preName +"_passElDenNoCross_","ht","elpt",maxOtherPT,ht_chs);

        //Mu side
        makeHTPlots(preName,"el_pt",ht_chs,maxOtherPT);
        if(passSEl )        makeHTPlots(preName+"_passSEl_"        ,"el_pt",ht_chs,maxOtherPT);
        if(passSEloHtEl)    makeHTPlots(preName+"_passSEloHtEl_"   ,"el_pt",ht_chs,maxOtherPT);
        if(passSEloHtEloBu) makeHTPlots(preName+"_passSEloHtEloBu_","el_pt",ht_chs,maxOtherPT);

        //2D
        if(maxOtherPT > 0){
            plotter.getOrMake2DPre(preName,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSEl )        plotter.getOrMake2DPre(preName+"_passSEl_"        ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSEloHtEl)    plotter.getOrMake2DPre(preName+"_passSEloHtEl_"   ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
            if(passSEloHtEloBu) plotter.getOrMake2DPre(preName+"_passSEloHtEloBu_","el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_chs,weight);
        }

    }


    void doMCLepton(const TString& prefix){
        float maxMu = probeMuons.size() ? probeMuons.front()->pt() : 0;
        float maxEl = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_MC";

        bool passBu = passTrig(HLT17_PFHT1050)||passTrig(HLT17_AK8PFJet500) || passTrig(HLT17_AK8PFHT850_TrimMass50)|| passTrig(HLT17_AK8PFJet400_TrimMass30)
        		|| passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight) || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) || passTrig(HLT17_PFMETNoMu140_PFMHTNoMu140_IDTight);
        bool passSMu = passTrig(HLT17_IsoMu27);
        bool passMCross = passTrig(FillerConstants::HLT17_Mu15_IsoVVVL_PFHT450);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passSMuoHtMuoBu = passSMuoHtMu || passBu || passTrig(HLT17_Mu50);
        bool passSEl = passTrig(HLT17_Ele35_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf) || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG);
        bool passECross = passTrig(FillerConstants::HLT17_Ele15_IsoVVVL_PFHT450);
        bool passSEloHtEl = passSEl  || passECross;
        bool passSEloHtEloBu = passSEloHtEl || passBu || passTrig(HLT17_Photon200) || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT);

        bool passMuDenNoCross = passBu || passSMu;
        bool passElDenNoCross = passBu || passSEl;


        //HT side
        makeLepPlots(preName,"ht","mupt",maxMu,ht_chs);
        if(passSMu ) makeLepPlots(preName +"_passSMu_","ht","mupt",maxMu,ht_chs);
        if(passSMuoHtMu) makeLepPlots(preName +"_passSMuoHtMu_","ht","mupt",maxMu,ht_chs);
        if(passSMuoHtMuoBu) makeLepPlots(preName +"_passSMuoHtMuoBu_","ht","mupt",maxMu,ht_chs);
        if(passMuDenNoCross) makeLepPlots(preName +"_passMuDenNoCross_","ht","mupt",maxMu,ht_chs);

        makeLepPlots(preName,"ht","elpt",maxEl,ht_chs);
        if(passSEl ) makeLepPlots(preName +"_passSEl_","ht","elpt",maxEl,ht_chs);
        if(passSEloHtEl) makeLepPlots(preName +"_passSEloHtEl_","ht","elpt",maxEl,ht_chs);
        if(passSEloHtEloBu) makeLepPlots(preName +"_passSEloHtEloBu_","ht","elpt",maxEl,ht_chs);
        if(passElDenNoCross) makeLepPlots(preName +"_passElDenNoCross_","ht","elpt",maxEl,ht_chs);


        //Mu side
        makeHTPlots(preName,"mu_pt",ht_chs,maxMu);
        if(passSMu )        makeHTPlots(preName+"_passSMu_"         ,"mu_pt",ht_chs,maxMu);
        if(passSMuoHtMu)    makeHTPlots(preName+"_passSMuoHtMu_"    ,"mu_pt",ht_chs,maxMu);
        if(passSMuoHtMuoBu) makeHTPlots(preName+"_passSMuoHtMuoBu_","mu_pt",ht_chs,maxMu);

        makeHTPlots(preName,"el_pt",ht_chs,maxEl);
        if(passSEl )        makeHTPlots(preName+"_passSEl_"        ,"el_pt",ht_chs,maxEl);
        if(passSEloHtEl)    makeHTPlots(preName+"_passSEloHtEl_"   ,"el_pt",ht_chs,maxEl);
        if(passSEloHtEloBu) makeHTPlots(preName+"_passSEloHtEloBu_","el_pt",ht_chs,maxEl);

        if(maxEl > 0){
            plotter.getOrMake2DPre(preName,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_chs,weight);
            if(passSMu )         plotter.getOrMake2DPre(preName+"_passSMu_"         ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_chs,weight);
            if(passSMuoHtMu)     plotter.getOrMake2DPre(preName+"_passSMuoHtMu_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_chs,weight);
            if(passSMuoHtMuoBu)  plotter.getOrMake2DPre(preName+"_passSMuoHtMuoBu_","mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_chs,weight);
        }

        if(maxMu > 0){
            plotter.getOrMake2DPre(preName,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_chs,weight);
            if(passSEl )        plotter.getOrMake2DPre(preName+"_passSEl_"        ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_chs,weight);
            if(passSEloHtEl)    plotter.getOrMake2DPre(preName+"_passSEloHtEl_"   ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_chs,weight);
            if(passSEloHtEloBu) plotter.getOrMake2DPre(preName+"_passSEloHtEloBu_","el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_chs,weight);
        }
    }


    void loadVariables() override {
        reader_event       =loadReader<EventReader>   ("event",isRealData());
        reader_jet_chs     =loadReader<JetReader>     ("ak4Jet",isRealData());
        reader_electron    =loadReader<ElectronReader>("electron");
        reader_muon        =loadReader<MuonReader>    ("muon");

        if(!isRealData()){
            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
        }

        checkConfig();
    }

    bool runEvent() override {

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(isRealData()) smpName = FillerConstants::DatasetNames[reader_event->dataset.val()];

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

        if(!passEventFilters) return false;

        triggerAccepts = reader_event->triggerAccepts.val();

        std::cout<<"event num muons (electrons) = "<<reader_muon->muons.size()<<" ("<<reader_electron->electrons.size()<<")"<<std::endl;

        if (reader_muon->muons.size()) {
        	for (const auto& mu : reader_muon->muons) {
        		printf("%i: pt = %.2f, eta = %.2f, mIso = %.2f, rIso = %.2f, ",mu.index(),mu.pt(),mu.eta(),mu.miniIso(),mu.pfIso());
        		std::cout<<"passTight = "<<mu.passTightID()<<std::endl;
        	}
        }
        if (reader_electron->electrons.size()) {
        	for (const auto& el : reader_electron->electrons) {
        		printf("%i: pt = %.2f, eta = %.2f, mIso = %.2f, rIso = %.2f, ",el.index(),el.pt(),el.eta(),el.miniIso(),el.pfIso());
        		std::cout<<"passTight = "<<el.passTightID_noIso()<<std::endl;
        	}
        }

        tagElectrons = LeptonProcessor::getElectrons(tagLeptonParam,*reader_electron);
        tagMuons     = LeptonProcessor::getMuons(tagLeptonParam,*reader_muon);
        std::cout<<"num tags for muons (electrons) = "<<tagMuons.size()<<" ("<<tagElectrons.size()<<")"<<std::endl;

        probeElectrons = LeptonProcessor::getElectrons(parameters.leptons,*reader_electron);
        probeMuons     = LeptonProcessor::getMuons(parameters.leptons,*reader_muon);
        std::cout<<"num probes for muons (electrons) = "<<probeMuons.size()<<" ("<<probeElectrons.size()<<")"<<std::endl;
        printf("\n");

        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleElectron){
        	testEachTriggerIndividually(smpName,false);
            doMuonLeg(smpName);
            doHTLegWithElDenom(smpName);
            doGrandLeptonWElDenom(smpName);
        }
        if(!isRealData() || reader_event->dataset.val() == FillerConstants::PD_SingleMuon){
        	testEachTriggerIndividually(smpName,true);
            doElectronLeg(smpName);
            doHTLegWithMuonDenom(smpName);
            doGrandLeptonWMuDenom(smpName);
        }
        if(!isRealData()) doMCLepton(smpName);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    size64 triggerAccepts =0;

    LeptonParameters tagLeptonParam;
//    std::unique_ptr<LeptonProcessor> tagLeptonProc ;
    std::vector<const Electron    *> tagElectrons;
    std::vector<const Muon        *> tagMuons;
    std::vector<const Electron    *> probeElectrons;
    std::vector<const Muon        *> probeMuons;

    static const int nHTBins = 16;
    const double htBins[nHTBins+1] = {0,50,100,150,200,250,300,350,400,450,500,550,600,800,1200,1600,2000};
    static const int nLepBins = 10;
    const double lepBins[nLepBins+1] = {5,10,15,20,25,30,35,50,75,100,500};

};

#endif

void getTriggerTurnons(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getTriggerTurnons(std::string fileName, int treeInt,  int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
