
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
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

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        tagLeptonProc .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*tagLeptonProc);
        tagLeptonProc->lepSelParams.mu_getID =  &Muon::passTightID;
        tagLeptonProc->lepSelParams.mu_minPT =  26;
        tagLeptonProc->lepSelParams.mu_getISO =&Muon::dbRelISO;
        tagLeptonProc->lepSelParams.mu_maxISO =  0.15;

        tagLeptonProc->lepSelParams.el_getISO =  &Electron::eaRelISO;
        tagLeptonProc->lepSelParams.el_maxISO =  0.15;
        tagLeptonProc->lepSelParams.el_minPT  =  30;
        tagLeptonProc->lepSelParams_dataABCDEF = tagLeptonProc->lepSelParams;

        leptonProc->lepSelParams.el_minPT = 5;
        leptonProc->lepSelParams.mu_minPT = 5;
        leptonProc->lepSelParams_dataABCDEF.el_minPT = 5;
        leptonProc->lepSelParams_dataABCDEF.mu_minPT = 5;

        turnOffCorr(CORR_TRIG);
    }


    bool passTrig(Triggers trig) {return doesPass(triggerAccepts,trig);}

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

    void doElectronLeg(const TString& prefix){
        if( !passTrig(HLT_IsoMu24) && !passTrig(HLT_IsoTkMu24)  ) return;
        if(!tagMuons.size()) return;
        float maxLepPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_passSMu";
        makeHTPlots(preName,"el_pt",ht_wlep,maxLepPT);
        bool passECross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);

        if(passECross)
            makeHTPlots(preName + "_passElHT_","el_pt",ht_wlep,maxLepPT);
    }
    void doMuonLeg(const TString& prefix){
        if( !passTrig(HLT_Ele27_WPTight_Gsf) ) return;
        if(!tagElectrons.size()) return;

        bool passOfflineE = false;
        float maxLepPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_passSE";
        makeHTPlots(preName,"mu_pt",ht_wlep,maxLepPT);
        bool passMCross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);

        if(passMCross)
            makeHTPlots(preName + "_passMuHT_","mu_pt",ht_wlep,maxLepPT);
    }

    void doHTLegWithMuonDenom(const TString& prefix){
        if( !passTrig(HLT_IsoMu24) && !passTrig(HLT_IsoTkMu24)  ) return;
        if(!tagMuons.size()) return;
        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_passSMu";
        makeLepPlots(preName,"ht","elpt",maxOtherPT,ht_wlep);
        makeLepPlots(preName,"ht","mupt",maxSamePT,ht_wlep);
        bool passMCross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        bool passECross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        if(passECross)
            makeLepPlots(preName + "_passElHT_","ht","elpt",maxOtherPT,ht_wlep);
        if(passMCross)
            makeLepPlots(preName + "_passMuHT_","ht","mupt",maxSamePT,ht_wlep);

        if(passECross|| passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30))
            makeLepPlots(preName + "_passElHToHad_","ht","elpt",maxOtherPT,ht_wlep);
        if(passMCross|| passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30))
            makeLepPlots(preName + "_passMuHToHad_","ht","mupt",maxSamePT,ht_wlep);
    }
    void doHTLegWithElDenom(const TString& prefix){
        if( !passTrig(HLT_Ele27_WPTight_Gsf) ) return;
        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_passSE";
        makeLepPlots(preName,"ht","mupt",maxOtherPT,ht_wlep);
        makeLepPlots(preName,"ht","elpt",maxSamePT,ht_wlep);
        bool passMCross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        bool passECross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);

        if(passECross)
            makeLepPlots(preName + "_passElHT_","ht","elpt",maxSamePT,ht_wlep);
        if(passMCross)
            makeLepPlots(preName + "_passMuHT_","ht","mupt",maxOtherPT,ht_wlep);

        if(passECross|| passTrig(HLT_PFHT800) ||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30))
            makeLepPlots(preName + "_passElHToHad_","ht","elpt",maxSamePT,ht_wlep);

        if(passMCross || passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30) )
            makeLepPlots(preName + "_passMuHToHad_","ht","mupt",maxOtherPT,ht_wlep);
    }

    void doGrandLeptonWElDenom(const TString& prefix){
        if( !passTrig(HLT_Ele27_WPTight_Gsf) ) return;
        if(!tagElectrons.size()) return;

        float maxSamePT = tagElectrons.front()->pt();
        float maxOtherPT = probeMuons.size() ? probeMuons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSE";

        bool passSMu = passTrig(HLT_IsoTkMu24) || passTrig(HLT_IsoTkMu24);
        bool passMCross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passBu = passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30)|| passTrig(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
        bool passSMuoHtMuoBu = passSMuoHtMu || passBu ||  passTrig(HLT_TkMu50)|| passTrig(HLT_Mu50);



        //HT side
        makeLepPlots(preName,"ht","mupt",maxOtherPT,ht_wlep);
        if(passSMu ) makeLepPlots(preName +"_passSMu_","ht","mupt",maxOtherPT,ht_wlep);
        if(passSMuoHtMu) makeLepPlots(preName +"_passSMuoHtMu_","ht","mupt",maxOtherPT,ht_wlep);
        if(passSMuoHtMuoBu) makeLepPlots(preName +"_passSMuoHtMuoBu_","ht","mupt",maxOtherPT,ht_wlep);

        //Mu side
        makeHTPlots(preName,"mu_pt",ht_wlep,maxOtherPT);
        if(passSMu )        makeHTPlots(preName+"_passSMu_"         ,"mu_pt",ht_wlep,maxOtherPT);
        if(passSMuoHtMu)    makeHTPlots(preName+"_passSMuoHtMu_"    ,"mu_pt",ht_wlep,maxOtherPT);
        if(passSMuoHtMuoBu) makeHTPlots(preName+"_passSMuoHtMuoBu_","mu_pt",ht_wlep,maxOtherPT);

        //2D
        if(maxOtherPT > 0){
            plotter.getOrMake2DPre(preName,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSMu )         plotter.getOrMake2DPre(preName+"_passSMu_"         ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSMuoHtMu)     plotter.getOrMake2DPre(preName+"_passSMuoHtMu_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSMuoHtMuoBu)  plotter.getOrMake2DPre(preName+"_passSMuoHtMuoBu_","mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
        }
    }

    void doGrandLeptonWMuDenom(const TString& prefix){
        if( !passTrig(HLT_IsoMu24) && !passTrig(HLT_IsoTkMu24)  ) return;
        if(!tagMuons.size()) return;
        float maxSamePT  = tagMuons.front()->pt();
        float maxOtherPT = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_GL_passSMu";

        bool passSEl = passTrig(HLT_Ele27_WPTight_Gsf);
        bool passECross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        bool passSEloHtEl = passSEl  || passECross;
        bool passBu = passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30)|| passTrig(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
        bool passSEloHtEloBu = passSEloHtEl || passBu ||  passTrig(HLT_Ele45_WPLoose_Gsf)|| passTrig(HLT_Ele115_CaloIdVT_GsfTrkIdT);



        //HT side
        makeLepPlots(preName,"ht","elpt",maxOtherPT,ht_wlep);
        if(passSEl ) makeLepPlots(preName +"_passSEl_","ht","elpt",maxOtherPT,ht_wlep);
        if(passSEloHtEl) makeLepPlots(preName +"_passSEloHtEl_","ht","elpt",maxOtherPT,ht_wlep);
        if(passSEloHtEloBu) makeLepPlots(preName +"_passSEloHtEloBu_","ht","elpt",maxOtherPT,ht_wlep);

        //Mu side
        makeHTPlots(preName,"el_pt",ht_wlep,maxOtherPT);
        if(passSEl )        makeHTPlots(preName+"_passSEl_"        ,"el_pt",ht_wlep,maxOtherPT);
        if(passSEloHtEl)    makeHTPlots(preName+"_passSEloHtEl_"   ,"el_pt",ht_wlep,maxOtherPT);
        if(passSEloHtEloBu) makeHTPlots(preName+"_passSEloHtEloBu_","el_pt",ht_wlep,maxOtherPT);

        //2D
        if(maxOtherPT > 0){
            plotter.getOrMake2DPre(preName,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSEl )        plotter.getOrMake2DPre(preName+"_passSEl_"        ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSEloHtEl)    plotter.getOrMake2DPre(preName+"_passSEloHtEl_"   ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
            if(passSEloHtEloBu) plotter.getOrMake2DPre(preName+"_passSEloHtEloBu_","el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxOtherPT,ht_wlep,weight);
        }

    }


    void doMCLepton(const TString& prefix){
        float maxMu = probeMuons.size() ? probeMuons.front()->pt() : 0;
        float maxEl = probeElectrons.size() ? probeElectrons.front()->pt() : 0;

        TString preName = prefix + "_MC";

        bool passBu = passTrig(HLT_PFHT800)||passTrig(HLT_PFHT900) || passTrig(HLT_AK8PFJet450)|| passTrig(HLT_AK8PFJet360_TrimMass30) || passTrig(HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);

        bool passSMu = passTrig(HLT_IsoTkMu24) || passTrig(HLT_IsoMu24);
        bool passMCross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        bool passSMuoHtMu = passSMu || passMCross;
        bool passSMuoHtMuoBu = passSMuoHtMu || passBu ||  passTrig(HLT_TkMu50)|| passTrig(HLT_Mu50);
        bool passSEl = passTrig(HLT_Ele27_WPTight_Gsf);
        bool passECross = (isRealData() && reader_event->run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        bool passSEloHtEl = passSEl  || passECross;
        bool passSEloHtEloBu = passSEloHtEl || passBu ||  passTrig(HLT_Ele45_WPLoose_Gsf)|| passTrig(HLT_Ele115_CaloIdVT_GsfTrkIdT);


        //HT side
        makeLepPlots(preName,"ht","mupt",maxMu,ht_wlep);
        if(passSMu ) makeLepPlots(preName +"_passSMu_","ht","mupt",maxMu,ht_wlep);
        if(passSMuoHtMu) makeLepPlots(preName +"_passSMuoHtMu_","ht","mupt",maxMu,ht_wlep);
        if(passSMuoHtMuoBu) makeLepPlots(preName +"_passSMuoHtMuoBu_","ht","mupt",maxMu,ht_wlep);

        makeLepPlots(preName,"ht","elpt",maxEl,ht_wlep);
        if(passSEl ) makeLepPlots(preName +"_passSEl_","ht","elpt",maxEl,ht_wlep);
        if(passSEloHtEl) makeLepPlots(preName +"_passSEloHtEl_","ht","elpt",maxEl,ht_wlep);
        if(passSEloHtEloBu) makeLepPlots(preName +"_passSEloHtEloBu_","ht","elpt",maxEl,ht_wlep);

        //Mu side
        makeHTPlots(preName,"mu_pt",ht_wlep,maxMu);
        if(passSMu )        makeHTPlots(preName+"_passSMu_"         ,"mu_pt",ht_wlep,maxMu);
        if(passSMuoHtMu)    makeHTPlots(preName+"_passSMuoHtMu_"    ,"mu_pt",ht_wlep,maxMu);
        if(passSMuoHtMuoBu) makeHTPlots(preName+"_passSMuoHtMuoBu_","mu_pt",ht_wlep,maxMu);

        makeHTPlots(preName,"el_pt",ht_wlep,maxEl);
        if(passSEl )        makeHTPlots(preName+"_passSEl_"        ,"el_pt",ht_wlep,maxEl);
        if(passSEloHtEl)    makeHTPlots(preName+"_passSEloHtEl_"   ,"el_pt",ht_wlep,maxEl);
        if(passSEloHtEloBu) makeHTPlots(preName+"_passSEloHtEloBu_","el_pt",ht_wlep,maxEl);

        if(maxEl > 0){
            plotter.getOrMake2DPre(preName,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_wlep,weight);
            if(passSMu )         plotter.getOrMake2DPre(preName+"_passSMu_"         ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_wlep,weight);
            if(passSMuoHtMu)     plotter.getOrMake2DPre(preName+"_passSMuoHtMu_"    ,"mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_wlep,weight);
            if(passSMuoHtMuoBu)  plotter.getOrMake2DPre(preName+"_passSMuoHtMuoBu_","mu_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxEl,ht_wlep,weight);
        }

        if(maxMu > 0){
            plotter.getOrMake2DPre(preName,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_wlep,weight);
            if(passSEl )        plotter.getOrMake2DPre(preName+"_passSEl_"        ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_wlep,weight);
            if(passSEloHtEl)    plotter.getOrMake2DPre(preName+"_passSEloHtEl_"   ,"el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_wlep,weight);
            if(passSEloHtEloBu) plotter.getOrMake2DPre(preName+"_passSEloHtEloBu_","el_pt_v_ht","lepton p_{T} [GeV]; #it{H}_{T} [GeV]; arbitrary units",nLepBins,lepBins,nHTBins,htBins )->Fill(maxMu,ht_wlep,weight);
        }
    }


    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;

        if(isRealData()) smpName = FillerConstants::DatasetNames[reader_event->dataset];

        triggerAccepts = reader_event->triggerAccepts;

        tagElectrons = tagLeptonProc->getElectrons(*reader_electron);
        tagMuons     = tagLeptonProc->getMuons(*reader_event,*reader_muon);

        probeElectrons = leptonProc->getElectrons(*reader_electron);
        probeMuons     = leptonProc->getMuons(*reader_event,*reader_muon);

        if(!isRealData() || reader_event->dataset == SINGLEE){
            doMuonLeg(smpName);
            doHTLegWithElDenom(smpName);
            doGrandLeptonWElDenom(smpName);
        }
        if(!isRealData() || reader_event->dataset == SINGLEMU){
            doElectronLeg(smpName);
            doHTLegWithMuonDenom(smpName);
            doGrandLeptonWMuDenom(smpName);
        }
        if(!isRealData()) doMCLepton(smpName);

        for(unsigned int iT = 0; iT <= triggerStrings.size(); ++iT ){
            if(FillerConstants::doesPass(reader_event->triggerAccepts, size64(1) << iT))
                plotter.getOrMake2DPre(smpName,"trigger_prescale",";trigger;isPrescaled",64,-0.5,63.5,2,-0.5,1.5)->Fill(iT,0.0);
            if(FillerConstants::doesPass(reader_event->triggerPrescales, size64(1) << iT))
                plotter.getOrMake2DPre(smpName,"trigger_prescale",";trigger;isPrescaled",64,-0.5,63.5,2,-0.5,1.5)->Fill(iT,1.0);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    size64 triggerAccepts =0;
    std::unique_ptr<LeptonProcessor> tagLeptonProc ;
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

void getTriggerTurnons(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getTriggerTurnons(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
