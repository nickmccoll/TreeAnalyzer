#include "Processors/Variables/interface/HiggsSolver.h"
#include "../interface/DefaultSearchRegionAnalyzer.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "Configuration/interface/FillerConstants.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/Hww2lSolver.h"
#include "Processors/Variables/interface/JetKinematics.h"

#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "TPRegexp.h"




namespace TAna {
//--------------------------------------------------------------------------------------------------
DefaultSearchRegionAnalyzer::DefaultSearchRegionAnalyzer(std::string fileName,
        std::string treeName, int treeInt, size randomSeed)
: BaseTreeAnalyzer(fileName,treeName,treeInt, randomSeed){
    TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
    auto match = r1.MatchS(fileName);
    const Int_t nrSubStr = match->GetLast()+1;
    if(nrSubStr>1){
        signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
    }
    fjProc      .reset(new FatJetProcessor ());
    trigSFProc  .reset(new TriggerScaleFactors (dataDirectory));
    puSFProc    .reset(new PUScaleFactors (dataDirectory));
    leptonSFProc.reset(new POGLeptonScaleFactors(dataDirectory));
    dileptonSFProc.reset(new POGLeptonScaleFactors(dataDirectory));
    ak4btagSFProc.reset(new JetBTagScaleFactors (dataDirectory));
    sjbtagSFProc.reset(new SubJetBTagScaleFactors (dataDirectory));
    hbbFJSFProc .reset(new HbbFatJetScaleFactors (dataDirectory));
    topPTProc   .reset(new TopPTWeighting (dataDirectory));

    JERProc   .reset(new JERCorrector (dataDirectory,randGen));;

    JESUncProc . reset(new JESUncShifter());
    METUncProc . reset(new METUncShifter());
    HEMIssueProc.reset(new HEM1516TestCorrector());
    hSolverChi . reset(new HSolverChi());
    hSolverLi . reset(new HSolverLi(dataDirectory));

    turnOnCorr(CORR_XSEC);
    turnOnCorr(CORR_TRIG);
    turnOnCorr(CORR_PU  );
    turnOnCorr(CORR_LEP );
    turnOnCorr(CORR_SJBTAG);
    turnOnCorr(CORR_AK4BTAG);
    turnOnCorr(CORR_SDMASS);
//    turnOnCorr(CORR_TOPPT);
    turnOnCorr(CORR_JER);
}
DefaultSearchRegionAnalyzer::~DefaultSearchRegionAnalyzer(){}

//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::resetCorr() {corrections = 0;}
bool DefaultSearchRegionAnalyzer::isCorrOn(Corrections corr) const
{return FillerConstants::doesPass(corrections,corr);}
void DefaultSearchRegionAnalyzer::turnOnCorr(Corrections corr)
{FillerConstants::addPass(corrections,corr);}
void DefaultSearchRegionAnalyzer::turnOffCorr(Corrections corr)
{FillerConstants::removePass(corrections,corr);}

//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::loadVariables()  {
    reader_event       =loadReader<EventReader>   ("event",isRealData());
    reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData(),true,true);
    reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,false,true);
    reader_jet         =loadReader<JetReader>     ("ak4Jet",isRealData());
    reader_electron    =loadReader<ElectronReader>("electron");
    reader_muon        =loadReader<MuonReader>    ("muon",isRealData());

    if(!isRealData()){
        reader_genpart =loadReader<GenParticleReader>   ("genParticle");
    }

    checkConfig();
}
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::checkConfig()  {
    auto mkErr = [&](const std::string& reader, const std::string corr) {
        throw std::invalid_argument(
                std::string("DefaultSearchRegionAnalyzer::checkConfig() -> Must load ")
        + reader + std::string(" if you want ") + corr);};
    if(!reader_event) mkErr("event","anything");

    if(isRealData()) return;

    if(isCorrOn(CORR_XSEC) && !reader_event) mkErr("event","CORR_XSEC");
    if(isCorrOn(CORR_TRIG) && !reader_jet) mkErr("ak4Jet","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_electron) mkErr("electron","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_muon) mkErr("muon","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_genpart) mkErr("genParticle","CORR_TRIG");
    if(isCorrOn(CORR_PU)   && !reader_event) mkErr("event","CORR_PU");
    if(isCorrOn(CORR_LEP) && !reader_electron) mkErr("electron","CORR_LEP");
    if(isCorrOn(CORR_LEP) && !reader_muon) mkErr("muon","CORR_LEP");
    if(isCorrOn(CORR_LEP) && !reader_genpart) mkErr("genParticle","CORR_LEP");
    if(isCorrOn(CORR_AK4BTAG) && !reader_jet) mkErr("ak4PuppiNoLepJet","CORR_AK4BTAG");
    if(isCorrOn(CORR_SJBTAG) && !reader_fatjet) mkErr("ak8PuppiNoLepJet","CORR_SJBTAG");
    if(isCorrOn(CORR_TOPPT) && !reader_event) mkErr("event","CORR_TOPPT");
    if(isCorrOn(CORR_TOPPT) && !reader_genpart) mkErr("genParticle","CORR_TOPPT");
    if(isCorrOn(CORR_JER) && !reader_fatjet) mkErr("fatjet","CORR_JER");
    if(isCorrOn(CORR_JER) && !reader_fatjet_noLep) mkErr("fatjet_noLep","CORR_JER");
    if(isCorrOn(CORR_JER) && !reader_jet) mkErr("jet","CORR_JER");

    if(isCorrOn(CORR_JES) && !reader_fatjet) mkErr("fatjet","CORR_JES");
    if(isCorrOn(CORR_JES) && !reader_fatjet_noLep) mkErr("fatjet_noLep","CORR_JES");
    if(isCorrOn(CORR_JES) && !reader_jet) mkErr("jet","CORR_JES");
}
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::setupParameters(){
    announceParameterChange();
    setupParametersFromEra();
    setupProcessorParameters();
}
void DefaultSearchRegionAnalyzer::announceParameterChange() {
    lastEra = FillerConstants::DataEra(*reader_event->dataEra);
    std::cout << " ++  Setting era: "<<FillerConstants::DataEraNames[lastEra]<<std::endl;
}
void DefaultSearchRegionAnalyzer::setupParametersFromEra(){
    switch(FillerConstants::DataEra(*reader_event->dataEra)){
        case FillerConstants::ERA_2018:
            parameters = ReaderConstants::set2018Parameters();
            break;
        case FillerConstants::ERA_2017:
            parameters = ReaderConstants::set2017Parameters();
            break;
        case FillerConstants::ERA_2016:
            parameters = ReaderConstants::set2016Parameters();
            break;
        default:
            throw std::invalid_argument(
                    "DefaultSearchRegionAnalyzer -> The era needs to be set to use this class");
    }
}
void DefaultSearchRegionAnalyzer::setupProcessorParameters(){
    if(isCorrOn(CORR_SJBTAG))  sjbtagSFProc->setParameters(parameters.jets);
    if(isCorrOn(CORR_AK4BTAG)) ak4btagSFProc->setParameters(parameters.jets);
    if(isCorrOn(CORR_TRIG))    trigSFProc->setParameters(parameters.event);
    if(isCorrOn(CORR_PU))      puSFProc->setParameters(parameters.event);
    if(isCorrOn(CORR_JER))     JERProc->setParameters(parameters.jets);
    if(isCorrOn(CORR_LEP))     leptonSFProc->setParameters(parameters.leptons);
    if(isCorrOn(CORR_LEP))     dileptonSFProc->setParameters(parameters.dileptons);

    hSolverLi->setParamters(parameters.hww);
}

bool DefaultSearchRegionAnalyzer::hasPromptLeptons() const {
    if(isRealData()) return false;
    return smDecayEvt.promptElectrons.size()
            || smDecayEvt.promptMuons.size()
            || *reader_event->process == FillerConstants::ZJETS;
}

//--------------------------------------------------------------------------------------------------
bool DefaultSearchRegionAnalyzer::runEvent() {
    fillEventLabels();

    if(*reader_event->dataEra != lastEra){
        setupParameters();
    }

    correctJetsAndMET();

    fillGenInfo();
    fillJetInfo();
    fillLeptonInfo();
    fillFilterResults();
    fillLeptonChannel();

    fillFatJetCandidates();
    fillHbbInfo();
    fillWjjInfo();
    fillHWWInfo();
    fillHHInfo();

    fillEventWeight();

    return true;
}

void DefaultSearchRegionAnalyzer::fillEventLabels() {
    if(isRealData()){
        mcProc = FillerConstants::NOPROCESS;
        smpName = "data";
        return;
    }
    mcProc = FillerConstants::MCProcess(*(reader_event->process));
    signal_mass = *reader_event->sampParam;
    if (mcProc == FillerConstants::SIGNAL){
        smpName = TString::Format("%s_m%i",
                FillerConstants::SignalTypeNames[*reader_event->signalType].data(),signal_mass);
    }
    else smpName = FillerConstants::MCProcessNames[mcProc];
}


void DefaultSearchRegionAnalyzer::correctJetsAndMET() {
    if(!isRealData()) return;

    if(isCorrOn(CORR_JES) ){
        Met dummyMET =reader_event->met;
        JESUncProc ->processJets(*reader_jet,reader_event->met);
        JESUncProc ->processFatJets(reader_fatjet_noLep->jets);
        JESUncProc ->processFatJets(reader_fatjet->jets);
    }
    if(isCorrOn(CORR_JER) ){
        JERProc   ->processJets(
                *reader_jet,reader_event->met,reader_jet->genJets,reader_event->rho.val());
        JERProc ->processFatJets(
                reader_fatjet_noLep->jets,std::vector<GenJet>(),reader_event->rho.val());
        JERProc ->processFatJets(
                reader_fatjet->jets,reader_fatjet->genJets,reader_event->rho.val());
    }
    if(isCorrOn(CORR_MET) ){
        METUncProc->process(reader_event->met,*reader_event);
    }
    if(isCorrOn(CORR_HEM1516) && *reader_event->dataEra == FillerConstants::ERA_2018) {
        HEMIssueProc->processJets(*reader_jet,reader_event->met);
        HEMIssueProc->processFatJets(reader_fatjet_noLep->jets);
        HEMIssueProc->processFatJets(reader_fatjet->jets);
    }
}

void DefaultSearchRegionAnalyzer::fillGenInfo() {
    if(!reader_genpart) return;
    if(mcProc == FillerConstants::SIGNAL)
        diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
    smDecayEvt.setDecayInfo(reader_genpart->genParticles,*reader_event->sampParam);
}

void DefaultSearchRegionAnalyzer::fillJetInfo() {
    if(!reader_jet) return;

    jets_ht = PhysicsUtilities::selObjsMom(reader_jet->jets,30);
    ht = JetKinematics::ht(jets_ht);

    jets = PhysicsUtilities::selObjsMom(reader_jet->jets,
            parameters.jets.minJetPT,parameters.jets.maxJetETA,
            [&](const Jet* j){return (j->*parameters.jets.passJetID)();} );
    nMedBTags = PhysicsUtilities::selObjsMomD(jets,
            parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
            [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();
}

void DefaultSearchRegionAnalyzer::fillLeptonInfo() {
    if(reader_muon == 0 || reader_electron == 0) return;

    resetLeptonInfo();
    selectLeptons();
    if (selectedDileptons.size() >= 2) fillTwoLepInfo();
}
void DefaultSearchRegionAnalyzer::resetLeptonInfo() {
    selectedLeptons.clear();
    selectedLepton = 0;
    selectedDileptons.clear();
    dilep1 = 0;
    dilep2 = 0;
    llMass = 9999;
    llDR = 9999;
    llMetDphi = 9999;

}
void DefaultSearchRegionAnalyzer::selectLeptons() {
    selectedLeptons = LeptonProcessor::getLeptons(parameters.leptons,
            *reader_muon,*reader_electron);
    selectedLepton = selectedLeptons.size() ? selectedLeptons.front() : 0;

    selectedDileptons = DileptonProcessor::getLeptons(
            parameters.dileptons,*reader_muon,*reader_electron);
}
void DefaultSearchRegionAnalyzer::fillTwoLepInfo() {
    dilep1 = selectedDileptons[0];
    dilep2 = selectedDileptons[1];
    llMass = (dilep1->p4()+dilep2->p4()).mass();
    llDR = PhysicsUtilities::deltaR(*dilep1,*dilep2);
    llMetDphi = PhysicsUtilities::absDeltaPhi(reader_event->met,(dilep1->p4()+dilep2->p4()));
}

void DefaultSearchRegionAnalyzer::fillFilterResults() {
    passEventFilters= EventSelection::passEventFilters(parameters.event,*reader_event);
    passTriggerPreselection= EventSelection::passTriggerPreselection(
            parameters.event,*reader_event,ht,selectedLeptons);
    passTriggerPreselection2l = EventSelection::passTriggerPreselection(
            parameters.event,*reader_event,ht,selectedDileptons);
}

void DefaultSearchRegionAnalyzer::fillLeptonChannel() {
    lepChan = NOCHANNEL;
    if (selectedDileptons.size() == 2 && dilep1->q() != dilep2->q() && passTriggerPreselection2l)
        lepChan = DILEP;
    else if (selectedDileptons.size() < 2 && selectedLepton && passTriggerPreselection)
        lepChan = SINGLELEP;
}

void DefaultSearchRegionAnalyzer::fillFatJetCandidates() {
    resetFatJetCandidates();
    if(lepChan == SINGLELEP) fillOneLepFatJetCandidates();
    else if(lepChan == DILEP) fillTwoLepFatJetCandidates();
}
void DefaultSearchRegionAnalyzer::resetFatJetCandidates() {
    hbbCand = 0;
    wjjCand = 0;
}
void DefaultSearchRegionAnalyzer::fillOneLepFatJetCandidates() {
    if(reader_fatjet && reader_fatjet_noLep) {
        fjProc->loadFatJets(parameters.fatJets,*reader_fatjet,*reader_fatjet_noLep,selectedLepton);
        hbbCand     = fjProc->getHBBCand();
        wjjCand     = fjProc->getWjjCand();
    }
}
void DefaultSearchRegionAnalyzer::fillTwoLepFatJetCandidates() {
    if(reader_fatjet) {
        fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,dilep1,dilep2);
        hbbCand     = fjProc->getDilepHbbCand();
    }
}

void DefaultSearchRegionAnalyzer::fillHbbInfo() {
    if(!hbbCand) {
        resetHbbInfo();
        return;
    }

    hbbCSVCat   = BTagging::getCSVSJCat(parameters.jets,hbbCand->subJets());
    hbbTag      = BTagging::getFatJetTagValue(parameters.jets,*hbbCand);
    hbbMass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbbCand) : hbbCand->sdMom().mass();

    jets_HbbV = PhysicsUtilities::selObjsD(jets,
            [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbbCand ) >= 1.2*1.2;});

    nMedBTags_HbbV = PhysicsUtilities::selObjsMomD(jets_HbbV,
            parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
            [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();
}
void DefaultSearchRegionAnalyzer::resetHbbInfo() {
    hbbCSVCat = BTagging::CSVSJ_INCL;
    hbbTag = 0;
    hbbMass = 0;
    jets_HbbV = jets;
    nMedBTags_HbbV = nMedBTags;
}

void DefaultSearchRegionAnalyzer::fillWjjInfo() {
    if(!wjjCand) {
        wjjCSVCat = BTagging::CSVSJ_INCL;
        return;
    }
    wjjCSVCat   = BTagging::getCSVSJCat(parameters.jets,wjjCand->subJets());
}

void DefaultSearchRegionAnalyzer::fillHWWInfo() {
    resetHWWInfo();
    if(lepChan == SINGLELEP) fillOneLepHWWInfo();
    else if(lepChan == DILEP) fillTwoLepHWWInfo();
}
void DefaultSearchRegionAnalyzer::resetHWWInfo() {
    hwwLi       =  1000;
    neutrino    =  MomentumF();
    wlnu        =  MomentumF();
    wqq         =  MomentumF();
    wwDM        =  0;
    hWW         =  MomentumF();
}
void DefaultSearchRegionAnalyzer::fillOneLepHWWInfo() {
    if(!wjjCand) return;

    HSolverLiInfo  hwwInfoLi;
    const double qqSDMass = wjjCand->sdMom().mass();
    hwwLi   = hSolverLi->minimize(selectedLepton->p4(),reader_event->met.p4(),
            wjjCand->p4(), qqSDMass, hwwInfoLi);

    neutrino = hwwInfoLi.neutrino;
    wlnu     = hwwInfoLi.wlnu;
    wqq      = hwwInfoLi.wqqjet;
    wwDM     = PhysicsUtilities::deltaR( wlnu,wqq) * hWW.pt()/2.0;
    hWW      = hwwInfoLi.hWW;
}
void DefaultSearchRegionAnalyzer::fillTwoLepHWWInfo() {
    hWW = Hww2lSolver::getSimpleHiggsMom(dilep1->p4()+dilep2->p4(),reader_event->met,
            parameters.hww.dilepInvMassGuess);
}

void DefaultSearchRegionAnalyzer::fillHHInfo() {
    resetHHInfo();

    bool hasSolution = (hbbCand && lepChan == DILEP) || (hbbCand && wjjCand);
    if(!hasSolution) return;

    hh = hbbCand->p4() + hWW.p4();
    if(lepChan == SINGLELEP) {
        auto oldN =  HSolverBasic::getInvisible(reader_event->met,
                (selectedLepton->p4() + wjjCand->p4()) );
        hh_basic =  hbbCand->p4() + oldN.p4() + selectedLepton->p4() + wjjCand->p4();
    }
}
void DefaultSearchRegionAnalyzer::resetHHInfo() {
    hh         =  MomentumF();
    hh_basic   =  MomentumF();
}

void DefaultSearchRegionAnalyzer::fillEventWeight() {
    weight = 1;
    if(isRealData()) return;

    if(isCorrOn(CORR_XSEC)) weight *= getXSecWeight();
    if(isCorrOn(CORR_TRIG)) weight *= getTriggerWeight();
    if(isCorrOn(CORR_PU)) weight *= getPUWeight();
    if(isCorrOn(CORR_LEP)) weight *= getLeptonWeight();
    if(isCorrOn(CORR_SJBTAG)) weight *= getSJBTagWeights();
    if(isCorrOn(CORR_AK4BTAG)) weight *= getAK4BTagWeights();
    if(isCorrOn(CORR_TOPPT)) weight *= getTopPTWeight();
}
float DefaultSearchRegionAnalyzer::getXSecWeight() {
    return EventWeights::getNormalizedEventWeight(
            *reader_event,xsec(),nSampEvt(),parameters.event,smDecayEvt.genMtt,smDecayEvt.nLepsTT);
}
float DefaultSearchRegionAnalyzer::getTriggerWeight() {
    double triggerWeight = 1;
    if(FillerConstants::DataEra(*reader_event->dataEra) == FillerConstants::ERA_2017) {
        triggerWeight *= EventSelection::get2017CrossTrigWeight(*reader_event);
    }

    if(!hasPromptLeptons()) return triggerWeight;

    if (lepChan == DILEP) {
        triggerWeight *= trigSFProc->getDileptonTriggerSF(ht,dilep2->pt(),dilep1->isMuon(),
                dilep2->isMuon());
    } else if(selectedLepton) {
        triggerWeight *= trigSFProc->getSingleLeptonTriggerSF(
                ht, selectedLepton->isMuon());
    }
    return triggerWeight;
}
float DefaultSearchRegionAnalyzer::getPUWeight() {
    return puSFProc->getCorrection(reader_event->nTruePUInts.val(),CorrHelp::NOMINAL);
}
float DefaultSearchRegionAnalyzer::getLeptonWeight() {
    LeptonScaleFactors * sfProc = 0;
    const std::vector<const Lepton*> * leptonCollection;

    if (lepChan == DILEP) {
        sfProc = dileptonSFProc.get();
        leptonCollection = &selectedDileptons;
    } else {
        sfProc = leptonSFProc.get();
        leptonCollection = &selectedLeptons;
    }

    if(*reader_event->process == FillerConstants::ZJETS) {
        sfProc->loadWithoutPromptCheck(*leptonCollection);
    } else {
        sfProc->load(smDecayEvt,*leptonCollection);
    }
    return sfProc->getSF();
}
float DefaultSearchRegionAnalyzer::getSJBTagWeights() {
    return sjbtagSFProc->getSF(parameters.jets,{hbbCand});
}
float DefaultSearchRegionAnalyzer::getAK4BTagWeights() {
    return ak4btagSFProc->getSF(jets_HbbV);
}
float DefaultSearchRegionAnalyzer::getTopPTWeight() {
    return topPTProc->getCorrection(mcProc,smDecayEvt);
}
}
//--------------------------------------------------------------------------------------------------

