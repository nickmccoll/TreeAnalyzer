#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {
//_____________________________________________________________________________
bool LepSelHelpers::isGoodMuon(const Muon* lep, const float minPT, const float maxETA,
        const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO,
        muFunBool getID, muFunFloat getISO) {
    if(lep->pt() < minPT) return false;
    if(lep->absEta() >= maxETA) return false;
    if(maxDZ > 0 && std::fabs(lep->dz()) >= maxDZ) return false;
    if(maxD0 > 0 && std::fabs(lep->d0()) >= maxD0) return false;
    if(maxSIP3D > 0 && lep->sip3D() >= maxSIP3D) return false;
    if(maxISO >0 && (lep->*getISO)()  >= maxISO) return false;
    if((lep->*getID)()  == false) return false;
    return true;
}
//_____________________________________________________________________________
bool LepSelHelpers::isGoodElectron(const Electron* lep, const float minPT, const float maxETA,
        const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO,
        elFunBool getID, elFunFloat getISO) {
    if(lep->pt() < minPT) return false;
    if(std::fabs(lep->scEta()) >= maxETA) return false;
    if(maxDZ > 0 && std::fabs(lep->dz()) >= maxDZ) return false;
    if(maxD0 > 0 && std::fabs(lep->d0()) >= maxD0) return false;
    if(maxSIP3D > 0 && lep->sip3D() >= maxSIP3D) return false;
    if(maxISO >0 && (lep->*getISO)()  >= maxISO) return false;
    if((lep->*getID)()  == false) return false;
    return true;
}
//_____________________________________________________________________________
bool LepSelHelpers::isGoodMuon(const Muon* lep, const LeptonParameters& params) {
    return isGoodMuon(lep,params.mu_minPT , params.mu_maxETA, params.mu_maxDZ ,
            params.mu_maxD0 ,params.mu_maxSip3D , params.mu_maxISO, params.mu_getID ,
            params.mu_getISO
            );
}
//_____________________________________________________________________________
bool LepSelHelpers::isGoodElectron(const Electron* lep, const LeptonParameters& params) {
    return isGoodElectron(lep,params.el_minPT , params.el_maxETA, params.el_maxDZ ,
            params.el_maxD0 ,params.el_maxSip3D,  params.el_maxISO, params.el_getID ,
            params.el_getISO
            );
}
//_____________________________________________________________________________
bool LeptonProcessor::isGoodMuon(const EventReader& reader_event, const Muon * lep) const {
   return   LepSelHelpers::isGoodMuon(lep, *params    );
}
//_____________________________________________________________________________
bool LeptonProcessor::isGoodElectron(const Electron * lep) const {
    return  LepSelHelpers::isGoodElectron(lep,*params);
}
//_____________________________________________________________________________
bool LeptonProcessor::isGoodLepton(const EventReader& reader_event, const Lepton * lep) const {
    return lep->isMuon() ? isGoodMuon(reader_event,(const Muon*)lep)
            : isGoodElectron((const Electron*)lep);
}
//_____________________________________________________________________________
std::vector<const Muon    *> LeptonProcessor::getMuons(const EventReader& reader_event,
        const MuonReader& reader_muon) const {
    std::vector<const Muon    *>  leps;
    for(const auto& lep :reader_muon.muons){
        if(LepSelHelpers::isGoodMuon(&lep,  *params    ))
            leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Muon>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Electron*> LeptonProcessor::getElectrons(
        const ElectronReader& reader_electron) const {
    std::vector<const Electron    *>  leps;
    for(const auto& lep :reader_electron.electrons){
        if(LepSelHelpers::isGoodElectron(&lep,*params))
            leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Electron>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Lepton*> LeptonProcessor::getLeptons(const EventReader& reader_event,
        const MuonReader& reader_muon, const ElectronReader& reader_electron) const {
    std::vector<const Lepton    *>  leps;

    auto electrons = getElectrons(reader_electron);
    auto muons     = getMuons    (reader_event,reader_muon);
    leps.reserve(electrons.size() + muons.size());
    for(const auto* l : electrons) leps.push_back(l);
    for(const auto* l : muons) leps.push_back(l);
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());
    return leps;
}
}
