#include "Processors/Variables/interface/DileptonSelection.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {
namespace DileptonProcessor {
//_____________________________________________________________________________
bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO, muFunBool getID1, muFunBool getID2, muFunFloat getISO) {
    if(lep->pt() < minPT) return false;
    if(lep->absEta() >= maxETA) return false;
    if(maxDZ > 0 && std::fabs(lep->dz()) >= maxDZ) return false;
    if(maxD0 > 0 && std::fabs(lep->d0()) >= maxD0) return false;
    if(maxSIP3D > 0 && lep->sip3D() >= maxSIP3D) return false;
    if(maxISO >0 && (lep->*getISO)()  >= maxISO) return false;
    if((lep->*getID1)()  == false && (lep->*getID2)() == false) return false;
    return true;
}
//_____________________________________________________________________________
bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO, elFunBool getID1, elFunBool getID2, elFunFloat getISO) {
    if(lep->pt() < minPT) return false;
    if(std::fabs(lep->scEta()) >= maxETA) return false;
    if(maxDZ > 0 && std::fabs(lep->dz()) >= maxDZ) return false;
    if(maxD0 > 0 && std::fabs(lep->d0()) >= maxD0) return false;
    if(maxSIP3D > 0 && lep->sip3D() >= maxSIP3D) return false;
    if(maxISO > 0 && (lep->*getISO)()  >= maxISO) return false;
    if((lep->*getID1)()  == false && (lep->*getID2)() == false) return false;
    return true;
}
//_____________________________________________________________________________
bool isGoodMuon(const Muon* lep, const DileptonParameters& params) {
    return isGoodMuon(lep, params.mu_minPT, params.mu_maxETA, params.mu_maxDZ, params.mu_maxD0,
            params.mu_maxSip3D, params.mu_maxISO, params.mu_getID1, params.mu_getID2, params.mu_getISO
            );
}
//_____________________________________________________________________________
bool isGoodElectron(const Electron* lep, const DileptonParameters& params) {
    return isGoodElectron(lep, params.el_minPT, params.el_maxETA, params.el_maxDZ, params.el_maxD0,
            params.el_maxSip3D, params.el_maxISO, params.el_getID1, params.el_getID2, params.el_getISO
            );
}
//_____________________________________________________________________________
bool isGoodLepton(const Lepton * lep, const DileptonParameters& params) {
    return lep->isMuon() ? isGoodMuon((const Muon*)lep,params) : isGoodElectron((const Electron*)lep,params);
}
//_____________________________________________________________________________
std::vector<const Muon*> getMuons(const DileptonParameters& params, const MuonReader& reader_muon) {
    std::vector<const Muon*> leps;
    for(const auto& lep : reader_muon.muons){
        if(isGoodMuon(&lep,params)) leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Muon>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Electron*> getElectrons(const DileptonParameters& params,
		const ElectronReader& reader_electron) {
    std::vector<const Electron*> leps;
    for(const auto& lep : reader_electron.electrons){
        if(isGoodElectron(&lep,params)) leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Electron>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Lepton*> getLeptons(const DileptonParameters& params,
		const MuonReader& reader_muon, const ElectronReader& reader_electron) {
    std::vector<const Lepton*>  leps;

    auto electrons = getElectrons(params,reader_electron);
    auto muons = getMuons(params, reader_muon);
    leps.reserve(electrons.size() + muons.size());
    for(const auto* l : electrons) leps.push_back(l);
    for(const auto* l : muons) leps.push_back(l);
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());
    return leps;
}

}
}
