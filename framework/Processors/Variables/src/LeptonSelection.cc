#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {
namespace LeptonProcessor{
//_____________________________________________________________________________
bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA,
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
bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA,
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
bool isGoodMuon(const LeptonParameters& params, const Muon* lep) {
    return isGoodMuon(lep,params.mu_minPT , params.mu_maxETA, params.mu_maxDZ ,
            params.mu_maxD0 ,params.mu_maxSip3D , params.mu_maxISO, params.mu_getID ,
            params.mu_getISO
    );
}
//_____________________________________________________________________________
bool isGoodElectron(const LeptonParameters& params, const Electron* lep) {
    return isGoodElectron(lep,params.el_minPT , params.el_maxETA, params.el_maxDZ ,
            params.el_maxD0 ,params.el_maxSip3D,  params.el_maxISO, params.el_getID ,
            params.el_getISO
    );
}
//_____________________________________________________________________________
bool isGoodLepton(const LeptonParameters& params, const Lepton * lep)  {
    return lep->isMuon() ? isGoodMuon(params,(const Muon*)lep)
            : isGoodElectron(params,(const Electron*)lep);
}
//_____________________________________________________________________________
std::vector<const Muon    *>  getMuons(const LeptonParameters& params,
        const MuonReader& reader_muon)  {
    std::vector<const Muon    *>  leps;
    for(const auto& lep :reader_muon.muons){
        if(isGoodMuon(params,&lep      ))
            leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Muon>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Electron*> getElectrons(const LeptonParameters& params,
        const ElectronReader& reader_electron) {
    std::vector<const Electron    *>  leps;
    for(const auto& lep :reader_electron.electrons){
        if(isGoodElectron(params,&lep))
            leps.push_back(&lep);
    }
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Electron>());
    return leps;
}
//_____________________________________________________________________________
std::vector<const Lepton*> getLeptons(const LeptonParameters& params,
        const MuonReader& reader_muon, const ElectronReader& reader_electron)  {
    std::vector<const Lepton    *>  leps;

    auto electrons = getElectrons(params,reader_electron);
    auto muons     = getMuons    (params,reader_muon);
    leps.reserve(electrons.size() + muons.size());
    for(const auto* l : electrons) leps.push_back(l);
    for(const auto* l : muons) leps.push_back(l);
    std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());
    return leps;
}
}
}

