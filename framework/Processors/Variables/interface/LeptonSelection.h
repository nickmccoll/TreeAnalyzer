
#ifndef PROCESSORS_VARIABLES_LEPTONSELECTION_H
#define PROCESSORS_VARIABLES_LEPTONSELECTION_H
#include <vector>
#include "Configuration/interface/ReaderConstants.h"

namespace TAna {
class Lepton;
class Electron;
class Muon;
class EventReader;
class MuonReader;
class ElectronReader;


namespace LeptonProcessor {
bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ,
        const float maxD0,const float maxSIP3D,  const float maxISO, muFunBool getID,
        muFunFloat getISO);
bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA,
        const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO,
        elFunBool getID, elFunFloat getISO);
bool isGoodMuon(const LeptonParameters& params, const Muon* lep);
bool isGoodElectron(const LeptonParameters& params, const Electron* lep);
bool isGoodLepton(const LeptonParameters& params, const Lepton * lep);
std::vector<const Muon    *> getMuons    (const LeptonParameters& params,
        const MuonReader& reader_muon);
std::vector<const Electron*> getElectrons(const LeptonParameters& params,
        const ElectronReader& reader_electron );
std::vector<const Lepton  *> getLeptons  (const LeptonParameters& params,
        const MuonReader& reader_muon, const ElectronReader& reader_electron);
}
}


#endif //ANALYSISTOOLS_KINEMATICVARIABLES_JETKINEMATICS_H

