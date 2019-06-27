
#ifndef PROCESSORS_VARIABLES_DILEPTONSELECTION_H
#define PROCESSORS_VARIABLES_DILEPTONSELECTION_H
#include <vector>
#include "TString.h"

#include "Configuration/interface/ReaderConstants.h"

namespace TAna {
class Lepton;
class Electron;
class Muon;
class EventReader;
class MuonReader;
class ElectronReader;

namespace DileptonProcessor {

enum DilepChan {LL_BAD, LL_EE, LL_MUMU, LL_EMU};
DilepChan getDilepChan(const Lepton* lep1, const Lepton* lep2);
TString   getDilepStr (const Lepton* lep1, const Lepton* lep2);

bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ,
        const float maxD0,const float maxSIP3D,  const float maxISO, muFunBool getID,
        muFunFloat getISO);
bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA,
        const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO,
        elFunBool getID, elFunFloat getISO);
bool isGoodMuon(const DileptonParameters& params, const Muon* lep);
bool isGoodElectron(const DileptonParameters& params, const Electron* lep);
bool isGoodLepton(const DileptonParameters& params, const Lepton * lep);
std::vector<const Muon    *> getMuons    (const DileptonParameters& params,
        const MuonReader& reader_muon);
std::vector<const Electron*> getElectrons(const DileptonParameters& params,
        const ElectronReader& reader_electron );
std::vector<const Lepton  *> getLeptons  (const DileptonParameters& params,
        const MuonReader& reader_muon, const ElectronReader& reader_electron);
}

}
#endif //ANALYSISTOOLS_KINEMATICVARIABLES_JETKINEMATICS_H

