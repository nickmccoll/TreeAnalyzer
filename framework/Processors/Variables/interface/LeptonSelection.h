
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

namespace LepSelHelpers {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
}



namespace LepSelHelpers {
    bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D,  const float maxISO, muFunBool getID, muFunFloat getISO);
    bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO, elFunBool getID, elFunFloat getISO);
    bool isGoodMuon(const Muon* lep, const LeptonParameters& params);
    bool isGoodElectron(const Electron* lep, const LeptonParameters& params);
}




class LeptonProcessor {
public:
    void setParameters(const LeptonParameters& inParams) {params = inParams;}

    bool isGoodMuon(const EventReader& reader_event, const Muon * lep)const;
    bool isGoodElectron(const Electron * lep) const;
    bool isGoodLepton(const EventReader& reader_event, const Lepton * lep) const;
    std::vector<const Muon    *> getMuons    (const EventReader& reader_event, const MuonReader& reader_muon) const;
    std::vector<const Electron*> getElectrons(const ElectronReader& reader_electron ) const ;
    std::vector<const Lepton  *> getLeptons  (const EventReader& reader_event, const MuonReader& reader_muon, const ElectronReader& reader_electron) const;
    LeptonParameters params;
};
}


#endif //ANALYSISTOOLS_KINEMATICVARIABLES_JETKINEMATICS_H

