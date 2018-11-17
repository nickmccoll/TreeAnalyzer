
#ifndef PROCESSORS_VARIABLES_LEPTONSELECTION_H
#define PROCESSORS_VARIABLES_LEPTONSELECTION_H
#include <vector>

using namespace std;

namespace TAna {
class Lepton;
class Electron;
class Muon;
class EventReader;
class MuonReader;
class ElectronReader;

namespace DilepSelHelpers {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
}

struct DilepSelParameters {
    float el_minPT  ;
    float el_maxETA ;
    float el_maxDZ  ;
    float el_maxD0  ;
    float el_maxSip3D  ;
    float el_maxISO ;
    DilepSelHelpers::elFunBool el_getID1 ;
    DilepSelHelpers::elFunBool el_getID2 ;
    DilepSelHelpers::elFunFloat el_getISO;

    float mu_minPT  ;
    float mu_maxETA ;
    float mu_maxDZ  ;
    float mu_maxD0  ;
    float mu_maxSip3D  ;
    float mu_maxISO ;
    DilepSelHelpers::muFunBool mu_getID1 ;
    DilepSelHelpers::muFunBool mu_getID2 ;
    DilepSelHelpers::muFunFloat mu_getISO;
};

namespace DilepSelHelpers {
    bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO, muFunBool getID1, muFunBool getID2, muFunFloat getISO);
    bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0,const float maxSIP3D, const float maxISO, elFunBool getID1, elFunBool getID2, elFunFloat getISO);
    bool isGoodMuon(const Muon* lep, const DilepSelParameters& params);
    bool isGoodElectron(const Electron* lep, const DilepSelParameters& params);
}

class DileptonProcessor {
public:
    bool isGoodMuon(const EventReader& reader_event, const Muon * lep)const;
    bool isGoodElectron(const Electron * lep) const;
    bool isGoodLepton(const EventReader& reader_event, const Lepton * lep) const;
    std::vector<const Muon    *> getMuons    (const EventReader& reader_event, const MuonReader& reader_muon) const;
    std::vector<const Electron*> getElectrons(const ElectronReader& reader_electron ) const ;
    std::vector<const Lepton  *> getLeptons  (const EventReader& reader_event, const MuonReader& reader_muon, const ElectronReader& reader_electron) const;
    DilepSelParameters lepSelParams;
    DilepSelParameters lepSelParams_dataABCDEF;
};

namespace DefaultDileptonSelections {
void setDefaultDilepSelParams(DilepSelParameters& par)       ;
void setDefaultDilepSelParams_dataAF(DilepSelParameters& par);
void setDefaultDileptonProcessor(DileptonProcessor& proc)     ;
}


}


#endif //ANALYSISTOOLS_KINEMATICVARIABLES_JETKINEMATICS_H

