
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

namespace LepSelHelpers {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
}

struct LepSelParameters {
    float el_minPT  ;
    float el_maxETA ;
    float el_maxDZ  ;
    float el_maxD0  ;
    float el_maxISO ;
    LepSelHelpers::elFunBool el_getID ;
    LepSelHelpers::elFunFloat el_getISO;

    float mu_minPT  ;
    float mu_maxETA ;
    float mu_maxDZ  ;
    float mu_maxD0  ;
    float mu_maxISO ;
    LepSelHelpers::muFunBool mu_getID ;
    LepSelHelpers::muFunFloat mu_getISO;
};

namespace LepSelHelpers {
    bool isGoodMuon(const Muon* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0, const float maxISO, muFunBool getID, muFunFloat getISO);
    bool isGoodElectron(const Electron* lep, const float minPT, const float maxETA, const float maxDZ, const float maxD0, const float maxISO, elFunBool getID, elFunFloat getISO);
    bool isGoodMuon(const Muon* lep, const LepSelParameters& params);
    bool isGoodElectron(const Electron* lep, const LepSelParameters& params);
}




class LeptonProcessor {
public:
    std::vector<const Muon    *> getMuons    (const EventReader& reader_event, const MuonReader& reader_muon);
    std::vector<const Electron*> getElectrons(const ElectronReader& reader_electron);
    std::vector<const Lepton  *> getLeptons  (const EventReader& reader_event, const MuonReader& reader_muon, const ElectronReader& reader_electron);
    LepSelParameters lepSelParams;
    LepSelParameters lepSelParams_dataABCDEF;
};

namespace DefaultLeptonSelections {
void setDefaultLepSelParams(LepSelParameters& par)       ;
void setDefaultLepSelParams_dataAF(LepSelParameters& par);
void setDefaultLeptonProcessor(LeptonProcessor& proc)     ;
}


}


#endif //ANALYSISTOOLS_KINEMATICVARIABLES_JETKINEMATICS_H

