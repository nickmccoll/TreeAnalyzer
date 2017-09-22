#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_SMDECAYEVENT_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_SMDECAYEVENT_H_

#include <vector>

namespace TAna {
class GenParticle;
typedef std::vector<GenParticle> GenParticleCollection;

struct BosonDecay {
    enum DECAYTYPE {BAD,Z_HAD,Z_DILEP, Z_INV, W_HAD, W_TAU_HAD,W_TAU_MU,W_TAU_E,W_MU,W_E};
    const GenParticle * boson = 0;
    const GenParticle * dau1  = 0; //sorted that in 1l lep is l
    const GenParticle * dau2  = 0; //sorted that in 1l lep is n
    DECAYTYPE type = BAD;
    void reset() {type = BAD; boson = 0; dau1 = 0; dau2 = 0;}
};
struct TopDecay {
    enum DECAYTYPE {BAD,HAD, TAU_HAD,TAU_MU,TAU_E,MU,E};
    const GenParticle * top   = 0;
    const GenParticle * b     = 0;
    BosonDecay W_decay  ;
    DECAYTYPE type = BAD;
    void reset() {type = BAD, top = 0; b = 0; W_decay.reset();}
};
class SMDecayEvent {
public:
    void setDecayInfo(const GenParticleCollection& genparts);
    void reset();
    std::vector<BosonDecay> bosonDecays;
    std::vector<TopDecay>   topDecays;
    std::vector<const GenParticle *>   promptElectrons   ;
    std::vector<const GenParticle *>   nonPromptElectrons;
    std::vector<const GenParticle *>   promptMuons       ;
    std::vector<const GenParticle *>   nonPromptMuons    ;
};


}



#endif
