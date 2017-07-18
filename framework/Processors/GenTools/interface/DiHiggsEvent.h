#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_

#include <vector>

namespace TAna {
class GenParticle;
typedef std::vector<GenParticle> GenParticleCollection;

class DiHiggsEvent {
public:
    enum DECAYTYPE {BAD,bbZZ,HAD,DILEP, TAU_HAD,TAU_MU,TAU_E, MU,E};
    void setDecayInfo(const GenParticleCollection& genparts);
    void reset();

    const GenParticle * hbb =0;
    const GenParticle * b1=0;
    const GenParticle * b2=0;

    const GenParticle * hww=0;
    const GenParticle * w1  =0;
    const GenParticle * w1_d1=0; //sorted that in lep is l
    const GenParticle * w1_d2=0; //sorted that in lep is n

    const GenParticle * w2  =0;
    const GenParticle * w2_d1=0;
    const GenParticle * w2_d2=0;
    DECAYTYPE type = BAD;


};


}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
