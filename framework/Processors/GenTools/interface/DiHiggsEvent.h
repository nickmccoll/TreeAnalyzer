#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_

#include <vector>
#include "DataFormats/interface/GenParticle.h"

namespace TAna {
class GenParticle;
typedef std::vector<GenParticle> GenParticleCollection;

struct WDecay {
    const GenParticle* id = 0;
    const GenParticle* d1 = 0;
    const GenParticle* d2 = 0;
    int decaytype = 0; // 0 is BAD, others defined above in enum decayidentifier
};

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
    int type = BAD;

private:
    int tau_search(CandidateRef<GenParticle> dau);
    bool isWpair(CandidateRef<GenParticle> f1, CandidateRef<GenParticle> f2);
    int isPair(CandidateRef<GenParticle> f1, CandidateRef<GenParticle> f2);
    int classify_W_pair(CandidateRef<GenParticle> p1, CandidateRef<GenParticle> p2);
    std::tuple<const GenParticle*, const GenParticle*> assign_gp(const GenParticle* p1, const GenParticle* p2);
    WDecay assign_W(const GenParticle* w);
    int classify_decaytype(std::vector<int> items);
    std::vector<int> search_4_daughters(CandidateRef<GenParticle> gp);
};


}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */