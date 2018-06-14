#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_GENTOOLS_INTERFACE_DIHIGGSEVENT_H_

#include <vector>
#include <tuple>

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

private:
    int tau_search(const GenParticle* dau);
    bool isWpair(const GenParticle* f1, const GenParticle* f2,const GenParticleCollection& g_parts);
    int isPair(const GenParticle* f1, const GenParticle* f2,const GenParticleCollection& g_parts);
    int classify_W_pair(const GenParticle* p1, const GenParticle* p2);
    std::tuple<const GenParticle*, const GenParticle*> assign_gp(const GenParticle* p1, const GenParticle* p2);
    struct WDecay {
        const GenParticle* id = 0;
        const GenParticle* d1 = 0;
        const GenParticle* d2 = 0;
        int decaytype = 0; // 0 is BAD, others defined above in enum decayidentifier
    };

    WDecay assign_W(const GenParticle* w, const GenParticleCollection& g_parts);
    DECAYTYPE classify_decaytype(const std::vector<int>& items);
    std::vector<int> search_4_daughters(const GenParticle* gp, const GenParticleCollection& g_parts);
};


}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
