
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
namespace TAna {

void DiHiggsEvent::reset() {
    hbb =0;
    b1=0;
    b2=0;

    hww=0;
    w1  =0;
    w1_d1=0; //sorted that in lep is l
    w1_d2=0; //sorted that in lep is n

    w2  =0;
    w2_d1=0;
    w2_d2=0;
    type = BAD;
}

void DiHiggsEvent::setDecayInfo(const GenParticleCollection& genparts) {
    reset();

    for(const auto& p : genparts){
        const auto* gp = &p;
        if(std::abs(gp->pdgId()) != ParticleInfo::p_h0) continue;
        if(!ParticleInfo::isLastInChain(gp)) continue;
        auto dau1 = ParticleInfo::getFinal(gp->daughterRef(0));
        auto dau2 = ParticleInfo::getFinal(gp->daughterRef(1));
        int dPDGID = std::abs(dau1->pdgId());
        if(dPDGID == ParticleInfo::p_b){
            hbb = gp;
            if(dau1->pt() > dau2->pt()){
                b1 = &(*dau1);
                b2 = &(*dau2);
            } else {
                b1 = &(*dau2);
                b2 = &(*dau1);
            }
        } else  {
            hww = gp;
            for(unsigned int iD = 0; iD <gp->numberOfDaughters(); ++iD ){
                auto dau = ParticleInfo::getFinal(gp->daughterRef(iD));
                if(dau->absPdgId() == ParticleInfo::p_Z0 ){
                    type = bbZZ;
                }
            }
            if(type == bbZZ) continue;


            struct WDecay {
                const GenParticle * dau1 = 0; //sorted that in lep is l
                const GenParticle * dau2=0; //sorted that in lep is n
                int type = 3;//0 lep //1 tau //2 had //3 bad
            };

            std::vector<WDecay> wDecays;

            auto getDaughters = [&](const GenParticle* mom) -> WDecay{ //0 lep //1 tau //2 had //3 bad
                WDecay wDecay;
                std::vector<const GenParticle*> quarks;
                std::vector<const GenParticle*> leps;
                for(unsigned int iD = 0; iD <mom->numberOfDaughters(); ++iD ){
                    auto dau = ParticleInfo::getFinal(mom->daughterRef(iD));
                    if(ParticleInfo::isLeptonOrNeutrino(dau->pdgId())){
                        leps.push_back(&*dau);
                    }
                    if(ParticleInfo::isQuark(dau->pdgId())) quarks.push_back(&*dau);
                }

                if(leps.size() == 2){
                    if(ParticleInfo::isANeutrino(leps[0]->pdgId()) && ParticleInfo::isLepton(leps[1]->pdgId()) ){
                        wDecay.dau1 = leps[1];
                        wDecay.dau2 = leps[0];
                        if(TMath::Abs(leps[1]->pdgId()) == ParticleInfo::p_tauminus){wDecay.type = 1;}
                        else wDecay.type = 0;
                    }
                    else if(ParticleInfo::isANeutrino(leps[1]->pdgId()) && ParticleInfo::isLepton(leps[0]->pdgId()) ){
                        wDecay.dau1 = leps[0];
                        wDecay.dau2  = leps[1];
                        if(TMath::Abs(leps[0]->pdgId()) == ParticleInfo::p_tauminus) {wDecay.type = 1;}
                        else wDecay.type = 0;
                    } else {
                        wDecay.dau1 = leps[0];
                        wDecay.dau2 = leps[1];
                        wDecay.type = 3;
                    }
                } else if(quarks.size() == 2){
                    wDecay.dau1 = quarks[1];
                    wDecay.dau2 = quarks[0];
                    wDecay.type = 2;
                } else {
                    wDecay.type = 3;
                }
                return wDecay;
            };


            for(unsigned int iD = 0; iD <gp->numberOfDaughters(); ++iD ){
                auto dau = ParticleInfo::getFinal(gp->daughterRef(iD));
                if(dau->absPdgId() == ParticleInfo::p_Wplus ){
                    wDecays.push_back(getDaughters(&(*dau)  ));
                }
            }

            //now for virtual
            auto virtW = getDaughters(gp);
            if(virtW.type != 3) wDecays.push_back(virtW);

            if(wDecays.size() != 2){ type = BAD;}
            else if(wDecays[0].type == 3 || wDecays[1].type == 3){
                type  = BAD;
                w1_d1 = wDecays[0].dau1;
                w1_d2 = wDecays[0].dau2;
                w2_d1 = wDecays[1].dau1;
                w2_d2 = wDecays[1].dau2;
            } else if(wDecays[0].type == 2 && wDecays[1].type == 2){
                type  = HAD;
                w1_d1 = wDecays[0].dau1;
                w1_d2 = wDecays[0].dau2;
                w2_d1 = wDecays[1].dau1;
                w2_d2 = wDecays[1].dau2;
            } else if(wDecays[0].type != 2 && wDecays[1].type != 2){
                type = DILEP;
                w1_d1 = wDecays[0].dau1;
                w1_d2 = wDecays[0].dau2;
                w2_d1 = wDecays[1].dau1;
                w2_d2 = wDecays[1].dau2;
            } else {
                if(wDecays[0].type == 1 || wDecays[1].type == 1) type = TAU_E;
                else type = MU;
                if(wDecays[0].type == 2){
                    w1_d1 = wDecays[1].dau1;
                    w1_d2 = wDecays[1].dau2;
                    w2_d1 = wDecays[0].dau1;
                    w2_d2 = wDecays[0].dau2;
                } else {
                    w1_d1 = wDecays[0].dau1;
                    w1_d2 = wDecays[0].dau2;
                    w2_d1 = wDecays[1].dau1;
                    w2_d2 = wDecays[1].dau2;
                }
            }
        }

    }

    if(hww== 0 || hbb == 0) type = BAD;
}


}



