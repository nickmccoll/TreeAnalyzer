
#include "../interface/SMDecayEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"

namespace TAna {

void SMDecayEvent::reset() {
    bosonDecays.clear();
    topDecays.clear();
    promptElectrons   .clear();
    nonPromptElectrons.clear();
    promptMuons       .clear();
    nonPromptMuons    .clear();

}

void SMDecayEvent::setDecayInfo(const GenParticleCollection& genparts, int ttSampParam) {
    reset();

    auto classifyLepW = [](const GenParticle* lep) -> BosonDecay::DECAYTYPE {
        switch(lep->absPdgId()){
        case ParticleInfo::p_eminus : return BosonDecay::W_E;
        case ParticleInfo::p_muminus : return BosonDecay::W_MU;
        case ParticleInfo::p_tauminus:
            for(unsigned int id = 0; id < lep->numberOfDaughters(); ++id){
                const int dPDG = lep->daughter(id)->absPdgId();
                if(dPDG == ParticleInfo::p_eminus) return BosonDecay::W_TAU_E;
                if(dPDG == ParticleInfo::p_muminus) return BosonDecay::W_TAU_MU;
            }
            return BosonDecay::W_TAU_HAD;
        default: return BosonDecay::BAD;
        }
    };

    auto getBosonDecay= [&](const GenParticleRef W) -> BosonDecay {
        BosonDecay bosonDecay;
        bosonDecay.boson = &*W;
        std::vector<const GenParticle*> quarks;
        std::vector<const GenParticle*> leps;
        for(size id = 0; id < W->numberOfDaughters(); ++id){
            auto d = ParticleInfo::getFinal(W->daughterRef(id));
            if(ParticleInfo::isLeptonOrNeutrino(d->pdgId())){
                leps.push_back(&*d);
            }
            else if(ParticleInfo::isQuark(d->pdgId())) quarks.push_back(&*d);
        }
        std::sort(quarks.begin(),quarks.end(),PhysicsUtilities::greaterPTDeref<GenParticle>());
        std::sort(leps.begin(),leps.end(),PhysicsUtilities::greaterPTDeref<GenParticle>());

        if(W->absPdgId() == ParticleInfo::p_Wplus){
            if(leps.size() == 2){
                if(ParticleInfo::isANeutrino(leps[0]->pdgId()) && ParticleInfo::isLepton(leps[1]->pdgId()) ){
                    bosonDecay.dau1 = leps[1];
                    bosonDecay.dau2 = leps[0];
                    bosonDecay.type = classifyLepW(bosonDecay.dau1);
                }
                else if(ParticleInfo::isANeutrino(leps[1]->pdgId()) && ParticleInfo::isLepton(leps[0]->pdgId()) ){
                    bosonDecay.dau1 = leps[0];
                    bosonDecay.dau2  = leps[1];
                    bosonDecay.type = classifyLepW(bosonDecay.dau1);
                } else {
                    bosonDecay.dau1 = leps[0];
                    bosonDecay.dau2  = leps[1];
                    bosonDecay.type = BosonDecay::BAD;
                }
            } else if(quarks.size() == 2){
                bosonDecay.dau1 = quarks[0];
                bosonDecay.dau2 = quarks[1];
                bosonDecay.type = BosonDecay::W_HAD;
            }
        } else {
            if(leps.size() == 2){
                if(ParticleInfo::isANeutrino(leps[0]->pdgId()) && ParticleInfo::isANeutrino(leps[1]->pdgId()) ){
                    bosonDecay.dau1 = leps[0];
                    bosonDecay.dau2 = leps[1];
                    bosonDecay.type = BosonDecay::Z_INV;
                }
                else if(ParticleInfo::isLepton(leps[0]->pdgId()) && ParticleInfo::isLepton(leps[1]->pdgId()) ){
                    bosonDecay.dau1 = leps[0];
                    bosonDecay.dau2  = leps[1];
                    bosonDecay.type = BosonDecay::Z_DILEP;
                } else {
                    bosonDecay.dau1 = leps[0];
                    bosonDecay.dau2  = leps[1];
                    bosonDecay.type = BosonDecay::BAD;
                }
            } else if(quarks.size() == 2){
                bosonDecay.dau1 = quarks[0];
                bosonDecay.dau2 = quarks[1];
                bosonDecay.type = BosonDecay::Z_HAD;
            }
        }

        return bosonDecay;
    };



    auto getTopDecay= [&](const GenParticleRef top) -> TopDecay {
        TopDecay topDecay;
        topDecay.top = &*top;
        for(size idt = 0; idt < top->numberOfDaughters(); ++idt){
            auto dt = ParticleInfo::getFinal(top->daughterRef(idt));
            if(dt->absPdgId() == ParticleInfo::p_b) topDecay.b = &*dt;
            else if(dt->absPdgId() == ParticleInfo::p_Wplus) topDecay.W_decay = getBosonDecay(dt);
        }
        if(topDecay.b && topDecay.W_decay.type >= BosonDecay::W_HAD){
            switch(topDecay.W_decay.type) {
            case BosonDecay::W_HAD    : topDecay.type = TopDecay::HAD    ; break;
            case BosonDecay::W_MU     : topDecay.type = TopDecay::MU     ; break;
            case BosonDecay::W_E      : topDecay.type = TopDecay::E      ; break;
            case BosonDecay::W_TAU_HAD: topDecay.type = TopDecay::TAU_HAD; break;
            case BosonDecay::W_TAU_MU : topDecay.type = TopDecay::TAU_MU ; break;
            case BosonDecay::W_TAU_E  : topDecay.type = TopDecay::TAU_E  ; break;
            default:
                topDecay.type = TopDecay::BAD;
            };
        }
        return topDecay;
    };

    for(unsigned int iP = 0;  iP < genparts.size(); ++iP){
        if(!ParticleInfo::isLastInChain(&genparts[iP])) continue;

        GenParticleRef gp(&genparts[iP],iP);
        const size absID = gp->absPdgId();

        if(absID == ParticleInfo::p_eminus || absID == ParticleInfo::p_muminus){
            auto fo = ParticleInfo::getOriginal(gp); //get original version to check if prompt

            auto isFromPromptEWK = [](const GenParticle* p) -> bool{
                return  ParticleInfo::isDoc(p->status() ) && ParticleInfo::isEWKBoson(p->pdgId()) && ParticleInfo::isLastInChain(p);
            };
            auto isFromPromptTau = [&](const GenParticle* p) -> bool{
                if( std::fabs(p->pdgId()) != ParticleInfo::p_tauminus ) return false;
                const auto* fo = ParticleInfo::getOriginal(p);
                return ParticleInfo::hasMother(&*fo,isFromPromptEWK);
            };

            if(ParticleInfo::hasMother(&*fo,
                    [&](const GenParticle* p) -> bool{ return isFromPromptEWK(p) || isFromPromptTau(p);  }))
            {
                (absID == ParticleInfo::p_eminus  ? promptElectrons :  promptMuons  ).push_back(&*gp);
            } else {
                (absID == ParticleInfo::p_eminus  ? nonPromptElectrons :  nonPromptMuons  ).push_back(&*gp);
            }
        } else if(absID == ParticleInfo::p_t){
            topDecays.push_back(getTopDecay(gp));
        } else if(absID == ParticleInfo::p_Wplus){
            auto fo = ParticleInfo::getOriginal(gp); //get original
            if(ParticleInfo::hasMother(&*fo,ParticleInfo::p_t,false))  continue;
            bosonDecays.push_back(getBosonDecay(gp));
        } else if(absID == ParticleInfo::p_Z0){
            bosonDecays.push_back(getBosonDecay(gp));
        }
    }
    if (topDecays.size() == 2) {
    	genMtt = (topDecays[0].top->p4() + topDecays[1].top->p4()).mass();
    	nLepsTT = 0;
    	if (ttSampParam > 2) {
    		if (topDecays[0].type > TopDecay::HAD) nLepsTT++;
    		if (topDecays[1].type > TopDecay::HAD) nLepsTT++;
    	} else {
    		nLepsTT = ttSampParam;
    	}
    }


}

}



