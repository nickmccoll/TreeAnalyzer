#include "Processors/Variables/interface/DileptonSelection.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "DataFormats/interface/GenParticle.h"
#include "Configuration/interface/FillerConstants.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "../interface/SignalHelper.h"

namespace TAna {
//_____________________________________________________________________________
//FatJet* getMatchedFJ(const MomentumF& genJet, const std::vector<FatJet>& fatjets, double maxDR) {
//	double nearestDR = 99;
//	if (!fatjets.size()) return 0;
//    int idx = PhysicsUtilities::findNearestDR(genJet,fatjets,nearestDR,maxDR);
//    if (idx < 0) return 0;
//    else return &fatjets[idx];
//}
//_____________________________________________________________________________
const GenParticle* SignalHelper::getGenLepFromTau(const GenParticle* TAU) {
	if (TAU->absPdgId() != ParticleInfo::p_tauminus) return 0;

//	std::cout<<"in gen tau func"<<std::endl;
	const GenParticle* tau=TAU;
	if (!ParticleInfo::isLastInChain(tau)) {
		std::cout<<"not last in chain"<<std::endl;
		tau = ParticleInfo::getFinal(TAU);
	} else {
		tau = TAU;
	}
//	std::cout<<"got final"<<std::endl;

	const GenParticle* gp=0;
	int nlepsfromtau = 0;
	for (unsigned int k = 0; k < tau->numberOfDaughters(); k++) {
		const auto* dau = tau->daughter(k);
		if (dau->absPdgId() == ParticleInfo::p_muminus || dau->absPdgId() == ParticleInfo::p_eminus) {

			if (nlepsfromtau >= 1) {
				continue;
			}
			gp = dau;
			nlepsfromtau++;
		}
	}
//	std::cout<<"after loop"<<std::endl;

	if (!gp) return 0;
	if (nlepsfromtau > 1) {
		return 0;
	}

	return gp;
}

SignalHelper::SignalHelper(DiHiggsEvent dhEvt, std::shared_ptr<MuonReader> reader_muon,
		std::shared_ptr<ElectronReader> reader_electron) {
	elReader = reader_electron;
	muReader = reader_muon;

	type = dhEvt.type;
	if (type == DiHiggsEvent::DILEP) {
		const GenParticle* glep1_0 = 0;
		const GenParticle* glep2_0 = 0;

		if(dhEvt.w1_d1->pt() > dhEvt.w2_d1->pt()) {
			glep1_0 = dhEvt.w1_d1;
			glep2_0 = dhEvt.w2_d1;
		} else {
			glep1_0 = dhEvt.w2_d1;
			glep2_0 = dhEvt.w1_d1;
		}

		if (glep1_0 && glep2_0 && (glep1_0->pdgId()<0 != glep2_0->pdgId()<0)) {

	    	if (glep1_0->absPdgId() == ParticleInfo::p_tauminus) genlep1 = getGenLepFromTau(glep1_0);
	    	else genlep1 = glep1_0;

	    	if (glep2_0->absPdgId() == ParticleInfo::p_tauminus) genlep2 = getGenLepFromTau(glep2_0);
	    	else genlep2 = glep2_0;
		} else {
			genlep1 = 0;
			genlep2 = 0;
		}
	} else if (type >= DiHiggsEvent::TAU_MU) {
//		std::cout<<"inside 1l constructor"<<std::endl;
		genlep2 = 0;

//		std::cout<<"is w1_d1? "<<dhEvt.w1_d1->absPdgId()<<std::endl;
//		std::cout<<"is w2_d1? "<<dhEvt.w2_d1->absPdgId()<<std::endl;

		if (ParticleInfo::isLepton(dhEvt.w1_d1->absPdgId())) {
//			std::cout<<"here"<<std::endl;
			genlep1 = (dhEvt.w1_d1->absPdgId() == ParticleInfo::p_tauminus) ?
					getGenLepFromTau(dhEvt.w1_d1) : dhEvt.w1_d1;
		} else if (ParticleInfo::isLepton(dhEvt.w2_d1->absPdgId())) {
			genlep1 = (dhEvt.w2_d1->absPdgId() == ParticleInfo::p_tauminus) ?
					getGenLepFromTau(dhEvt.w2_d1) : dhEvt.w2_d1;
		} else {
			throw std::invalid_argument("one of w1_d1 or w2_d1 should be a lepton");
		}
//		std::cout<<"ending 1l constructor"<<std::endl;

	}
}


Lepton* SignalHelper::getMatchedLepton(const GenParticle* genLep,double maxDR) {
    if (!genLep) return 0;

	double nearestDR =10;
    int idx = -1;
	if(genLep->absPdgId() == ParticleInfo::p_muminus){
	    if(!candMuons.size()) return 0;
   	    for (const auto& mu : candMuons) {
   		    if ((mu->q() > 0) == (genLep->pdgId() > 0)) continue;
   	        double dr = PhysicsUtilities::deltaR(*genLep,*mu);
   	        if (dr < nearestDR) {
   	            nearestDR = dr;
       	        idx = mu->index();
   	        }
   	    }
   	    if (nearestDR > maxDR) return 0;
   	    if (idx < 0) return 0;
   	    return &muReader->muons[idx];
    } else {
	    if(!candElectrons.size()) return 0;
   	    for (const auto& el : candElectrons) {
   		    if ((el->q() > 0) == (genLep->pdgId() > 0)) continue;
   	        double dr = PhysicsUtilities::deltaR(*genLep,*el);
   	        if (dr < nearestDR) {
   	            nearestDR = dr;
       	        idx = el->index();
   	        }
   	    }
   	    if (nearestDR > maxDR) return 0;
   	    if (idx < 0) return 0;
   	    return &elReader->electrons[idx];
    }
}

void SignalHelper::setRecoLeptons(double matchDR) {
	candMuons = PhysicsUtilities::selObjsMom(muReader->muons,minMuRecoPt,maxMuRecoEta);
	candElectrons = PhysicsUtilities::selObjsMom(elReader->electrons,minElRecoPt,maxElRecoEta);

	if (!candMuons.size() && !candElectrons.size()) {
		recolep1 = 0;
		recolep2 = 0;
	} else {
		recolep1 = getMatchedLepton(genlep1,matchDR);
		recolep2 = getMatchedLepton(genlep2,matchDR);
	}

}

bool SignalHelper::hasMatchedSingleLep() {
	if (!genlep1 || !recolep1) return false;
	return true;
}

bool SignalHelper::hasMatchedDileps() {
	if (!genlep1 || !genlep2 || !recolep1 || !recolep2) return false;
	return true;
}

}
