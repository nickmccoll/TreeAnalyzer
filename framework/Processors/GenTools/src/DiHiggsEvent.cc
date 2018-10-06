#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"

namespace TAna {

enum decayidentifier {Z = 1, HADRON, LEPTON, ELECTRON, MUON, TAU_H, TAU_EL, TAU_MUON};

// Function to search and classify the decay channel of a Tau lepton
int DiHiggsEvent::tau_search(const GenParticle* dau) {
    // get final state of the tau
    auto p = ParticleInfo::getFinal(dau);

    int n = p->numberOfDaughters();
    for (int i = 0; i<n; i++) {
        int pp = p->daughter(i)->absPdgId();
        if (pp == ParticleInfo::p_eminus || pp == ParticleInfo::p_nu_e) {return TAU_EL;}
        if (pp == ParticleInfo::p_muminus || pp == ParticleInfo::p_nu_mu) {return TAU_MUON;}
    }
    return TAU_H;
}

// Function that determines if a pair of GenParticles were birthed by a W boson -> 
bool DiHiggsEvent::isWpair(const GenParticle* f1, const GenParticle* f2, const GenParticleCollection& g_parts) {
    int f1id = f1->pdgId();
    int f2id = f2->pdgId();

    // return false if one of the two particles is a photon
    if (f1id==22 || f2id==22) {
    	std::cout << "Found photon, not a W pair"<< std::endl;
    	return false;
    }
    // if first particle is a lepton and the second is a neutrino of the same flavor, then they are products of a W
    else if (ParticleInfo::isLepton(f1id)) {
        if (ParticleInfo::isANeutrino(f2id) && (std::abs(f2id) == std::abs(f1id) + 1) && ((f1id < 0) != (f2id < 0))) {
            return true;
        }
        else {return false;}
    }
    // if the first is a neutrino and the second is a charged lepton of the same flavour, they come from a W
    else if (ParticleInfo::isANeutrino(f1id)) {
        if (ParticleInfo::isLepton(f2id) && (std::abs(f2id) == std::abs(f1id) - 1) && ((f1id < 0) != (f2id < 0))) {
            return true;
        }
        else {return false;}
    }
    // if both particles comprise a quark-antiquark pair, then they come from a W 
    // (IMPORTANT: products of a Z will never enter this function in the code)
    else if (ParticleInfo::isQuark(f1id)) {
        if ((ParticleInfo::isQuark(f2id)) && ((f1id<0) != (f2id<0))) {
            return true;
        }
        else {return false;}
    }
    else {
        std::cout << "Function isWpair has failed" << std::endl;
	ParticleInfo::printGenInfo(g_parts,-1);
        return false;
    }
}

// Function that returns 1 if a pair of particles comes from a Z (particle-antiparticle pair), 2 if from a W, and 0 otherwise
int DiHiggsEvent::isPair(const GenParticle* f1, const GenParticle* f2, const GenParticleCollection& g_parts) {
    // Z
    if (f1->pdgId() == (-1)*f2->pdgId()) {return 1;}
    // W
    else if (isWpair(f1,f2,g_parts)) {return 2;}
    // otherwise return 0
    else {return 0;}
}

// Function that intakes a pair of candidate W daughters (and so far does nothing if it fails) and classifies the W decay channel
int DiHiggsEvent::classify_W_pair(const GenParticle* p1, const GenParticle* p2) {
    int p1id = p1->absPdgId();

    if (ParticleInfo::isQuark(p1id)) {return HADRON;}
    else if ((p1id == ParticleInfo::p_eminus) || (p1id == ParticleInfo::p_nu_e)) {return ELECTRON;}
    else if ((p1id == ParticleInfo::p_muminus) || (p1id == ParticleInfo::p_nu_mu)) {return MUON;}

    else if (p1id == ParticleInfo::p_tauminus) {
        int tau_daughter = tau_search(p1);
        return tau_daughter;
    }
    else if (p1id == ParticleInfo::p_nu_tau) {
        int tau_daughter = tau_search(p2);
        return tau_daughter;
    }
    else {return 0;}
}

// Function that takes two GenParticles and organizes such that if they are leptonic, d1 is the Lepton, d2 is the Neutrino
std::tuple<const GenParticle*, const GenParticle*> DiHiggsEvent::assign_gp(const GenParticle* p1, const GenParticle* p2) {
    const GenParticle* d1;
    const GenParticle* d2;
    if (ParticleInfo::isLepton(p2->absPdgId()) && ParticleInfo::isANeutrino(p1->absPdgId())) {
        d1 = p2;
        d2 = p1;
    } else {
        d1 = p1;
        d2 = p2;
    }
    return std::make_tuple(d1,d2);
}

// Function that fills the data variables of the WDecay structure 
DiHiggsEvent::WDecay DiHiggsEvent::assign_W(const GenParticle* w, const GenParticleCollection& gen_parts) {
    WDecay w_obj;
    auto w_final = ParticleInfo::getFinal(w);

    w_obj.id = w_final;
    for (int j=0; j<w_final->numberOfDaughters()-1; j++) {
	for (int k=1; k<w_final->numberOfDaughters(); k++) {
	    if (isWpair(w_final->daughter(j), w_final->daughter(k), gen_parts)) {
	       auto dtup = assign_gp(w_final->daughter(j), w_final->daughter(k));
	       w_obj.d1 = std::get<0>(dtup);
	       w_obj.d2 = std::get<1>(dtup);
	       break;
	    }
	}
    }
    w_obj.decaytype = classify_W_pair(w_obj.d1, w_obj.d2);
    return w_obj;
}

// Function that takes two identifiers for how the Higgs daughters decay, and then classifies the total decay channel of X -> HH
DiHiggsEvent::DECAYTYPE DiHiggsEvent::classify_decaytype(const std::vector<int>& items) {
    if (items.size() != 2) {
        type = BAD;
    }
    else if (items.size() == 2) {
        // if items is holding two Z's, then fill the event as a ZZ decay into the histogram
        if ((items[0] == Z) && (items[1] == Z)) {
            type = bbZZ;
        }
        // if items is holding two Hadrons, then fill the event as an HH decay
        if ((items[0] == HADRON) && (items[1] == HADRON)) {
            type = HAD;
        }

        int electron_num = 0;
        int muon_num = 0;
        int tauH_num = 0;
        int tauE_num = 0;
        int tauMU_num = 0;
        int hadron_num = 0;

        // for each value in items:
        for (int k = 0 ; k < 2 ; k++) {
            if (items[k] == HADRON) {hadron_num++;}
            else if (items[k] == ELECTRON) {electron_num++;}
            else if (items[k] == MUON) {muon_num++;}
            else if (items[k] == TAU_H) {tauH_num++;}
            else if (items[k] == TAU_EL) {tauE_num++;}
            else if (items[k] == TAU_MUON) {tauMU_num++;}
        }
        // if there are two kinds of Leptons stored in items, then fill the event as a LL decay
        if (electron_num + muon_num + tauH_num + tauE_num + tauMU_num == 2) {type = DILEP;}
        // if there is one Lepton type and one Hadron in items:
        if ((electron_num + muon_num + tauH_num + tauE_num + tauMU_num == 1) && (hadron_num == 1)) {
            if (electron_num == 1) {type = E;}
            else if (muon_num == 1) {type = MU;}
            else if (tauH_num == 1) {type = TAU_HAD;}
            else if (tauE_num == 1) {type = TAU_E;}
            else if (tauMU_num == 1) {type = TAU_MU;}
            else {type = BAD;}
        }
    }
    return type;
}

// Function to classify the decay mode when the Higgs has 4 daughters in the Event. Returns a vector of ints (enums), items
std::vector<int> DiHiggsEvent::search_4_daughters(const GenParticle* gp, const GenParticleCollection& g_parts) {

    std::vector<int> items;

    auto mkW = [&] (const GenParticle* d11, const GenParticle* d12, const GenParticle* d21, const GenParticle* d22) {
        items.push_back(classify_W_pair(d11,d12));
        auto op12 = assign_gp(d11,d12);
        w1_d1 = std::get<0>(op12);
        w1_d2 = std::get<1>(op12);

        items.push_back(classify_W_pair(d21,d22));
        auto op34 = assign_gp(d21,d22);
        w2_d1 = std::get<0>(op34);
        w2_d2 = std::get<1>(op34);
    };

    auto d1 = gp->daughter(0);
    auto d2 = gp->daughter(1);
    auto d3 = gp->daughter(2);
    auto d4 = gp->daughter(3);

    int pair12 = isPair(d1,d2,g_parts);
    int pair13 = isPair(d1,d3,g_parts);
    int pair14 = isPair(d1,d4,g_parts);
    int pair23 = isPair(d2,d3,g_parts);
    int pair24 = isPair(d2,d4,g_parts);
    int pair34 = isPair(d3,d4,g_parts);

    // The vector bosons are not explicitly defined in these decays, so they will not be filled

    if(pair12==1 && pair34==1) {items.push_back(Z);items.push_back(Z);}
    else if(pair12==2 && pair34==2) mkW(d1,d2,d3,d4);
    else if(pair13==1 && pair24==1) {items.push_back(Z);items.push_back(Z);}
    else if(pair13==2 && pair24==2) mkW(d1,d3,d2,d4);
    else if(pair14==1 && pair23==1) {items.push_back(Z);items.push_back(Z);}
    else if(pair14==2 && pair23==2) mkW(d1,d4,d2,d3);
    
    return items;
}

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

  // For each event, iterate down the list of genparticles to classify the decay channel
    for(const auto& p : genparts) {
        const auto* gp = &p;
        // if not a Higgs in final state, continue
        if (gp->absPdgId() != ParticleInfo::p_h0) {continue;}
        if (!ParticleInfo::isLastInChain(gp)) {continue;}
        // if this Higgs decays to b-bbar, set the appropriate GenParticles and then continue to the next iteration
        if ((gp->numberOfDaughters() == 2) && (gp->daughter(0)->absPdgId() == ParticleInfo::p_b) && (gp->daughter(1)->absPdgId() == ParticleInfo::p_b)) {
            hbb = gp;
            if(hbb->daughter(0)->pt() > hbb->daughter(1)->pt()) {
                b1 = hbb->daughter(0);
                b2 = hbb->daughter(1);
            } else { 
                b1 = hbb->daughter(1);
                b2 = hbb->daughter(0);
            }
        }
        else {
            hww = gp;
    // Here is where I will organize the daughters of the Hww appropriately
            std::vector<int> zbosons;
            std::vector<const GenParticle*> quarleps;
            std::vector<WDecay> wDecays;

            // iterate over the daughters of this Higgs
            for (int k = 0; k < hww->numberOfDaughters(); k++) {
		auto dau = hww->daughter(k);
                if (dau->absPdgId() == ParticleInfo::p_Z0) {
                    zbosons.push_back(Z);
                } else if (dau->absPdgId() == ParticleInfo::p_Wplus) {
                    wDecays.push_back(assign_W(dau, genparts));
                } else if (ParticleInfo::isQuark(dau->pdgId()) || ParticleInfo::isLeptonOrNeutrino(dau->pdgId())){
                    quarleps.push_back(dau);
                } else {continue;}
            }

    // Now analyze the combinations of sizes of each of these vectors as different cases
            std::vector<int> fermIDs;
            // if zbosons is of size 2 already, then it must have two Z bosons and we can immediately write down the class variable type
            if (zbosons.size() == 2) {type = bbZZ;}

            // if there is one Z and two quarleps, check to ensure the two quarleps can form a Z. Otherwise, it is a BAD decay
            else if ((zbosons.size() == 1) && (quarleps.size() == 2)) {
                if (isPair(quarleps[0], quarleps[1],genparts) == 1) {
                    type = bbZZ;
                } else {type = BAD;}
            }
            // if quarleps is of size 4, and if the W and Z vectors are empty, the protocol will be same as search_4_daughters
            else if ((quarleps.size() == 4) && (zbosons.size() + wDecays.size() == 0)) {
                fermIDs = search_4_daughters(hww,genparts); // assignment of GenParticle class variables executed within this function
                type = classify_decaytype(fermIDs);
            }
	    // if there is one W and more than 2 quarleps, then need to identify a combination that is a valid W pair
	    else if (wDecays.size() == 1 && quarleps.size() > 2) {
                std::vector<std::pair <const GenParticle*, const GenParticle*>> pair_cands;
		for (unsigned int j=0; j<(quarleps.size()-1); j++) {
		    for (unsigned int k=0; k<quarleps.size(); k++) {
			if (isWpair(quarleps[j], quarleps[k], genparts)) {
			   pair_cands.push_back(std::make_pair(quarleps[j], quarleps[k]));
			}
		    }
		}
		std::vector<const GenParticle*> valid_pair;
		valid_pair.push_back(pair_cands[0].first);
		valid_pair.push_back(pair_cands[0].second);
		if (valid_pair.size() != 2) {std::cout << "Error when one W and > 2 quarleps of the Higgs daughters" << std::endl; ParticleInfo::printGenInfo(genparts,-1);}
		else {
		  if (ParticleInfo::isLepton(valid_pair[0]->absPdgId()) || ParticleInfo::isANeutrino(valid_pair[0]->absPdgId())) {
                    auto fs = assign_gp(valid_pair[0], valid_pair[1]);
                    w1_d1 = std::get<0>(fs);
                    w1_d2 = std::get<1>(fs);

                    w2 = wDecays[0].id;
                    w2_d1 = wDecays[0].d1;
                    w2_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(valid_pair[0], valid_pair[1]));

                    if (dtypes[0] == 0) {std::cout<<"Decaytype unfilled"<<std::endl; ParticleInfo::printGenInfo(genparts,-1);}
                    type = classify_decaytype(dtypes);

                  } else {
                    auto fs = assign_gp(valid_pair[0], valid_pair[1]);
                    w2_d1 = std::get<0>(fs);
                    w2_d2 = std::get<1>(fs);

                    w1 = wDecays[0].id;
                    w1_d1 = wDecays[0].d1;
                    w1_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(valid_pair[0], valid_pair[1]));
                    type = classify_decaytype(dtypes);
                  } 
		}
	    }
            // if there is one W and two quarleps
            else if (wDecays.size() == 1 && quarleps.size() == 2){
                if (ParticleInfo::isLepton(quarleps[0]->absPdgId()) || ParticleInfo::isANeutrino(quarleps[0]->absPdgId())) {
                    auto fs = assign_gp(quarleps[0], quarleps[1]);
                    w1_d1 = std::get<0>(fs);
                    w1_d2 = std::get<1>(fs);

                    w2 = wDecays[0].id;
                    w2_d1 = wDecays[0].d1;
                    w2_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(quarleps[0], quarleps[1]));

                    if (dtypes[0] == 0) {std::cout<<"Decaytype unfilled"<<std::endl; ParticleInfo::printGenInfo(genparts,-1);}
                    type = classify_decaytype(dtypes);

                } else {
                    auto fs = assign_gp(quarleps[0], quarleps[1]);
                    w2_d1 = std::get<0>(fs);
                    w2_d2 = std::get<1>(fs);

                    w1 = wDecays[0].id;
                    w1_d1 = wDecays[0].d1;
                    w1_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(quarleps[0], quarleps[1]));
                    type = classify_decaytype(dtypes);
                }
            }
            // if there are two W objects
            else if (wDecays.size() == 2) {
                std::vector<int> dtypes;
                dtypes.push_back(wDecays[0].decaytype);
                dtypes.push_back(wDecays[1].decaytype);
                type = classify_decaytype(dtypes);

                // we need to classify the leptonic W as w1, and we need to classify the daughters now
                if (ParticleInfo::isLepton(wDecays[1].d1->absPdgId()) || ParticleInfo::isANeutrino(wDecays[1].d1->absPdgId())) {
                    w1 = wDecays[1].id;
                    w1_d1 = wDecays[1].d1;
                    w1_d2 = wDecays[1].d2;

                    w2 = wDecays[0].id;
                    w2_d1 = wDecays[0].d1;
                    w2_d2 = wDecays[0].d2;
                } else {
                    w1 = wDecays[0].id;
                    w1_d1 = wDecays[0].d1;
                    w1_d2 = wDecays[0].d2;

                    w2 = wDecays[1].id;
                    w2_d1 = wDecays[1].d1;
                    w2_d2 = wDecays[1].d2;
                }
            } else {type = BAD;}

            if (hbb && hww) {break;}
        }
    } 
}

}

