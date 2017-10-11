#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"

namespace TAna {

enum decayidentifier {Z = 1, HADRON, LEPTON, ELECTRON, MUON, TAU_H, TAU_EL, TAU_MUON};

// Function to search and classify the decay channel of a Tau lepton
int DiHiggsEvent::tau_search(CandidateRef<GenParticle> dau) {
    // get final state of the tau
    auto p = ParticleInfo::getFinal(dau);

    int n = p->numberOfDaughters();
    int e_num = 0;
    int mu_num = 0;
    for (int i = 0; i<n; i++) {
        int pp = p->daughter(i)->absPdgId();
        if (pp == ParticleInfo::p_eminus || pp == ParticleInfo::p_nu_e) {e_num++;}
        if (pp == ParticleInfo::p_muminus || pp == ParticleInfo::p_nu_mu) {mu_num++;}
    }

    if (e_num > 0) {return TAU_EL;}
    if (mu_num > 0) {return TAU_MUON;}
    return TAU_H;
}

// Function that determines if a pair of GenParticles were birthed by a W boson -> 
bool DiHiggsEvent::isWpair(CandidateRef<GenParticle> f1, CandidateRef<GenParticle> f2) {
    // if first particle is a lepton and the second is a neutrino of the same flavor, then they are products of a W
    if (ParticleInfo::isLepton(f1->absPdgId())) {
        if (ParticleInfo::isANeutrino(f2->absPdgId()) && (f2->pdgId() == f1->absPdgId() + 1) && ((f1->pdgId()/f1->absPdgId()) == (-1)*(f2->pdgId()/f2->absPdgId()))) {
            return true;
        }
        else {return false;}
    }
    // if the first is a neutrino and the second is a charged lepton of the same flavour, they come from a W
    else if (ParticleInfo::isANeutrino(f1->absPdgId())) {
        if (ParticleInfo::isLepton(f2->absPdgId()) && (f2->absPdgId() == f1->absPdgId() - 1) && ((f1->pdgId()/f1->absPdgId()) == (-1)*(f2->pdgId()/f2->absPdgId()))) {
            return true;
        }
        else {return false;}
    }
    // if both particles comprise a quark-antiquark pair, then they come from a W 
    // (IMPORTANT: products of a Z will never enter this function in the code)
    else if (ParticleInfo::isQuark(f1->absPdgId())) {
        if ((ParticleInfo::isQuark(f2->absPdgId())) && ((f1->pdgId()/f1->absPdgId()) == (-1)*(f2->pdgId()/f2->absPdgId()))) {
            return true;
        }
        else {return false;}
    }
    else {
        std::cout << "Function isWpair has failed" << std::endl;
        return false;
    }
}

// Function that returns 1 if a pair of particles comes from a Z (particle-antiparticle pair), 2 if from a W, and 0 otherwise
int DiHiggsEvent::isPair(CandidateRef<GenParticle> f1, CandidateRef<GenParticle> f2) {
    // Z
    if (f1->pdgId() == (-1)*f2->pdgId()) {return 1;}
    // W
    else if (isWpair(f1,f2)) {return 2;}
    // otherwise return 0
    else {return 0;}
}

// Function that intakes a pair of candidate W daughters (and so far does nothing if it fails) and classifies the W decay channel
int DiHiggsEvent::classify_W_pair(CandidateRef<GenParticle> p1, CandidateRef<GenParticle> p2) {
    if (ParticleInfo::isQuark(p1->absPdgId())) {return HADRON;}
    else if ((p1->absPdgId() == ParticleInfo::p_eminus) || (p1->absPdgId() == ParticleInfo::p_nu_e)) {return ELECTRON;}
    else if ((p1->absPdgId() == ParticleInfo::p_muminus) || (p1->absPdgId() == ParticleInfo::p_nu_mu)) {return MUON;}

    else if (p1->absPdgId() == ParticleInfo::p_tauminus) {
        int tau_daughter = tau_search(p1);
        return tau_daughter;
    }
    else if (p1->absPdgId() == ParticleInfo::p_nu_tau) {
        int tau_daughter = tau_search(p2);
        return tau_daughter;
    }
    else {return 0;}
}

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

WDecay DiHiggsEvent::assign_W(const GenParticle* w) {
    WDecay w_obj;
    CandidateRef<GenParticle> w_ref(w,0);
    auto w_final = ParticleInfo::getFinal(w_ref);
    w_obj.id = &*w_final;
    if (w_final->numberOfDaughters() == 2) {
        auto dtup = assign_gp(w_final->daughter(0), w_final->daughter(1));
        w_obj.d1 = std::get<0>(dtup);
        w_obj.d2 = std::get<1>(dtup);

        const CandidateRef<GenParticle> dau1(w_obj.d1,0);
        const CandidateRef<GenParticle> dau2(w_obj.d2,0);

        w_obj.decaytype = classify_W_pair(dau1, dau2);
    }
    return w_obj;
}

int DiHiggsEvent::classify_decaytype(std::vector<int> items) {
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
std::vector<int> DiHiggsEvent::search_4_daughters(CandidateRef<GenParticle> gp) {

    std::vector<int> items;

    auto d1 = gp->daughterRef(0);
    auto d2 = gp->daughterRef(1);
    auto d3 = gp->daughterRef(2);
    auto d4 = gp->daughterRef(3);

    auto dd1 = gp->daughter(0);
    auto dd2 = gp->daughter(1);
    auto dd3 = gp->daughter(2);
    auto dd4 = gp->daughter(3);

    int pair12 = isPair(d1,d2);
    int pair13 = isPair(d1,d3);
    int pair14 = isPair(d1,d4);
    int pair23 = isPair(d2,d3);
    int pair24 = isPair(d2,d4);
    int pair34 = isPair(d3,d4);

    // The vector bosons are not explicitly defined in these decays, so they will not be filled

    // If the pair is from a Z, then we have a ZZ decay
    if ((pair12 == 1) && (pair34 == 1)) {
        items.push_back(Z);
        items.push_back(Z);
    }

    // If the pair is from a W
    else if (pair12 == 2) {
        items.push_back(classify_W_pair(d1,d2));
        
        // assignment of pair12
        auto op12 = assign_gp(dd1, dd2);
        w1_d1 = std::get<0>(op12);
        w1_d2 = std::get<1>(op12);

        // Now look at the d3,d4 pair, which should be another W if there are no errors
        if (pair34 == 2) {
            items.push_back(classify_W_pair(d3,d4));

            // assignment of pair34
            auto op34 = assign_gp(dd3, dd4);
            w2_d1 = std::get<0>(op34);
            w2_d2 = std::get<1>(op34);
        }
    }

    // If a valid pair is not found, search between f1 and f3, and then if that fails as well, search between f1 and f4
    else if (pair12 == 0) {
        if (pair13 == 1) {
            items.push_back(Z);
            items.push_back(Z);
        }
        else if (pair13 == 2) {
            items.push_back(classify_W_pair(d1,d3));

            auto op13 = assign_gp(dd1, dd3);
            w1_d1 = std::get<0>(op13);
            w1_d2 = std::get<1>(op13);

            // Now look at the d2,d4 pair, which should be another W if there are no errors
            if (pair24 == 2) {
                items.push_back(classify_W_pair(d2,d4));

                auto op24 = assign_gp(dd2, dd4);
                w2_d1 = std::get<0>(op24);
                w2_d2 = std::get<1>(op24);
            }
        }
        else {
            if (pair14 == 1) {
                items.push_back(Z);
                items.push_back(Z);
            } 
            else if (pair14 == 2) {
                items.push_back(classify_W_pair(d1,d4));

                auto op14 = assign_gp(dd1, dd4);
                w1_d1 = std::get<0>(op14);
                w1_d2 = std::get<1>(op14);

                if (pair23 == 2) {
                    items.push_back(classify_W_pair(d2,d3));

                    auto op23 = assign_gp(dd2, dd3);
                    w2_d1 = std::get<0>(op23);
                    w2_d2 = std::get<1>(op23);
                }
            }
        }
    }
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
        if ((gp->daughter(0)->absPdgId() == ParticleInfo::p_b) && (gp->daughter(1)->absPdgId() == ParticleInfo::p_b) && (gp->numberOfDaughters() == 2)) {
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
            std::vector<CandidateRef<GenParticle>> fermions;
            std::vector<WDecay> wDecays;

            // iterate over the daughters of this Higgs
            for (int k = 0; k < hww->numberOfDaughters(); k++) {
                if (hww->daughter(k)->absPdgId() == ParticleInfo::p_Z0) {
                    zbosons.push_back(Z);
                } else if (hww->daughter(k)->absPdgId() == ParticleInfo::p_Wplus) {
                    wDecays.push_back(assign_W(hww->daughter(k)));
                } else {
                    fermions.push_back(hww->daughterRef(k));
                }
            }

    // Now analyze the combinations of sizes of each of these vectors as different cases
            std::vector<int> fermIDs;
            // if items is of size 2 already, then it must have two Z bosons and we can immediately write down the class variable type
            if (zbosons.size() == 2) {type = bbZZ;}

            // if there is one Z and two fermions, check to ensure the two fermions can form a Z. Otherwise, it is a BAD decay
            else if ((zbosons.size() == 1) && (fermions.size() == 2)) {
                if (isPair(fermions[0], fermions[1]) == 1) {
                    type = bbZZ;
                } else {type = BAD;}
            }
            // if fermions is of size 4, the protocol will be same as search_4_daughters
            else if (fermions.size() == 4) {
                CandidateRef<GenParticle> href(hww,0);
                fermIDs = search_4_daughters(href); // assignment of GenParticle class variables executed within this function
                type = classify_decaytype(fermIDs);
            }
            // if there is one W and two fermions
            else if (wDecays.size() == 1 && fermions.size() == 2){
                if (ParticleInfo::isLepton(fermions[0]->absPdgId()) || ParticleInfo::isANeutrino(fermions[0]->absPdgId())) {
                    auto fs = assign_gp(&*fermions[0], &*fermions[1]);
                    w1_d1 = std::get<0>(fs);
                    w1_d2 = std::get<1>(fs);

                    w2 = wDecays[0].id;
                    w2_d1 = wDecays[0].d1;
                    w2_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(fermions[0], fermions[1]));

                    if (dtypes[0] == 0) {ParticleInfo::printGenInfo(genparts,-1);}
                    type = classify_decaytype(dtypes);
                    if (type == BAD) {std::cout << "yo ova here dickwad"<<std::endl;}

                } else {
                    auto fs = assign_gp(&*fermions[0], &*fermions[1]);
                    w2_d1 = std::get<0>(fs);
                    w2_d2 = std::get<1>(fs);

                    w1 = wDecays[0].id;
                    w1_d1 = wDecays[0].d1;
                    w1_d2 = wDecays[0].d2;

                    std::vector<int> dtypes;
                    dtypes.push_back(wDecays[0].decaytype);
                    dtypes.push_back(classify_W_pair(fermions[0], fermions[1]));
                    type = classify_decaytype(dtypes);
                    if (type == BAD) {std::cout << "yo ova here damlingus"<<std::endl;}

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

