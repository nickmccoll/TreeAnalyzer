#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"

namespace TAna {

enum decayidentifier {Z = 1, HADRON, LEPTON, ELECTRON, MUON, TAU_H, TAU_EL, TAU_MUON};

// Function to search and classify the decay channel of a Tau lepton
int DiHiggsEvent::Tau_search(CandidateRef<GenParticle> dau) {
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
    if ((p1->absPdgId() == ParticleInfo::p_eminus) || (p1->absPdgId() == ParticleInfo::p_nu_e)) {return ELECTRON;}
    if ((p1->absPdgId() == ParticleInfo::p_muminus) || (p1->absPdgId() == ParticleInfo::p_nu_mu)) {return MUON;}

    if (p1->absPdgId() == ParticleInfo::p_tauminus) {
        int tau_daughter = Tau_search(p1);
        return tau_daughter;
    }
    if (p1->absPdgId() == ParticleInfo::p_nu_tau) {
        int tau_daughter = Tau_search(p2);
        return tau_daughter;
    }
    return 0;
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

//-------------------------------------------------------------------------------------------------------------------------------------//
                                                // FUNCTION SEARCH_4_DAUGHTERS //
//-------------------------------------------------------------------------------------------------------------------------------------//

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

    // const GenParticle* w1;
    // const GenParticle* w2;
    // const GenParticle* w1_d1;
    // const GenParticle* w1_d2;
    // const GenParticle* w2_d1;
    // const GenParticle* w2_d2;

    // If the pair is from a Z, then we have a ZZ decay
    if ((pair12 == 1) && (pair34 == 1)) {
        items.push_back(Z);
        items.push_back(Z);

        // Pair12 assignments
        w1_d1 = dd1;
        w1_d2 = dd2;

        // Pair34 assignments 
        w2_d1 = dd3;
        w2_d2 = dd4;
    }

    // If the pair is from a W
    if (pair12 == 2) {
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
    if (pair12 == 0) {
        if (pair13 == 1) {
            items.push_back(Z);
            items.push_back(Z);

            w1_d1 = dd1;
            w1_d2 = dd3;

            w2_d1 = dd2;
            w2_d2 = dd4;
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

                w1_d1 = dd1;
                w1_d2 = dd4;

                w2_d1 = dd2;
                w2_d2 = dd3;

            } else if (pair14 == 2) {
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
//-------------------------------------------------------------------------------------------------------------------------------------//
                                                // FUNCTION SEARCH_3_DAUGHTERS //
//-------------------------------------------------------------------------------------------------------------------------------------//

std::vector<int> DiHiggsEvent::search_3_daughters(CandidateRef<GenParticle> gp) {

    std::vector<int> items;
    std::vector<CandidateRef<GenParticle>> fermions;

    // const GenParticle* w1;
    // const GenParticle* w2;
    // const GenParticle* w1_d1;
    // const GenParticle* w1_d2;
    // const GenParticle* w2_d1;
    // const GenParticle* w2_d2;

    for (int i=0 ; i < 3 ; i++) {
        auto dau = gp->daughterRef(i);

        // if the first daughter particle is a Z, then the other decay product of note is a Z. Push two Zs to items and break the loop
        if (dau->absPdgId() == ParticleInfo::p_Z0) {
            items.push_back(Z);
            items.push_back(Z);

            // assign the Z to w1 and set w2 = 0
            w1 = &*ParticleInfo::getFinal(dau);

            // Assign the daughters of the explicit Z to daughters of w1 (assign the remaining daughters after the loop)
            if (w1->numberOfDaughters() == 2) {
                w1_d1 = w1->daughter(0);
                w1_d2 = w1->daughter(1);
            } else {std::cout << "Warning: this Z does not have two daughters" << std::endl;}
            

        } else if (dau->absPdgId() == ParticleInfo::p_Wplus) {
        // if the daughter particle is a W, then we need to probe the daughters of this W but can directly determine one item from here

            // Assign this W and its two direct daughters
            w1 = &*ParticleInfo::getFinal(dau);

            // Conditional statements to classify W daughters (want d1 to be lepton if possible)
            if (w1->numberOfDaughters() == 2) {
                auto wdaus = assign_gp(w1->daughter(0),w1->daughter(1));
                w1_d1 = std::get<0>(wdaus);
                w1_d2 = std::get<1>(wdaus);
            } else {std::cout << "Warning: this W does not have two daughters" << std::endl;}

            // Make these GenParticles a Ref so I can plug them into classify_W_pair
            const CandidateRef<GenParticle> d1(w1_d1,0);
            const CandidateRef<GenParticle> d2(w1_d2,0);

            int identifier = classify_W_pair(d1,d2);
            items.push_back(identifier);
        } else {
        // otherwise, this is a fermion; I will classify this later
            fermions.push_back(dau);
        }
    }
    // At this point, due to the structure of the event tree, I expect two Genparticles from a boson to populate the fermions vector
    // Thus, I only need to run them through function classify_W_pair and apend the output to the items vector
    if (fermions.size() == 2) {
        auto ferms = assign_gp(&*fermions[0], &*fermions[1]);
        w2_d1 = std::get<0>(ferms);
        w2_d2 = std::get<1>(ferms);

        if (items.size() != 2) {items.push_back(classify_W_pair(fermions[0],fermions[1]));}
    } else {
        std::cout << "Error encountered: fermions vector not of size 2" << std::endl;
    }
    return items;
}
//-------------------------------------------------------------------------------------------------------------------------------------//
                                                // FUNCTION SEARCH_2_DAUGHTERS //
//-------------------------------------------------------------------------------------------------------------------------------------//

std::vector<int> DiHiggsEvent::search_2_daughters(CandidateRef<GenParticle> gp) {
    std::vector<int> items;

    // const GenParticle* w1;
    // const GenParticle* w2;
    // const GenParticle* w1_d1;
    // const GenParticle* w1_d2;
    // const GenParticle* w2_d1;
    // const GenParticle* w2_d2;

    // if both daughters are Z bosons, fill the items vector with Z twice
    if ((gp->daughter(0)->absPdgId() == ParticleInfo::p_Z0) && (gp->daughter(1)->absPdgId() == ParticleInfo::p_Z0)) {
        items.push_back(Z);
        items.push_back(Z);

        // Assign both Z's (their final forms) to w1 and w2
        w1 = &*ParticleInfo::getFinal(gp->daughterRef(0));
        w2 = &*ParticleInfo::getFinal(gp->daughterRef(1));

        if ((w1->numberOfDaughters() == 2) && (w2->numberOfDaughters() == 2)) {
            auto zdaus = assign_gp(w1->daughter(0), w1->daughter(1));
            auto z2daus = assign_gp(w2->daughter(0), w2->daughter(1));

            w1_d1 = std::get<0>(zdaus);
            w1_d2 = std::get<1>(zdaus);
            w2_d1 = std::get<0>(z2daus);
            w2_d2 = std::get<1>(z2daus);

        } else {std::cout << "Warning: At least one Z does not have two daughters here" << std::endl;}
    }

        // otherwise, if both daughters are a W+ W- pair, then probe further and call classify_W_pair twice to fill the items vector
    if ((gp->daughter(0)->absPdgId() == ParticleInfo::p_Wplus) && (gp->daughter(1)->absPdgId() == ParticleInfo::p_Wplus)) {

        // Assign both W's to their appropriate daughter GenParticles
        w1 = &*ParticleInfo::getFinal(gp->daughterRef(0));
        w2 = &*ParticleInfo::getFinal(gp->daughterRef(1));

        //std::cout << "the Higgs->WW loop has been initiated" << std::endl;
        if ((w1->numberOfDaughters() == 2) && (w2->numberOfDaughters() == 2)) {
            items.push_back(classify_W_pair(w1->daughterRef(0), w1->daughterRef(1)));
            items.push_back(classify_W_pair(w2->daughterRef(0), w2->daughterRef(1)));

            auto w1daus = assign_gp(w1->daughter(0), w1->daughter(1));
            auto w2daus = assign_gp(w2->daughter(0), w2->daughter(1));

            w1_d1 = std::get<0>(w1daus);
            w1_d2 = std::get<1>(w1daus);
            w2_d1 = std::get<0>(w2daus);
            w2_d2 = std::get<1>(w2daus);
        } else {std::cout << "Warning: At least one W does not have two daughters here" << std::endl;}
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

    std::vector<int> items;
  // For each event, iterate down the list of genparticles to classify the decay channel
    for(unsigned int i = 0; i < genparts.size(); i++) {
        const auto& gp = genparts[i];

        // Start search when I hit a Higgs
        if (gp.absPdgId() == ParticleInfo::p_h0) {
            const auto* gpp = &gp;
            CandidateRef<GenParticle> curr_Ref(gpp,i); 

            auto gpfinal = ParticleInfo::getFinal(curr_Ref);

            // if this Higgs decays to b-bbar, set the appropriate GenParticles and then continue to the next iteration
            if ((gpfinal->daughter(0)->absPdgId() == ParticleInfo::p_b) && (gpfinal->numberOfDaughters() != 4)) {
                hbb = &*gpfinal;
                if(gpfinal->daughter(0)->pt() > gpfinal->daughter(1)->pt()) {
                    b1 = &(*gpfinal->daughter(0));
                    b2 = &(*gpfinal->daughter(1));
                } else { 
                    b1 = &(*gpfinal->daughter(1));
                    b2 = &(*gpfinal->daughter(0));
                }
                continue;
            }
            else {
                hww = &*gpfinal;
// Approach the classification by separately analyzing cases where the Higgs has different daughter amounts

            // 2 daughters
                if (gpfinal->numberOfDaughters() == 2){
                    items = search_2_daughters(gpfinal);
                    //std::cout << "items size = " << items.size() << std::endl;
                }
            // 3 daughters
                if (gpfinal->numberOfDaughters() == 3) {
                    items = search_3_daughters(gpfinal);
                    //std::cout << "items size = " << items.size() << std::endl;
                }
            // 4 daughters
                if (gpfinal->numberOfDaughters() == 4) {
                    items = search_4_daughters(gpfinal);
                    //std::cout << "items size = " << items.size() << std::endl;
                }

            //printf("Second Higgs has number of daughters: %i\n", gp.numberOfDaughters());
/*
            if ((gp.numberOfDaughters() == 3) && (gp.daughter(0)->absPdgId() == ParticleInfo::p_Z0))
            {
                ParticleInfo::printGenInfo(genparticles,-1);
            } 

                std::cout << "Elements of items: ";
                for (unsigned int i = 0; i < items.size() ; i++) {
                    if (i == items.size()-1) {
                        std::cout << items[i] << std::endl;
                        continue;
                    }
                    std::cout << items[i] << ", ";
                }
                */
                break;
            }
        } 
    }
// At this point, we should have completed filling the vector items. Now next bit of code will analyze this vector and 
// appropriately fill a histogram

    // if items does not have two elements, then fill the event as an error/miscellaneous
    if (items.size() != 2) {
        std::cout << "Filled as Error" << std::endl;
        ParticleInfo::printGenInfo(genparts,-1);
        // ParticleInfo::printGenInfo(genparticles,-1);
    }
    else if (items.size() == 2) {
        // if items is holding two Z's, then fill the event as a ZZ decay into the histogram
        if ((items[0] == Z) && (items[1] == Z)) {
            //std::cout << "Filled as ZZ" << std::endl;
            type = bbZZ;
        }

        // if items is holding two Hadrons, then fill the event as an HH decay
        if ((items[0] == HADRON) && (items[1] == HADRON)) {
            //std::cout << "Filled as Hadron-Hadron" << std::endl;
            type = HAD;
        }

// If the histogram has not been filled for this event at this point, then we need to do some extra work to fill it appropriately:

        // some variables to keep track of what values are stored in the vector item
        int electron_num = 0;
        int muon_num = 0;
        int tauH_num = 0;
        int tauE_num = 0;
        int tauMU_num = 0;
        int hadron_num = 0;

        // for each value in items:
        for (int k = 0 ; k < 2 ; k++) {
            // if the value is a Hadron
            if (items[k] == HADRON) {
                hadron_num++;
            }
            // if the value is an Electron
            if (items[k] == ELECTRON) {
                electron_num++;
            }
            // if the value is a Muon
            if (items[k] == MUON) {
                muon_num++;
            }
            // if the value is a Tau
            if (items[k] == TAU_H) {
                tauH_num++;
            }
            if (items[k] == TAU_EL) {
                tauE_num++;
            }
            if (items[k] == TAU_MUON) {
                tauMU_num++;
            }
        }

        // if there are two kinds of Leptons stored in items, then fill the event as a LL decay
        if (electron_num + muon_num + tauH_num + tauE_num + tauMU_num == 2) {
            //std::cout << "Filled as Lepton-Lepton" << std::endl;
            type = DILEP;
        }

        // if there is one Lepton type and one Hadron in items:
        if ((electron_num + muon_num + tauH_num + tauE_num + tauMU_num == 1) && (hadron_num == 1)) {
            // Fill as an EH decay if the lepton is electron flavor
            if (electron_num == 1) {
                //std::cout << "Filled as Electron-Hadron" << std::endl;
                type = E;
            }
            // Fill as an MH decay if the lepton is muon flavor
            if (muon_num == 1) {
                //std::cout << "Filled as Muon-Hadron" << std::endl;
                type = MU;
            }
            // Fill as a THH decay if the lepton is tau flavored
            if (tauH_num == 1) {
                //ParticleInfo::printGenInfo(genparticles,-1);
                //std::cout << "Filled as Tau(Hadron)-Hadron" << std::endl;
                type = TAU_HAD;
            }
            // Fill as TEH decay
            if (tauE_num == 1) {
                //ParticleInfo::printGenInfo(genparticles,-1);
                //std::cout << "Filled as Tau(Electron)-Hadron" << std::endl;
                type = TAU_E;
            }
            // Fill as T Mu H decay
            if (tauMU_num == 1) {
                //ParticleInfo::printGenInfo(genparticles,-1);
                //std::cout << "Filled as Tau(Muon)-Hadron" << std::endl;
                type = TAU_MU;
            }
        }
    }
}

}

