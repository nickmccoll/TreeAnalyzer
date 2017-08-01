
#ifndef PROCESSORS_VARIABLES_FATJETSELECTION_H
#define PROCESSORS_VARIABLES_FATJETSELECTION_H


#include <vector>
#include <utility>

#include "Processors/Variables/interface/BTagging.h"
#include "AnalysisSupport/Utilities/interface/Types.h"

namespace TAna {
template <class CoordSystem> class Momentum;
typedef Momentum<ASTypes::CylLorentzCoordF> MomentumF;
class FatJet;
class FatJetReader;


namespace FatJetSelHelpers {
    typedef  bool (FatJet::*fjFunBool)() const;
    std::vector<const FatJet*> selectFatJets(const FatJetReader* reader_fatjet, float minPT, float maxETA, fjFunBool jetID = 0 );
    const FatJet* getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, float minPT, float maxLepDR);
    const FatJet* getHbbCand(const FatJet* wjjCand, const MomentumF* lepton, const std::vector<const FatJet*>& jets,float minPT, float hbb_minLepDPhi);
    bool passHbbSelection(const FatJet* fj, double maxTau2oTau1 = -1,BTagging::CSVWP firMinCSVWP = BTagging::CSV_INCL, BTagging::CSVWP secMinCSVWP = BTagging::CSV_INCL, double minMass=-1, double maxMass=-1);
    bool passWjjSelection(const FatJet* fj, double maxTau2oTau1 = -1,BTagging::CSVWP maxCSVWP    = BTagging::CSV_INCL, float minMass = -1, float maxMass = -1 );
}

class FatJetProcessor {
public:

    std::vector<const FatJet *> loadFatJets(const MomentumF* lepton, const FatJetReader* reader_fatjet);

    const FatJet * getHBBCand() const;
    const FatJet * getWjjCand() const;
    bool passWjjSel() const;
    bool passHbbSel() const;
    bool passHbbSelTightBTag() const;

    float cand_minPT     =-1;
    float cand_maxETA    =-1;
    FatJetSelHelpers::fjFunBool fjJetID =0 ;


    float wjj_maxLepDR   =-1;
    float wjj_minPT      =-1;
    float wjj_maxT2oT1   =-1;
    float wjj_minMass    =-1;
    float wjj_maxMass    =-1;
    BTagging::CSVWP wjj_maxCSVWP   = BTagging::CSV_INCL;

    float hbb_minLepDPhi =-1;
    float hbb_minPT      =-1;
    float hbb_maxT2oT1   =-1;
    float hbb_minMass    =-1;
    float hbb_maxMass    =-1;
    BTagging::CSVWP hbb_l_firMinCSVWP= BTagging::CSV_INCL;
    BTagging::CSVWP hbb_l_secMinCSVWP= BTagging::CSV_INCL;
    BTagging::CSVWP hbb_t_firMinCSVWP= BTagging::CSV_INCL;
    BTagging::CSVWP hbb_t_secMinCSVWP= BTagging::CSV_INCL;

private:
    const FatJet* hbbCand = 0;
    const FatJet* wjjCand = 0;

};

namespace DefaultFatJetSelections {
FatJetProcessor  getDefaultFatJetProcessor()     ;
}


}


#endif

