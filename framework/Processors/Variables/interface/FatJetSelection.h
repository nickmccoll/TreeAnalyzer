
#ifndef PROCESSORS_VARIABLES_FATJETSELECTION_H
#define PROCESSORS_VARIABLES_FATJETSELECTION_H


#include <vector>
#include <utility>

#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Configuration/interface/ReaderConstants.h"



namespace TAna {
template <class CoordSystem> class Momentum;
typedef Momentum<ASTypes::CylLorentzCoordF> MomentumF;
class FatJet;
class FatJetReader;

namespace FatJetSelHelpers {
    std::vector<const FatJet*> selectFatJets(const FatJetParameters& params,
            const FatJetReader& reader_fatjet);
    const FatJet* getWjjCand(const FatJetParameters& params, const MomentumF* lepton,
            const std::vector<const FatJet*>& jets);
    const FatJet* getHbbCand(const FatJetParameters& params, const FatJet* wjjCand,
            const MomentumF* lepton, const std::vector<const FatJet*>& jets);
    const FatJet* getDilepHbbCand(const FatJetParameters& params, const MomentumF* lep1,
            const MomentumF* lep2, const std::vector<const FatJet*>& fatjets);
}

class FatJetProcessor {
public:

    //uses built in FatJetParameters
    void loadFatJets(const FatJetParameters& params,
            const FatJetReader& reader_fatjet, const FatJetReader& reader_fatjet_noLep,
            const MomentumF* lepton);
    void loadDilepFatJet (const FatJetParameters& params, const FatJetReader& reader_fatjet,
            const MomentumF* lep1,const MomentumF* lep2);

    const FatJet * getHBBCand() const;
    const FatJet * getWjjCand() const;
    const FatJet * getDilepHbbCand() const;



private:
    const FatJet* hbbCand = 0;
    const FatJet* wjjCand = 0;
    const FatJet* dilepHbbCand = 0;
};
}


#endif

