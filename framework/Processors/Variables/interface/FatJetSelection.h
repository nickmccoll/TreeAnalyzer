
#ifndef PROCESSORS_VARIABLES_FATJETSELECTION_H
#define PROCESSORS_VARIABLES_FATJETSELECTION_H


#include <vector>
#include <utility>

#include "Processors/Variables/interface/BTagging.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Configuration/interface/ReaderConstants.h"



namespace TAna {
template <class CoordSystem> class Momentum;
typedef Momentum<ASTypes::CylLorentzCoordF> MomentumF;
class FatJet;
class FatJetReader;

namespace FatJetSelHelpers {
    typedef  bool (FatJet::*fjFunBool)() const;
    std::vector<const FatJet*> selectFatJets(const FatJetReader& reader_fatjet, const FatJetParameters& param );
    const FatJet* getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat);
    const FatJet* getHbbCand(const FatJet* wjjCand, const MomentumF* lepton, const std::vector<const FatJet*>& jets,const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat);
    const FatJet* getDilepHbbCand(const MomentumF* lep1, const MomentumF* lep2, const std::vector<const FatJet*>& fatjets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat);
}

class FatJetProcessor {
public:
    void setParameters(const FatJetParameters& inParams) {param = inParams;}

    //uses built in FatJetParameters
    void loadFatJets( const FatJetReader& reader_fatjet,const FatJetReader& reader_fatjet_noLep, const MomentumF* lepton);
    void loadDilepFatJet (const FatJetReader& reader_fatjet, const MomentumF* lep1, const MomentumF* lep2);

    const FatJet * getHBBCand() const;
    const FatJet * getWjjCand() const;
    const FatJet * getDilepHbbCand() const;
    BTagging::CSVSJ_CAT getHbbCSVCat() const;
    BTagging::CSVSJ_CAT getWjjCSVCat() const;
    BTagging::CSVSJ_CAT getDilepHbbCSVCat() const;

    FatJetParameters param;


private:
    const FatJet* hbbCand = 0;
    const FatJet* wjjCand = 0;
    const FatJet* dilepHbbCand = 0;
    BTagging::CSVSJ_CAT hbbCSVCat = BTagging::CSVSJ_INCL;
    BTagging::CSVSJ_CAT wjjCSVCat = BTagging::CSVSJ_INCL;
    BTagging::CSVSJ_CAT dilepHbbCSVCat = BTagging::CSVSJ_INCL;

};

namespace DefaultFatJetSelections {
void setDefaultFatJetProcessor(FatJetParameters& proc);
void setDefaultFatJetProcessor(FatJetProcessor& proc);
}


}


#endif

