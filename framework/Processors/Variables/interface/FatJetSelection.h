
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


struct FatJetParameters;

namespace FatJetSelHelpers {
    typedef  bool (FatJet::*fjFunBool)() const;
    std::vector<const FatJet*> selectFatJets(const FatJetReader& reader_fatjet, const FatJetParameters& param );
    const FatJet* getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat);
    const FatJet* getHbbCand(const FatJet* wjjCand, const MomentumF* lepton, const std::vector<const FatJet*>& jets,const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat);
}


struct FatJetParameters{
    //parameters to select jets om the loadFatJets func
    float cand_minPT     =-1;
    float cand_maxETA    =-1;
    FatJetSelHelpers::fjFunBool fjJetID =0 ;

    //parameters applied to both wjj and hbb
    float sj_minPT       =-1; //cand sel
    float sj_maxETA      =-1; //cand sel
    float sj_minBTagPT       =-1;  //For counting btags in pass sel
    float sj_maxBTagETA      =-1;  //For counting btags in pass sel

    //parameters applied to wjj candidate sel
    float wjj_maxLepDR   =-1;
    float wjj_minPT      =-1;
    float wjj_minSJs     =-1;

    //parameters applied to hbb candidate sel
    float hbb_minLepDPhi =-1;
    float hbb_minPT      =-1;
    float hbb_minSJs     =-1;


};

class FatJetProcessor {
public:

    //uses built in FatJetParameters
    void loadFatJets( const FatJetReader& reader_fatjet,const FatJetReader& reader_fatjet_noLep, const MomentumF* lepton);

    const FatJet * getHBBCand() const;
    const FatJet * getWjjCand() const;
    BTagging::CSVSJ_CAT getHbbCSVCat() const;
    BTagging::CSVSJ_CAT getWjjCSVCat() const;

    FatJetParameters param;


private:
    const FatJet* hbbCand = 0;
    const FatJet* wjjCand = 0;
    BTagging::CSVSJ_CAT hbbCSVCat = BTagging::CSVSJ_INCL;
    BTagging::CSVSJ_CAT wjjCSVCat = BTagging::CSVSJ_INCL;

};

namespace DefaultFatJetSelections {
void setDefaultFatJetProcessor(FatJetParameters& proc);
void setDefaultFatJetProcessor(FatJetProcessor& proc);
}


}


#endif

