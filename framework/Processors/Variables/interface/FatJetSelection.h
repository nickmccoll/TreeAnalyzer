
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
    bool passHbbSelection(const FatJet* fj, const BTagging::CSVSJ_CAT csvCat, const FatJetParameters& param);
    bool passWjjSelection(const FatJet* fj, const BTagging::CSVSJ_CAT csvCat,const FatJetParameters& param);
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
    //parameters applied to wjj pass sel
    float wjj_maxT2oT1   =-1;
    float wjj_minMass    =-1;
    float wjj_maxMass    =-1;
    BTagging::CSVSJ_CAT wjj_min_CSVSJCat   = BTagging::CSVSJ_INCL; //max means <, min means >=
    BTagging::CSVSJ_CAT wjj_max_CSVSJCat   = BTagging::CSVSJ_INCL; //max means <, min means >=

    //parameters applied to hbb candidate sel
    float hbb_minLepDPhi =-1;
    float hbb_minPT      =-1;
    //parameters applied to hbb pass sel
    float hbb_maxT2oT1   =-1;
    BTagging::CSVSJ_CAT hbb_min_CSVSJCat = BTagging::CSVSJ_INCL; //max means <, min means >=
    BTagging::CSVSJ_CAT hbb_max_CSVSJCat = BTagging::CSVSJ_INCL; //max means <, min means >=

};

class FatJetProcessor {
public:

    //uses built in FatJetParameters
    std::vector<const FatJet *> loadFatJets( const FatJetReader& reader_fatjet, const MomentumF* lepton);

    const FatJet * getHBBCand() const;
    const FatJet * getWjjCand() const;
    BTagging::CSVSJ_CAT getHbbCSVCat() const;
    BTagging::CSVSJ_CAT getWjjCSVCat() const;
    //use built in FatJetParameters
    bool passWjjSel() const;
    bool passHbbSel() const;
    bool passHbbSelTightBTag() const;

    //use custom in FatJetParameters
    bool passWjjSel(const FatJetParameters& param ) const;
    bool passHbbSel(const FatJetParameters& param ) const;



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

void setTTBarCRFatJetProcessor(FatJetParameters& proc);
void setTTBarCRFatJetProcessor(FatJetProcessor& proc);

void setNonTTBarCRFatJetProcessor(FatJetParameters& proc);
void setNonTTBarCRFatJetProcessor(FatJetProcessor& proc);
}


}


#endif

