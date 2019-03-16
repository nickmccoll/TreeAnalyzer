#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/FatJet.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {

//_____________________________________________________________________________
std::vector<const FatJet*> FatJetSelHelpers::selectFatJets(const FatJetParameters& params,
        const FatJetReader& reader_fatjet){
    auto outJ = params.fjJetID
                      ? PhysicsUtilities::selObjsMom(reader_fatjet.jets,
                              params.cand_minPT,
                              params.cand_maxETA < 0 ? 999. : params.cand_maxETA, //bounds
                              [&](const FatJet* j) {return (j->*params.fjJetID)(); })
                      : PhysicsUtilities::selObjsMom(reader_fatjet.jets,
                              params.cand_minPT,
                              params.cand_maxETA < 0 ? 999. : params.cand_maxETA);
    std::sort(outJ.begin(),outJ.end(),PhysicsUtilities::greaterPTDeref<FatJet>());
    return outJ;
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getWjjCand(const FatJetParameters& params, const MomentumF* lepton,
        const std::vector<const FatJet*>& jets, BTagging::CSVSJ_CAT& bCat){
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,params.wjj_minPT);
    if(fjIDX < 0 )return 0;
    if(jets[fjIDX]->pt() < params.wjj_minPT) return 0;
    if(params.wjj_maxLepDR > 0 && minDR > params.wjj_maxLepDR) return 0;

    nSJs = PhysicsUtilities::selObjsMom(jets[fjIDX]->subJets(),
            params.sj_minPT, params.sj_maxETA < 0 ? 999.0 : params.sj_maxETA).size();
    if(nSJs < params.wjj_minSJs) return 0;

    bCat = BTagging::getCSVSJCat(jets[fjIDX]->subJets(), params.sj_minBTagPT, params.sj_maxBTagETA);
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJetParameters& params,
        const FatJet* wjjCand, const MomentumF* lepton, const std::vector<const FatJet*>& jets,
        BTagging::CSVSJ_CAT& bCat){
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    //assuming the jets collection is ordered by pT
    const FatJet * fj = 0;
    for(unsigned int iJ = 0; iJ < 2 && iJ < jets.size(); ++iJ){
        if(jets[iJ] ->pt() < params.hbb_minPT) break;
        if(wjjCand &&  PhysicsUtilities::deltaR2(*wjjCand,*jets[iJ]  ) < 2.56) continue; //1.6^2
        if(params.hbb_minLepDPhi > 0
                && PhysicsUtilities::absDeltaPhi(*lepton,*jets[iJ]) < params.hbb_minLepDPhi
                ) continue;
        fj = jets[iJ];
        break;
    }
    if(fj==0) return 0;
    nSJs = PhysicsUtilities::selObjsMom(fj->subJets(),
            params.sj_minPT, params.sj_maxETA < 0 ? 999.0 : params.sj_maxETA).size();
    if(nSJs < params.hbb_minSJs) return 0;
    bCat = BTagging::getCSVSJCat(fj->subJets(), params.sj_minBTagPT, params.sj_maxBTagETA);
    return fj;
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getDilepHbbCand(const FatJetParameters& params,
        const MomentumF* lep1, const MomentumF* lep2, const std::vector<const FatJet*>& jets,
        BTagging::CSVSJ_CAT& bCat) {
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    const FatJet* selFJ = 0;
    if (!(lep1 && lep2)) return 0;
    const MomentumF dilepmom = lep1->p4() + lep2->p4();
    // assuming the fatjet collection is ordered by pt
    for (unsigned int k = 0; k < 2 && k < jets.size(); ++k) {
    	if (jets[k]->pt() < params.hbbLL_minPT) break;
    	if (params.hbbLL_minDphiBBLL > 0
    	        && PhysicsUtilities::absDeltaPhi(*jets[k],dilepmom) < params.hbbLL_minDphiBBLL)
    	    continue;
    	if (params.hbbLL_minDRbbLL > 0
    	        && PhysicsUtilities::deltaR2(*jets[k],dilepmom)
    	        < params.hbbLL_minDRbbLL*params.hbbLL_minDRbbLL
    	        ) continue;
    	if (PhysicsUtilities::deltaR2(*jets[k],*lep1) < 0.8*0.8
    	        || PhysicsUtilities::deltaR2(*jets[k],*lep2) < 0.8*0.8) continue;
    	selFJ = jets[k];
    	break;
    }
    if (selFJ==0) return 0;
    nSJs = PhysicsUtilities::selObjsMom(selFJ->subJets(),
            params.sj_minPT, params.sj_maxETA < 0 ? 999 : params.sj_maxETA
                    ).size();
    if (nSJs < params.hbbLL_minSJs) return 0;
    bCat = BTagging::getCSVSJCat(selFJ->subJets(), params.sj_minBTagPT, params.sj_maxBTagETA);
    return selFJ;
}
//_____________________________________________________________________________
void FatJetProcessor::loadFatJets(const FatJetParameters& params,
        const FatJetReader& reader_fatjet,const FatJetReader& reader_fatjet_noLep,
        const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(params,reader_fatjet);
    auto fjs_noLep = FatJetSelHelpers::selectFatJets(params,reader_fatjet_noLep);
    wjjCand = FatJetSelHelpers::getWjjCand(params,lepton,fjs_noLep,wjjCSVCat);
    hbbCand = FatJetSelHelpers::getHbbCand(params,wjjCand,lepton,fjs,hbbCSVCat);
}
//_____________________________________________________________________________
void FatJetProcessor::loadDilepFatJet(const FatJetParameters& params,
        const FatJetReader& reader_fatjet, const MomentumF* lep1,
        const MomentumF* lep2) {
	auto fjs = FatJetSelHelpers::selectFatJets(params,reader_fatjet);
	dilepHbbCand = FatJetSelHelpers::getDilepHbbCand(params,lep1,lep2,fjs,dilepHbbCSVCat);
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
const FatJet * FatJetProcessor::getDilepHbbCand() const {return dilepHbbCand;}
//_____________________________________________________________________________
BTagging::CSVSJ_CAT FatJetProcessor::getHbbCSVCat() const {return hbbCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getWjjCSVCat() const {return wjjCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getDilepHbbCSVCat() const {return dilepHbbCSVCat;}

}
