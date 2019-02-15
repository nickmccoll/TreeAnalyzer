#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/FatJet.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {

//_____________________________________________________________________________
std::vector<const FatJet*> FatJetSelHelpers::selectFatJets(const FatJetReader& reader_fatjet, const FatJetParameters& param  ){
    auto outJ = param.fjJetID ?
            PhysicsUtilities::selObjsMom(reader_fatjet.jets,param.cand_minPT,param.cand_maxETA < 0 ? 999. : param.cand_maxETA, //bounds
                    [&](const FatJet* j) {return (j->*param.fjJetID)(); })
                      : PhysicsUtilities::selObjsMom(reader_fatjet.jets,param.cand_minPT,param.cand_maxETA < 0 ? 999. : param.cand_maxETA);
    std::sort(outJ.begin(),outJ.end(),PhysicsUtilities::greaterPTDeref<FatJet>());
    return outJ;
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat){
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,param.wjj_minPT);
    if(fjIDX < 0 )return 0;
    if(jets[fjIDX]->pt() < param.wjj_minPT) return 0;
    if(param.wjj_maxLepDR > 0 && minDR > param.wjj_maxLepDR) return 0;

    nSJs = PhysicsUtilities::selObjsMom(jets[fjIDX]->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size();
    if(nSJs < param.wjj_minSJs) return 0;

    bCat = BTagging::getCSVSJCat(jets[fjIDX]->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJet* wjjCand, const MomentumF* lepton,  const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat){
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    //assuming the jets collection is ordered by pT
    const FatJet * fj = 0;
    for(unsigned int iJ = 0; iJ < 2 && iJ < jets.size(); ++iJ){
        if(jets[iJ] ->pt() < param.hbb_minPT) break;
        if(wjjCand &&  PhysicsUtilities::deltaR2(*wjjCand,*jets[iJ]  ) < 2.56) continue; //1.6^2
        if(param.hbb_minLepDPhi > 0 && PhysicsUtilities::absDeltaPhi(*lepton,*jets[iJ]) < param.hbb_minLepDPhi) continue;
        fj = jets[iJ];
        break;
    }
    if(fj==0) return 0;
    nSJs = PhysicsUtilities::selObjsMom(fj->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size();
    if(nSJs < param.hbb_minSJs) return 0;
    bCat = BTagging::getCSVSJCat(fj->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return fj;
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getDilepHbbCand(const MomentumF* lep1, const MomentumF* lep2, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat) {
    bCat = BTagging::CSVSJ_INCL;
    int nSJs = 0;
    const FatJet* selFJ = 0;
    if (!(lep1 && lep2)) return 0;
    const MomentumF dilepmom = lep1->p4() + lep2->p4();
    // assuming the fatjet collection is ordered by pt
    for (unsigned int k = 0; k < 2 && k < jets.size(); ++k) {
    	if (jets[k]->pt() < param.hbbLL_minPT) break;
    	if (param.hbbLL_minDphiBBLL > 0 && PhysicsUtilities::absDeltaPhi(*jets[k],dilepmom) < param.hbbLL_minDphiBBLL) continue;
    	if (param.hbbLL_minDRbbLL > 0 && PhysicsUtilities::deltaR2(*jets[k],dilepmom) < param.hbbLL_minDRbbLL*param.hbbLL_minDRbbLL) continue;
    	if (PhysicsUtilities::deltaR2(*jets[k],*lep1) < 0.8*0.8 || PhysicsUtilities::deltaR2(*jets[k],*lep2) < 0.8*0.8) continue;
    	selFJ = jets[k];
    	break;
    }
    if (selFJ==0) return 0;
    nSJs = PhysicsUtilities::selObjsMom(selFJ->subJets(), param.sj_minPT, param.sj_maxETA < 0 ? 999 : param.sj_maxETA).size();
    if (nSJs < param.hbbLL_minSJs) return 0;
    bCat = BTagging::getCSVSJCat(selFJ->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return selFJ;
}
//_____________________________________________________________________________
void FatJetProcessor::loadFatJets(const FatJetReader& reader_fatjet,const FatJetReader& reader_fatjet_noLep, const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,param);
    auto fjs_noLep = FatJetSelHelpers::selectFatJets(reader_fatjet_noLep,param);
    wjjCand = FatJetSelHelpers::getWjjCand(lepton,fjs_noLep,param,wjjCSVCat);
    hbbCand = FatJetSelHelpers::getHbbCand(wjjCand,lepton,fjs,param,hbbCSVCat);
}
//_____________________________________________________________________________
void FatJetProcessor::loadDilepFatJet(const FatJetReader& reader_fatjet, const MomentumF* lep1, const MomentumF* lep2) {
	auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,param);
	dilepHbbCand = FatJetSelHelpers::getDilepHbbCand(lep1,lep2,fjs,param,dilepHbbCSVCat);
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
const FatJet * FatJetProcessor::getDilepHbbCand() const {return dilepHbbCand;}
//_____________________________________________________________________________
BTagging::CSVSJ_CAT FatJetProcessor::getHbbCSVCat() const {return hbbCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getWjjCSVCat() const {return wjjCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getDilepHbbCSVCat() const {return dilepHbbCSVCat;}
//_____________________________________________________________________________
void DefaultFatJetSelections::setDefaultFatJetProcessor(FatJetParameters& proc) {
    proc.cand_minPT     = 50                   ;
    proc.cand_maxETA    = 2.4                  ;
    proc.fjJetID        = &FatJet::passTightID ;

    proc.sj_minPT       = 20 ;
    proc.sj_maxETA      = 2.4;
    proc.sj_minBTagPT    = 30 ;
    proc.sj_maxBTagETA  = 2.4;

    proc.wjj_maxLepDR   = 1.2    ;
    proc.wjj_minPT      = 50     ;
    proc.wjj_minSJs     = 2      ;

    proc.hbb_minLepDPhi = 2.0    ;
    proc.hbb_minPT      = 200    ;
    proc.hbb_minSJs     = 2      ;

    proc.hbbLL_minDphiBBLL = 2.0;
    proc.hbbLL_minPT       = 200;
    proc.hbbLL_minSJs      = 2;
    proc.hbbLL_minDRbbLL   = 2.0;
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setDefaultFatJetProcessor(FatJetProcessor& proc) {
    setDefaultFatJetProcessor(proc.param);
}
}
