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
const FatJet* FatJetSelHelpers::getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat ){
    bCat = BTagging::CSVSJ_INCL;
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,param.wjj_minPT);
    if(fjIDX < 0 )return 0;
    if(param.wjj_maxLepDR > 0 && minDR > param.wjj_maxLepDR) return 0;
    if(  PhysicsUtilities::selObjsMom(jets[fjIDX]->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size() != 2
            ) return 0;
    bCat = BTagging::getCSVSJCat(jets[fjIDX]->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJet* wjjCand, const MomentumF* lepton,  const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat){
    bCat = BTagging::CSVSJ_INCL;

    //assuming the jets collection is ordered by pT
    const FatJet * fj = 0;
    if(jets.size() == 0){
        return 0;
    } else if(jets.size() == 1){
        if(wjjCand == jets[0]) return 0;
        fj = jets[0];
    } else {
        fj = (wjjCand == jets[0] ? jets[1] : jets[0]);
    }
    if(param.hbb_minPT >= 0 && fj->pt() < param.hbb_minPT) return 0;
    if(param.hbb_minLepDPhi > 0 && PhysicsUtilities::absDeltaPhi(*lepton,*fj) < param.hbb_minLepDPhi) return 0;
    if(  PhysicsUtilities::selObjsMom(fj->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size() != 2
            ) return 0;
    bCat = BTagging::getCSVSJCat(fj->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return fj;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passHbbSelection(const FatJet* fj, const BTagging::CSVSJ_CAT csvCat, const FatJetParameters& param){
    if(param.hbb_maxT2oT1 > 0 && fj->tau2otau1() >= param.hbb_maxT2oT1 ) return false;
    if(param.hbb_min_CSVSJCat > BTagging::CSVSJ_INCL && csvCat < param.hbb_min_CSVSJCat) return false;
    if(param.hbb_max_CSVSJCat > BTagging::CSVSJ_INCL && csvCat >= param.hbb_max_CSVSJCat) return false;
    return true;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passWjjSelection(const FatJet* fj,const BTagging::CSVSJ_CAT csvCat, const FatJetParameters& param){
    const float mass = fj->sdMom().mass();
    if(param.wjj_minMass > 0 && mass < param.wjj_minMass  ) return false;
    if(param.wjj_maxMass > 0 && mass >= param.wjj_maxMass ) return false;
    if(param.wjj_maxT2oT1 > 0 && fj->tau2otau1() >= param.wjj_maxT2oT1) return false;
    if(param.wjj_min_CSVSJCat > BTagging::CSVSJ_INCL && csvCat < param.wjj_min_CSVSJCat) return false;
    if(param.wjj_max_CSVSJCat > BTagging::CSVSJ_INCL && csvCat >= param.wjj_max_CSVSJCat) return false;
    return true;
}
//_____________________________________________________________________________
std::vector<const FatJet *> FatJetProcessor::loadFatJets(const FatJetReader& reader_fatjet, const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,param);
    wjjCand = FatJetSelHelpers::getWjjCand(lepton,fjs,param,wjjCSVCat);
    hbbCand = FatJetSelHelpers::getHbbCand(wjjCand,lepton,fjs,param,hbbCSVCat);
    return fjs;
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
//_____________________________________________________________________________
BTagging::CSVSJ_CAT FatJetProcessor::getHbbCSVCat() const {return hbbCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getWjjCSVCat() const {return wjjCSVCat;}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel(const FatJetParameters& param ) const {
    if(wjjCand == 0) return 0;
    return FatJetSelHelpers::passWjjSelection(wjjCand,wjjCSVCat,param);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel(const FatJetParameters& param ) const {
    if(hbbCand == 0) return 0;
    return FatJetSelHelpers::passHbbSelection(hbbCand,hbbCSVCat,param);
}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel() const { return passWjjSel(param);}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel() const { return  passHbbSel(param);}
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
    proc.wjj_maxT2oT1   = 0.55   ;
    proc.wjj_minMass    = 10     ;
    proc.wjj_maxMass    = -1     ;
    proc.wjj_min_CSVSJCat    = BTagging::CSVSJ_INCL;
    proc.wjj_max_CSVSJCat    = BTagging::CSVSJ_MF   ;

    proc.hbb_minLepDPhi = 2.0    ;
    proc.hbb_minPT      = 200    ;
    proc.hbb_maxT2oT1   = -1     ;
    proc.hbb_min_CSVSJCat= BTagging::CSVSJ_MF;
    proc.hbb_max_CSVSJCat= BTagging::CSVSJ_INCL;
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setDefaultFatJetProcessor(FatJetProcessor& proc) {
    setDefaultFatJetProcessor(proc.param);
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setTTBarCRFatJetProcessor(FatJetParameters& proc){
    setDefaultFatJetProcessor(proc);
    proc.hbb_min_CSVSJCat= BTagging::CSVSJ_INCL;
    proc.hbb_max_CSVSJCat= BTagging::CSVSJ_INCL;

    proc.wjj_min_CSVSJCat   = BTagging::CSVSJ_MF   ;
    proc.wjj_max_CSVSJCat   = BTagging::CSVSJ_INCL   ;
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setTTBarCRFatJetProcessor(FatJetProcessor& proc) {
    setTTBarCRFatJetProcessor(proc.param);
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setNonTTBarCRFatJetProcessor(FatJetParameters& proc){
    setDefaultFatJetProcessor(proc);
    proc.hbb_min_CSVSJCat= BTagging::CSVSJ_INCL;
    proc.hbb_max_CSVSJCat= BTagging::CSVSJ_LL  ;
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setNonTTBarCRFatJetProcessor(FatJetProcessor& proc) {
    setNonTTBarCRFatJetProcessor(proc.param);
}

}
