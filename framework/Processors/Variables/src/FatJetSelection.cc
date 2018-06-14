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
const FatJet* FatJetSelHelpers::getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat, int& nSJs ){
    bCat = BTagging::CSVSJ_INCL;
    nSJs = 0;
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,param.wjj_minPT);
    if(fjIDX < 0 )return 0;
    if(param.wjj_maxLepDR > 0 && minDR > param.wjj_maxLepDR) return 0;
    nSJs = PhysicsUtilities::selObjsMom(jets[fjIDX]->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size();
    bCat = BTagging::getCSVSJCat(jets[fjIDX]->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJet* wjjCand, const MomentumF* lepton,  const std::vector<const FatJet*>& jets, const FatJetParameters& param, BTagging::CSVSJ_CAT& bCat, int& nSJs){
    bCat = BTagging::CSVSJ_INCL;
    nSJs = 0;
    //assuming the jets collection is ordered by pT
    const FatJet * fj = 0;
    for(unsigned int iJ = 0; iJ < 2 && iJ < jets.size(); ++iJ){
        if(param.hbb_minPT >= 0 && jets[iJ] ->pt() < param.hbb_minPT) break;
        if(wjjCand &&  PhysicsUtilities::deltaR2(*wjjCand,*jets[iJ]  ) < 2.56) continue; //1.6^2
        if(param.hbb_minLepDPhi > 0 && PhysicsUtilities::absDeltaPhi(*lepton,*jets[iJ]) < param.hbb_minLepDPhi) continue;
        fj = jets[iJ];
        break;
    }
    if(fj==0) return 0;
    nSJs = PhysicsUtilities::selObjsMom(fj->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size();
    bCat = BTagging::getCSVSJCat(fj->subJets(), param.sj_minBTagPT, param.sj_maxBTagETA);
    return fj;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passHbbSelection(const FatJet* fj, const BTagging::CSVSJ_CAT csvCat,const int nSJs, const FatJetParameters& param){
    if(param.hbb_minSJs > 0 && nSJs < param.hbb_minSJs) return false;
    if(param.hbb_maxT2oT1 > 0 && fj->tau2otau1() >= param.hbb_maxT2oT1 ) return false;
    if(param.hbb_min_CSVSJCat > BTagging::CSVSJ_INCL && csvCat < param.hbb_min_CSVSJCat) return false;
    if(param.hbb_max_CSVSJCat > BTagging::CSVSJ_INCL && csvCat >= param.hbb_max_CSVSJCat) return false;
    return true;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passWjjSelection(const FatJet* fj,const BTagging::CSVSJ_CAT csvCat,const int nSJs,  const FatJetParameters& param){
    const float mass = fj->sdMom().mass();
    if(param.wjj_minSJs > 0 && nSJs < param.wjj_minSJs) return false;
    if(param.wjj_minMass > 0 && mass < param.wjj_minMass  ) return false;
    if(param.wjj_maxMass > 0 && mass >= param.wjj_maxMass ) return false;
    if(param.wjj_maxT2oT1 > 0 && fj->tau2otau1() >= param.wjj_maxT2oT1) return false;
    if(param.wjj_min_CSVSJCat > BTagging::CSVSJ_INCL && csvCat < param.wjj_min_CSVSJCat) return false;
    if(param.wjj_max_CSVSJCat > BTagging::CSVSJ_INCL && csvCat >= param.wjj_max_CSVSJCat) return false;
    return true;
}
//_____________________________________________________________________________
void FatJetProcessor::loadFatJets(const FatJetReader& reader_fatjet,const FatJetReader& reader_fatjet_noLep, const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,param);
    auto fjs_noLep = FatJetSelHelpers::selectFatJets(reader_fatjet_noLep,param);
    wjjCand = FatJetSelHelpers::getWjjCand(lepton,fjs_noLep,param,wjjCSVCat,wjjNSJs);
    hbbCand = FatJetSelHelpers::getHbbCand(wjjCand,lepton,fjs,param,hbbCSVCat,hbbNSJs);
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
//_____________________________________________________________________________
BTagging::CSVSJ_CAT FatJetProcessor::getHbbCSVCat() const {return hbbCSVCat;}
BTagging::CSVSJ_CAT FatJetProcessor::getWjjCSVCat() const {return wjjCSVCat;}
//_____________________________________________________________________________
int FatJetProcessor::getHbbNSJs() const {return hbbNSJs;}
int FatJetProcessor::getWjjNSJs() const {return wjjNSJs;}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel(const FatJetParameters& param ) const {
    if(wjjCand == 0) return 0;
    return FatJetSelHelpers::passWjjSelection(wjjCand,wjjCSVCat,wjjNSJs,param);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel(const FatJetParameters& param ) const {
    if(hbbCand == 0) return 0;
    return FatJetSelHelpers::passHbbSelection(hbbCand,hbbCSVCat,hbbNSJs,param);
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
    proc.wjj_minSJs     = 2      ;
    proc.wjj_min_CSVSJCat    = BTagging::CSVSJ_INCL;
    proc.wjj_max_CSVSJCat    = BTagging::CSVSJ_INCL   ;

    proc.hbb_minLepDPhi = 2.0    ;
    proc.hbb_minPT      = 200    ;
    proc.hbb_minSJs     = 2      ;
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
