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
const FatJet* FatJetSelHelpers::getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, const FatJetParameters& param ){
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,param.wjj_minPT);
    if(fjIDX < 0 )return 0;
    if(param.wjj_maxLepDR > 0 && minDR > param.wjj_maxLepDR) return 0;
    if(  PhysicsUtilities::selObjsMom(jets[fjIDX]->subJets(),
            param.sj_minPT, param.sj_maxETA < 0 ? 999.0 : param.sj_maxETA).size() != 2
            ) return 0;
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJet* wjjCand, const MomentumF* lepton,  const std::vector<const FatJet*>& jets, const FatJetParameters& param){
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
    return fj;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passHbbSelection(const FatJet* fj, const FatJetParameters& param, const bool tight){
    if(param.hbb_minMass > 0 && fj->sdMom().mass() < param.hbb_minMass  ) return false;
    if(param.hbb_maxMass > 0 && fj->sdMom().mass() >= param.hbb_maxMass ) return false;
    if(param.hbb_maxT2oT1 > 0 && fj->tau2otau1() >= param.hbb_maxT2oT1 ) return false;

    auto btagSJs = PhysicsUtilities::selObjsMom(fj->subJets(),
            param.sj_minBTagPT, param.sj_maxBTagETA < 0 ? 999.0 : param.sj_maxBTagETA);

    std::sort(btagSJs.begin(),btagSJs.end(), [](const SubJet* a,const SubJet* b) {return a->csv() > b->csv();} );

    auto passBSelBTag = [&](const BTagging::CSVWP firMinCSVWP, const BTagging::CSVWP secMinCSVWP  ) -> bool{
        if(firMinCSVWP > BTagging::CSV_INCL)
            if(btagSJs.size() == 0 || btagSJs[0]->csv() < BTagging::CSVWP_VALS[firMinCSVWP] ) return false;
        if(secMinCSVWP > BTagging::CSV_INCL)
            if(btagSJs.size() < 2 || btagSJs[1]->csv() < BTagging::CSVWP_VALS[secMinCSVWP] ) return false;
        return true;
    };
    return tight ? passBSelBTag(param.hbb_t_firMinCSVWP,param.hbb_t_secMinCSVWP) : passBSelBTag(param.hbb_l_firMinCSVWP,param.hbb_l_secMinCSVWP);
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passWjjSelection(const FatJet* fj, const FatJetParameters& param){
    if(param.wjj_minMass > 0 && fj->sdMom().mass() < param.wjj_minMass  ) return false;
    if(param.wjj_maxMass > 0 && fj->sdMom().mass() >= param.wjj_maxMass ) return false;
    if(param.wjj_maxT2oT1 > 0 && fj->tau2otau1() >= param.wjj_maxT2oT1) return false;

    auto btagSJs = PhysicsUtilities::selObjsMom(fj->subJets(),
            param.sj_minBTagPT, param.sj_maxBTagETA < 0 ? 999.0 : param.sj_maxBTagETA);
    std::sort(btagSJs.begin(),btagSJs.end(), [](const SubJet* a,const SubJet* b) {return a->csv() > b->csv();} );

    if(param.wjj_maxCSVWP > BTagging::CSV_INCL && btagSJs.size())
        if(btagSJs.front()->csv() >= BTagging::CSVWP_VALS[param.wjj_maxCSVWP]) return false;
    return true;
}
//_____________________________________________________________________________
std::vector<const FatJet *> FatJetProcessor::loadFatJets(const FatJetReader& reader_fatjet, const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,param);
    wjjCand = FatJetSelHelpers::getWjjCand(lepton,fjs,param);
    hbbCand = FatJetSelHelpers::getHbbCand(wjjCand,lepton,fjs,param);
    return fjs;
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel(const FatJetParameters& param ) const {
    if(wjjCand == 0) return 0;
    return FatJetSelHelpers::passWjjSelection(wjjCand,param);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel(const FatJetParameters& param ) const {
    if(hbbCand == 0) return 0;
    return FatJetSelHelpers::passHbbSelection(hbbCand,param,false);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSelTightBTag(const FatJetParameters& param ) const{
    if(hbbCand == 0) return 0;
    return FatJetSelHelpers::passHbbSelection(hbbCand,param,true);
}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel() const { return passWjjSel(param);}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel() const { return  passHbbSel(param);}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSelTightBTag() const{ return passHbbSelTightBTag(param);}
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
    proc.wjj_maxCSVWP   = BTagging::CSV_M  ;

    proc.hbb_minLepDPhi = 2.0    ;
    proc.hbb_minPT      = 50     ;
    proc.hbb_maxT2oT1   = -1     ;
    proc.hbb_minMass    = 10     ;
    proc.hbb_maxMass    = -1     ;
    proc.hbb_l_firMinCSVWP= BTagging::CSV_M;
    proc.hbb_l_secMinCSVWP= BTagging::CSV_L;
    proc.hbb_t_firMinCSVWP= BTagging::CSV_M;
    proc.hbb_t_secMinCSVWP= BTagging::CSV_M;
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setDefaultFatJetProcessor(FatJetProcessor& proc) {
    setDefaultFatJetProcessor(proc.param);
}
}
