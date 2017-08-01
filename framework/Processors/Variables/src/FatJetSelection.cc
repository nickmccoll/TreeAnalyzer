#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "DataFormats/interface/FatJet.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {

//_____________________________________________________________________________
std::vector<const FatJet*> FatJetSelHelpers::selectFatJets(const FatJetReader& reader_fatjet, float minPT, float maxETA, fjFunBool jetID ){
    const auto& inJ = reader_fatjet.jets;
    std::vector<const FatJet*> outJ;
    outJ.reserve(inJ.size());
    for(const auto& j : inJ){
        if(jetID && (j.*jetID)()  == false) continue;
        if(minPT > 0 && j.pt() < minPT ) continue;
        if(maxETA > 0 && j.absEta() >= maxETA ) continue;
        outJ.push_back(&j);
    }
    std::sort(outJ.begin(),outJ.end(),PhysicsUtilities::greaterPTDeref<FatJet>());
    return outJ;
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getWjjCand(const MomentumF* lepton, const std::vector<const FatJet*>& jets, float minPT, float maxLepDR){
    double minDR = 100000;
    int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,jets,minDR,minPT);
    if(fjIDX < 0 )return 0;
    if(maxLepDR > 0 && minDR > maxLepDR) return 0;
    return jets[fjIDX];
}
//_____________________________________________________________________________
const FatJet* FatJetSelHelpers::getHbbCand(const FatJet* wjjCand, const MomentumF* lepton,  const std::vector<const FatJet*>& jets, float minPT, float hbb_minLepDPhi){
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
    if(minPT >= 0 && fj->pt() < minPT) return 0;
    if(hbb_minLepDPhi > 0 && PhysicsUtilities::absDeltaPhi(*lepton,*fj) < hbb_minLepDPhi) return 0;
    return fj;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passHbbSelection(const FatJet* fj, double maxTau2oTau1,  BTagging::CSVWP firMinCSVWP, BTagging::CSVWP secMinCSVWP, double minMass, double maxMass){
    if(minMass > 0 && fj->sdMom().mass() < minMass  ) return false;
    if(maxMass > 0 && fj->sdMom().mass() >= maxMass ) return false;
    if(maxTau2oTau1 > 0 && fj->tau2otau1() >= 0.55 ) return false;
    if(firMinCSVWP > BTagging::CSV_INCL && fj->maxSJCSV() < BTagging::CSVWP_VALS[firMinCSVWP] ) return false;
    if(secMinCSVWP > BTagging::CSV_INCL && fj->minSJCSV() < BTagging::CSVWP_VALS[secMinCSVWP] ) return false;
    return true;
}
//_____________________________________________________________________________
bool FatJetSelHelpers::passWjjSelection(const FatJet* fj, double maxTau2oTau1,  BTagging::CSVWP maxCSVWP, float minMass, float maxMass){
    if(minMass > 0 && fj->sdMom().mass() < minMass  ) return false;
    if(maxMass > 0 && fj->sdMom().mass() >= maxMass ) return false;
    if(maxTau2oTau1 > 0 && fj->tau2otau1() >= 0.55 ) return false;
    if(maxCSVWP > BTagging::CSV_INCL && fj->maxSJCSV() >= BTagging::CSVWP_VALS[maxCSVWP] ) return false;
    return true;
}
//_____________________________________________________________________________
std::vector<const FatJet *> FatJetProcessor::loadFatJets(const FatJetReader& reader_fatjet, const MomentumF* lepton) {
    auto fjs = FatJetSelHelpers::selectFatJets(reader_fatjet,cand_minPT, cand_maxETA, fjJetID);
    wjjCand = FatJetSelHelpers::getWjjCand(lepton,fjs,wjj_minPT,wjj_maxLepDR);
    hbbCand = FatJetSelHelpers::getHbbCand(wjjCand,lepton,fjs,hbb_minPT,hbb_minLepDPhi);
    return fjs;
}
//_____________________________________________________________________________
const FatJet * FatJetProcessor::getHBBCand() const {return hbbCand;}
const FatJet * FatJetProcessor::getWjjCand() const {return wjjCand;}
//_____________________________________________________________________________
bool FatJetProcessor::passWjjSel() const {
    if(wjjCand == 0) return 0;
    return FatJetSelHelpers::passWjjSelection(wjjCand,wjj_maxT2oT1,wjj_maxCSVWP,wjj_minMass,wjj_maxMass);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSel() const {
    if(hbbCand == 0) return 0;
    return FatJetSelHelpers::passHbbSelection(hbbCand,hbb_maxT2oT1,hbb_l_firMinCSVWP,hbb_l_secMinCSVWP,hbb_minMass,hbb_maxMass);
}
//_____________________________________________________________________________
bool FatJetProcessor::passHbbSelTightBTag() const{
    return passHbbSel() &&
            FatJetSelHelpers::passHbbSelection(hbbCand,hbb_maxT2oT1,hbb_t_firMinCSVWP,hbb_t_secMinCSVWP,hbb_minMass,hbb_maxMass);
}
//_____________________________________________________________________________
void DefaultFatJetSelections::setDefaultFatJetProcessor(FatJetProcessor& proc) {
    proc.cand_minPT     = 50                   ;
    proc.cand_maxETA    = 2.4                  ;
    proc.fjJetID        = &FatJet::passTightID ;

    proc.wjj_maxLepDR   = 1.2    ;
    proc.wjj_minPT      = 50     ;
    proc.wjj_maxT2oT1   = 0.55   ;
    proc.wjj_minMass    = 10     ;
    proc.wjj_maxMass    = -1     ;
    proc.wjj_maxCSVWP   = BTagging::CSV_M  ;

    proc.hbb_minLepDPhi = 2.0    ;
    proc.hbb_minPT      = 50     ;
    proc.hbb_maxT2oT1   = 0.55   ;
    proc.hbb_minMass    = 10     ;
    proc.hbb_maxMass    = -1     ;
    proc.hbb_l_firMinCSVWP= BTagging::CSV_M;
    proc.hbb_l_secMinCSVWP= BTagging::CSV_L;
    proc.hbb_t_firMinCSVWP= BTagging::CSV_M;
    proc.hbb_t_secMinCSVWP= BTagging::CSV_M;
}
}
