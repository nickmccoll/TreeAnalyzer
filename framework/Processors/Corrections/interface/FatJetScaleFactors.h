#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_FATJETSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_FATJETSCALEFACTORS_H_
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
namespace TAna {

template<class T> class Momentum;
typedef Momentum<ASTypes::CylLorentzCoordF> MomentumF;
class FatJet;

class HbbFatJetScaleFactors {
public:
    HbbFatJetScaleFactors(const std::string& dataDir, const std::string& sfFile = "corrections/hbbFatJetScaleFactors.root", bool verbose = false);
    float getMassScaleFactor(const float corrJetPT, const float jetAbsETA) const;
    float getJEC(const FatJet* fj) const;
    float getCorrSDMass(const FatJet* fj) const;
private:


    std::unique_ptr<TObjectHelper::GraphEContainer> sdMassCorr_eta_0p0_0p8;
    std::unique_ptr<TObjectHelper::GraphEContainer> sdMassCorr_eta_0p8_1p6;
    std::unique_ptr<TObjectHelper::GraphEContainer> sdMassCorr_eta_1p6_2p4;

    std::unique_ptr<TObjectHelper::TF1Container   > jec_eta_0p0_0p6    ;
    std::unique_ptr<TObjectHelper::TF1Container   > jec_eta_0p6_1p4    ;
    std::unique_ptr<TObjectHelper::TF1Container   > jec_eta_1p4_1p8    ;
    std::unique_ptr<TObjectHelper::TF1Container   > jec_eta_1p8_2p0    ;
    std::unique_ptr<TObjectHelper::TF1Container   > jec_eta_2p0_2p4    ;

};

class SoftDropMassScaleFactorsFromW {
public:
    SoftDropMassScaleFactorsFromW(const std::string& dataDir, const std::string& sfFile = "corrections/std_sdMass_puppiCorr.root", bool verbose = false);
    float getScaleFactor(const MomentumF* mom) const;
    float getCorrSDMass(const FatJet* fj) const;

private:

    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrGEN     ;
    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrRECO_cen;
    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrRECO_for;


};
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
