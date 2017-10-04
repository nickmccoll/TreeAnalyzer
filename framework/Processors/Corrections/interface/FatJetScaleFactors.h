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

class SoftDropMassScaleFactors {
public:
    SoftDropMassScaleFactors(const std::string& dataDir, const std::string& sfFile = "corrections/std_sdMass_puppiCorr.root", bool verbose = false);
    float getScaleFactor(const MomentumF* mom) const;
    float getCorrSDMass(const FatJet* fj) const;

private:
    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrGEN     ;
    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrRECO_cen;
    std::unique_ptr<TObjectHelper::TF1Container> puppisd_corrRECO_for;

};
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
