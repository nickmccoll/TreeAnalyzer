#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_BTAGGINGSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_BTAGGINGSCALEFACTORS_H_
#include <string>
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Processors/Variables/interface/BTagging.h"
#include "CorrectionHelper.h"

class BTagCalibrationReader;
class BTagCalibration;
namespace TAna {

class Jet;

class BTagScaleFactors {
public:
    BTagScaleFactors(const std::string& dataDir, const std::string& sfFile = "corrections/CSVv2_Moriond17_B_H.csv",
            const std::string& effFile = "corrections/ak4_csvEff.root", bool verbose=false);
    virtual ~BTagScaleFactors();
    float getJetCorr(const Jet* jet,  CorrHelp::CORRTYPE lightT = CorrHelp::NOMINAL,  CorrHelp::CORRTYPE heavyT = CorrHelp::NOMINAL ) const;
    float getJetEff(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp) const;
    float getJetSF(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp, CorrHelp::CORRTYPE corrT) const;

    float getSF(const std::vector<const Jet*>& jets, CorrHelp::CORRTYPE lightT = CorrHelp::NOMINAL,  CorrHelp::CORRTYPE heavyT = CorrHelp::NOMINAL) const;


protected:
std::unique_ptr<BTagCalibration> calib;
std::vector< std::unique_ptr<BTagCalibrationReader> > calibReaders;
std::vector< std::vector< std::unique_ptr<TObjectHelper::Hist2DContainer> >> efficiencies;
const std::vector<std::string> systNames = {"none","down","central","up"};

};
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
