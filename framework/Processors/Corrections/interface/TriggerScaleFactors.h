#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_TRIGGERSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_TRIGGERSCALEFACTORS_H_
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "Configuration/interface/ReaderConstants.h"

namespace TAna {
class TriggerScaleFactors {
public:
    TriggerScaleFactors(const std::string& dataDir);
    void setParameters(const EventParameters& evtParam, bool verbose=false);
    float getElectronTriggerSF(const float ht) const;
    float getMuonTriggerSF(const float ht) const;
    float getLeptonTriggerSF(const float ht, const bool leadingLepIsMuon) const;

private:
    std::unique_ptr<TObjectHelper::Hist1DContainer> electronSFs;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonSFs;
    std::string dataDirectory;
};
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
