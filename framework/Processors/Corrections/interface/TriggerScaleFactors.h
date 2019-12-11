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
    float getElectronTriggerSF_1l(const float ht) const;
    float getMuonTriggerSF_1l(const float ht) const;
    float getSingleLeptonTriggerSF(const float ht, const bool leadingLepIsMuon) const;

    float getDileptonTriggerSF(const float ht, const float pt2, const bool lep1IsMuon, const bool lep2IsMuon) const;
    float getElectronTriggerEff_2l(const float val, const bool getHT, const bool isData) const;
    float getMuonTriggerEff_2l(const float value, const bool getHT, const bool isData) const;


private:
    std::unique_ptr<TObjectHelper::Hist1DContainer> electronSFs_1l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonSFs_1l;

    std::unique_ptr<TObjectHelper::Hist1DContainer> electronEffs_mcHT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonEffs_mcHT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> electronEffs_mcPT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonEffs_mcPT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> electronEffs_dataHT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonEffs_dataHT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> electronEffs_dataPT_2l;
    std::unique_ptr<TObjectHelper::Hist1DContainer> muonEffs_dataPT_2l;

    std::string dataDirectory;
};
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
