#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_

#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "CorrectionHelper.h"
#include "Configuration/interface/ReaderConstants.h"

namespace TAna {
class EventReader;
class SMDecayEvent;
namespace EventWeights {
    //function to calculate weight normalized to some lumi
    float calcNormalizedEventWeight(const EventReader& reader_event, const float cross,
            const float numE, const float lumi);
    //Returns event weight normalized to some lumi
    //If this is a data event, set cross section to -1
    //If it is not, will calculate it on the fly
    float getNormalizedEventWeight(const EventReader& reader_event, const float cross,
            const float numE, const EventParameters& evtParam, const float genMtt=-1,
            const int nLepsTT=-1);
}

class PUScaleFactors {
public:
    PUScaleFactors(const std::string& dataDir);
    void setParameters(const EventParameters& evtParam, bool verbose=false);
    float getCorrection(const unsigned int trueNumInteractions,
            const CorrHelp::CORRTYPE corrType = CorrHelp::NOMINAL) const;

private:
    std::unique_ptr<TObjectHelper::Hist1DContainer> nominalSF;
    std::unique_ptr<TObjectHelper::Hist1DContainer> upSF;
    std::unique_ptr<TObjectHelper::Hist1DContainer> downSF;
    std::string dataDirectory;

};
class TopPTWeighting {
public:
    TopPTWeighting(const std::string& dataDir,
            const std::string& sfFile = "corrections/topPTWeight.root", bool verbose = false);
    float getCorrection(const ASTypes::size8 process,const SMDecayEvent& decayEvent) const;
    float getCorrectionNoNorm(const ASTypes::size8 process,const SMDecayEvent& decayEvent) const;
    float getAvgPT(const SMDecayEvent& decayEvent) const;

private:
    std::unique_ptr<TObjectHelper::Hist1DContainer> weightConsts;
    float a;
    float b;
    float nf;
};
}








#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
