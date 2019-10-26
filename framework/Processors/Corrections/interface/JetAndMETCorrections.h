#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_JETANDMETCORRECTIONS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_JETANDMETCORRECTIONS_H_
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"

#include "AnalysisSupport/Utilities/interface/Types.h"
#include "CorrectionHelper.h"
#include "JetResolutionObject.h"
#include "DataFormats/interface/Met.h"
#include "DataFormats/interface/FatJet.h"
#include "Configuration/interface/ReaderConstants.h"

namespace TAna {
class JetReader;
class EventReader;

class JERCorrector {
public:
    JERCorrector (const std::string& dataDir, std::shared_ptr<TRandom3> rndGen, const CorrHelp::CORRTYPE cT = CorrHelp::NOMINAL );
    ~JERCorrector(){}

    void setParameters(const JetParameters& param);

    void setCorrType(const CorrHelp::CORRTYPE newCT) {cT=newCT;}
    int getSFCount(const CorrHelp::CORRTYPE c) ;
    void processFatJets(FatJetCollection& jets,const GenFatJetCollection& genjets, const float rho);
    void processJets(JetReader& jetreader,Met& met,const GenJetCollection& genjets, const float rho);

    float getSF(const CorrHelp::CORRTYPE cT, const float absETA) const ;
private:
    float correctJet(Jet* jet,const std::vector<const GenJet*> genjets, const float jres,
            const float jSF,const float coneR);
    const std::string dataDir;
    CorrHelp::CORRTYPE cT;
    std::unique_ptr<JMEStand::JetResolutionObject> ak8Puppi_resObj;
    std::unique_ptr<JMEStand::JetResolutionObject> ak8Puppi_sfObj ;
    std::unique_ptr<JMEStand::JetResolutionObject> ak4CHS_resObj  ;
    std::unique_ptr<JMEStand::JetResolutionObject> ak4CHS_sfObj   ;
    std::shared_ptr<TRandom3> rndGen;
};


class JESUncShifter { //only does something for down or up
public:
    JESUncShifter (const CorrHelp::CORRTYPE cT = CorrHelp::NONE );
    ~JESUncShifter(){}
    void setCorrType(const CorrHelp::CORRTYPE newCT) {cT=newCT;}
    void processFatJets(FatJetCollection& jets);
    void processJets(JetReader& jetreader,Met& met);

private:
    float correctJet(BaseRecoJet* jet) const;
    CorrHelp::CORRTYPE cT;
};

class METUncShifter{ //only does something for up or down
public:
    METUncShifter (const CorrHelp::CORRTYPE cT = CorrHelp::NONE );
    ~METUncShifter(){}
    void setCorrType(const CorrHelp::CORRTYPE newCT) {cT=newCT;}
    void process(Met& met, const EventReader& eventreader) const;
private:
    CorrHelp::CORRTYPE cT;
};

}
#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
