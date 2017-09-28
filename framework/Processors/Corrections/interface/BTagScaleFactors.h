#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_BTAGINGSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_BTAGINGSCALEFACTORS_H_
#include <string>
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Processors/Variables/interface/BTagging.h"
#include "CorrectionHelper.h"

class BTagCalibrationReader;
class BTagCalibration;
namespace TAna {
class BaseRecoJet;
class Jet;
class FatJet;

class BTagScaleFactors {
public:
    BTagScaleFactors(BTagging::CSVWP maxWP,std::string ht,std::string lt, const std::string& dataDir, const std::string& sfFile,
            const std::string& effFile, bool verbose=false);
    virtual ~BTagScaleFactors();
    virtual float getJetEff(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp) const;
    virtual float getJetSF(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp, CorrHelp::CORRTYPE corrT) const;
    float getJetCorr(const BaseRecoJet* jet,  CorrHelp::CORRTYPE lightT = CorrHelp::NOMINAL,  CorrHelp::CORRTYPE heavyT = CorrHelp::NOMINAL ) const;


protected:
    virtual void assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv,   const CorrHelp::CORRTYPE corrT,
            float& lE,float& hE,float& lSF,float& hSF) const =0;


std::unique_ptr<BTagCalibration> calib;
std::vector< std::unique_ptr<BTagCalibrationReader> > calibReaders;
std::vector< std::vector< std::unique_ptr<TObjectHelper::Hist2DContainer> >> efficiencies;
const std::vector<std::string> systNames = {"none","down","central","up"};

};


class JetBTagScaleFactors : public BTagScaleFactors {
public:
    JetBTagScaleFactors(const std::string& dataDir, const std::string& sfFile = "corrections/CSVv2_Moriond17_B_H.csv",
            const std::string& effFile = "corrections/ak4_csvEff.root", bool verbose=false);
    virtual ~JetBTagScaleFactors();

    float getSF(const std::vector<const Jet*>& jets, CorrHelp::CORRTYPE lightT = CorrHelp::NOMINAL,  CorrHelp::CORRTYPE heavyT = CorrHelp::NOMINAL) const;

protected:
    void assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv, const CorrHelp::CORRTYPE corrT,
                float& lE,float& hE,float& lSF,float& hSF) const override;
};

class SubJetBTagScaleFactors : public BTagScaleFactors {
public:
    SubJetBTagScaleFactors(const std::string& dataDir, const std::string& sfFile = "corrections/subjet_CSVv2_Moriond17_B_H.csv",
            const std::string& effFile = "corrections/sj_csvEff.root", bool verbose=false);
    virtual ~SubJetBTagScaleFactors();

    float getSF(const std::vector<const FatJet*>& fatJets, CorrHelp::CORRTYPE lightT = CorrHelp::NOMINAL,  CorrHelp::CORRTYPE heavyT = CorrHelp::NOMINAL) const;

protected:
    void assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv, const CorrHelp::CORRTYPE corrT,
                float& lE,float& hE,float& lSF,float& hSF) const override;
};


}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
