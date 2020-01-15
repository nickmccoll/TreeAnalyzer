#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_LEPTONSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_LEPTONSCALEFACTORS_H_
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "DataFormats/interface/GenParticle.h"
#include "Configuration/interface/ReaderConstants.h"
#include "CorrectionHelper.h"

namespace TAna {
class SMDecayEvent;
class Lepton;
class Electron;
class Muon;
class Jet;

class LeptonScaleFactors {
public:
    LeptonScaleFactors(const std::string& dataDir);
    void setParameters(const LeptonParameters& param);
    void setParameters(const DileptonParameters& param);
    virtual void setParameters(const std::string& el_fileName,const std::string& mu_fileName) =0;


    virtual ~LeptonScaleFactors() {}

    virtual float getElectronSF(
            const CorrHelp::CORRTYPE recoT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE isoT = CorrHelp::NONE
    ) const = 0;
    virtual float getMuonSF(
            const CorrHelp::CORRTYPE recoT = CorrHelp::NONE,
            const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE isoT = CorrHelp::NONE
    ) const = 0;
    virtual float getSF(
            const CorrHelp::CORRTYPE elReco = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE elID = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE elISO = CorrHelp::NONE,
            const CorrHelp::CORRTYPE muReco = CorrHelp::NONE,
            const CorrHelp::CORRTYPE muID = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE muISO = CorrHelp::NONE
    ) const;

    virtual void load(const SMDecayEvent& genDecays,
            const std::vector<const Lepton*>& selectedLeptons);
    virtual void loadWithoutPromptCheck(const std::vector<const Lepton*>& selectedLeptons);

    const std::vector<const Muon*>& getPromptMuons() const {return promptMuons;}
    const std::vector<const Electron*>& getPromptElectrons() const {return promptElectrons;}

protected:
    template<typename RecoL>
    void loadPromptLeptons(const std::vector<const GenParticle*>& genLeptons,
            const std::vector<const RecoL*>& selectedLeptons,
            std::vector<const RecoL*>& promptLeptons) {
        const ASTypes::size nSL = selectedLeptons.size();
        std::vector<int> matched(nSL,0);
        for(const auto * gl : genLeptons){
            std::vector<std::pair<float, ASTypes::size> > matchedRLeps;
            for(ASTypes::size iR = 0; iR < nSL; ++iR){
                if(matched[iR]) continue;
                const auto * rl = selectedLeptons[iR];
                if(std::fabs(rl->pt() - gl->pt() )/gl->pt() >= 1.0 ) continue;
                const float dr2 = PhysicsUtilities::deltaR2(*gl,*rl);
                if(dr2 >= 0.3*0.3) continue;
                matchedRLeps.emplace_back(dr2,iR);
            }
            if(!matchedRLeps.size()) continue;
            ASTypes::size rLIDX = matchedRLeps[0].second;
            if(matchedRLeps.size() > 1){
                for(const auto& m : matchedRLeps){
                    //negative charge == positive pdgid for e/mu
                    if((selectedLeptons[m.second]->q() >0)  != (gl->pdgId() > 0) ){
                        rLIDX = m.second;
                        break;
                    }
                }
            }
            promptLeptons.push_back(selectedLeptons[rLIDX]);
            matched[rLIDX] = 1;
        }
    }

    const std::string            dataDir        ;
    std::vector<const Muon*>     promptMuons    ;
    std::vector<const Electron*> promptElectrons;
};

//--------------------------------------------------------------------------------------------------
class POGLeptonScaleFactors : public LeptonScaleFactors {
public:
    POGLeptonScaleFactors(const std::string& dataDir);

    float getElectronSF(
            const CorrHelp::CORRTYPE recoT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE isoT = CorrHelp::NOMINAL
    ) const;
    float getMuonSF(
            const CorrHelp::CORRTYPE recoT = CorrHelp::NONE,
            const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL,
            const CorrHelp::CORRTYPE isoT = CorrHelp::NOMINAL
    ) const;

    float flatSFUNC_m_reco = 0.00;
    float flatSFUnc_m_id   = 0.02;
    float flatSFUnc_m_iso  = 0.00;
    float flatSFUNC_e_reco = 0.00;
    float flatSFUnc_e_id   = 0.02;
    float flatSFUnc_e_iso  = 0.00;

private:
    void setParameters(const std::string& el_fileName,const std::string& mu_fileName);

    std::unique_ptr<TObjectHelper::Hist2DContainer> electronRecoSFs;
    std::unique_ptr<TObjectHelper::Hist2DContainer> electronIDSFs;
    std::unique_ptr<TObjectHelper::Hist2DContainer> electronISOSFs;

    std::unique_ptr<TObjectHelper::GraphAEContainer> muonRecoSFs;
    std::unique_ptr<TObjectHelper::Hist2DContainer> muonIDSFs;
    std::unique_ptr<TObjectHelper::Hist2DContainer> muonISOSFs;
};

//
//class ActParamScaleFactors : public LeptonScaleFactors {
//public:
//    ActParamScaleFactors(const std::string& dataDir, const std::string& electronSFFile = "corrections/electronSF_tightID_mini0p1ISO.root",
//            const std::string& muonSFFile = "corrections/muonSF_medID_mini0p2ISO.root",
//            bool verbose = false);
//
//    float getElectronSF(
//            const CorrHelp::CORRTYPE recoT = CorrHelp::NOMINAL, const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL, const CorrHelp::CORRTYPE isoT = CorrHelp::NOMINAL
//    ) const;
//    float getMuonSF(
//            const CorrHelp::CORRTYPE recoT = CorrHelp::NONE, const CorrHelp::CORRTYPE idT = CorrHelp::NOMINAL, const CorrHelp::CORRTYPE isoT = CorrHelp::NOMINAL
//    ) const;
//
//    float flatSFUNC_m_reco = 0.00;
//    float flatSFUnc_m_id   = 0.01;
//    float flatSFUnc_m_iso  = 0.05;
//    float flatSFUNC_e_reco = 0.05;
//    float flatSFUNC_e_reco_ex = 0.10;
//    float flatSFUnc_e_id   = 0.00;
//    float flatSFUnc_e_iso  = 0.00;
//
//private:
//    std::unique_ptr<TObjectHelper::Hist2DContainer> electronRecoSFs;
//    std::unique_ptr<TObjectHelper::Hist2DContainer> electronIDSFs;
//    std::unique_ptr<TObjectHelper::Hist2DContainer> electronISOSFs;
//
//    std::unique_ptr<TObjectHelper::GraphAEContainer> muonRecoSFs;
//    std::unique_ptr<TObjectHelper::Hist2DContainer> muonIDSFs;
//    std::unique_ptr<TObjectHelper::Hist2DContainer> muonISOSFs;
//    std::unique_ptr<TObjectHelper::Hist2DContainer> muonISOActSFs;
//};

}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
