
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "DataFormats/interface/Muon.h"
#include "DataFormats/interface/Electron.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"



namespace TAna {
using namespace CorrHelp;
//--------------------------------------------------------------------------------------------------
// LeptonScaleFactors
//--------------------------------------------------------------------------------------------------
LeptonScaleFactors::LeptonScaleFactors(const std::string& dataDir): dataDir(dataDir){}
//--------------------------------------------------------------------------------------------------
void LeptonScaleFactors::setParameters(const LeptonParameters& param){
    setParameters(param.el_SFFile,param.mu_SFFile);
}
//--------------------------------------------------------------------------------------------------
void LeptonScaleFactors::setParameters(const DileptonParameters& param){
    if(param.el_getID1 != param.el_getID2)
        throw std::invalid_argument(
                std::string("LeptonScaleFactors does not support two different electron IDs"));
    if(param.mu_getID1 != param.mu_getID2)
        throw std::invalid_argument(
                std::string("LeptonScaleFactors does not support two different muon IDs"));
    setParameters(param.el_SFFile,param.mu_SFFile);
}
//--------------------------------------------------------------------------------------------------
float LeptonScaleFactors::getSF(
        const CORRTYPE elReco, const CORRTYPE elID, const CORRTYPE elISO,
        const CORRTYPE muReco, const CORRTYPE muID, const CORRTYPE muISO
) const {
    return getElectronSF(elReco,elID,elISO) * getMuonSF(muReco,muID,muISO);
}
//--------------------------------------------------------------------------------------------------
void LeptonScaleFactors::load(
        const SMDecayEvent& genDecays, const std::vector<const Lepton*>& selectedLeptons) {
    promptMuons    .clear();
    promptElectrons.clear();

    std::vector<const Muon*>     selectedMuons    ;
    std::vector<const Electron*> selectedElectrons;
    for(const auto* sl : selectedLeptons){
        if(sl->isMuon()) selectedMuons.push_back(static_cast<const Muon*>(sl));
        else selectedElectrons.push_back(static_cast<const Electron*>(sl));
    }

    loadPromptLeptons(genDecays.promptElectrons,selectedElectrons,promptElectrons);
    loadPromptLeptons(genDecays.promptMuons,selectedMuons,promptMuons);

}
//--------------------------------------------------------------------------------------------------
void LeptonScaleFactors::loadWithoutPromptCheck(const std::vector<const Lepton*>& selectedLeptons) {
    promptMuons    .clear();
    promptElectrons.clear();

    std::vector<const Muon*>     selectedMuons    ;
    std::vector<const Electron*> selectedElectrons;
    for(const auto* sl : selectedLeptons){
        if(sl->isMuon()) promptMuons.push_back(static_cast<const Muon*>(sl));
        else promptElectrons.push_back(static_cast<const Electron*>(sl));
    }
}


//--------------------------------------------------------------------------------------------------
// POGLeptonScaleFactors
//--------------------------------------------------------------------------------------------------
POGLeptonScaleFactors::POGLeptonScaleFactors(const std::string& dataDir)
    : LeptonScaleFactors(dataDir) {}
//--------------------------------------------------------------------------------------------------
void POGLeptonScaleFactors::setParameters(const std::string& el_fileName,
        const std::string& mu_fileName){
    TFile * efile = TObjectHelper::getFile(dataDir+el_fileName,"read");
    electronRecoSFs.reset(new  TObjectHelper::Hist2DContainer(efile,"reco") );
    electronIDSFs  .reset(new  TObjectHelper::Hist2DContainer(efile,"id"  ) );
//    electronISOSFs .reset(new  TObjectHelper::Hist2DContainer(efile,"iso" ) );
    delete efile;
    TFile * mfile = TObjectHelper::getFile(dataDir+mu_fileName,"read");
//    muonRecoSFs.reset(new  TObjectHelper::GraphAEContainer(mfile,"reco") );
    muonIDSFs  .reset(new  TObjectHelper::Hist2DContainer(mfile,"id"  ) );
//    muonISOSFs .reset(new  TObjectHelper::Hist2DContainer(mfile,"iso") );
    delete mfile;

}
//--------------------------------------------------------------------------------------------------
float POGLeptonScaleFactors::getElectronSF(
        const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{

    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr)->float {
        if (type == NOMINAL)
            return sf.val();
        else if(type == UP)
            return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
        else if(type == DOWN)
            return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr),float(0.0));
        else return 1;
    };
    float sf = 1;
    for(const auto* l : promptElectrons){
        const float eta = l->scEta();
        const float pt = l->pt();
        if(recoT >= DOWN){
            sf *= getSF(recoT,electronRecoSFs->getBinContentByValue(eta,pt),flatSFUNC_e_reco);}
        if(idT   >= DOWN){
            sf *= getSF(idT  ,electronIDSFs->getBinContentByValue(eta,pt)  ,flatSFUnc_e_id  );}
        if(isoT  >= DOWN){
            sf *= getSF(isoT ,electronISOSFs->getBinContentByValue(eta,pt) ,flatSFUnc_e_iso );}
    }
    return sf;
}
//--------------------------------------------------------------------------------------------------
float POGLeptonScaleFactors::getMuonSF(
        const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{

    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr)->float {
        if (type == NOMINAL)
            return sf.val();
        else if(type == UP)
            return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
        else if(type == DOWN)
            return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr),float(0.0));
        else return 1;
    };

    auto getSFNoStat = [](const CORRTYPE type, const float sf, const float systErr  ) -> float {
        if (type == NOMINAL) return sf;
        else if(type == UP)  return   sf + systErr  ;
        else if(type == DOWN) return   std::max(sf - systErr,float(0.0));
        else return 1;
    };

    float sf = 1;
    for(const auto* l : promptMuons){
        const float eta = l->eta();
        const float pt = l->pt();
        if(recoT >= DOWN){
            sf *= getSFNoStat(recoT,muonRecoSFs->eval(eta),flatSFUNC_m_reco);}
        if(idT   >= DOWN){
            sf *= getSF(idT  ,muonIDSFs->getBinContentByValue(eta,pt) ,flatSFUnc_m_id  );}
        if(isoT  >= DOWN){
            sf *= getSF(isoT ,muonISOSFs->getBinContentByValue(eta,pt),flatSFUnc_m_iso );}
    }
    return sf;
}

////--------------------------------------------------------------------------------------------------
//
////--------------------------------------------------------------------------------------------------
//ActParamScaleFactors::ActParamScaleFactors(const std::string& dataDir, const std::string& electronSFFile,
//        const std::string& muonSFFile, bool verbose ) {
//    TFile * efile = TObjectHelper::getFile(dataDir+electronSFFile,"read",verbose);
//    electronRecoSFs.reset(new  TObjectHelper::Hist2DContainer(efile,"reco",verbose) );
//    electronIDSFs  .reset(new  TObjectHelper::Hist2DContainer(efile,"id"  ,verbose) );
//    electronISOSFs .reset(new  TObjectHelper::Hist2DContainer(efile,"iso" ,verbose) );
//    delete efile;
//    TFile * mfile = TObjectHelper::getFile(dataDir+muonSFFile,"read",verbose);
//    muonRecoSFs.reset(new  TObjectHelper::GraphAEContainer(mfile,"reco",verbose) );
//    muonIDSFs  .reset(new  TObjectHelper::Hist2DContainer(mfile,"id"  ,verbose) );
//    muonISOSFs .reset(new  TObjectHelper::Hist2DContainer(mfile,"iso" ,verbose) );
//    muonISOActSFs .reset(new  TObjectHelper::Hist2DContainer(mfile,"iso_act" ,verbose) );
//    delete mfile;
//}
////--------------------------------------------------------------------------------------------------
//float ActParamScaleFactors::getElectronSF( const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{
//
//    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr  ) -> float {
//        if (type == NOMINAL) return sf.val();
//        else if(type == UP)  return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
//        else if(type == DOWN) return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr  ),float(0.0));
//        else return 1;
//    };
//    float sf = 1;
//    for(const auto* l : promptElectrons){
//        const float eta = l->scEta();
//        const float pt = l->pt();
//        if(recoT >= DOWN){
//            const float systError = (((l->sc_act_o_pt() -1.0 ) > 0.6) && (l->sc_dr_act() < 0.08)) ?  flatSFUNC_e_reco_ex :flatSFUNC_e_reco;
//            sf *= getSF(recoT,electronRecoSFs->getBinContentByValue(eta,pt),systError);
//        }
//        if(idT   >= DOWN){ sf *= getSF(idT  ,electronIDSFs->getBinContentByValue(eta,pt)  ,flatSFUnc_e_id  );}
//        if(isoT  >= DOWN){ sf *= getSF(isoT ,electronISOSFs->getBinContentByValue(eta,pt) ,flatSFUnc_e_iso );}
//    }
//    return sf;
//}
////--------------------------------------------------------------------------------------------------
//float ActParamScaleFactors::getMuonSF( const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{
//
//    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr  ) -> float {
//        if (type == NOMINAL) return sf.val();
//        else if(type == UP)  return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
//        else if(type == DOWN) return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr  ),float(0.0));
//        else return 1;
//    };
//
//    auto getSFNoStat = [](const CORRTYPE type, const float sf, const float systErr  ) -> float {
//        if (type == NOMINAL) return sf;
//        else if(type == UP)  return   sf + systErr  ;
//        else if(type == DOWN) return   std::max(sf - systErr,float(0.0));
//        else return 1;
//    };
//
//    float sf = 1;
//    for(const auto* l : promptMuons){
//        const float eta = l->absEta();
//        const float pt = l->pt();
//        if(recoT >= DOWN){ sf *= getSFNoStat(recoT,muonRecoSFs->eval(l->eta()),flatSFUNC_m_reco);}
//        if(idT   >= DOWN){ sf *= getSF(idT  ,muonIDSFs->getBinContentByValue(eta,pt) ,flatSFUnc_m_id  );}
//        if(isoT  >= DOWN){
//            if(l->lepAct_o_pt() < 0.2 || l->dRnorm() > 1.0 )
//                sf *= getSF(isoT ,muonISOSFs->getBinContentByValue(eta,pt),flatSFUnc_m_iso );
//            else
//                sf *= getSF(isoT ,muonISOActSFs->getBinContentByValue(l->dRnorm(),l->lepAct_o_pt()),flatSFUnc_m_iso );
//
//        }
//    }
//    return sf;
//}
}



