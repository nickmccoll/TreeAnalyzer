
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "DataFormats/interface/Muon.h"
#include "DataFormats/interface/Electron.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"



namespace TAna {
using namespace CorrHelp;

//--------------------------------------------------------------------------------------------------
float LeptonScaleFactors::getSF(
        const CORRTYPE elReco, const CORRTYPE elID, const CORRTYPE elISO,
        const CORRTYPE muReco, const CORRTYPE muID, const CORRTYPE muISO
) const {
    return getElectronSF(elReco,elID,elISO) * getMuonSF(muReco,muID,muISO);
}
//--------------------------------------------------------------------------------------------------
void LeptonScaleFactors::load(const SMDecayEvent& genDecays, const std::vector<const Lepton*>& selectedLeptons,
        const std::vector<const Jet*>* jets) {
    jets = 0; //warnign removal :D
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

//--------------------------------------------------------------------------------------------------
POGLeptonScaleFactors::POGLeptonScaleFactors(const std::string& dataDir, const std::string& electronSFFile,
        const std::string& muonSFFile, bool verbose ) {
    TFile * efile = TObjectHelper::getFile(dataDir+electronSFFile,"read",verbose);
    electronRecoSFs.reset(new  TObjectHelper::Hist2DContainer(efile,"reco",verbose) );
    electronIDSFs  .reset(new  TObjectHelper::Hist2DContainer(efile,"id"  ,verbose) );
    electronISOSFs .reset(new  TObjectHelper::Hist2DContainer(efile,"iso" ,verbose) );
    delete efile;
    TFile * mfile = TObjectHelper::getFile(dataDir+muonSFFile,"read",verbose);
    muonRecoSFs.reset(new  TObjectHelper::GraphAEContainer(mfile,"reco",verbose) );
    muonIDSFs  .reset(new  TObjectHelper::Hist2DContainer(mfile,"id"  ,verbose) );
    muonISOSFs .reset(new  TObjectHelper::Hist2DContainer(mfile,"iso" ,verbose) );
    delete mfile;
}
//--------------------------------------------------------------------------------------------------
float POGLeptonScaleFactors::getElectronSF( const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{

    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr  ) -> float {
        if (type == NOMINAL) return sf.val();
        else if(type == UP)  return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
        else if(type == DOWN) return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr  ),float(0.0));
        else return 1;
    };
    float sf = 1;
    for(const auto* l : promptElectrons){
        const float eta = l->scEta();
        const float pt = l->pt();
        if(recoT >= DOWN){ sf *= getSF(recoT,electronRecoSFs->getBinContentByValue(eta,pt),(pt >= 80 || pt < 20 ? flatSFUNC_e_reco : 0.0));}
        if(idT   >= DOWN){ sf *= getSF(idT  ,electronIDSFs->getBinContentByValue(eta,pt)  ,flatSFUnc_e_id  );}
        if(isoT  >= DOWN){ sf *= getSF(isoT ,electronISOSFs->getBinContentByValue(pt,std::fabs(eta)) ,flatSFUnc_e_iso );}
    }
    return sf;
}
//--------------------------------------------------------------------------------------------------
float POGLeptonScaleFactors::getMuonSF( const CORRTYPE recoT, const CORRTYPE idT, const CORRTYPE isoT) const{

    auto getSF = [](const CORRTYPE type, const ASTypes::ValAndErrF sf, const float systErr  ) -> float {
        if (type == NOMINAL) return sf.val();
        else if(type == UP)  return   sf.val() + std::sqrt(sf.err()*sf.err() + systErr*systErr  );
        else if(type == DOWN) return   std::max(sf.val() - std::sqrt(sf.err()*sf.err() + systErr*systErr  ),float(0.0));
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
        const float eta = l->absEta();
        const float pt = l->pt();
        if(recoT >= DOWN){ sf *= getSFNoStat(recoT,muonRecoSFs->eval(l->eta()),flatSFUNC_m_reco);}
        if(idT   >= DOWN){ sf *= getSF(idT  ,muonIDSFs->getBinContentByValue(eta,pt) ,flatSFUnc_m_id  );}
        if(isoT  >= DOWN){ sf *= getSF(isoT ,muonISOSFs->getBinContentByValue(pt,eta),flatSFUnc_m_iso );}
    }
    return sf;
}
}



