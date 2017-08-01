#include "../interface/DefaultSearchRegionAnalyzer.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Variables/interface/JetKinematics.h"

#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"

#include "TPRegexp.h"


namespace TAna {
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::setupProcessors(std::string fileName) {
    TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
    auto match = r1.MatchS(fileName);
    const Int_t nrSubStr = match->GetLast()+1;
    if(nrSubStr>1){
        signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
    }
    fjProc     .reset(new FatJetProcessor ()); DefaultFatJetSelections::setDefaultFatJetProcessor(*fjProc);
    leptonProc .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProc);
}
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::loadVariables()  {
    reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
    reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
    reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
    reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData());  load(reader_jet     );
    reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
    reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

    if(!isRealData()){
        reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
    }
}
//--------------------------------------------------------------------------------------------------
bool DefaultSearchRegionAnalyzer::runEvent() {
    if(reader_jetwlep){
        auto jets = JetKinematics::selectObjects(reader_jetwlep->jets,30);
        ht_wlep = JetKinematics::ht(jets);
    }
    weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
    passEventFilters= EventWeights::passEventFilters(*reader_event);

    if(reader_genpart && reader_event->process == FillerConstants::SIGNAL)
        diHiggsEvt.setDecayInfo(reader_genpart->genParticles);

    if(reader_electron && reader_muon){
        selectedLeptons = leptonProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
        selectedLepton = selectedLeptons.size() ? selectedLeptons.front() : 0;
    }

    if(reader_fatjet && selectedLepton){
        fatjetCands = fjProc->loadFatJets(*reader_fatjet,selectedLepton);
        hbbCand     = fjProc->getHBBCand();
        wjjCand     = fjProc->getWjjCand();
        passHbbSel  = fjProc->passHbbSel();
        passHbbTSel = fjProc->passHbbSelTightBTag();
        passWjjSel  = fjProc->passWjjSel();

        neutrino = wjjCand ? HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + wjjCand->sdMom().p4()) ) : MomentumF();
        hh =  wjjCand && hbbCand ? (selectedLepton->p4() + neutrino.p4()+ wjjCand->sdMom().p4() + hbbCand->sdMom().p4()) :  MomentumF();

    } else {
        fatjetCands.clear();
        wjjCand    =  0;
        hbbCand    =  0;
        passWjjSel =  false;
        passHbbSel =  false;
        passHbbTSel=  false;
        neutrino   =  MomentumF();
        hh         =  MomentumF();
    }

    return true;
}
}
