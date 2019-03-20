
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/EventSelection/interface/EventSelection.h"


using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){

        leptonProcNoISO .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProcNoISO);
        leptonProcNoISO->lepSelParams.el_maxISO =-1;
        leptonProcNoISO->lepSelParams.mu_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.el_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.mu_maxISO =-1;
    }


    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }
    }


    bool isRealLepton(const Lepton* lepton){
        int l_pdgid =  lepton->isMuon() ?  ParticleInfo::p_muminus  : ParticleInfo::p_eminus;
        if(lepton->q() > 0) l_pdgid *= -1;

        for(const auto& p : reader_genpart->genParticles ){
            if(p.pdgId() != l_pdgid ) continue;
            if(!ParticleInfo::isLastInChain(&p)) continue;
            if(PhysicsUtilities::deltaR2(p,*lepton) > 0.2*0.2 ) continue;
            return true;
        }
        return false;
    }


    void doPlots(const TString& prefix, const Lepton* l, float r, float rN, float ptRat){
        plotter.getOrMake1DPre(prefix,"drN"           ,";#DeltaR_{N}",200,0,4)->Fill(rN,weight);
        plotter.getOrMake1DPre(prefix,"dr"            ,";#DeltaR",200,0,4)->Fill(r,weight);
        plotter.getOrMake1DPre(prefix,"ptRat"         ,";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill(ptRat,weight);

        plotter.getOrMake1DPre(prefix,"pt"            ,";lepton p_{T}",500,0,500)->Fill(l->pt(),weight);
        plotter.getOrMake1DPre(prefix,"eta"           ,";lepton |#eta|",25,0,2.5)->Fill(l->absEta(),weight);
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(!isRealData()){
            switch(reader_event->process){
            case SIGNAL:
            case TTBAR:
                smpName = FillerConstants::MCProcessNames[reader_event->process];
                break;
            case SINGLET:
            case TTX:
                smpName = "singlet";
                break;
            default:
                "other";
            }
        }
        if(reader_event->process >= FillerConstants::SINGLET && reader_event->process <= FillerConstants::TTX )
            smpName = "other";

        if(!passEventFilters) return false;
        if(!EventSelection::passElectronTriggerSuite(*reader_event)) return false;

        std::vector<const Electron*> tagLeptons = leptonProc->getElectrons(*reader_electron);
        std::vector<const Muon*>     probeLeptons = leptonProcNoISO->getMuons(*reader_event,*reader_muon);
        const std::vector<const Jet*>jets      = PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20,10);

        if(tagLeptons.size() != 1) return false;
        const auto* tagLepton = tagLeptons.front();
        if(tagLepton->pt() < 30) return false;

        std::vector<const Jet*>  filteredJets = PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20.0,2.4);
        std::vector<const Jet*> bjets;
        for(const auto* j : filteredJets) if(BTagging::isMediumCSVTagged(*j)) bjets.push_back(j);
        const size nBjs = bjets.size();
        if(nBjs == 0) return false;

        std::sort(probeLeptons.begin(), probeLeptons.end(),     [&](const Lepton * a, const Lepton * b) -> bool
                {
                    return PhysicsUtilities::deltaR2(*a,*tagLepton) > PhysicsUtilities::deltaR2(*b,*tagLepton);
                });

        plotter.getOrMake1DPre(smpName,"nCands"         ,";number of lepton candidates",20,-0.5,19.5)->Fill(probeLeptons.size(),weight);

        int nCandPio2 = 0;
        for(unsigned int iR = 0; iR < probeLeptons.size(); ++iR){
            if(PhysicsUtilities::deltaR2(*probeLeptons[iR],*tagLepton) < TMath::PiOver2()*TMath::PiOver2() ) break;
            nCandPio2++;
        }
        plotter.getOrMake1DPre(smpName,"nHighETACands"         ,";number of lepton candidates",20,-0.5,19.5)->Fill(nCandPio2,weight);

        if(!nCandPio2) return false;
        const auto* probeLepton = probeLeptons.front();
        bool realLep = false;
        if(!isRealData()) realLep = isRealLepton(probeLepton);


        TString prefix = "data";
        if(!isRealData()) prefix = realLep ? "real" : "fake";

        double nearDR = 10;
        int iJ = PhysicsUtilities::findNearestDRDeref(*probeLepton,jets,nearDR);
        if(nearDR > 0.4) return false;
        const MomentumF jetMLep = jets[iJ]->p4() - probeLepton->p4();
        const float jetML_frac = jetMLep.pt()/probeLepton->pt();

        const float rC  = std::max(0.05,std::min(0.2, 10.0/probeLepton->pt()));
        const float dR = PhysicsUtilities::deltaR(probeLepton->p4(),jetMLep.p4());
        const float ratN = dR/rC;

        bool passISO = leptonProc->isGoodMuon(*reader_event,probeLepton);
        doPlots(prefix,probeLepton,dR,ratN,jetML_frac);
        if(passISO) doPlots(prefix + "_iso",probeLepton,dR,ratN,jetML_frac);
        if(jetML_frac >= 0.5){
            doPlots(prefix+"_jml0p5",probeLepton,dR,ratN,jetML_frac);
            if(passISO) doPlots(prefix + "_jml0p5_iso",probeLepton,dR,ratN,jetML_frac);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;
};

#endif

void getHighPTLepSF(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getHighPTLepSF(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000,100000);
    a.write(outFileName);
}
