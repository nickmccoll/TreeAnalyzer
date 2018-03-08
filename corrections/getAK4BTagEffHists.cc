
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"

using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_LEP  );
    }


    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

        checkConfig();
    }



    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(!hbbCand || !wjjCand) return false;
        if(hbbNSJs < 2 || wjjNSJs < 2) return false;
        if(hbbCand->sdMom().mass() < 10 || wjjCand->sdMom().mass() < 10) return false;

        auto jets = PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4);
        for(const auto* j : jets){
            if(!j->passTightID()) continue;
            auto flvI = BTagging::jetFlavor(*j);
            TString flvS = "l";
            if(flvI == BTagging::FLV_B) flvS = "b";
            else if(flvI == BTagging::FLV_C) flvS = "c";

            const float pt = j->pt();
            const float absETA = j->absEta();

            auto fill = [&](const TString& label) {
                plotter.getOrMake2DPre(flvS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(pt,absETA,weight);
            };

            fill("incl");
            if(BTagging::isLooseCSVTagged(*j)) fill("loose");
            if(BTagging::isMediumCSVTagged(*j)) fill("med");
            if(BTagging::isTightCSVTagged(*j)) fill("tight");
        }
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getAK4BTagEffHists(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getAK4BTagEffHists(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
