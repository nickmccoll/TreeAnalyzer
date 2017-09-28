
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
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/FatJetSelection.h"

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
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

        checkConfig();
    }



    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";

        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(!hbbCand || !wjjCand) return false;
        if(hbbCand->sdMom().mass() < 10 || wjjCand->sdMom().mass() < 10) return false;
        if(hbbCand->nSubJets() < 2 || wjjCand->nSubJets() < 2) return false;
        std::vector<const SubJet*> subjets;
        subjets.push_back(&hbbCand->subJet(0));
        subjets.push_back(&hbbCand->subJet(1));
        subjets.push_back(&wjjCand->subJet(0));
        subjets.push_back(&wjjCand->subJet(1));

        const float xsecW  = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());


        for(const auto* j : subjets){
            auto flvI = BTagging::jetFlavor(*j);
            TString flvS = "l";
            if(flvI == BTagging::FLV_B) flvS = "b";
            else if(flvI == BTagging::FLV_C) flvS = "c";

            const float pt = j->pt();
            const float absETA = std::min(j->absEta(),float(2.8));

            auto fill = [&](const TString& label) {
                plotter.getOrMake2DPre(flvS, label,";jet p_{T}[GeV];jet |#eta|",200,0,1000,14,0,2.8)->Fill(pt,absETA,weight);
                plotter.getOrMake2DPre(flvS, label+"_noxsec",";jet p_{T}[GeV];jet |#eta|",200,0,1000,14,0,2.8)->Fill(pt,absETA,weight/xsecW);
            };

            fill("incl");
            if(BTagging::isLooseCSVTagged(*j)) fill("loose");
            if(BTagging::isMediumCSVTagged(*j)) fill("med");
            if(BTagging::isTightCSVTagged(*j)) fill("tight");
        }

        auto mkPlt = [&](const TString& name, float pt, float eta){
            plotter.getOrMake1DPre(name, "minSJPT",";jet p_{T}[GeV]",100,0,500)->Fill(pt,weight);
            plotter.getOrMake1DPre(name, "maxSJETA",";jet |#eta|",50,0,5.0)->Fill(eta,weight);
            plotter.getOrMake2DPre(name, "minSJPT_maxSJETA",";jet p_{T}[GeV];jet |#eta|",100,0,500,50,0,5.0)->Fill(pt,eta,weight);

            if(hh.mass() > 1000){
                plotter.getOrMake1DPre(name+"_m1000", "minSJPT",";jet p_{T}[GeV]",100,0,500)->Fill(pt,weight);
                plotter.getOrMake1DPre(name+"_m1000", "maxSJETA",";jet |#eta|",50,0,5.0)->Fill(eta,weight);
                plotter.getOrMake2DPre(name+"_m1000", "minSJPT_maxSJETA",";jet p_{T}[GeV];jet |#eta|",100,0,500,50,0,5.0)->Fill(pt,eta,weight);
            }

        };

        if(passHbbSel && passWjjSel ) mkPlt(smpName + "_hbb",
                std::min(hbbCand->subJet(0).pt(),hbbCand->subJet(1).pt()),
                std::max(hbbCand->subJet(0).absEta(),hbbCand->subJet(1).absEta()));

        if(passHbbSel && passWjjSel) mkPlt(smpName + "_wjj",
                std::min(wjjCand->subJet(0).pt(),wjjCand->subJet(1).pt()),
                std::max(wjjCand->subJet(0).absEta(),wjjCand->subJet(1).absEta()));

        const bool passWJJI =  FatJetSelHelpers::passWjjSelection(wjjCand,fjProc->wjj_maxT2oT1,BTagging::CSV_INCL,fjProc->wjj_minMass,fjProc->wjj_maxMass);

        if(passHbbSel && passWJJI) mkPlt(smpName + "_wjjI",
                std::min(wjjCand->subJet(0).pt(),wjjCand->subJet(1).pt()),
                std::max(wjjCand->subJet(0).absEta(),wjjCand->subJet(1).absEta()));

        std::vector<const SubJet*> btaggedSJs = PhysicsUtilities::selObjs(wjjCand->subJets(), [](const SubJet* sj){return BTagging::isMediumCSVTagged(*sj); } );
        if(btaggedSJs.size()){
            float minbpt = btaggedSJs.front()->pt();
            float maxbeta = btaggedSJs.front()->absEta();
            if(btaggedSJs.size() > 1){
                minbpt  = std::min(btaggedSJs[1]->pt(),minbpt);
                maxbeta =  std::max(btaggedSJs[1]->absEta(),maxbeta);
            }
            if(passHbbSel && passWJJI) mkPlt(smpName + "_wjjB",
                    minbpt,
                    maxbeta);
        } else {
            if(passHbbSel && passWJJI) mkPlt(smpName + "_wjjB",
                    -1,
                    -1);
        }


        return true;


    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getSJTagEffHists(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getSJTagEffHists(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
