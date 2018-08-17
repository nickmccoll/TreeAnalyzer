
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/PhotonReader.h"

#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());            load(reader_event   );
        reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        reader_electron=std::make_shared<ElectronReader>("electron",false,true);           load(reader_electron);
        reader_photon=std::make_shared<PhotonReader>("photon");                            load(reader_photon);
    }

    bool runEvent() override {
        const float weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        DiHiggsEvent diHiggsEvt; diHiggsEvt.setDecayInfo(reader_genpart->genParticles);


        if(diHiggsEvt.type != DiHiggsEvent::TAU_E && diHiggsEvt.type != DiHiggsEvent::E ) return false;

        bool passGenAcc = diHiggsEvt.w1_d1->pt() > 30 && diHiggsEvt.w1_d1->absEta()<2.5;
        if(!passGenAcc) return false;


        auto makePlots =[&]( int ecIDX) {
            plotter.getOrMake1D(TString::Format("event_count"),";event counts; arbitrary units",20,-0.5,19.5 )->Fill(ecIDX,weight);

        };
        makePlots(0);

        double nearestDR;
        const Electron* electron = 0;
        const Photon* photon = 0;
        int idxE = PhysicsUtilities::findNearestDR(*diHiggsEvt.w1_d1,reader_electron->electrons,nearestDR,0.2);
        int idxP = PhysicsUtilities::findNearestDR(*diHiggsEvt.w1_d1,reader_photon->photons,nearestDR,0.2);

        if(idxE>=0) electron = &reader_electron->electrons[idxE];
        if(idxP>=0) photon = &reader_photon->photons[idxP];

        if(photon) makePlots(1);
        bool goodPhoton = photon && (*reader_photon->hadOvEm)[photon->index()] < 0.15;
        if(goodPhoton) makePlots(2);
        bool goodelectron = goodPhoton && electron;
        if(goodelectron)  makePlots(3);

        if(goodelectron){
            plotter.getOrMake1D(TString::Format("energryRatio"),";electron energy resolution; arbitrary units",400,0,4 )->Fill(electron->energy()/diHiggsEvt.w1_d1->energy(),weight);
        }

        if(goodelectron){
            bool passTrk = FillerConstants::doesPass((*reader_electron->reco_flag)[electron->index()],FillerConstants::ELRECO_TrckDrv);
            bool passEcal =FillerConstants::doesPass((*reader_electron->reco_flag)[electron->index()],FillerConstants::ELRECO_ECALDrv);
            if(passEcal&&passTrk) makePlots(4);
            if(passEcal&&!passTrk ) makePlots(5);
            if(!passEcal&&passTrk) makePlots(6);
            if(!passEcal && !passTrk) makePlots(7);
        }


        if(goodPhoton && (*reader_photon->hadOvEm)[photon->index()]*1.1 < 0.15)  makePlots(8);


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<PhotonReader       > reader_photon     ;
    HistGetter plotter;

};

#endif

void understandElectronReco(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void understandElectronReco(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
