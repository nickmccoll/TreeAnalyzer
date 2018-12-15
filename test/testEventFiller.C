
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Corrections/interface/EventWeights.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : BaseTreeAnalyzer(fileName,treeName,treeInt, randSeed){

    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
    }

    bool runEvent() override {
        float eventWeight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());

        plotter.getOrMake1D("met",";#slash{E}_{T}",200,0,1000)->Fill(reader_event->met.pt(),eventWeight);
        plotter.getOrMake1D("met_phi",";#slash{E}_{T} #phi",200,-7,7)->Fill(reader_event->met.phi(),eventWeight);
        plotter.getOrMake1D("npv",";npv",200,-0.5,99.5)->Fill(*reader_event->npv,eventWeight);
        plotter.getOrMake1D("rho",";rho",200,0,100)->Fill(*reader_event->rho,eventWeight);
        if(!isRealData()){
            plotter.getOrMake1D("nTruePUInts",";nTruePUInts",200,-0.5,99.5)->Fill(*reader_event->nTruePUInts,eventWeight);
            plotter.getOrMake2D("process_weight",";weight;process",1000,0,500,20,-0.5,19.5)->Fill(eventWeight,*reader_event->process);
        } else {
            plotter.getOrMake2D("dataset_dataRun",";dataset;dataRun",20,-0.5,19.5,20,-0.5,19.5)->Fill(*reader_event->dataset,*reader_event->dataRun);
        }

        for(unsigned int iT = 0; iT <= FillerConstants::HLT17_NTrig; ++iT ){
            if(FillerConstants::doesPass(*reader_event->triggerAccepts, iT)){
                plotter.getOrMake1D("trigger_accept",";trigger",64,-0.5,63.5)->Fill(iT,eventWeight);
                for(unsigned int iT2 = iT+1; iT2 <= FillerConstants::HLT17_NTrig; ++iT2 ){
                    if(FillerConstants::doesPass(*reader_event->triggerAccepts, iT2)){
                        plotter.getOrMake2D("trigger_accept2D",";trigger;trigger",64,-0.5,63.5,64,-0.5,63.5)->Fill(iT,iT2,eventWeight);
                    }
                }
            }
        }

        for(unsigned int iT = 0; iT <= 25; ++iT ){
            if(FillerConstants::doesPass(*reader_event->metFilters,iT))
                plotter.getOrMake1D("metfilters",";metfilters",33,-1.5,31.5)->Fill(iT);
        }
        if(*reader_event->goodVtx)
            plotter.getOrMake1D("metfilters",";metfilters",33,-1.5,31.5)->Fill(-1);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader      > reader_event    ;
    HistGetter plotter;

};

#endif
void testEventFiller(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
