
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, bool realData) : BaseTreeAnalyzer(fileName,treeName), realData(realData){

    }
    void loadVariables() override {
        reader_event = (EventReader*)load(new EventReader("event",realData));
    }

    bool runEvent() override {
        plotter.getOrMake1D("met",";#slash{E}_{T}",200,0,1000)->Fill(reader_event->met.pt(),reader_event->weight);
        plotter.getOrMake1D("met_o_rawmet",";#slash{E}_{T}/raw #slash{E}_{T}",200,0,10)->Fill(reader_event->met.pt()/reader_event->rawMet.pt(),reader_event->weight);
        plotter.getOrMake1D("met_phi",";#slash{E}_{T} #phi",200,-7,7)->Fill(reader_event->met.phi(),reader_event->weight);
        plotter.getOrMake1D("npv",";npv",200,-0.5,99.5)->Fill(reader_event->npv,reader_event->weight);
        plotter.getOrMake1D("rho",";rho",200,0,100)->Fill(reader_event->rho,reader_event->weight);
        if(!realData){
            plotter.getOrMake1D("nTruePUInts",";nTruePUInts",200,-0.5,99.5)->Fill(reader_event->nTruePUInts,reader_event->weight);
            plotter.getOrMake2D("process_weight",";weight;process",1000,0,500,20,-0.5,19.5)->Fill(reader_event->weight,reader_event->process);
        } else {
            plotter.getOrMake2D("dataset_dataRun",";dataset;dataRun",20,-0.5,19.5,20,-0.5,19.5)->Fill(reader_event->dataset,reader_event->dataRun);
        }

        for(unsigned int iT = 0; iT <= FillerConstants::HLT_PFMET120_PFMHT120_IDTight; ++iT ){
            if(FillerConstants::doesPass(reader_event->triggerAccepts, iT))
                plotter.getOrMake2D("trigger_prescale",";trigger;isPrescaled",64,-0.5,63.5,2,-0.5,1.5)->Fill(iT,0.0);
            if(FillerConstants::doesPass(reader_event->triggerPrescales, iT))
                plotter.getOrMake2D("trigger_prescale",";trigger;isPrescaled",64,-0.5,63.5,2,-0.5,1.5)->Fill(iT,1.0);
        }

        for(unsigned int iT = 0; iT <= FillerConstants::AnaTM_hitsNotReplaced; ++iT ){
            if(FillerConstants::doesPass(reader_event->metFilters, iT))
                plotter.getOrMake1D("metfilters",";metfilters",33,-1.5,31.5)->Fill(iT);
        }
        if(reader_event->goodVtx)
            plotter.getOrMake1D("metfilters",";metfilters",33,-1.5,31.5)->Fill(-1);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    EventReader * reader_event = 0;
    HistGetter plotter;

};

#endif

void testEventFiller(std::string fileName ="output.root",std::string outFileName = "plots.root"){
    TString path = fileName;
    Analyzer a(fileName,"treeMaker/Events",path.Contains("data",TString::kIgnoreCase));
    a.analyze();
    a.write(outFileName);
}
