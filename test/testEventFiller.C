
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
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
