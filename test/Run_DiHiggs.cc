#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

//#include "/Users/brentstone01/HistoPlotting/include/Plotter.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

#include <cmath>
#include <vector>

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt) {

    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        if(!isRealData()){
            reader_genParticles =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genParticles   );
        }
    }

// Function that runs the analysis on each event and will called iteratively for each event
    bool runEvent() override {
    	DiHiggsEvent event;
	//std::cout << reader_event->run << " : " << reader_event->lumi << " : " << reader_event->event << std::endl;
        event.setDecayInfo(reader_genParticles->genParticles);

        plotter.getOrMake1D("h_decaymode","Distribution of Decay Modes",9,-0.5,8.5)->Fill(event.type);
        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader> reader_event;
    std::shared_ptr<GenParticleReader> reader_genParticles;
    HistGetter plotter;
};

void Run_DiHiggs(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();

    auto *h_decaymode = a.plotter.getOrMake1D("h_decaymode","Distribution of Decay Modes",9,-0.5,8.5);
    for (int iX = 1; iX <= h_decaymode->GetNbinsX(); ++iX) {
     	std::cout << h_decaymode->GetBinContent(iX) << std::endl;
    }
    a.write(outFileName);
}

void Run_DiHiggs(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();

    auto *h_decaymode = a.plotter.getOrMake1D("h_decaymode","Distribution of Decay Modes",9,-0.5,8.5);
    for (int iX = 1; iX <= h_decaymode->GetNbinsX(); ++iX) {
        std::cout << h_decaymode->GetBinContent(iX) << std::endl;
    }
    a.write(outFileName);
}

