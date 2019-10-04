
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;
//--------------------------------------------------------------------------------------------------
// testJetFiller
//--------------------------------------------------------------------------------------------------
class Analyzer : public BaseTreeAnalyzer {
public:
    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
    : BaseTreeAnalyzer(fileName,treeName,treeInt, randSeed){
    }
    void loadVariables() override {
        reader_event = loadReader<EventReader>("event",isRealData());
        reader_jet         =loadReader<JetReader>     ("ak4Jet",isRealData());
    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& jets = reader_jet->jets;
        const auto& genjets = reader_jet->genJets;

        auto getJS =[] (float pt)-> TString {
            TString jetstr = "pt_";
            if(pt < 20) jetstr += "lt20";
            else if(pt < 30) jetstr += "20to30";
            else if(pt < 50) jetstr += "30to50";
            else if(pt < 100) jetstr += "50to100";
            else jetstr += "geq100";
            return jetstr;
        };
        int n20 = 0;
        for(const auto& j : jets){
            auto jetstr = getJS(j.pt()).Data();
            plotter.getOrMake1D(TString::Format("%s_eta",jetstr),";#eta", 10, -5,5)->Fill(j.eta(),weight);
            if(j.absEta() > 2.4) continue;

            if(j.genJet())
                plotter.getOrMake1D(TString::Format("%s_res",getJS(j.genJet()->pt()).Data()),";jet res",100,0,5)->
                Fill(j.pt()/j.genJet()->pt(), weight);

            plotter.getOrMake2D(TString::Format("%s_hadronflv_partonflv",jetstr),";hadron flv ;parton flv",10,-0.5,9.5,10,-0.5,9.5)->
            Fill(std::abs(j.hadronFlv()),std::abs(j.partonFlv()), weight);
            plotter.getOrMake1D("incl_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passPUID()) plotter.getOrMake1D("passPU_pt"      ,";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passTightID()) plotter.getOrMake1D("passTight_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passTightNoLepID()) plotter.getOrMake1D("passTightNoLep_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.pt() > 30) n20++;

            plotter.getOrMake2D(TString::Format("%s_csv_hadronflv",jetstr),";csv ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.csv(), j.hadronFlv(), weight);

            plotter.getOrMake1D("csv",";csv", 80, -2,2)->Fill(j.csv(),weight);
            plotter.getOrMake1D("deep_csv",";deep csv", 80, -2,2)->Fill(j.deep_csv(),weight);
            plotter.getOrMake1D("deep_flavor",";deep csv", 80, -2,2)->Fill(j.deep_flavor(),weight);


        }
        plotter.getOrMake1D("njet_20",";N. jets", 10, -0.5,9.5)->Fill(n20,weight);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<JetReader        > reader_jet        ;
    HistGetter plotter;

};

#endif

void testJetFiller(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
