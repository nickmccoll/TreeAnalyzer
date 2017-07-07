
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, bool realData) : BaseTreeAnalyzer(fileName,treeName), realData(realData){

    }
    void loadVariables() override {
        reader_event = (EventReader*)load(new EventReader("event",realData));
        reader_jets = (JetReader*)load(new JetReader("ak4Jet",realData));
    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& jets = reader_jets->jets;
        const auto& genjets = reader_jets->genJets;

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
            if(j.passLooseID()) plotter.getOrMake1D("passLoose_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passTightID()) plotter.getOrMake1D("passTight_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.pt() > 30) n20++;

            plotter.getOrMake2D(TString::Format("%s_csv_hadronflv",jetstr),";csv ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.csv(), j.hadronFlv(), weight);


        }

        plotter.getOrMake1D("njet_20",";N. jets", 10, -0.5,9.5)->Fill(n20,weight);



        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    EventReader * reader_event = 0;
    JetReader * reader_jets = 0;
    HistGetter plotter;

};

#endif

void testJetFiller(std::string fileName ="output.root",std::string outFileName = "plots.root"){
    TString path = fileName;
    Analyzer a(fileName,"treeMaker/Events",path.Contains("data",TString::kIgnoreCase));
    a.analyze();
    a.write(outFileName);
}
