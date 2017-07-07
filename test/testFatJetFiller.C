
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, bool realData) : BaseTreeAnalyzer(fileName,treeName), realData(realData){

    }
    void loadVariables() override {
        reader_event = (EventReader*)load(new EventReader("event",realData));
        reader_fatjets = (FatJetReader*)load(new FatJetReader("ak8PuppiNoLepJet",realData));
    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& jets = reader_fatjets->jets;
        const auto& genjets = reader_fatjets->genJets;

        auto getJS =[] (float pt)-> TString {
            TString jetstr = "pt_";
            if(pt < 50) jetstr += "lt50";
            else if(pt < 100) jetstr += "50to100";
            else if(pt < 200) jetstr += "100to200";
            else jetstr += "geq200";
            return jetstr;
        };
        int n20 = 0;
        for(const auto& j : jets){
            auto jetstr = getJS(j.pt()).Data();
            plotter.getOrMake1D(TString::Format("%s_eta",jetstr),";#eta", 10, -5,5)->Fill(j.eta(),weight);
            if(j.absEta() > 2.4) continue;

            if(j.genJet()){
                plotter.getOrMake1D(TString::Format("%s_res",getJS(j.genJet()->pt()).Data()),";jet res",100,0,5)->
                Fill(j.pt()/j.genJet()->pt(), weight);

                plotter.getOrMake1D(TString::Format("%s_sdres",getJS(j.genJet()->pt()).Data()),";jet sdres",100,0,5)->
                Fill(j.sdMom().pt()/j.genJet()->pt(), weight);

                plotter.getOrMake1D(TString::Format("%s_rawsdres",getJS(j.genJet()->pt()).Data()),";jet raw sdres",100,0,5)->
                Fill(j.rawSdMom().pt()/j.genJet()->pt(), weight);
            }


            plotter.getOrMake2D(TString::Format("%s_hadronflv_partonflv",jetstr),";hadron flv ;parton flv",10,-0.5,9.5,10,-0.5,9.5)->
            Fill(std::abs(j.hadronFlv()),std::abs(j.partonFlv()), weight);
            plotter.getOrMake1D("incl_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passPUID()) plotter.getOrMake1D("passPU_pt"      ,";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passLooseID()) plotter.getOrMake1D("passLoose_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.passTightID()) plotter.getOrMake1D("passTight_pt",";p_{T}", 20, 0,200)->Fill(j.pt(),weight);
            if(j.pt() > 30) n20++;

            plotter.getOrMake2D(TString::Format("%s_csv_hadronflv",jetstr),";csv ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.csv(), j.hadronFlv(), weight);

            plotter.getOrMake2D(TString::Format("%s_mass_sdmass",jetstr),";mass; sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.sdMom().mass(),weight);
            plotter.getOrMake2D(TString::Format("%s_mass_rawsdmass",jetstr),";mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.rawSdMom().mass(),weight);
            plotter.getOrMake2D(TString::Format("%s_sdmass_rawsdmass",jetstr),";sd mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.sdMom().mass(),j.rawSdMom().mass(),weight);

            plotter.getOrMake1D(TString::Format("%s_tau2otau1",jetstr),";tau2otau1", 100,0,1)->Fill(j.tau2otau1(),weight);
            plotter.getOrMake1D(TString::Format("%s_tau3otau1",jetstr),";tau3otau1", 100,0,1)->Fill(j.tau3otau1(),weight);
            plotter.getOrMake1D(TString::Format("%s_tau3otau2",jetstr),";tau3otau2", 100,0,1)->Fill(j.tau3otau2(),weight);

            plotter.getOrMake2D(TString::Format("%s_bbt_hadronflv",jetstr),";bbt ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.bbt(), j.hadronFlv(), weight);

            plotter.getOrMake2D(TString::Format("%s_minsjcsv_hadronflv",jetstr),";minsjcsv ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.minSJCSV(), j.hadronFlv(), weight);

            plotter.getOrMake2D(TString::Format("%s_maxsjcsv_hadronflv",jetstr),";maxsjcsv ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.maxSJCSV(), j.hadronFlv(), weight);

            plotter.getOrMake1D("nsd",";N. sub jets", 3, -0.5,2.5)->Fill(j.nSubJets(),weight);
            int nSJ20=0;
            for(unsigned int iSJ = 0; iSJ < j.nSubJets(); ++iSJ){
                if(j.subJet(iSJ).pt() > 20) nSJ20++;
            }
            plotter.getOrMake1D("nsd_20",";N. sub jets", 3, -0.5,2.5)->Fill(nSJ20,weight);




        }

        plotter.getOrMake1D("njet_20",";N. jets", 10, -0.5,9.5)->Fill(n20,weight);



        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    EventReader * reader_event = 0;
    FatJetReader * reader_fatjets = 0;
    HistGetter plotter;

};

#endif

void testFatJetFiller(std::string fileName ="output.root",std::string outFileName = "plots.root"){
    TString path = fileName;
    Analyzer a(fileName,"treeMaker/Events",path.Contains("data",TString::kIgnoreCase));
    a.analyze();
    a.write(outFileName);
}
