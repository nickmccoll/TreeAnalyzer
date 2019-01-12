
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;
//--------------------------------------------------------------------------------------------------
// testFatJetFiller
//--------------------------------------------------------------------------------------------------
class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: BaseTreeAnalyzer(fileName,treeName,treeInt, randSeed){
    }
    void loadVariables() override {
        reader_event = loadReader<EventReader>("event",isRealData());
        //        reader_fatjets = loadReader<FatJetReader>("ak8PuppiJet",isRealData(),true,true);
        reader_fatjets = loadReader<FatJetReader>("ak8PuppiNoLepJet",isRealData(),false,false);
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
            auto jetstr = getJS(j.pt());
            plotter.getOrMake1DPre(jetstr,"eta",";#eta", 10, -5,5)->Fill(j.eta(),weight);
            if(j.absEta() > 2.4) continue;

            if(j.genJet()){
                auto genJetstr = getJS(j.genJet()->pt());
                plotter.getOrMake1DPre(genJetstr,"res",";jet res",100,0,5)->
                        Fill(j.pt()/j.genJet()->pt(), weight);

                plotter.getOrMake1DPre(genJetstr,"sdres",";jet sdres",100,0,5)->
                        Fill(j.sdMom().pt()/j.genJet()->pt(), weight);

                plotter.getOrMake1DPre(genJetstr,"rawsdres",";jet raw sdres",100,0,5)->
                        Fill(j.rawSdMom().pt()/j.genJet()->pt(), weight);
            }


            plotter.getOrMake2DPre(jetstr,"hadronflv_partonflv",";hadron flv ;parton flv",10,-0.5,9.5,10,-0.5,9.5)->
                    Fill(std::abs(j.hadronFlv()),std::abs(j.partonFlv()), weight);
            plotter.getOrMake1D("incl_pt",";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
            if(j.passPUID()) plotter.getOrMake1D("passPU_pt"      ,";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
            if(j.passLooseID()) plotter.getOrMake1D("passLoose_pt",";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
            if(j.passTightID()) plotter.getOrMake1D("passTight_pt",";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
            if(j.pt() > 30) n20++;

            plotter.getOrMake2DPre(jetstr,"mass_sdmass",";mass; sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.sdMom().mass(),weight);
            plotter.getOrMake2DPre(jetstr,"mass_rawsdmass",";mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.rawSdMom().mass(),weight);
            plotter.getOrMake2DPre(jetstr,"sdmass_rawsdmass",";sd mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.sdMom().mass(),j.rawSdMom().mass(),weight);

            //
            plotter.getOrMake1DPre(jetstr,"tau2otau1",";tau2otau1", 100,0,1)->Fill(j.tau2otau1(),weight);
            plotter.getOrMake1DPre(jetstr,"tau3otau1",";tau3otau1", 100,0,1)->Fill(j.tau3otau1(),weight);
            plotter.getOrMake1DPre(jetstr,"tau3otau2",";tau3otau2", 100,0,1)->Fill(j.tau3otau2(),weight);

            plotter.getOrMake1DPre(jetstr,"ecfb1",";ecfb1", 100,0,1)->Fill(j.ecfb1(),weight);
            plotter.getOrMake1DPre(jetstr,"ecfb2",";ecfb2", 100,0,1)->Fill(j.ecfb2(),weight);

            plotter.getOrMake2DPre(jetstr,"bbt_hadronflv",";bbt ;hadron flv",50,0,1,10,-0.5,9.5)->
                    Fill(j.bbt(), j.hadronFlv(), weight);

            plotter.getOrMake1DPre(jetstr,"bbt",";bbt ;hadron flv",40,-1,1)->Fill(j.bbt(), weight);

            double minSJCSV=0,maxSJCSV=0,minSJDCSV=0,maxSJDCSV=0;

            std::vector<float> sjCSV;
            std::vector<float> sjDCSV;
            for(unsigned int iSJ = 0; iSJ < j.nSubJets(); ++iSJ){
                sjCSV.push_back(j.subJet(iSJ).csv());
                sjDCSV.push_back(j.subJet(iSJ).deep_csv());
            }
            std::sort(sjCSV.begin(), sjCSV.end(), [](float a, float b){return a>b;});
            std::sort(sjDCSV.begin(), sjDCSV.end(), [](float a, float b){return a>b;});

            if(sjCSV.size()){
                plotter.getOrMake1DPre(jetstr,"minSJCSV" ,";minSJCSV", 100,0,1)->Fill(sjCSV.back(),weight);
                plotter.getOrMake1DPre(jetstr,"maxSJCSV" ,";maxSJCSV", 100,0,1)->Fill(sjCSV.front(),weight);
                plotter.getOrMake1DPre(jetstr,"minSJDCSV",";minSJDCSV", 100,0,1)->Fill(sjDCSV.back(),weight);
                plotter.getOrMake1DPre(jetstr,"maxSJDCSV",";maxSJDCSV", 100,0,1)->Fill(sjDCSV.front(),weight);
            }



            //            plotter.getOrMake2D(TString::Format("%s_minsjcsv_hadronflv",jetstr),";minsjcsv ;hadron flv",50,0,1,10,-0.5,9.5)->
            //                    Fill(j.minSJCSV(), j.hadronFlv(), weight);
            //
            //            plotter.getOrMake2D(TString::Format("%s_maxsjcsv_hadronflv",jetstr),";maxsjcsv ;hadron flv",50,0,1,10,-0.5,9.5)->
            //                    Fill(j.maxSJCSV(), j.hadronFlv(), weight);

            plotter.getOrMake1DPre(jetstr,"nsd",";N. sub jets", 3, -0.5,2.5)->Fill(j.nSubJets(),weight);
            int nSJ20=0;
            for(unsigned int iSJ = 0; iSJ < j.nSubJets(); ++iSJ){
                if(j.subJet(iSJ).pt() > 20) nSJ20++;
                plotter.getOrMake1DPre(jetstr,"sd_pt",";sub jet pt", 20, 0,200)->Fill(j.subJet(iSJ).pt(),weight);
            }
            plotter.getOrMake1DPre(jetstr,"nsd_20",";N. sub jets", 3, -0.5,2.5)->Fill(nSJ20,weight);




        }

        plotter.getOrMake1D("njet_20",";N. jets", 10, -0.5,9.5)->Fill(n20,weight);



        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    std::shared_ptr<EventReader> reader_event = 0;
    std::shared_ptr<FatJetReader> reader_fatjets = 0;
    HistGetter plotter;

};

#endif

void testFatJetFiller(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
