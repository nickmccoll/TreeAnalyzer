
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
        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData(),true,true);
        reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,false,true);

    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& jetsWLep = reader_fatjet->jets;
        const auto& jetsNoLeps = reader_fatjet_noLep->jets;
        const auto& genjets = reader_fatjet->genJets;

        auto getJS =[] (float pt)-> TString {
            TString jetstr = "pt_";
            if(pt < 50) jetstr += "lt50";
            else if(pt < 100) jetstr += "50to100";
            else if(pt < 200) jetstr += "100to200";
            else jetstr += "geq200";
            return jetstr;
        };

        auto procJets=[&](TString prefix, const std::vector<FatJet>& jets){
            int n20 = 0;
            for(const auto& j : jets){
                auto jetstr = prefix+"_"+getJS(j.pt());
                plotter.getOrMake1DPre(jetstr,"eta",";#eta", 10, -5,5)->Fill(j.eta(),weight);
                if(j.absEta() > 2.4) continue;

                if(j.genJet()){
                    auto genJetstr = prefix+"_"+getJS(j.genJet()->pt());
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
                if(j.passPUID()) plotter.getOrMake1DPre(prefix,"passPU_pt"      ,";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
                if(j.passTightID()) plotter.getOrMake1DPre(prefix,"passTight_pt",";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
                if(j.passTightNoLepID()) plotter.getOrMake1DPre(prefix,"passTightNoLep_pt",";p_{T}", 20, 0,2000)->Fill(j.pt(),weight);
                if(j.pt() > 30) n20++;

                plotter.getOrMake2DPre(jetstr,"mass_sdmass",";mass; sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.sdMom().mass(),weight);
                plotter.getOrMake2DPre(jetstr,"mass_rawsdmass",";mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.mass(),j.rawSdMom().mass(),weight);
                plotter.getOrMake2DPre(jetstr,"sdmass_rawsdmass",";sd mass; raw sd mass", 100,0,1000,100,0,1000)->Fill(j.sdMom().mass(),j.rawSdMom().mass(),weight);


                plotter.getOrMake1DPre(jetstr,"tau2otau1",";tau2otau1", 100,0,1)->Fill(j.tau2otau1(),weight);


                plotter.getOrMake2DPre(jetstr,"bbt_hadronflv",";bbt ;hadron flv",50,0,1,10,-0.5,9.5)->
                        Fill(j.bbt(), j.hadronFlv(), weight);

                plotter.getOrMake1DPre(jetstr,"bbt",";bbt ;hadron flv",40,-1,1)->Fill(j.bbt(), weight);

                plotter.getOrMake1DPre(jetstr,"deep_MDZHbb",";deep_MDZHbb",40,-1,1)->Fill(j.deep_MDZHbb(), weight);
                plotter.getOrMake1DPre(jetstr,"deep_MDHbb",";deep_MDHbb",40,-1,1)->Fill(j.deep_MDHbb(), weight);
                plotter.getOrMake1DPre(jetstr,"deep_Hbb",";deep_Hbb",40,-1,1)->Fill(j.deep_Hbb(), weight);
                plotter.getOrMake1DPre(jetstr,"deep_W",";deep_MDZHbb",40,-1,1)->Fill(j.deep_W(), weight);

                plotter.getOrMake1DPre(jetstr,"raw_SDNom",";raw_SDNom",100,0,300)->Fill(j.rawMom().mass(), weight);
                plotter.getOrMake1DPre(jetstr,"rawSDMassUp",";rawSDMassUp",100,0,300)->Fill(j.rawSDMassUp(), weight);
                plotter.getOrMake1DPre(jetstr,"rawSDMassDown",";rawSDMassDown",100,0,300)->Fill(j.rawSDMassDown(), weight);

                double minSJDCSV=0,maxSJDCSV=0;

                std::vector<float> sjDCSV;
                for(unsigned int iSJ = 0; iSJ < j.nSubJets(); ++iSJ){
                    sjDCSV.push_back(j.subJet(iSJ).deep_csv());
                }
                std::sort(sjDCSV.begin(), sjDCSV.end(), [](float a, float b){return a>b;});

                if(sjDCSV.size()){
                    plotter.getOrMake1DPre(jetstr,"minSJDCSV",";minSJDCSV", 100,0,1)->Fill(sjDCSV.back(),weight);
                    plotter.getOrMake1DPre(jetstr,"maxSJDCSV",";maxSJDCSV", 100,0,1)->Fill(sjDCSV.front(),weight);
                }


                plotter.getOrMake1DPre(jetstr,"nsd",";N. sub jets", 3, -0.5,2.5)->Fill(j.nSubJets(),weight);
                int nSJ20=0;
                for(unsigned int iSJ = 0; iSJ < j.nSubJets(); ++iSJ){
                    if(j.subJet(iSJ).pt() > 20) nSJ20++;
                    plotter.getOrMake1DPre(jetstr,"sd_pt",";sub jet pt", 20, 0,200)->Fill(j.subJet(iSJ).pt(),weight);
                }
                plotter.getOrMake1DPre(jetstr,"nsd_20",";N. sub jets", 3, -0.5,2.5)->Fill(nSJ20,weight);




            }

            plotter.getOrMake1DPre(prefix,"njet_20",";N. jets", 10, -0.5,9.5)->Fill(n20,weight);



        };

        procJets("jetsWLep",jetsWLep);
        procJets("jetsNoLeps",jetsNoLeps);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    std::shared_ptr<EventReader> reader_event = 0;
    std::shared_ptr<FatJetReader> reader_fatjet = 0;
    std::shared_ptr<FatJetReader> reader_fatjet_noLep = 0;
    HistGetter plotter;

};

#endif

void testFatJetFiller(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
