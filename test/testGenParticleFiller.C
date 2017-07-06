
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, bool realData) : BaseTreeAnalyzer(fileName,treeName), realData(realData){

    }
    void loadVariables() override {
        reader_event = (EventReader*)load(new EventReader("event",realData));
        if(!realData)
            reader_genParticles = (GenParticleReader*)load(new GenParticleReader("genParticle"));
    }

    bool runEvent() override {
        const auto& genparticles = reader_genParticles->genParticles;
//        if(getEventNumber() < 10)
            std::cout <<getEventNumber()<<std::endl;
            ParticleInfo::printGenInfo(genparticles,-1);

        std::vector<const GenParticle *> genMuons;
        for(const auto& gp : genparticles){
            if(gp.absPdgId() != ParticleInfo::p_Z0 ) continue;
            if(!ParticleInfo::isLastInChain(&gp)) continue;
//            if(gp.numberOfDaughters()!= 2) continue;
            for(unsigned int iD = 0; iD < gp.numberOfDaughters(); ++iD ){
                if(gp.daughter(iD)->absPdgId() != ParticleInfo::p_muminus) continue;
                genMuons.push_back(&*gp.daughter(iD));
            }
        }
        if(genMuons.size() != 2) return false;

        MomentumF diMuon(genMuons[0]->p4() + genMuons[1]->p4());

        int nMuonsETA2p8 =0;
        int nMuonsETA2p4 =0;
        int nMuonsETA2p4_25 =0;
        int nMuonsETA2p4_16 =0;

        for(const auto* m : genMuons){
            if(m->pt() < 5) continue;
            if(m->absEta() < 2.4) nMuonsETA2p4++;
            if(m->absEta() < 2.8) nMuonsETA2p8++;
            if(m->absEta() < 2.4 && m->pt() >= 16) nMuonsETA2p4_16++;
            if(m->absEta() < 2.4 && m->pt() >= 25) nMuonsETA2p4_25++;
        }

        bool passScen1 = nMuonsETA2p4==2 && nMuonsETA2p4_25 >= 1;
        bool passScen2 = nMuonsETA2p8==2 && nMuonsETA2p4_25 >= 1;
        bool passScen3 = nMuonsETA2p8==2 && nMuonsETA2p4_16 >= 1;



        static const double diMuonETAs[] = {0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.0,5.0};
        static const int nBins = 9;
        if(diMuon.pt() > 50){
            plotter.getOrMake1D("pt_gt50_norm",";norm",1,0,2)->Fill(1,reader_event->weight);
                plotter.getOrMake1D("pt_gt50_all_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
            if(passScen1)
                plotter.getOrMake1D("pt_gt50_scen1_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
            if(passScen2)
                plotter.getOrMake1D("pt_gt50_scen2_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
            if(passScen3)
                plotter.getOrMake1D("pt_gt50_scen3_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);

        }
        plotter.getOrMake1D("incl_norm",";norm",1,0,2)->Fill(1,reader_event->weight);

        plotter.getOrMake1D("incl_all_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
    if(passScen1)
        plotter.getOrMake1D("incl_scen1_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
    if(passScen2)
        plotter.getOrMake1D("incl_scen2_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);
    if(passScen3)
        plotter.getOrMake1D("incl_scen3_zeta",";#eta(ll)",nBins,diMuonETAs)->Fill(diMuon.absEta(),reader_event->weight);

    plotter.getOrMake1D("all_zpt",";p_{T}(ll)",500,0,500)->Fill(diMuon.pt(),reader_event->weight);
    plotter.getOrMake1D("all_min_mu_pt",";min p_{T}(l)",500,0,500)->Fill(std::min(genMuons[0]->pt(),genMuons[1]->pt()),reader_event->weight);

//    if(diMuon.absEta() > 3.0 && passScen1)
//        ParticleInfo::printGenInfo(genparticles,-1);


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    EventReader * reader_event = 0;
    GenParticleReader * reader_genParticles = 0;
    HistGetter plotter;

};

#endif

void testGenParticleFiller(std::string fileName ="output.root",std::string outFileName = "plots.root"){
    TString path = fileName;
    Analyzer a(fileName,"treeMaker/Events",path.Contains("data",TString::kIgnoreCase));
    a.analyze();
    a.write(outFileName);
}
