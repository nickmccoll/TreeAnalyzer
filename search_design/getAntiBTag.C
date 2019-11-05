
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TPRegexp.h"


using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) :
        DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed)
{
}

    void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.event.doTTBarStitching = false;
    }

    std::vector<float> DeepJET_WP    = {-100,0.0521,0.3033,0.7489};
    std::vector<float> CSV_WP        = {-100,0.5803,0.8838,0.9693};
    std::vector<float> DeepCSV_WP    = {-100,0.1522,0.4941,0.8001};


    bool runEvent() override {
        if(! DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(lepChan != SINGLELEP) return false;
        if(isSignal() && diHiggsEvt.type < DiHiggsEvent::TAU_MU ) return false;

        weight = EventWeights::getNormalizedEventWeight(
                *reader_event,xsec(),nSampEvt(),parameters.event,smDecayEvt.genMtt,smDecayEvt.nLepsTT);

        plotter.getOrMake1DPre(smpName,"nEvts",";# of events",1,0,2)->Fill(1,weight);


        auto passWP = [&](int bType, const Jet* j, const BTagging::BTAGWP wp) -> bool{
            if(bType==0){
                return j->csv() >= CSV_WP[wp];
            } else if(bType==1){
                return j->deep_csv() >= DeepCSV_WP[wp];
            }
            return j->deep_flavor() >= DeepJET_WP[wp];
        };

        auto getDisc = [&](int bType, const Jet* j) -> double{
            if(bType==0){
                return j->csv() ;
            } else if(bType==1){
                return j->deep_csv() ;
            }
            return j->deep_flavor() ;
        };

        auto getDiscName = [&](int bType) -> TString{
            if(bType==0){
                return "csv" ;
            } else if(bType==1){
                return "deep_csv" ;
            }
            return "deep_flavor" ;
        };


        auto jetPlts = [&](const TString& prefix, int bType,
                const std::vector<const Jet*>& jets){
            int nJ = 0;
            int nL = 0;
            int nM = 0;
            int nT = 0;

            int nJ3 = 0;
            int nL3 = 0;
            int nM3 = 0;
            int nT3 = 0;

            for(const auto * j: jets){
                plotter.getOrMake1DPre(prefix,"bIncl_jetPt",";jet #it{p}_T",60,0,300)->Fill(j->pt(),weight);
                if(passWP(bType,j,BTagging::BTAG_L))
                    plotter.getOrMake1DPre(prefix,"bL_jetPt",";jet #it{p}_T",60,0,300)->Fill(j->pt(),weight);
                if(passWP(bType,j,BTagging::BTAG_M))
                    plotter.getOrMake1DPre(prefix,"bM_jetPt",";jet #it{p}_T",60,0,300)->Fill(j->pt(),weight);
                if(passWP(bType,j,BTagging::BTAG_T))
                    plotter.getOrMake1DPre(prefix,"bT_jetPt",";jet #it{p}_T",60,0,300)->Fill(j->pt(),weight);

                if(j->pt()>20){
                    plotter.getOrMake1DPre(prefix,"pt20_disc",";"+getDiscName(bType),100,0,1)->Fill(getDisc(bType,j),weight);
                    nJ++;
                    if(passWP(bType,j,BTagging::BTAG_L)) nL++;
                    if(passWP(bType,j,BTagging::BTAG_M)) nM++;
                    if(passWP(bType,j,BTagging::BTAG_T)) nT++;
                }

                if(j->pt()>30){
                    plotter.getOrMake1DPre(prefix,"pt30_disc",";"+getDiscName(bType),100,0,1)->Fill(getDisc(bType,j),weight);
                    nJ3++;
                    if(passWP(bType,j,BTagging::BTAG_L)) nL3++;
                    if(passWP(bType,j,BTagging::BTAG_M)) nM3++;
                    if(passWP(bType,j,BTagging::BTAG_T)) nT3++;
                }
            }
            plotter.getOrMake1DPre(prefix,"pt20_nJ",";# of Jets",10,-0.5,9.5)->Fill(nJ,weight);
            plotter.getOrMake1DPre(prefix,"pt20_nL",";# of Jets",10,-0.5,9.5)->Fill(nL,weight);
            plotter.getOrMake1DPre(prefix,"pt20_nM",";# of Jets",10,-0.5,9.5)->Fill(nM,weight);
            plotter.getOrMake1DPre(prefix,"pt20_nT",";# of Jets",10,-0.5,9.5)->Fill(nT,weight);

            plotter.getOrMake1DPre(prefix,"pt30_nJ",";# of Jets",10,-0.5,9.5)->Fill(nJ3,weight);
            plotter.getOrMake1DPre(prefix,"pt30_nL",";# of Jets",10,-0.5,9.5)->Fill(nL3,weight);
            plotter.getOrMake1DPre(prefix,"pt30_nM",";# of Jets",10,-0.5,9.5)->Fill(nM3,weight);
            plotter.getOrMake1DPre(prefix,"pt30_nT",";# of Jets",10,-0.5,9.5)->Fill(nT3,weight);


        };

        auto doFlvs  = [&](const TString& prefix,const std::vector<const Jet*>& jets){

            jetPlts(prefix +"_csv", 0,jets);
            jetPlts(prefix +"_deep_csv", 1,jets);
            jetPlts(prefix +"_deep_flavor", 2,jets);
        };


        doFlvs(smpName, PhysicsUtilities::selObjsMom(reader_jet->jets,0.0,2.4));
        doFlvs(smpName+"_bF", PhysicsUtilities::selObjsMom(reader_jet->jets,0.0,2.4,
                [&](const Jet* j){ return BTagging::jetFlavor(*j) == BTagging::FLV_B;}));
        doFlvs(smpName+"_cF", PhysicsUtilities::selObjsMom(reader_jet->jets,0.0,2.4,
                [&](const Jet* j){ return BTagging::jetFlavor(*j) == BTagging::FLV_C;}));
        doFlvs(smpName+"_lF", PhysicsUtilities::selObjsMom(reader_jet->jets,0.0,2.4,
                [&](const Jet* j){ return BTagging::jetFlavor(*j) != BTagging::FLV_B && BTagging::jetFlavor(*j) != BTagging::FLV_C ;}));



        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;



};

#endif

void getAntiBTag(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
