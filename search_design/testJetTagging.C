
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "DataFormats/interface/FatJet.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_JER);

        fjProcNoSJ      .reset(new FatJetProcessor ());
    }

    void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.jets.minJetPT = 10;
        parameters.event.passTrigger = &EventSelection::alwaysTrue;


    };
    void testSubjets(const std::string& prefix){
        auto noSJFatJetParams = parameters.fatJets;
        noSJFatJetParams.wjj_minSJs =0;
        noSJFatJetParams.hbb_minSJs =0;

        fjProcNoSJ->loadFatJets(noSJFatJetParams,
                *reader_fatjet,*reader_fatjet_noLep,selectedLepton);
        auto hbbCand     = fjProcNoSJ->getHBBCand();
        auto wjjCand     = fjProcNoSJ->getWjjCand();
        if(!hbbCand || !wjjCand) return;

        auto neutrino    = HiggsSolver::getInvisible(reader_event->met,
                (selectedLepton->p4() + wjjCand->p4()) );
        auto wlnu        =  neutrino.p4() + selectedLepton->p4();
        auto hWW         = wlnu + wjjCand->p4();
        auto hh =  (hWW + hbbCand->p4()) ;
        const std::string hhMassString
        = hh.mass() > 2000 ? "hhMass_gt2000" :
                (hh.mass() > 1000 ? "hhMass_1000to2000" : "hhMass_lt1000");


        auto getNSubjets=[](const FatJet* fj, float minPT) -> int {
            return PhysicsUtilities::selObjsMom(fj->subJets(),
                    minPT, 2.4).size();
        };

        auto mkPlots = [&](const std::string& prefix, const float minPT){
            int nHbbSJ = getNSubjets(hbbCand,minPT);
            int nWjjSJ = getNSubjets(wjjCand,minPT);

            plotter.getOrMake1DPre((prefix + "_"+ "hh_incl").c_str()
                    ,"nHbbSJ",";#it{N.} of b#bar{b} subjets",4,-0.5,3.5)->Fill(nHbbSJ,weight);
            plotter.getOrMake1DPre((prefix + "_"+ hhMassString).c_str()
                    ,"nHbbSJ",";#it{N.} of b#bar{b} subjets",4,-0.5,3.5)->Fill(nHbbSJ,weight);
            plotter.getOrMake1DPre((prefix + "_"+ "hh_incl").c_str()
                    ,"nWjjSJ",";#it{N.} of q#bar{q} subjets",4,-0.5,3.5)->Fill(nWjjSJ,weight);
            plotter.getOrMake1DPre((prefix + "_"+ hhMassString).c_str()
                    ,"nWjjSJ",";#it{N.} of q#bar{q} subjets",4,-0.5,3.5)->Fill(nWjjSJ,weight);
        };
        mkPlots(prefix +"_minSJPT_0",0);
        mkPlots(prefix +"_minSJPT_10",10);
        mkPlots(prefix +"_minSJPT_20",20);
        mkPlots(prefix +"_minSJPT_30",30);
    }

    void testWjjDisc(const std::string& prefix){
        plotter.getOrMake1DPre((prefix).c_str()
                ,"tau2otau1",";#tau_{2}/#tau_{1}",100,0,1)->Fill(wjjCand->tau2otau1(),weight);
        plotter.getOrMake1DPre((prefix).c_str()
                ,"N2_b1",";#it{N}_{2} (#beta=1)",101,-0.004,0.4)->Fill(wjjCand->ecfb1(),weight);
        plotter.getOrMake1DPre((prefix).c_str()
                ,"N2_b2",";#it{N}_{2} (#beta=2)",101,-0.004,0.4)->Fill(wjjCand->ecfb2(),weight);
    }

    void testAK4Btagging(const std::string& prefix){
        for(const auto* j : jets_HbbV){
            const std::string flvStr = BTagging::jetFlavorString(BTagging::jetFlavor(*j));
            std::string ptStr = "pt_gt200";
            if(j->pt() < 20){
                ptStr = "pt_10to20";
            } else if(j->pt() < 30){
                ptStr = "pt_20to30";
            } else if(j->pt() < 50){
                ptStr = "pt_30to50";
            } else if(j->pt() < 100){
                ptStr = "pt_50to100";
            } else if(j->pt() < 200){
                ptStr = "pt_1000to200";
            }
            plotter.getOrMake1DPre((prefix + "_"+flvStr+"_"+ ptStr).c_str()
                    ,"csv",";AK4 jet CSV",100,0,1)->Fill(j->csv(),weight);
            plotter.getOrMake1DPre((prefix + "_"+flvStr+"_"+ ptStr).c_str()
                    ,"deepcsv",";AK4 jet deep CSV",100,0,1)->Fill(j->deep_csv(),weight);
        }
        const int nMCSV20 = PhysicsUtilities::selObjsMomD(jets_HbbV,20,2.4,
                [&](const Jet* j){return j->csv() >= parameters.jets.CSV_WP[BTagging::BTAG_M];}
        ).size();
        const int nMCSV30 = PhysicsUtilities::selObjsMomD(jets_HbbV,30,2.4,
                [&](const Jet* j){return j->csv() >= parameters.jets.CSV_WP[BTagging::BTAG_M];}
        ).size();
        const int nLCSV20 = PhysicsUtilities::selObjsMomD(jets_HbbV,20,2.4,
                [&](const Jet* j){return j->csv() >= parameters.jets.CSV_WP[BTagging::BTAG_L];}
        ).size();
        const int nLCSV30 = PhysicsUtilities::selObjsMomD(jets_HbbV,30,2.4,
                [&](const Jet* j){return j->csv() >= parameters.jets.CSV_WP[BTagging::BTAG_L];}
        ).size();

        const int nMDCSV20 = PhysicsUtilities::selObjsMomD(jets_HbbV,20,2.4,
                [&](const Jet* j){return j->deep_csv() >= parameters.jets.DeepCSV_WP[BTagging::BTAG_M];}
        ).size();
        const int nMDCSV30 = PhysicsUtilities::selObjsMomD(jets_HbbV,30,2.4,
                [&](const Jet* j){return j->deep_csv() >= parameters.jets.DeepCSV_WP[BTagging::BTAG_M];}
        ).size();
        const int nLDCSV20 = PhysicsUtilities::selObjsMomD(jets_HbbV,20,2.4,
                [&](const Jet* j){return j->deep_csv() >= parameters.jets.DeepCSV_WP[BTagging::BTAG_L];}
        ).size();
        const int nKDCSV30 = PhysicsUtilities::selObjsMomD(jets_HbbV,30,2.4,
                [&](const Jet* j){return j->deep_csv() >= parameters.jets.DeepCSV_WP[BTagging::BTAG_L];}
        ).size();

        auto mkPlt = [&](const std::string prefix, const int n20, const int n30){
            plotter.getOrMake1DPre((prefix+"pt_20to30").c_str(),"nAK4Btags",
                    ";#it{N.} of AK4 b tags (20-30 GeV)",4,-0.5,3.5)->Fill(n20-n30,weight);
            plotter.getOrMake1DPre((prefix+"_pt_gt20").c_str(),"nAK4Btags",
                    ";#it{N.} of AK4 b tags (>20 GeV)",4,-0.5,3.5)->Fill(n20,weight);
            plotter.getOrMake1DPre((prefix+"_pt_gt30").c_str(),"nAK4Btags",
                    ";#it{N.} of AK4 b tags (>30 GeV)",4,-0.5,3.5)->Fill(n30,weight);
        };
        mkPlt(prefix+"_"+"mCSV",nMCSV20,nMCSV30);
        mkPlt(prefix+"_"+"lCSV",nLCSV20,nLCSV30);
        mkPlt(prefix+"_"+"lnmCSV",nLCSV20-nMCSV20,nLCSV30-nMCSV30);
    };

    void testAK8Btagging(const std::string& prefix){
        auto dCSVParams30 = parameters.jets;
        dCSVParams30.minBtagJetPT = 30;
        dCSVParams30.getSJBTagVal =&BaseRecoJet::deep_csv;
        dCSVParams30.sjBTagLWP = parameters.jets.DeepCSV_WP [BTagging::BTAG_L];
        dCSVParams30.sjBTagMWP = parameters.jets.DeepCSV_WP [BTagging::BTAG_M];
        auto dCSVParams20 = dCSVParams30;
        dCSVParams20.minBtagJetPT = 20;

        auto csvParams30 = parameters.jets;
        csvParams30.minBtagJetPT = 30;
        csvParams30.getSJBTagVal =&BaseRecoJet::csv;
        csvParams30.sjBTagLWP = parameters.jets.CSV_WP [BTagging::BTAG_L];
        csvParams30.sjBTagMWP = parameters.jets.CSV_WP [BTagging::BTAG_M];
        auto csvParams20 = csvParams30;
        csvParams20.minBtagJetPT = 20;

        auto csv20Cat  = BTagging::getCSVSJCat(csvParams20,hbbCand->subJets());
        auto csv30Cat  = BTagging::getCSVSJCat(csvParams30,hbbCand->subJets());
        auto dCSV20Cat = BTagging::getCSVSJCat(dCSVParams20,hbbCand->subJets());
        auto dCSV30Cat = BTagging::getCSVSJCat(dCSVParams30,hbbCand->subJets());

        //        if(hbbCand->subJets().size() > 1){
        //            std::cout <<  hbbCand->subJets()[0].csv() <<" - " <<  hbbCand->subJets()[0].deep_csv()
        //                    <<" - "<<hbbCand->subJets()[1].csv() <<" - " <<  hbbCand->subJets()[1].deep_csv()
        //                    << " -> "<< BTagging::passSubjetBTagLWP(dCSVParams30,hbbCand->subJets()[0])
        //                    <<","<< BTagging::passSubjetBTagMWP(dCSVParams30,hbbCand->subJets()[0])
        //            << ","<< BTagging::passSubjetBTagLWP(dCSVParams30,hbbCand->subJets()[1])
        //            <<","<< BTagging::passSubjetBTagMWP(dCSVParams30,hbbCand->subJets()[1])
        //                    <<","<< dCSV30Cat<<std::endl;
        //        }


        plotter.getOrMake1DPre((prefix).c_str(),"csv20Cat",
                ";hbbCat",6,-0.5,5.5)->Fill(csv20Cat-1,weight);
        plotter.getOrMake1DPre((prefix).c_str(),"csv30Cat",
                ";hbbCat",6,-0.5,5.5)->Fill(csv30Cat-1,weight);
        plotter.getOrMake1DPre((prefix).c_str(),"dCSV20Cat",
                ";hbbCat",6,-0.5,5.5)->Fill(dCSV20Cat-1,weight);
        plotter.getOrMake1DPre((prefix).c_str(),"dCSV30Cat",
                ";hbbCat",6,-0.5,5.5)->Fill(dCSV30Cat-1,weight);

        int bbtCat = 4;
        if(hbbCand->bbt() < 0.3 ) bbtCat = 0;
        else if(hbbCand->bbt() < 0.6 ) bbtCat = 1;
        else if(hbbCand->bbt() < 0.8 ) bbtCat = 2;
        else if(hbbCand->bbt() < 0.9 ) bbtCat = 3;

        plotter.getOrMake1DPre((prefix).c_str(),"bbtCat",
                ";hbbCat",5,-0.5,4.5)->Fill(bbtCat,weight);
    };

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;

        LeptonParameters paramsNoID = parameters.leptons;
        paramsNoID.el_getID = &Electron::passInclID;
        auto electrons = LeptonProcessor::getElectrons(paramsNoID,*reader_electron);

        for(const auto* e : electrons){
            bool passTight = e->passTightID_noIso();
            bool passMed   = e->passMedID_noIso();
            bool passMVA90 =e->passMVA90ID_noIso();

            auto fillBin = [&](const float bin){
                plotter.getOrMake1DPre(smpName,"leptonSel",";leptonSel",9,-0.5,8.5)->Fill(bin,weight);
            };

            if(passMed)fillBin(0);
            if(passTight)fillBin(1);
            if(passMVA90)fillBin(2);

            if(passMed&&passMVA90)fillBin(3);
            if(!passMed&&passMVA90)fillBin(4);
            if(passMed&&!passMVA90)fillBin(5);
            if(passTight&&passMVA90)fillBin(6);
            if(!passTight&&passMVA90)fillBin(7);
            if(passTight&&!passMVA90)fillBin(8);

        }





      if(selectedLepton && selectedLepton->isElectron()){
          const Electron * e = (const Electron*)selectedLepton;
          std::cout << e->p4() <<" -> "<< e->passMedID_noIso() <<" "<< e->passTightID_noIso()<< " "<< e->passMVA90ID_noIso()
                 <<" :: "<<e->miniIso()<<" "<<e->miniIsoFP()   <<" :: ";
                 if(wjjCand) std::cout << wjjCand->p4()<<std::endl;
                 else std::cout <<std::endl;
          for(const auto& j : reader_fatjet->jets){
              std::cout << j.p4()<<std::endl;
          }
      }

        if(!passTriggerPreselection) return false;
        if(!selectedLepton) return false;
        std::string prefix = smpName.Data();

        testSubjets(prefix);

        if(!hbbCand || !wjjCand) return false;
        const std::string hhMassString
        = hh.mass() > 2000 ? "hhMass_gt2000" :
                (hh.mass() > 1000 ? "hhMass_1000to2000" : "hhMass_lt1000");

        if(hh.mass() > 700){
            testWjjDisc(prefix+"_"+"hhMass_gt700");
            testAK4Btagging(prefix+"_"+"hhMass_gt700");
            testAK8Btagging(prefix+"_"+"hhMass_gt700");
        }
        testWjjDisc(prefix+"_"+hhMassString);
        testAK4Btagging(prefix+"_"+hhMassString);
        testAK8Btagging(prefix+"_"+hhMassString);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<FatJetProcessor>        fjProcNoSJ     ;


};

#endif

void testJetTagging(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
