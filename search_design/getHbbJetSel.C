
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        TPRegexp r1(".*(m\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
//        std::cout <<" ->>"<<match->GetSize()  <<std::endl;
//        std::cout <<" ->>"<<((TObjString *)match->At(0))->GetString()  <<std::endl;
//        std::cout <<" ->>"<<((TObjString *)match->At(1))->GetString()  <<std::endl;
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1)
            prefix = ((TObjString *)match->At(1))->GetString();
        else std::cout <<" No pre!"<<std::endl;

        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc);

    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        if(treeType == TREE_OTHER)
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false);  load(reader_fatjet  );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );

        setBranchAddress("skim" ,"ht"         ,   &ht                  ,true);
        setBranchAddress("skim" ,"selLep_pt"  ,   &selLep_pt           ,true);
        setBranchAddress("skim" ,"selLep_eta" ,   &selLep_eta          ,true);
        setBranchAddress("skim" ,"selLep_phi" ,   &selLep_phi          ,true);
        setBranchAddress("skim" ,"selLep_muon",   &selLep_muon         ,true);
    }

    const FatJet * getSignalFJ(const DiHiggsEvent& diHiggsEvt, const MomentumF* lepton, const std::vector<const FatJet*>& fjs){

        double genLepDR = PhysicsUtilities::deltaR(*lepton,*diHiggsEvt.w1_d1);
        if(genLepDR > 0.2) return 0;
        const MomentumF hbb = diHiggsEvt.b1->p4() + diHiggsEvt.b2->p4();
        if(hbb.pt() < 50) return 0;
        if(hbb.absEta() > 2.4) return 0;
        if(PhysicsUtilities::deltaR(*diHiggsEvt.b1,*diHiggsEvt.b2) > 0.8 ) return 0;
        double minDR = 100000;
        int fjIDX = PhysicsUtilities::findNearestDRDeref(hbb,fjs,minDR);
        if(fjIDX < 0) return 0;
        if(minDR > 0.5) return 0;
        return fjs[fjIDX];
    }

    const std::pair<const Jet *, const Jet*> getSignaljets(const DiHiggsEvent& diHiggsEvt, const MomentumF* lepton, const std::vector<Jet*>& jets){
        auto nullResult = std::make_pair((const Jet*)(0),(const Jet*)(0));
        double genLepDR = PhysicsUtilities::deltaR(*lepton,*diHiggsEvt.w1_d1);
        if(genLepDR > 0.2) return nullResult;
        if(PhysicsUtilities::deltaR(*diHiggsEvt.b1,*diHiggsEvt.b2) <= 0.8 ) return nullResult;
        if(diHiggsEvt.b1->pt() < 20 || diHiggsEvt.b2->pt() < 20 || diHiggsEvt.b1->absEta() > 2.4 ||diHiggsEvt.b1->absEta() > 2.4)
            return nullResult;
        double minDR1 = 0;
        int j1IDX = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.b1,jets,minDR1);
        if(j1IDX < 0) return nullResult;
        std::vector<bool> vetoed(jets.size(),false); vetoed[j1IDX] = true;
        double minDR2 = 0;
        int j2IDX = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.b2,jets,minDR2,99999,20,&vetoed);
        if(j2IDX < 0) return nullResult;
        if(minDR2 > 0.4 || minDR1 > 0.4) return nullResult;
        const auto* jet1 = jets[j1IDX];
        const auto* jet2 = jets[j2IDX];
        return std::make_pair(jet1,jet2);
    }
    const std::pair<const Jet *, const Jet*> getJetPair(const MomentumF* lepton, const std::vector<Jet*>& jets){
        std::vector<std::pair<const Jet*,const Jet*> > jetPairs;
        for(unsigned int iJ1 = 0; iJ1 < jets.size(); ++iJ1){
            if(PhysicsUtilities::absDeltaPhi(*jets[iJ1],*lepton) < TMath::PiOver2()) continue;
            for(unsigned int iJ2 = iJ1+1; iJ2 < jets.size(); ++iJ2){
                if(PhysicsUtilities::absDeltaPhi(*jets[iJ2],*lepton) < TMath::PiOver2()) continue;
                if(PhysicsUtilities::deltaR(*jets[iJ1],*jets[iJ2]) > 1.5) continue;
                jetPairs.emplace_back(std::make_pair(jets[iJ1],jets[iJ2]));
            }
        }
        std::vector<std::pair<float,size>> rankedPairPTS(jetPairs.size());
        for(unsigned int iJ = 0; iJ < jetPairs.size(); ++iJ){
            rankedPairPTS[iJ] = std::make_pair((jetPairs[iJ].first->p4()+jetPairs[iJ].second->p4()).pt(),iJ);
        }
        std::sort(rankedPairPTS.begin(), rankedPairPTS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        if(rankedPairPTS.size())
            return jetPairs[rankedPairPTS.front().second];
        else
            return  std::make_pair((const Jet*)(0),(const Jet*)(0));
    }


    void makeHbbPairPlots(TString pre, const std::pair<const Jet*,const Jet*>& jetPair){

        const float pairMass = (jetPair.first->p4()+jetPair.second->p4() ).mass();
        const float minCSV = std::min(jetPair.first->csv(),jetPair.second->csv());
        const float maxCSV = std::max(jetPair.first->csv(),jetPair.second->csv());

        plotter.getOrMake1DPre(pre,"hbb_pair_mass"             ,";hbb pair mass [GeV]; a.u."       ,250,0,500)->Fill(pairMass,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_minsdcsv"         ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(minCSV,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxsdcsv"         ,";hbb pair bb max sj csv; a.u."    ,100,0,1  )->Fill(maxCSV,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxMed_minsdcsv"  ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M] ? minCSV :0.0,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxTight_minsdcsv",";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_T] ? minCSV :0.0,weight );

        bool passMassL = (pairMass > 20);
        bool passMass = (pairMass > 90 && pairMass< 140);
        bool passCSV = (maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M] && minCSV >= BTagging::CSVWP_VALS[BTagging::CSV_L]  );
        bool passCSVT = (minCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M]  );

        if(passMass)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(12.0,weight);
        if(passMass&&passCSV)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(13.0,weight);
        if(passMass && passCSVT)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(14.0,weight);


        if(passMass){
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_minsdcsv"         ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(minCSV,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxsdcsv"         ,";hbb pair bb max sj csv; a.u."    ,100,0,1  )->Fill(maxCSV,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxMed_minsdcsv"  ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M] ? minCSV :0.0,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxTight_minsdcsv",";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_T] ? minCSV :0.0,weight );
        }

        if(passCSV&&passMassL){
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_mass"             ,";hbb pair mass [GeV]; a.u."       ,250,0,500)->Fill(pairMass,weight );
        }
    }




    void makeHbbPlots(TString pre, const FatJet* fj, bool goodSignal){

        plotter.getOrMake1DPre(pre,"hbb_fj_mass"             ,";hbb fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_sd_mass"          ,";hbb fj sd mass [GeV]; a.u."    ,250,0,500)->Fill(fj->sdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_rawsd_mass"       ,";hbb fj raw sd mass [GeV]; a.u.",250,0,500)->Fill(fj->rawSdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_tau2otau1"        ,";hbb fj #tau_{2}/#tau_{1}; a.u.",100,0,1  )->Fill(fj->tau2otau1(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_csv"              ,";hbb fj csv; a.u."              ,100,0,1  )->Fill(fj->csv(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_bbcsv"            ,";hbb fj bb csv; a.u."           ,100,0,1  )->Fill(fj->bbt(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_minsdcsv"         ,";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->minSJCSV(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxsdcsv"         ,";hbb fj bb max sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxMed_minsdcsv"  ,";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_M] ? fj->minSJCSV() :0.0,weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxTight_minsdcsv",";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_T]? fj->minSJCSV() :0.0,weight );

        bool passMassL = (fj->sdMom().mass() > 10);
        bool passMass = (fj->sdMom().mass() > 90 && fj->sdMom().mass() < 140);
        bool passTau = fj->tau2otau1() < 0.55;
        bool passCSV = (fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_M] && fj->minSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_L]  );
        bool passCSVT = (fj->minSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_M]  );

        if(passMassL){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(2.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(7.0,weight);
        }
        if(passTau&&passMassL){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(3.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(8.0,weight);
        }
        if(passTau&&passCSV&&passMassL){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(4.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(9.0,weight);

        }
        if(passTau&& passCSVT&&passMassL){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(5.0,weight);
             if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(10.0,weight);
        }




        if(passMass && passTau){
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_csv"              ,";hbb fj csv; a.u."          ,100,0,1)->Fill(fj->csv(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_bbcsv"            ,";hbb fj bb csv; a.u."       ,100,0,1)->Fill(fj->bbt(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_minsdcsv"         ,";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->minSJCSV(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxsdcsv"         ,";hbb fj bb max sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxMed_minsdcsv"  ,";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_M] ? fj->minSJCSV() :0.0,weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxTight_minsdcsv",";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_T]? fj->minSJCSV() :0.0,weight );
        }
        if(passMass&&passCSV){
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_tau2otau1",";hbb fj #tau_{2}/#tau_{1}; a.u.",100,0,1)->Fill(fj->tau2otau1(),weight );
        }
        if(passTau&&passCSV&&passMassL){
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_mass"             ,";hbb fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_sd_mass"          ,";hbb fj sd mass [GeV]; a.u."    ,250,0,500)->Fill(fj->sdMom().mass(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_rawsd_mass"       ,";hbb fj raw sd mass [GeV]; a.u.",250,0,500)->Fill(fj->rawSdMom().mass(),weight );
        }
        if(passTau&&passCSV&&passMass){
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_rawmass"  ,";hbb fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass()     ,weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_tau3otau2",";hbb fj #tau_{3}/#tau_{2}; a.u.",100,0,1)  ->Fill(fj->tau3otau2(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_tau3otau1",";hbb fj #tau_{3}/#tau_{1}; a.u.",100,0,1)  ->Fill(fj->tau3otau1(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxSJMass",";hbb fj max sj mass [GeV; a.u." ,250,0,500)  ->Fill(std::max(fj->subJet(0).mass(),fj->subJet(1).mass()),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_lowCSVSJMass" ,";hbb fj lowest CSV sj mass[GeV; a.u." ,250,0,500)  ->Fill(fj->subJet(0).csv() < fj->subJet(1).csv()
                    ? fj->subJet(0).mass() :fj->subJet(1).mass()  ,weight );

        }
    }


    bool runEvent() override {
        DiHiggsEvent diHiggsEvt;
        if(!EventSelection::passEventFilters(*reader_event)) return false;

        if(reader_event->process !=
                FillerConstants::SIGNAL){
            prefix = (reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)  ? "rare" :
                    FillerConstants::MCProcessNames[reader_event->process];
        } else {
            diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
            if(diHiggsEvt.type < DiHiggsEvent::MU) return 0;
        }

        MomentumF lepton(ASTypes::CylLorentzVectorF(selLep_pt,selLep_eta,selLep_phi,0));
        weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        if(reader_event->process == FillerConstants::QCD) weight /= 20;

        auto fjs = fjProc.loadFatJets(*reader_fatjet,&lepton);

        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(0.0,weight);

        //Do fat jets
        const auto* fj = fjProc.getHBBCand();
        if(fj){
            bool goodSignal = true;
            if(reader_event->process == FillerConstants::SIGNAL) {
                const auto* sigfj = getSignalFJ(diHiggsEvt,&lepton,fjs);
                if(sigfj != fj) goodSignal = false;
            }

            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(1.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(6.0,weight);

            makeHbbPlots(prefix,fj,goodSignal);
        }
        //do unmerged jets
        auto jets = JetKinematics::selectObjects(reader_jet->jets,20,2.4);
        const auto jetPair = getJetPair(&lepton,jets);
        if(jetPair.first && jetPair.second){
            bool process = true;
            if(reader_event->process == FillerConstants::SIGNAL) {
                const auto sigfj = getSignaljets(diHiggsEvt,&lepton,jets);
                if(sigfj.first != jetPair.first && sigfj.first != jetPair.second ) process = false;
                if(sigfj.second != jetPair.first && sigfj.second != jetPair.second ) process = false;
            }
            if(process){
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(11.0,weight);
                makeHbbPairPlots(prefix,jetPair);
            }
        }




        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    std::shared_ptr<JetReader        > reader_jet      ;
    float   ht         =0;
    float   selLep_pt  =0;
    float   selLep_eta =0;
    float   selLep_phi =0;
    size8   selLep_muon=0;
    float weight  = 0;
    HistGetter plotter;
    TString prefix;

    FatJetProcessor fjProc;


};

#endif

void getHbbJetSel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getHbbJetSel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
