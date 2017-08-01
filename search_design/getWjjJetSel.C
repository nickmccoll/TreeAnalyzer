
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
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/FatJetSelection.h"

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        TPRegexp r1(".*(m\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1)
            prefix = ((TObjString *)match->At(1))->GetString();
        else std::cout <<" No pre!"<<std::endl;

        fjProc     .reset(new FatJetProcessor ()); DefaultFatJetSelections::setDefaultFatJetProcessor(*fjProc);

    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        if(treeType == TREE_OTHER)
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData());  load(reader_jet     );

        setBranchAddress("skim" ,"ht"         ,   &ht                  ,true);
        setBranchAddress("skim" ,"selLep_pt"  ,   &selLep_pt           ,true);
        setBranchAddress("skim" ,"selLep_eta" ,   &selLep_eta          ,true);
        setBranchAddress("skim" ,"selLep_phi" ,   &selLep_phi          ,true);
        setBranchAddress("skim" ,"selLep_muon",   &selLep_muon         ,true);
    }

    const FatJet * getSignalFJ(const DiHiggsEvent& diHiggsEvt, const MomentumF* lepton, const std::vector<const FatJet*>& fjs){
        double genLepDR = PhysicsUtilities::deltaR(*lepton,*diHiggsEvt.w1_d1);
        if(genLepDR > 0.2) return 0;
        const MomentumF wjj = diHiggsEvt.w2_d1->p4() + diHiggsEvt.w2_d2->p4();
        if(wjj.pt() < 50) return 0;
        if(wjj.absEta() > 2.4) return 0;
        if(PhysicsUtilities::deltaR(*diHiggsEvt.w2_d1,*diHiggsEvt.w2_d2) > 0.8 ) return 0;
        double minDR = 100000;
        int fjIDX = PhysicsUtilities::findNearestDRDeref(wjj,fjs,minDR);
        if(fjIDX < 0) return 0;
        if(minDR > 0.5) return 0;
        return fjs[fjIDX];
    }


    void makeWjjPlots(TString pre, const FatJet* fj, bool goodSignal){

        plotter.getOrMake1DPre(pre,"wjj_fj_mass"             ,";wjj fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_sd_mass"          ,";wjj fj sd mass [GeV]; a.u."    ,250,0,500)->Fill(fj->sdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_rawsd_mass"       ,";wjj fj raw sd mass [GeV]; a.u.",250,0,500)->Fill(fj->rawSdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_tau2otau1"        ,";wjj fj #tau_{2}/#tau_{1}; a.u.",100,0,1  )->Fill(fj->tau2otau1(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_csv"              ,";wjj fj csv; a.u."              ,100,0,1  )->Fill(fj->csv(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_minsdcsv"         ,";wjj fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->minSJCSV(),weight );
        plotter.getOrMake1DPre(pre,"wjj_fj_maxsdcsv"         ,";wjj fj bb max sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV(),weight );


        bool passMass = (fj->sdMom().mass() > 10);
        bool passTau = fj->tau2otau1() < 0.55;
        bool passCSV = (fj->maxSJCSV() < BTagging::CSVWP_VALS[BTagging::CSV_M]);

        if(passMass){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(2.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(6.0,weight);

        }
        if(passMass&&passTau){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(3.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(7.0,weight);

        }
        if(passMass&&passTau&&passCSV){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(4.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(8.0,weight);

        }

        if(passMass && passTau){
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_csv"              ,";wjj fj csv; a.u."          ,100,0,1)->Fill(fj->csv(),weight );
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_minsdcsv"         ,";wjj fj bb min sj csv; a.u.",100,0,1)->Fill(fj->minSJCSV(),weight );
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_maxsdcsv"         ,";wjj fj bb max sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV(),weight );
        }
        if(passMass&&passCSV){
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_tau2otau1",";wjj fj #tau_{2}/#tau_{1}; a.u.",100,0,1)->Fill(fj->tau2otau1(),weight );
        }
        if(passTau&&passCSV){
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_mass"             ,";wjj fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass(),weight );
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_sd_mass"          ,";wjj fj sd mass [GeV]; a.u."    ,250,0,500)->Fill(fj->sdMom().mass(),weight );
            plotter.getOrMake1DPre(pre,"wjj_fj_oM_rawsd_mass"       ,";wjj fj raw sd mass [GeV]; a.u.",250,0,500)->Fill(fj->rawSdMom().mass(),weight );
        }
    }


    bool runEvent() override {
        if(!EventWeights::passEventFilters(*reader_event)) return false;

        DiHiggsEvent diHiggsEvt;

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


        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(0.0,weight);

        //Do fat jets

        auto fjs = fjProc->loadFatJets(*reader_fatjet,&lepton);
        const auto* wjjfj = fjProc->getWjjCand();
        if(wjjfj){
            bool goodSignal = true;
            if(reader_event->process == FillerConstants::SIGNAL) {
                const auto* sigfj = getSignalFJ(diHiggsEvt,&lepton,fjs);
                if(sigfj != wjjfj) goodSignal = false;
            }
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(1.0,weight);
            if(goodSignal) plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(5.0,weight);
            makeWjjPlots(prefix,wjjfj, goodSignal);
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

    std::unique_ptr<FatJetProcessor> fjProc     ;


};

#endif

void getWjjJetSel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getWjjJetSel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
