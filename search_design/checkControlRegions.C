
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
#include "Processors/Variables/interface/HiggsSolver.h"

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1){
            mass = (((TObjString *)match->At(1))->GetString()).Atoi();
            prefix = "m" + ((TObjString *)match->At(1))->GetString();
        }
        else std::cout <<" No pre!"<<std::endl;

        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc);
        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc_inclWBtag);
        fjProc_inclWBtag.wjj_maxCSVWP =BTagging::CSV_INCL;
        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc_inclHbbBtag);
        fjProc_inclHbbBtag.hbb_l_firMinCSVWP =BTagging::CSV_INCL;
        fjProc_inclHbbBtag.hbb_l_secMinCSVWP =BTagging::CSV_INCL;
    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        if(treeType == TREE_OTHER){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");                load(reader_genpart   );
        }
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false);  load(reader_fatjet  );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );

        setBranchAddress("skim" ,"ht"         ,   &ht                  ,true);
        setBranchAddress("skim" ,"selLep_pt"  ,   &selLep_pt           ,true);
        setBranchAddress("skim" ,"selLep_eta" ,   &selLep_eta          ,true);
        setBranchAddress("skim" ,"selLep_phi" ,   &selLep_phi          ,true);
        setBranchAddress("skim" ,"selLep_muon",   &selLep_muon         ,true);
    }

    void makeCRPlots(TString pre, const MomentumF& hWW, const MomentumF& hbb){
        auto hh = hWW.p4() + hbb.p4();
        plotter.getOrMake1DPre(pre, "hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        if(hbb.mass() < 90)
            plotter.getOrMake1DPre(pre, "hbblt90_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        else if(hbb.mass() < 140)
            plotter.getOrMake1DPre(pre, "hbb90to140_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        else
            plotter.getOrMake1DPre(pre, "hbbgeq140_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);


        plotter.getOrMake1DPre(pre, "hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 700 && hh.mass() < 900)
            plotter.getOrMake1DPre(pre, "hh700to900_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 900 && hh.mass() < 1100)
            plotter.getOrMake1DPre(pre, "hh900to1100_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 1400 && hh.mass() < 1800)
            plotter.getOrMake1DPre(pre, "hh1400to1800_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
    }


    bool runEvent() override {
        if(reader_event->process !=
                FillerConstants::SIGNAL){
            prefix = (reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)  ? "rare" :
                    FillerConstants::MCProcessNames[reader_event->process];
        }
        if(ht < 500) return false;
        if(!EventWeights::passEventFilters(*reader_event)) return false;


        MomentumF lepton(ASTypes::CylLorentzVectorF(selLep_pt,selLep_eta,selLep_phi,0));
        weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
//        if(reader_event->process == FillerConstants::SIGNAL) weight *=  (0.8241887906 * EventWeights::get4bXSecLimit(mass)/1000.0);
        if(reader_event->process == FillerConstants::SIGNAL) weight *=  20.0/1000.0;

        fjProc.loadFatJets(*reader_fatjet,&lepton);
        fjProc_inclHbbBtag.loadFatJets(*reader_fatjet,&lepton);
        fjProc_inclWBtag.loadFatJets(*reader_fatjet,&lepton);
        const auto* hbbfj = fjProc.getHBBCand();
        const auto* wjjfj = fjProc.getWjjCand();
        const bool goodHBBFJ  = fjProc.passHbbSel();
        const bool goodHBBFJT = fjProc.passHbbSelTightBTag();
        const bool goodWJJFJ  = fjProc.passWjjSel();

        if(!wjjfj || !hbbfj) return false;
        auto neutrino = HiggsSolver::getInvisible(reader_event->met,(lepton.p4() + wjjfj->sdMom().p4()) );
        MomentumF Wlnu = lepton.p4() + neutrino.p4();
        MomentumF hWW = Wlnu.p4() + wjjfj->sdMom().p4();
        if(PhysicsUtilities::deltaR(Wlnu, wjjfj->sdMom().p4()) > 0.5) return false;

        if(goodWJJFJ && goodHBBFJ) makeCRPlots(prefix+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
        if(goodWJJFJ && goodHBBFJT) makeCRPlots(prefix+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());
        if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() >= BTagging::BBT_VALS[BTagging::CSV_M] && hbbfj->minSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L]  )
            makeCRPlots(prefix+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
        if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L])
            makeCRPlots(prefix+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
        if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(prefix+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
        if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(prefix+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());

        if(selLep_muon==1){
            TString tempP = prefix + "_mu";
            if(goodWJJFJ && goodHBBFJ) makeCRPlots(tempP+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ && goodHBBFJT) makeCRPlots(tempP+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() >= BTagging::BBT_VALS[BTagging::CSV_M] && hbbfj->minSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L]  )
                makeCRPlots(tempP+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L])
                makeCRPlots(tempP+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(tempP+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(tempP+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());
        } else {
            TString tempP = prefix + "_el";
            if(goodWJJFJ && goodHBBFJ) makeCRPlots(tempP+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ && goodHBBFJT) makeCRPlots(tempP+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() >= BTagging::BBT_VALS[BTagging::CSV_M] && hbbfj->minSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L]  )
                makeCRPlots(tempP+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && hbbfj->maxSJCSV() < BTagging::BBT_VALS[BTagging::CSV_L])
                makeCRPlots(tempP+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(tempP+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && wjjfj->maxSJCSV()>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(tempP+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());
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
    int mass=0;

    FatJetProcessor fjProc           ;
    FatJetProcessor fjProc_inclWBtag      ;
    FatJetProcessor fjProc_inclHbbBtag ;



};

#endif

void checkControlRegions(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setLumi(36);
    a.analyze();
    a.write(outFileName);
}
void checkControlRegions(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.setLumi(36);
    a.analyze();
    a.write(outFileName);
}
