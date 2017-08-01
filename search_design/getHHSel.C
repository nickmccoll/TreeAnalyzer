
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
        fjProc.hbb_minMass = 90;
        fjProc.hbb_maxMass = 140;
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

    void makeHHPlots(TString pre, const float ht,  const MomentumF& hbb, const bool passTightCSV, const MomentumF& wjj,const MomentumF& wjjNoSD, const MomentumF& lepton, const MomentumF& met){
        auto neutrino =HiggsSolver::getInvisible(met,lepton.p4() + wjj.p4() );
        MomentumF Wlnu = lepton.p4()+ neutrino.p4();
        MomentumF hWW = Wlnu.p4()+ wjj.p4();
        MomentumF hh = hWW.p4() + hbb.p4();
        if(PhysicsUtilities::deltaR(Wlnu, wjj) > 0.5) return;
        auto makeMPlots =[&](const TString& pre) {
            plotter.getOrMake1DPre(pre, "lepW_pt"  ,";lep W pT [GeV]; arbitrary units",200,0,1000)->Fill((met.p4()+lepton.p4() ).pt(),weight);
            plotter.getOrMake1DPre(pre, "met"      ,";met [GeV]; arbitrary units",200,0,1000)->Fill(met.pt(),weight);
            plotter.getOrMake1DPre(pre, "met_o_fjlep" ,";met/fj+lep [GeV]; arbitrary units",200,0,10)->Fill(met.pt()/(wjj.p4()+lepton.p4()).pt(),weight);
            plotter.getOrMake1DPre(pre, "met_o_h" ,";met/h; arbitrary units",200,0,10)->Fill(met.pt()/hWW.pt(),weight);
            plotter.getOrMake1DPre(pre, "met_o_fj" ,";met/fj [GeV]; arbitrary units",200,0,10)->Fill(met.pt()/wjj.pt(),weight);
            plotter.getOrMake1DPre(pre, "met_o_fjNoSD" ,";met/fj [GeV]; arbitrary units",200,0,10)->Fill(met.pt()/wjjNoSD.pt(),weight);
            plotter.getOrMake1DPre(pre, "met_dPhifj" ,";#Delta#phi(met,fj); arbitrary units",160,0,4)->Fill(PhysicsUtilities::absDeltaPhi(met,wjj),weight);
            plotter.getOrMake1DPre(pre, "met_dPhifjNoSD" ,";#Delta#phi(met,fj); arbitrary units",160,0,4)->Fill(PhysicsUtilities::absDeltaPhi(met,wjjNoSD),weight);
            plotter.getOrMake1DPre(pre, "met_dPhihbb" ,";#Delta#phi(met,hbb); arbitrary units",160,0,4)->Fill(PhysicsUtilities::absDeltaPhi(met,hbb),weight);
            if(PhysicsUtilities::absDeltaPhi(met,hbb) > 0.8)
                plotter.getOrMake1DPre(pre, "highDPhi_lepW_pt"  ,";lep W pT [GeV]; arbitrary units",200,0,1000)->Fill((met.p4()+lepton.p4() ).pt(),weight);
            plotter.getOrMake1DPre(pre, "lepPT"    ,";lep pt [GeV]; arbitrary units",200,0,1000)->Fill(lepton.pt(),weight);
            plotter.getOrMake1DPre(pre, "hWW_mass",";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(hWW.mass(),weight);
            plotter.getOrMake1DPre(pre, "hWW_pt"  ,";hWW #it{p}_{T} [GeV]; arbitrary units",250,0,2500)->Fill(hWW.pt(),weight);
            plotter.getOrMake1DPre(pre, "hh_mass" ,";hh mass [GeV]; arbitrary units",500,0,5000)->Fill(hh.mass(),weight);
            plotter.getOrMake1DPre(pre, "W_W_dR"  ,";#DeltaR(W,W); arbitrary units",50,0,5.0)->Fill( PhysicsUtilities::deltaR(Wlnu, wjj),weight);
        };

        makeMPlots(pre);
        if(passTightCSV)
            makeMPlots(pre+"_tCSV");
        else
            makeMPlots(pre+"_lCSV");
        if(ht >= 500)
            makeMPlots(pre+"_hHT");
        else if( ht < 500 && lepton.pt() >= 30)
            makeMPlots(pre+"_lHT");

        if(passTightCSV && ht >= 500 )
            makeMPlots(pre+"_hHT_tCSV");
        if(passTightCSV &&  ht < 500 && lepton.pt() >= 30 )
            makeMPlots(pre+"_lHT_tCSV");

        if(!passTightCSV && ht >= 500 )
            makeMPlots(pre+"_hHT_lCSV");
        if(!passTightCSV &&  ht < 500 && lepton.pt() >= 30 )
            makeMPlots(pre+"_lHT_lCSV");
    }


    bool runEvent() override {
//        DiHiggsEvent diHiggsEvt;

        if(reader_event->process !=
                FillerConstants::SIGNAL){
            prefix = (reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)  ? "rare" :
                    FillerConstants::MCProcessNames[reader_event->process];
        }
//        else {
//            diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
//            if(diHiggsEvt.type < DiHiggsEvent::MU) return 0;
//        }

        MomentumF lepton(ASTypes::CylLorentzVectorF(selLep_pt,selLep_eta,selLep_phi,0));
        weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        if(reader_event->process == FillerConstants::SIGNAL) weight *=  (0.8241887906 * EventWeights::get4bXSecLimit(mass)/1000.0);
        fjProc.loadFatJets(*reader_fatjet,&lepton);
        const auto* hbbfj = fjProc.getHBBCand();
        const auto* wjjfj = fjProc.getWjjCand();
        const bool goodHBBFJ  = fjProc.passHbbSel();
        const bool goodHBBFJT = fjProc.passHbbSelTightBTag();
        const bool goodWJJFJ  = fjProc.passWjjSel();

        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(0.0,weight);
        if(!EventWeights::passEventFilters(*reader_event)) return false;

        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(1.0,weight);


        auto jets = JetKinematics::selectObjects(reader_jet->jets,20,2.4);
        const auto jetPair = getJetPair(&lepton,jets);
        const float ht = JetKinematics::ht(jets,30,2.4);


        bool passHbbPair = false;
        bool passHbbPairT = false;
        if(jetPair.first && jetPair.second){
            const float pairMass = (jetPair.first->p4()+jetPair.second->p4() ).mass();
            const float minCSV = std::min(jetPair.first->csv(),jetPair.second->csv());
            const float maxCSV = std::max(jetPair.first->csv(),jetPair.second->csv());
            passHbbPair = pairMass > 90 && pairMass< 140 && maxCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M] && minCSV >= BTagging::CSVWP_VALS[BTagging::CSV_L];
            passHbbPairT = passHbbPair && minCSV >= BTagging::CSVWP_VALS[BTagging::CSV_M];
        }

        if(!goodWJJFJ) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(2.0,weight);


        if(goodHBBFJ){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(3.0,weight);
            makeHHPlots(prefix + "_hbb",ht,hbbfj->sdMom(),goodHBBFJT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(goodHBBFJT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(4.0,weight);
        }
        if(passHbbPair){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(5.0,weight);
            makeHHPlots(prefix + "_hbbpair",ht,jetPair.first->p4()+jetPair.second->p4(),passHbbPairT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbPairT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(6.0,weight);
        }
        if(passHbbPair && !goodHBBFJ){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(7.0,weight);
            makeHHPlots(prefix + "_hbbpairNoHbb",ht,jetPair.first->p4()+jetPair.second->p4(),passHbbPairT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbPairT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(8.0,weight);
        }
        if(goodHBBFJ && !passHbbPair){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(9.0,weight);
            makeHHPlots(prefix + "_hbbNoPair",ht,hbbfj->sdMom(),goodHBBFJT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(goodHBBFJT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(10.0,weight);
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

    FatJetProcessor fjProc;



};

#endif

void getHHSel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setLumi(36);
    a.analyze();
    a.write(outFileName);
}
void getHHSel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.setLumi(36);
    a.analyze();
    a.write(outFileName);
}
