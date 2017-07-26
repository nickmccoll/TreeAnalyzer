
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
    }
    void loadVariables() override {
        reader_event    = (EventReader*)load(new EventReader("event",isRealData()));
        if(treeType == TREE_OTHER)
            reader_genpart  = (GenParticleReader*)load(new GenParticleReader("genParticle"));
        reader_fatjet   = (FatJetReader*)load(new FatJetReader("ak8PuppiNoLepJet",isRealData(),false));
        reader_jet = (JetReader*)load(new JetReader("ak4PuppiNoLepJet",isRealData(),false));

        setBranchAddress("skim" ,"ht"         ,   &ht                  ,true);
        setBranchAddress("skim" ,"selLep_pt"  ,   &selLep_pt           ,true);
        setBranchAddress("skim" ,"selLep_eta" ,   &selLep_eta          ,true);
        setBranchAddress("skim" ,"selLep_phi" ,   &selLep_phi          ,true);
        setBranchAddress("skim" ,"selLep_muon",   &selLep_muon         ,true);
    }

    MomentumF getnz(const MomentumF& met, const MomentumF& vis) {
        const double mH = 125;
        const double a = mH*mH - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
        const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
        const double B = -4*a* vis.z();
        const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
        const double delta = B*B -4*A*C;

        double metZ = 0;
        if(delta < 0) {
            metZ= -B/(2*A);
        } else {
            const double pos = (-B + std::sqrt(delta))/(2*A);
            const double neg = (-B - std::sqrt(delta))/(2*A);
            if(std::fabs(pos) < std::fabs(neg)) metZ = pos;
            else metZ = neg;
        }
        ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
        return MomentumF(neutrino);
    };


    const FatJet* getWjjFJ(const MomentumF* lepton, const std::vector<FatJet*>& fjs) {
        double minDR = 100000;
        int fjIDX = PhysicsUtilities::findNearestDRDeref(*lepton,fjs,minDR);
        if(fjIDX < 0 || minDR > 1.2) return 0;
        return fjs[fjIDX];
    }
    const FatJet* getHbbFJ(const MomentumF* lepton, const std::vector<FatJet*>& fjs) {
        std::vector<FatJet*> sortedFJ = fjs;
        std::sort(sortedFJ.begin(), sortedFJ.end(),PhysicsUtilities::greaterPTDeref<FatJet>());

        const size nFJ = sortedFJ.size();
        const FatJet * fj = 0;
        if (nFJ >= 2) fj =
                PhysicsUtilities::absDeltaPhi(*lepton,*sortedFJ[0]) > PhysicsUtilities::absDeltaPhi(*lepton,*sortedFJ[1])
                ?  sortedFJ[0] : sortedFJ[1];
        else if(nFJ == 1) fj = sortedFJ[0];
        else return 0;

        if(PhysicsUtilities::absDeltaPhi(*lepton,*fj) < 2.0) return 0;
        return fj;
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

    const float CSVT = 0.9535;
    const float CSVM = 0.8484;
    const float CSVL = 0.5426;
    const float DoubleBL  = 0.3;
    const float DoubleBM1 = 0.6;
    const float DoubleBM2 = 0.8;
    const float DoubleBT  = 0.9;





    void makeHHPlots(TString pre, const float ht,  const MomentumF& hbb, const bool passTightCSV, const MomentumF& wjj,const MomentumF& wjjNoSD, const MomentumF& lepton, const MomentumF& met){
        auto neutrino =getnz(met,lepton.p4() + wjj.p4() );
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
        if(ht >= 400)
            makeMPlots(pre+"_hHT");
        else if( ht < 400 && lepton.pt() >= 30)
            makeMPlots(pre+"_lHT");

        if(passTightCSV && ht >= 400 )
            makeMPlots(pre+"_hHT_tCSV");
        if(passTightCSV &&  ht < 400 && lepton.pt() >= 30 )
            makeMPlots(pre+"_lHT_tCSV");

        if(!passTightCSV && ht >= 400 )
            makeMPlots(pre+"_hHT_lCSV");
        if(!passTightCSV &&  ht < 400 && lepton.pt() >= 30 )
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
        weight = EventWeights::getNormalizedEventWeight(reader_event,xsec(),nSampEvt(),lumi());
        if(reader_event->process == FillerConstants::SIGNAL) weight *=  (0.8241887906 * EventWeights::get4bXSecLimit(mass)/1000.0);
        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(0.0,weight);


        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_goodVertices) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_globalTightHalo2016Filter) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_HBHENoiseFilter) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_HBHENoiseIsoFilter) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::Flag_eeBadScFilter) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::AnaTM_badMuons) ) return false;
        if(!FillerConstants::doesPass(reader_event->metFilters,FillerConstants::AnaTM_badChargedHadrons) ) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(1.0,weight);


        auto fjs = JetKinematics::selectObjects(reader_fatjet->jets,50,2.4);
        const auto* hbbfj = getHbbFJ(&lepton,fjs);

        auto jets = JetKinematics::selectObjects(reader_jet->jets,20,2.4);
        const auto jetPair = getJetPair(&lepton,jets);
        const float ht = JetKinematics::ht(jets,30,2.4);

        const auto* wjjfj = getWjjFJ(&lepton,fjs);
        bool passWjj = wjjfj && wjjfj->sdMom().mass() > 10 && wjjfj->tau2otau1() < 0.55 && wjjfj->maxSJCSV() < CSVM;
        bool passHbb = hbbfj && hbbfj->sdMom().mass() > 90 && hbbfj->sdMom().mass() < 140 && hbbfj->tau2otau1() < 0.55 && hbbfj->maxSJCSV() >= CSVM && hbbfj->minSJCSV() >= CSVL;
        bool passHbbT = passHbb && hbbfj->minSJCSV() >= CSVM;

        bool passHbbPair = false;
        bool passHbbPairT = false;
        if(jetPair.first && jetPair.second){
            const float pairMass = (jetPair.first->p4()+jetPair.second->p4() ).mass();
            const float minCSV = std::min(jetPair.first->csv(),jetPair.second->csv());
            const float maxCSV = std::max(jetPair.first->csv(),jetPair.second->csv());
            passHbbPair = pairMass > 90 && pairMass< 140 && maxCSV >= CSVM && minCSV >= CSVL;
            passHbbPairT = passHbbPair && minCSV >= CSVM;
        }

        if(!passWjj || !wjjfj->passTightID()) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(2.0,weight);


        if(passHbb && hbbfj->passTightID() ){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(3.0,weight);
            makeHHPlots(prefix + "_hbb",ht,hbbfj->sdMom(),passHbbT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(4.0,weight);
        }
        if(passHbbPair){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(5.0,weight);
            makeHHPlots(prefix + "_hbbpair",ht,jetPair.first->p4()+jetPair.second->p4(),passHbbPairT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbPairT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(6.0,weight);
        }
        if(passHbbPair && !passHbb){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(7.0,weight);
            makeHHPlots(prefix + "_hbbpairNoHbb",ht,jetPair.first->p4()+jetPair.second->p4(),passHbbPairT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbPairT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(8.0,weight);
        }
        if(passHbb && !passHbbPair){
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(9.0,weight);
            makeHHPlots(prefix + "_hbbNoPair",ht,hbbfj->sdMom(),passHbbT,wjjfj->sdMom(),wjjfj->p4(),lepton,reader_event->met.p4());
            if(passHbbT)
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(10.0,weight);
        }


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    EventReader       * reader_event = 0;
    GenParticleReader * reader_genpart  = 0;
    ElectronReader    * reader_electron = 0;
    MuonReader        * reader_muon     = 0;
    JetReader         * reader_jet      = 0;
    FatJetReader      * reader_fatjet   = 0;
    float   ht         =0;
    float   selLep_pt  =0;
    float   selLep_eta =0;
    float   selLep_phi =0;
    size8   selLep_muon=0;
    float weight  = 0;
    HistGetter plotter;
    TString prefix;
    int mass=0;



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
