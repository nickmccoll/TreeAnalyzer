
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
        TPRegexp r1(".*(m\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
//        std::cout <<" ->>"<<match->GetSize()  <<std::endl;
//        std::cout <<" ->>"<<((TObjString *)match->At(0))->GetString()  <<std::endl;
//        std::cout <<" ->>"<<((TObjString *)match->At(1))->GetString()  <<std::endl;
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1)
            prefix = ((TObjString *)match->At(1))->GetString();
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

    const FatJet * getSignalFJ(const DiHiggsEvent& diHiggsEvt, const MomentumF* lepton, const std::vector<FatJet*>& fjs){

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

    const FatJet* getFJ(const MomentumF* lepton, const std::vector<FatJet*>& fjs) {
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

    const float CSVT = 0.9535;
    const float CSVM = 0.8484;
    const float CSVL = 0.5426;
    const float DoubleBL  = 0.3;
    const float DoubleBM1 = 0.6;
    const float DoubleBM2 = 0.8;
    const float DoubleBT  = 0.9;


    void makeHbbPairPlots(TString pre, const std::pair<const Jet*,const Jet*>& jetPair){

        const float pairMass = (jetPair.first->p4()+jetPair.second->p4() ).mass();
        const float minCSV = std::min(jetPair.first->csv(),jetPair.second->csv());
        const float maxCSV = std::max(jetPair.first->csv(),jetPair.second->csv());

        plotter.getOrMake1DPre(pre,"hbb_pair_mass"             ,";hbb pair mass [GeV]; a.u."       ,250,0,500)->Fill(pairMass,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_minsdcsv"         ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(minCSV,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxsdcsv"         ,";hbb pair bb max sj csv; a.u."    ,100,0,1  )->Fill(maxCSV,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxMed_minsdcsv"  ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= CSVM ? minCSV :0.0,weight );
        plotter.getOrMake1DPre(pre,"hbb_pair_maxTight_minsdcsv",";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= CSVT ? minCSV :0.0,weight );

        bool passMassL = (pairMass > 20);
        bool passMass = (pairMass > 90 && pairMass< 140);
        bool passCSV = (maxCSV >= CSVM && minCSV >= CSVL  );
        bool passCSVT = (minCSV >= CSVM  );

        if(passMass)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(7.0,weight);
        if(passMass&&passCSV)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(8.0,weight);
        if(passMass && passCSVT)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(9.0,weight);


        if(passMass){
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_minsdcsv"         ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(minCSV,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxsdcsv"         ,";hbb pair bb max sj csv; a.u."    ,100,0,1  )->Fill(maxCSV,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxMed_minsdcsv"  ,";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= CSVM ? minCSV :0.0,weight );
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_maxTight_minsdcsv",";hbb pair bb min sj csv; a.u."    ,100,0,1  )->Fill(maxCSV >= CSVT ? minCSV :0.0,weight );
        }

        if(passCSV&&passMassL){
            plotter.getOrMake1DPre(pre,"hbb_pair_oM_mass"             ,";hbb pair mass [GeV]; a.u."       ,250,0,500)->Fill(pairMass,weight );
        }
    }




    void makeHbbPlots(TString pre, const FatJet* fj){

        plotter.getOrMake1DPre(pre,"hbb_fj_mass"             ,";hbb fj mass [GeV]; a.u."       ,250,0,500)->Fill(fj->mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_sd_mass"          ,";hbb fj sd mass [GeV]; a.u."    ,250,0,500)->Fill(fj->sdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_rawsd_mass"       ,";hbb fj raw sd mass [GeV]; a.u.",250,0,500)->Fill(fj->rawSdMom().mass(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_tau2otau1"        ,";hbb fj #tau_{2}/#tau_{1}; a.u.",100,0,1  )->Fill(fj->tau2otau1(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_csv"              ,";hbb fj csv; a.u."              ,100,0,1  )->Fill(fj->csv(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_bbcsv"            ,";hbb fj bb csv; a.u."           ,100,0,1  )->Fill(fj->bbt(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_minsdcsv"         ,";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->minSJCSV(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxsdcsv"         ,";hbb fj bb max sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV(),weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxMed_minsdcsv"  ,";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV() >= CSVM ? fj->minSJCSV() :0.0,weight );
        plotter.getOrMake1DPre(pre,"hbb_fj_maxTight_minsdcsv",";hbb fj bb min sj csv; a.u."    ,100,0,1  )->Fill(fj->maxSJCSV() >= CSVT? fj->minSJCSV() :0.0,weight );

        bool passMassL = (fj->rawSdMom().mass() > 20);
        bool passMass = (fj->sdMom().mass() > 90 && fj->sdMom().mass() < 140);
        bool passTau = fj->tau2otau1() < 0.55;
        bool passCSV = (fj->maxSJCSV() >= CSVM && fj->minSJCSV() >= CSVL  );
        bool passCSVT = (fj->minSJCSV() >= CSVM  );

        if(passMass)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(2.0,weight);
        if(passMass&&passTau)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(3.0,weight);
        if(passMass&&passTau&&passCSV)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(4.0,weight);
        if(passMass && passTau&& passCSVT)
            plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(5.0,weight);


        if(passMass && passTau){
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_csv"              ,";hbb fj csv; a.u."          ,100,0,1)->Fill(fj->csv(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_bbcsv"            ,";hbb fj bb csv; a.u."       ,100,0,1)->Fill(fj->bbt(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_minsdcsv"         ,";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->minSJCSV(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxsdcsv"         ,";hbb fj bb max sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV(),weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxMed_minsdcsv"  ,";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV() >= CSVM ? fj->minSJCSV() :0.0,weight );
            plotter.getOrMake1DPre(pre,"hbb_fj_oM_maxTight_minsdcsv",";hbb fj bb min sj csv; a.u.",100,0,1)->Fill(fj->maxSJCSV() >= CSVT? fj->minSJCSV() :0.0,weight );
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

        if(reader_event->process !=
                FillerConstants::SIGNAL){
            prefix = reader_event->process >= FillerConstants::ZJETS ? "rare" :
                    FillerConstants::MCProcessNames[reader_event->process];
        } else {
            diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
            if(diHiggsEvt.type < DiHiggsEvent::MU) return 0;
        }

        MomentumF lepton(ASTypes::CylLorentzVectorF(selLep_pt,selLep_eta,selLep_phi,0));
        weight = EventWeights::getNormalizedEventWeight(reader_event,xsec(),nSampEvt(),lumi());
        auto fjs = JetKinematics::selectObjects(reader_fatjet->jets,50,2.4);

        plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(0.0,weight);

        //Do fat jets
        const auto* fj = getFJ(&lepton,fjs);
        if(fj){
            bool process = true;
            if(reader_event->process == FillerConstants::SIGNAL) {
                const auto* sigfj = getSignalFJ(diHiggsEvt,&lepton,fjs);
                if(sigfj != fj) process = false;
            }
            if(process){
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(1.0,weight);
                makeHbbPlots(prefix,fj);
            }
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
                plotter.getOrMake1DPre(prefix,"selection",";selection; a.u.",20,-0.5,19.5 )->Fill(6.0,weight);
                makeHbbPairPlots(prefix,jetPair);
            }
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
