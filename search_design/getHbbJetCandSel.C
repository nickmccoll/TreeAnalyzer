
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
        TPRegexp r1(".*(m\\d+).*");
        auto match = r1.MatchS(fileName);
        if(match->GetSize())
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



    void getAK4Sel(const DiHiggsEvent& diHiggsEvt, const MomentumF& lepton, const float weight ) {
        const float selOff = 10;
        if(diHiggsEvt.b1->pt() < 20 || diHiggsEvt.b2->pt() < 20 || diHiggsEvt.b1->absEta() > 2.4 ||diHiggsEvt.b1->absEta() > 2.4)
            return;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(0.0 +  selOff,weight);
        auto jets = JetKinematics::selectObjects(reader_jet->jets,20,2.4);
        double minDR1 = 0;
        int j1IDX = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.b1,jets,minDR1);
        if(j1IDX < 0) return;
        std::vector<bool> vetoed(jets.size(),false); vetoed[j1IDX] = true;
        double minDR2 = 0;
        int j2IDX = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.b2,jets,minDR2,99999,20,&vetoed);
        if(j2IDX < 0) return;
        plotter.getOrMake1DPre(prefix,"maxDRgenAK4Jet",";max gen #DeltaR(jet); arbitrary units",500,0,5.0  )->Fill(std::max(minDR2,minDR1),weight);
        if(std::max(minDR2,minDR1) > 0.4) return;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(1.0 +  selOff,weight);
        const auto* jet1 = jets[j1IDX];
        const auto* jet2 = jets[j2IDX];
        plotter.getOrMake1DPre(prefix,"drAK4Jets",";#DeltaR(jet,jet); arbitrary units",500,0,5.0  )->Fill(PhysicsUtilities::deltaR(*jet1,*jet2),weight);
        plotter.getOrMake1DPre(prefix,"min_dphiAK4Jetslep",";min #Delta#phi(jet,lepton); arbitrary units",500,0,5.0  )->Fill(
                std::min(PhysicsUtilities::absDeltaPhi(*jet1,lepton),PhysicsUtilities::absDeltaPhi(*jet2,lepton)),weight);
        plotter.getOrMake1DPre(prefix,"min_dRAK4Jetslep",";min #Delta#R(jet,lepton); arbitrary units",500,0,5.0  )->Fill(
                std::min(PhysicsUtilities::deltaR(*jet1,lepton),PhysicsUtilities::deltaR(*jet2,lepton)),weight);

        if(std::min(PhysicsUtilities::absDeltaPhi(*jet1,lepton),PhysicsUtilities::absDeltaPhi(*jet2,lepton)) < TMath::PiOver2()) return;
        if(PhysicsUtilities::deltaR(*jet1,*jet2) > 1.5) return;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(2.0 +  selOff,weight);


        int ptRank1 = 0;
        int ptRank2 = 0;
        std::vector<std::pair<float,int>> rankedPTS(jets.size());
        for(unsigned int iJ = 0; iJ < jets.size(); ++iJ){
            if(PhysicsUtilities::absDeltaPhi(*jets[iJ],lepton) < TMath::PiOver2()) continue;
            rankedPTS[iJ] = std::make_pair(jets[iJ]->pt(),iJ);
        }
        std::sort(rankedPTS.begin(), rankedPTS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < jets.size(); ++iJ) {
            if(rankedPTS[iJ].second == j1IDX) {ptRank1 = iJ;}
            if(rankedPTS[iJ].second == j2IDX) {ptRank2 = iJ;}
        }

        plotter.getOrMake1DPre(prefix,"ak4jetptRank",";higgs max(j) #it{p}_{T} rank; arbitrary units",10,-0.5,9.5 )->Fill(std::max(ptRank1,ptRank2),weight);


        std::vector<std::pair<const Jet*,const Jet*> > jetPairs;
        for(unsigned int iJ1 = 0; iJ1 < jets.size(); ++iJ1){
            if(PhysicsUtilities::absDeltaPhi(*jets[iJ1],lepton) < TMath::PiOver2()) continue;
            for(unsigned int iJ2 = iJ1+1; iJ2 < jets.size(); ++iJ2){
                if(PhysicsUtilities::absDeltaPhi(*jets[iJ2],lepton) < TMath::PiOver2()) continue;
                if(PhysicsUtilities::deltaR(*jets[iJ1],*jets[iJ2]) > 1.5) continue;
                jetPairs.emplace_back(std::make_pair(jets[iJ1],jets[iJ2]));
            }
        }
        int ptPairRank1 = 0;
        unsigned int goodPairIDX= 0;
        std::vector<std::pair<float,size>> rankedPairPTS(jetPairs.size());
        for(unsigned int iJ = 0; iJ < jetPairs.size(); ++iJ){
            if((jetPairs[iJ].first == jet1 || jetPairs[iJ].first == jet2) &&(jetPairs[iJ].second == jet1 || jetPairs[iJ].second == jet2))
                goodPairIDX = iJ;
            rankedPairPTS[iJ] = std::make_pair((jetPairs[iJ].first->p4()+jetPairs[iJ].second->p4()).pt(),iJ);
        }
        std::sort(rankedPairPTS.begin(), rankedPairPTS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < rankedPairPTS.size(); ++iJ) {
            if(rankedPairPTS[iJ].second == goodPairIDX) {ptPairRank1 = iJ;}
        }
        plotter.getOrMake1DPre(prefix,"ak4jetptPairRank",";higgs pair #it{p}_{T} rank; arbitrary units",10,-0.5,9.5 )->Fill(ptPairRank1,weight);

        for(unsigned int iJ = 0; iJ < jetPairs.size(); ++iJ){
            auto pairMom = jetPairs[iJ].first->p4() + jetPairs[iJ].second->p4();
            float pairDPHI = PhysicsUtilities::absDeltaPhi(pairMom,lepton);
            float pairDR = PhysicsUtilities::deltaR(*jetPairs[iJ].first,*jetPairs[iJ].second);
            if(iJ ==goodPairIDX ){
                plotter.getOrMake1DPre(prefix,"dphipairlep",";min #Delta#phi(jet pair,lepton); arbitrary units",500,0,5.0  )->Fill(pairDPHI,weight);
                plotter.getOrMake1DPre(prefix,"drpairjj",";min #Delta#R(jet,jet); arbitrary units",500,0,5.0  )->Fill(pairDR,weight);
            }else {
                plotter.getOrMake1DPre(prefix,"bkg_dphipairlep",";min #Delta#phi(jet pair,lepton); arbitrary units",500,0,5.0  )->Fill(pairDPHI,weight);
                plotter.getOrMake1DPre(prefix,"bkg_drpairjj",";min #Delta#R(jet,jet); arbitrary units",500,0,5.0  )->Fill(pairDR,weight);

            }

        }

        int ptPairBYDRRank = -1;
        std::vector<std::pair<float,size>> rankedPairDRS(jetPairs.size());
        for(unsigned int iJ = 0; iJ < jetPairs.size(); ++iJ){
            rankedPairDRS[iJ] = std::make_pair(PhysicsUtilities::deltaR(*jetPairs[iJ].first,*jetPairs[iJ].second),iJ);
        }
        std::sort(rankedPairPTS.begin(), rankedPairPTS.end(),PhysicsUtilities::lesserAbsFirst<float,int>());
        std::vector<int> usedJets;
        std::vector<std::pair<float,size>> rankedPairPTSByDRS;
        for(unsigned int iJ = 0; iJ < rankedPairDRS.size(); ++iJ){
            const auto* pj1 = jetPairs[rankedPairDRS[iJ].second].first;
            const auto* pj2 = jetPairs[rankedPairDRS[iJ].second].second;
            bool used = false;
            for(const auto& uj : usedJets){
                if(pj1->index() == uj) used = true;
                if(pj2->index() == uj) used = true;
            }
            if(used) continue;
            rankedPairPTSByDRS.emplace_back(std::make_pair((pj1->p4() + pj2->p4()).pt(),rankedPairDRS[iJ].second));
            usedJets.push_back(pj1->index());
            usedJets.push_back(pj2->index());
        }
        std::sort(rankedPairPTSByDRS.begin(), rankedPairPTSByDRS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < rankedPairPTSByDRS.size(); ++iJ) {
            if(rankedPairPTSByDRS[iJ].second == goodPairIDX) {ptPairBYDRRank = iJ;}
        }
        plotter.getOrMake1DPre(prefix,"ak4jetptPairRank_byDR",";higgs pair #it{p}_{T} rank (DR match); arbitrary units",11,-1.5,9.5 )->Fill(ptPairBYDRRank,weight);

        if(ptPairRank1 == 0)
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(3.0 +  selOff,weight);


    }

    bool runEvent() override {
        if(reader_event->process != FillerConstants::SIGNAL) return false;

        MomentumF lepton(ASTypes::CylLorentzVectorF(selLep_pt,selLep_eta,selLep_phi,0));
        MomentumF recoW = lepton.p4() + reader_event->met.p4();

        const float weight = EventWeights::getNormalizedEventWeight(reader_event,xsec(),nSampEvt(),lumi());
        DiHiggsEvent diHiggsEvt; diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;

        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(0.0,weight);

        double genLepDR = PhysicsUtilities::deltaR(lepton,*diHiggsEvt.w1_d1);
        if(genLepDR > 0.2) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(1.0,weight);

        if(PhysicsUtilities::deltaR(*diHiggsEvt.b1,*diHiggsEvt.b2) > 0.8 )
            getAK4Sel(diHiggsEvt,lepton,weight);



        const MomentumF hbb = diHiggsEvt.b1->p4() + diHiggsEvt.b2->p4();
        if(hbb.pt() < 50) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(2.0,weight);
        if(hbb.absEta() > 2.4) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(3.0,weight);
        if(PhysicsUtilities::deltaR(*diHiggsEvt.b1,*diHiggsEvt.b2) > 0.8 ) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(4.0,weight);

        auto fjs = JetKinematics::selectObjects(reader_fatjet->jets,50,2.4);
        double minDR = 0;
        int fjIDX = PhysicsUtilities::findNearestDRDeref(hbb,fjs,minDR);
        const auto* fj = fjs[fjIDX];
        if(fjIDX < 0) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(5.0,weight);

        plotter.getOrMake1DPre(prefix,"nearestDR",";nearest fj #DeltaR; arbitrary units",500,0,5.0 )->Fill(minDR,weight);

        if(minDR > 0.5) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(6.0,weight);

        int ptRank = 0;
        std::vector<std::pair<float,int>> rankedPTS(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedPTS[iJ] = std::make_pair(fjs[iJ]->pt(),iJ);
        std::sort(rankedPTS.begin(), rankedPTS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedPTS[iJ].second == fjIDX) {ptRank = iJ; break;}

        int drlepRank = 0;
        std::vector<std::pair<float,int>> rankedDRLeps(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedDRLeps[iJ] = std::make_pair(PhysicsUtilities::deltaR2(lepton,*fjs[iJ]),iJ);
        std::sort(rankedDRLeps.begin(), rankedDRLeps.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDRLeps[iJ].second == fjIDX) {drlepRank = iJ; break;}

        int dphiWRank = 0;
        std::vector<std::pair<float,int>> rankedDphiW(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedDphiW[iJ] = std::make_pair(PhysicsUtilities::absDeltaPhi(recoW,*fjs[iJ]),iJ);
        std::sort(rankedDphiW.begin(), rankedDphiW.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDphiW[iJ].second == fjIDX) {dphiWRank = iJ; break;}

        int dphilepRank = 0;
        std::vector<std::pair<float,int>> rankedDphilep(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedDphilep[iJ] = std::make_pair(PhysicsUtilities::absDeltaPhi(lepton,*fjs[iJ]),iJ);
        std::sort(rankedDphilep.begin(), rankedDphilep.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDphilep[iJ].second == fjIDX) {dphilepRank = iJ; break;}



        plotter.getOrMake1DPre(prefix,"ptRank",";higgs fj #it{p}_{T} rank; arbitrary units",10,-0.5,9.5 )->Fill(ptRank,weight);
        plotter.getOrMake1DPre(prefix,"drlepRank",";higgs fj #DeltaR(lep) rank; arbitrary units",10,-0.5,9.5 )->Fill(drlepRank,weight);
        plotter.getOrMake1DPre(prefix,"dphiWRank",";higgs fj #Delta#phi(W) rank; arbitrary units",10,-0.5,9.5 )->Fill(dphiWRank,weight);
        plotter.getOrMake1DPre(prefix,"dphilepRank",";higgs fj #Delta#phi(lep) rank; arbitrary units",10,-0.5,9.5 )->Fill(dphilepRank,weight);

        int leadpt_drlepRank = -1;
        int leadpt_dphiWRank = -1;
        int leadpt_dphilepRank = -1;
        if(fjs.size() == 1){
            leadpt_drlepRank =0;
            leadpt_dphiWRank =0;
            leadpt_dphilepRank =0;
        }
        else if(ptRank < 2){
            const auto * oj = fjs[ptRank==0 ? 1 : 0 ];
            leadpt_drlepRank  = PhysicsUtilities::deltaR2(lepton,*fj) >  PhysicsUtilities::deltaR2(lepton,*oj) ? 0 : 1;
            leadpt_dphiWRank  = PhysicsUtilities::absDeltaPhi(recoW,*fj) >  PhysicsUtilities::absDeltaPhi(recoW,*oj) ? 0 : 1;
            leadpt_dphilepRank  = PhysicsUtilities::absDeltaPhi(lepton,*fj) >  PhysicsUtilities::absDeltaPhi(lepton,*oj) ? 0 : 1;

        }

        plotter.getOrMake1DPre(prefix,"leadpt_drlepRank",";higgs fj #DeltaR(lep) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_drlepRank,weight);
        plotter.getOrMake1DPre(prefix,"leadpt_dphiWRank",";higgs fj #Delta#phi(W) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_dphiWRank,weight);
        plotter.getOrMake1DPre(prefix,"leadpt_dphilepRank",";higgs fj #Delta#phi(lep) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_dphilepRank,weight);

        std::vector<FatJet*> highDPhiJets;
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ){
            if(PhysicsUtilities::absDeltaPhi(lepton,*fjs[iJ]) < 2) continue;
            highDPhiJets.push_back(fjs[iJ]);
        }
        int dPhi_leadPTRank = -1;
        std::vector<std::pair<float,int>> dPhi_rankedPTS(fjs.size());
        for(unsigned int iJ = 0; iJ < highDPhiJets.size(); ++iJ) dPhi_rankedPTS[iJ] = std::make_pair(highDPhiJets[iJ]->pt(),iJ);
        std::sort(dPhi_rankedPTS.begin(), dPhi_rankedPTS.end(),PhysicsUtilities::greaterAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < highDPhiJets.size(); ++iJ) if(dPhi_rankedPTS[iJ].second == fjIDX) {dPhi_leadPTRank = iJ; break;}
        plotter.getOrMake1DPre(prefix,"dPhi_leadPTRank",";higgs fj #it{p}_{T} rank; arbitrary units",11,-1.5,9.5 )->Fill(dPhi_leadPTRank,weight);


        plotter.getOrMake1DPre(prefix,"pt",";higgs fj #it{p}_{T}; arbitrary units",600,0,3000 )->Fill(fj->pt(),weight);
        plotter.getOrMake1DPre(prefix,"drlep",";higgs #DeltaR(lep); arbitrary units",500,0,5.0  )->Fill(PhysicsUtilities::deltaR(lepton,*fj),weight);
        plotter.getOrMake1DPre(prefix,"dphiW",";higgs fj #Delta#phi(W); arbitrary units",500,0,5.0 )->Fill(PhysicsUtilities::absDeltaPhi(recoW,*fj),weight);
        plotter.getOrMake1DPre(prefix,"dphilep",";higgs fj #Delta#phi(lep); arbitrary units",500,0,5.0 )->Fill(PhysicsUtilities::absDeltaPhi(lepton,*fj),weight);

        if(ptRank > 1) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(7.0,weight);
        if(leadpt_dphilepRank != 0) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(8.0,weight);
        if(PhysicsUtilities::absDeltaPhi(lepton,*fj) < 2.0) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(9.0,weight);



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

    HistGetter plotter;
    TString prefix;



};

#endif

void getHbbJetCandSel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getHbbJetCandSel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
