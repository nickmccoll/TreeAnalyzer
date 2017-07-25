
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


        const MomentumF wjj = diHiggsEvt.w2_d1->p4() + diHiggsEvt.w2_d2->p4();

        plotter.getOrMake1DPre(prefix,"wjj_pt",";wjj p_{T}; arbitrary units",400,0,2000.0 )->Fill(wjj.pt(),weight);


        if(wjj.pt() < 50) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(2.0,weight);
        if(wjj.absEta() > 2.4) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(3.0,weight);
        if(PhysicsUtilities::deltaR(*diHiggsEvt.w2_d1,*diHiggsEvt.w2_d2) > 0.8 ) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(4.0,weight);

        auto fjs = JetKinematics::selectObjects(reader_fatjet->jets,50,2.4);
        double minDR = 0;
        int fjIDX = PhysicsUtilities::findNearestDRDeref(wjj,fjs,minDR);
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
        std::sort(rankedDRLeps.begin(), rankedDRLeps.end(),PhysicsUtilities::lesserAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDRLeps[iJ].second == fjIDX) {drlepRank = iJ; break;}

        int dphiWRank = 0;
        std::vector<std::pair<float,int>> rankedDphiW(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedDphiW[iJ] = std::make_pair(PhysicsUtilities::absDeltaPhi(recoW,*fjs[iJ]),iJ);
        std::sort(rankedDphiW.begin(), rankedDphiW.end(),PhysicsUtilities::lesserAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDphiW[iJ].second == fjIDX) {dphiWRank = iJ; break;}

        int dphilepRank = 0;
        std::vector<std::pair<float,int>> rankedDphilep(fjs.size());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) rankedDphilep[iJ] = std::make_pair(PhysicsUtilities::absDeltaPhi(lepton,*fjs[iJ]),iJ);
        std::sort(rankedDphilep.begin(), rankedDphilep.end(),PhysicsUtilities::lesserAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ) if(rankedDphilep[iJ].second == fjIDX) {dphilepRank = iJ; break;}



        plotter.getOrMake1DPre(prefix,"ptRank",";wjj fj #it{p}_{T} rank; arbitrary units",10,-0.5,9.5 )->Fill(ptRank,weight);
        plotter.getOrMake1DPre(prefix,"drlepRank",";wjj fj #DeltaR(lep) rank; arbitrary units",10,-0.5,9.5 )->Fill(drlepRank,weight);
        plotter.getOrMake1DPre(prefix,"dphiWRank",";wjj fj #Delta#phi(W) rank; arbitrary units",10,-0.5,9.5 )->Fill(dphiWRank,weight);
        plotter.getOrMake1DPre(prefix,"dphilepRank",";wjj fj #Delta#phi(lep) rank; arbitrary units",10,-0.5,9.5 )->Fill(dphilepRank,weight);

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

        plotter.getOrMake1DPre(prefix,"leadpt_drlepRank",";wjj fj #DeltaR(lep) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_drlepRank,weight);
        plotter.getOrMake1DPre(prefix,"leadpt_dphiWRank",";wjj fj #Delta#phi(W) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_dphiWRank,weight);
        plotter.getOrMake1DPre(prefix,"leadpt_dphilepRank",";wjj fj #Delta#phi(lep) rank; arbitrary units",3,-1.5,1.5 )->Fill(leadpt_dphilepRank,weight);

        std::vector<FatJet*> highDPhiJets;
        for(unsigned int iJ = 0; iJ < fjs.size(); ++iJ){
            if(PhysicsUtilities::absDeltaPhi(lepton,*fjs[iJ]) > 2) continue;
            highDPhiJets.push_back(fjs[iJ]);
        }
        int dPhi_leadPTRank = -1;
        std::vector<std::pair<float,int>> dPhi_rankedPTS(fjs.size());
        for(unsigned int iJ = 0; iJ < highDPhiJets.size(); ++iJ) dPhi_rankedPTS[iJ] = std::make_pair(highDPhiJets[iJ]->pt(),iJ);
        std::sort(dPhi_rankedPTS.begin(), dPhi_rankedPTS.end(),PhysicsUtilities::lesserAbsFirst<float,int>());
        for(unsigned int iJ = 0; iJ < highDPhiJets.size(); ++iJ) if(dPhi_rankedPTS[iJ].second == fjIDX) {dPhi_leadPTRank = iJ; break;}
        plotter.getOrMake1DPre(prefix,"dPhi_leadPTRank",";higgs fj #it{p}_{T} rank; arbitrary units",11,-1.5,9.5 )->Fill(dPhi_leadPTRank,weight);


        plotter.getOrMake1DPre(prefix,"pt",";higgs fj #it{p}_{T}; arbitrary units",600,0,3000 )->Fill(fj->pt(),weight);
        plotter.getOrMake1DPre(prefix,"drlep",";higgs #DeltaR(lep); arbitrary units",500,0,5.0  )->Fill(PhysicsUtilities::deltaR(lepton,*fj),weight);
        plotter.getOrMake1DPre(prefix,"dphiW",";higgs fj #Delta#phi(W); arbitrary units",500,0,5.0 )->Fill(PhysicsUtilities::absDeltaPhi(recoW,*fj),weight);
        plotter.getOrMake1DPre(prefix,"dphilep",";higgs fj #Delta#phi(lep); arbitrary units",500,0,5.0 )->Fill(PhysicsUtilities::absDeltaPhi(lepton,*fj),weight);

        if(drlepRank > 0) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(7.0,weight);
        if(PhysicsUtilities::deltaR(lepton,*fj) > 1.2) return false;
        plotter.getOrMake1DPre(prefix,"selection",";selection; arbitrary units",20,-0.5,19.5 )->Fill(8.0,weight);



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

void getWjjJetCandSel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getWjjJetCandSel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
