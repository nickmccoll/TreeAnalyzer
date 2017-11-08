
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
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
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }


    void doSignal(){
        if(diHiggsEvt.type < DiHiggsEvent::MU) return;
        bool passPre = passTriggerPreselection && passEventFilters && selectedLeptons.size() == 1;
        if( wjjCand && passHbbSel && passWlnuDR && passWWDM && nMedBTags_HbbV == 0){
            if(wjjCand->tau2otau1() < 0.55)
                plotter.getOrMake1DPre(smpName, "wjj_mass" ,";Wjj mass [GeV]; arbitrary units",24,0,120)->Fill(wjjCand->sdMom().mass(),weight);
        }
    }



    bool runEvent() override {
        const bool passBase = DefaultSearchRegionAnalyzer::runEvent();
        const float wjjSDMass = wjjCand ? wjjCand->sdMom().mass() : 0;
        const bool passLepton = selectedLeptons.size() ==1;
        if(isSignal() && diHiggsEvt.type < DiHiggsEvent::MU) return false;
        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";



        //WJJ mass
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passHbbSel && passWWDM && passWlnuDR && !nMedBTags_HbbV  ){
            if(wjjCand && wjjCand->tau2otau1() < 0.55)
                plotter.getOrMake1DPre(smpName, "wjj_mass" ,";Wjj soft-drop mass [GeV]; arbitrary units",24,0,120)->Fill(wjjSDMass,weight);

            if(wjjCand && wjjSDMass >= 10)
                plotter.getOrMake1DPre(smpName, "wjj_tau" ,";Wjj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(wjjCand->tau2otau1(),weight);

            if(passWjjSel){
                plotter.getOrMake1DPre(smpName, "wlnu_mass" ,";Wl#nu mass [GeV]; arbitrary units",25,0,250)->Fill(wlnu.mass(),weight);
                plotter.getOrMake1DPre(smpName, "hww_mass" ,";HWW* mass [GeV]; arbitrary units",30,0,300)->Fill(hWW.mass(),weight);
            }
        }

        //MND
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passHbbSel && passWjjSel && passWlnuDR && !nMedBTags_HbbV  ){
                plotter.getOrMake1DPre(smpName, "mnd" ,";#it{m}_{ND}; arbitrary units",50,0,5)->Fill(wwDM,weight);
        }

        //SJ BTagging categories
        if(isSignal()){
            const float signalV = float(signal_mass)/1000.0;
            const int cat = int(hbbCSVCat) -1;
            if(passBase && passTriggerPreselection && passEventFilters && passLepton && hbbCand && passWjjSel && passWWDM && passWlnuDR && !nMedBTags_HbbV  )
            plotter.getOrMake2DPre("signal", "mass_hbbCat" ,";#it{m}(X) [TeV];Hbb sj b-tagging",40,0.550,4.550, 6,-0.5,5.5)->Fill(signalV,cat,weight);

        }

        //Hbb HH mass
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passHbbSel && passWjjSel && passWWDM && passWlnuDR && !nMedBTags_HbbV  ){
            plotter.getOrMake1DPre(smpName,"hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
            plotter.getOrMake1DPre(smpName,"hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
        }

        //QCD reduction
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passHbbSel && passWjjSel && passWWDM  && !nMedBTags_HbbV  ){
            plotter.getOrMake1DPre(smpName, "wlnudr" ,";#DeltaR[lepton,#nu]; arbitrary units",40,0,4)->Fill(wlnuDR,weight);
        }
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passHbbSel && wjjCand && passWWDM  && !nMedBTags_HbbV  ){
            plotter.getOrMake1DPre(smpName +"_loose", "wlnudr" ,";#DeltaR[lepton,#nu]; arbitrary units",40,0,4)->Fill(wlnuDR,weight);
        }

        //Signal selection
        if(isSignal()){
            auto * tots = plotter.getOrMake1DPre(smpName,"signal_evtCounts","N. of events",60,-0.5,59.5);
            tots->Fill(0.0,weight);
            if(passBase && passTriggerPreselection && passEventFilters && passLepton){
                tots->Fill(1.0,weight);
                if(passWjjSel) tots->Fill(2.0,weight);
                if(passWjjSel && passHbbSel)tots->Fill(3.0,weight);
                if(passWjjSel && passHbbSel && passWWDM)tots->Fill(4.0,weight);
                if(passWjjSel && passHbbSel && passWWDM && passWlnuDR )tots->Fill(5.0,weight);
                if(passWjjSel && passHbbSel && passWWDM && passWlnuDR  &&  !nMedBTags_HbbV ){
                    tots->Fill(6.0,weight);
                    if(hbbCSVCat == BTagging::CSVSJ_MF) tots->Fill(7.0,weight);
                    if(hbbCSVCat == BTagging::CSVSJ_ML) tots->Fill(8.0,weight);
                    if(hbbCSVCat == BTagging::CSVSJ_MM) tots->Fill(9.0,weight);
                }
            }
        }

        //Background composition
        if(!isSignal() && passBase && passTriggerPreselection && passEventFilters){
            plotter.getOrMake1D("bkg_trigger_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton )
                plotter.getOrMake1D("bkg_lepVeto_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton && passWjjSel )
                plotter.getOrMake1D("bkg_wjj_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton && passWjjSel && passHbbSel )
                plotter.getOrMake1D("bkg_hbb_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton && passWjjSel && passHbbSel && passWWDM )
                plotter.getOrMake1D("bkg_mnd_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton && passWjjSel && passHbbSel && passWWDM  && passWlnuDR )
                plotter.getOrMake1D("bkg_wlnu_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            if(passLepton && passWjjSel && passHbbSel && passWWDM  && passWlnuDR && !nMedBTags_HbbV  ){
                plotter.getOrMake1D("bkg_abtag_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
                if(hbbCSVCat == BTagging::CSVSJ_MF) plotter.getOrMake1D("bkg_L_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
                if(hbbCSVCat == BTagging::CSVSJ_ML) plotter.getOrMake1D("bkg_M_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
                if(hbbCSVCat == BTagging::CSVSJ_MM) plotter.getOrMake1D("bkg_T_evtCounts","; process",7,1.5,8.5)->Fill(reader_event->process,weight);
            }

        }

        //fit region
        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passWjjSel && passHbbSel && passWWDM && passWlnuDR  &&  !nMedBTags_HbbV){
            plotter.getOrMake2DPre(smpName, "LI_hhmass_hbbmass" ,"; HH mass [TeV];Hbb SD mass",45,0.5,5,25,0,250)->Fill(hh.mass() / 1000.,hbbMass,weight);
            if(hbbCSVCat == BTagging::CSVSJ_MF) plotter.getOrMake2DPre(smpName, "L_hhmass_hbbmass" ,"; HH mass [TeV];Hbb SD mass",45,0.5,5,25,0,250)->Fill(hh.mass() / 1000.,hbbMass,weight);
            if(hbbCSVCat == BTagging::CSVSJ_ML) plotter.getOrMake2DPre(smpName, "M_hhmass_hbbmass" ,"; HH mass [TeV];Hbb SD mass",45,0.5,5,25,0,250)->Fill(hh.mass() / 1000.,hbbMass,weight);
            if(hbbCSVCat == BTagging::CSVSJ_MM) plotter.getOrMake2DPre(smpName, "T_hhmass_hbbmass" ,"; HH mass [TeV];Hbb SD mass",45,0.5,5,25,0,250)->Fill(hh.mass() / 1000.,hbbMass,weight);
            bool q_in = false;
            const float matchR = 0.8*0.8;

            for(const auto& d : smDecayEvt.bosonDecays  ){
                if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
                int nQIn = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < matchR);
                if(nQIn > 0) q_in = true;
            }

            for(const auto& d : smDecayEvt.topDecays  ){
                if(d.type != TopDecay::HAD ) continue;
                bool passB = PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchR;
                int nWIn = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < matchR);
                int nQIn = nWIn + passB;
                if(nQIn > 0) q_in = true;

            }

            if(!isSignal() && reader_event->process != FillerConstants::QCD){
                TString bkgName = q_in ? "res" : "nres";
                plotter.getOrMake1DPre(bkgName+"_LI","hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
                plotter.getOrMake1DPre(bkgName+"_LI","hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
                if(hbbCSVCat == BTagging::CSVSJ_MF) {
                    plotter.getOrMake1DPre(bkgName+"_L","hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
                    plotter.getOrMake1DPre(bkgName+"_L","hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
                }
                if(hbbCSVCat == BTagging::CSVSJ_ML) {
                    plotter.getOrMake1DPre(bkgName+"_M","hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
                    plotter.getOrMake1DPre(bkgName+"_M","hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
                }
                if(hbbCSVCat == BTagging::CSVSJ_MM) {
                    plotter.getOrMake1DPre(bkgName+"_T","hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
                    plotter.getOrMake1DPre(bkgName+"_T","hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
                }
            }
        }

        if(passBase && passTriggerPreselection && passEventFilters && passLepton && passWjjSel && hbbCand && passWWDM && passWlnuDR  &&  !nMedBTags_HbbV && hbbCSVCat < BTagging::CSVSJ_MF){
            plotter.getOrMake2DPre(smpName, "wjj_hhmass_hbbmass" ,"; HH mass [TeV];Hbb SD mass",45,0.5,5,25,0,250)->Fill(hh.mass() / 1000.,hbbMass,weight);

            bool q_in = false;
            const float matchR = 0.8*0.8;

            for(const auto& d : smDecayEvt.bosonDecays  ){
                if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
                int nQIn = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < matchR);
                if(nQIn > 0) q_in = true;
            }

            for(const auto& d : smDecayEvt.topDecays  ){
                if(d.type != TopDecay::HAD ) continue;
                bool passB = PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchR;
                int nWIn = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < matchR);
                int nQIn = nWIn + passB;
                if(nQIn > 0) q_in = true;

            }
            if(!isSignal() && reader_event->process != FillerConstants::QCD){
                TString bkgName = q_in ? "res" : "nres";
                plotter.getOrMake1DPre(bkgName+"_wJJ","hh_mass",";HH mass [TeV]; arbitrary units",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
                plotter.getOrMake1DPre(bkgName+"_wJJ","hbb_mass",";Hbb soft-drop mass [TeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            }
        }
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;




};

#endif

void drawVariables(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void drawVariables(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
