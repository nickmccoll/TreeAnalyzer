
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

        leptonNoISOProc  .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonNoISOProc);
        leptonNoISOProc->lepSelParams.mu_maxISO  = -1;
        leptonNoISOProc->lepSelParams.el_maxISO  = -1;
        leptonNoISOProc->lepSelParams_dataABCDEF.mu_maxISO  = -1;
        leptonNoISOProc->lepSelParams_dataABCDEF.el_maxISO  = -1;
    }



    void doLeptonEff() {
        if(diHiggsEvt.type < DiHiggsEvent::MU) return;

        const float flMass = float(signal_mass)/1000.;
        auto slLepNoISO_mu = leptonNoISOProc->getMuons(*reader_event,*reader_muon);
        auto slLepNoISO_e  = leptonNoISOProc->getElectrons(*reader_electron);
        auto slLep_mu = leptonProc->getMuons(*reader_event,*reader_muon);
        auto slLep_e  = leptonProc->getElectrons(*reader_electron);

        plotter.getOrMake1DPre("signal_incl", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
        if(ht_chs >= 500)
            plotter.getOrMake1DPre("signal_ht", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);


        if(diHiggsEvt.w1_d1->absPdgId() == ParticleInfo::p_eminus){
            plotter.getOrMake1DPre("signal_e_incl", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
            if(slLepNoISO_e.size() && slLepNoISO_e.front()->pt() >= 30  )
                plotter.getOrMake1DPre("signal_e_id", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
            if(slLep_e.size() && slLep_e.front()->pt() >= 30 )
                plotter.getOrMake1DPre("signal_e_iso", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);


            if(slLepNoISO_e.size() && slLepNoISO_e.front()->pt() >= 20  )
                plotter.getOrMake1DPre("signal_e_id2", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
            if(slLep_e.size() && slLep_e.front()->pt() >= 20 )
                plotter.getOrMake1DPre("signal_e_iso2", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);


        } else if(diHiggsEvt.w1_d1->absPdgId() == ParticleInfo::p_muminus){
            plotter.getOrMake1DPre("signal_mu_incl", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);

            if(slLepNoISO_mu.size() && slLepNoISO_mu.front()->pt() >= 26)
                plotter.getOrMake1DPre("signal_mu_id", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
            if(slLep_mu.size() && slLep_mu.front()->pt() >= 26 )
                plotter.getOrMake1DPre("signal_mu_iso", "mass" ,";#it{m}(X) [TeV]; arbitrary units",40,0.550,4.550)->Fill(flMass,weight);
        }
    };


    bool runEvent() override {
        const bool passBase = DefaultSearchRegionAnalyzer::runEvent();

        doLeptonEff();


        float wjjSDMASS = 0;
        float passWjjAB = nMedBTags_HbbV == 0;
        const bool hBBTM = (hbbMass > 105 && hbbMass < 145);

        if(wjjCand){
            wjjSDMASS = wjjCand->sdMom().mass();
        }

        auto mkStatSect = [&](TH1* tots, int startN){
            if(passHbbSel) tots->Fill(startN,weight);
            if(passHbbSel && passWWDM) tots->Fill(startN+1,weight);
            if(passHbbSel && passWWDM && passWlnuDR){
                tots->Fill(startN+2,weight);
                if(hbbCSVCat == BTagging::CSVSJ_MF) tots->Fill(startN+3,weight);
                if(hbbCSVCat == BTagging::CSVSJ_ML) tots->Fill(startN+4,weight);
                if(hbbCSVCat == BTagging::CSVSJ_MM) tots->Fill(startN+5,weight);
        }
        };

        auto makeStats=[&](const TString& prefix){
            auto * tots = plotter.getOrMake1DPre(prefix,"evtCounts","N. of events",60,-0.5,59.5);
            tots->Fill(0.0,weight);
            if(ht_chs < 500) return;
            tots->Fill(1.0,weight);
            if(!passTriggerPreselection) return;
            tots->Fill(2.0,weight);
            if(!passEventFilters) return;
            tots->Fill(3.0,weight);
            if(selectedLeptons.size() != 1) return;
            tots->Fill(4.0,weight);
            if(wjjCand == 0) return;
            tots->Fill(5.0,weight);
            if(wjjSDMASS < 10) return;
            tots->Fill(6.0,weight);
            if(wjjCand->tau2otau1() >= 0.55) return;
            tots->Fill(7.0,weight);
            if(!passWjjAB) return;
            tots->Fill(8.0,weight);
            if(hbbCand == 0) return;
            tots->Fill(9.0,weight);

            mkStatSect(tots,10);

            if(hBBTM){
            tots->Fill(15.0,weight);
            mkStatSect(tots,17);
            }

            if(hh.mass() >= 900 && hh.mass() < 1100){
                tots->Fill(23.0,weight);
                mkStatSect(tots,24);
                if(hBBTM){
                tots->Fill(30.0,weight);
                mkStatSect(tots,31);
                }
             }

            if(hh.mass() >= 1400 && hh.mass() < 1800){
                tots->Fill(37.0,weight);
                mkStatSect(tots,38);
                if(hBBTM){
                tots->Fill(44.0,weight);
                mkStatSect(tots,45);
                }
             }

        };



        makeStats(smpName+"_emu");
        if(selectedLepton && selectedLepton->isMuon()) makeStats(smpName+"_mu");
        if(selectedLepton && !selectedLepton->isMuon()) makeStats(smpName+"_e");

        if(isSignal() && diHiggsEvt.type >= DiHiggsEvent::MU){
            makeStats(smpName+"_genemu_emu");
            if( diHiggsEvt.type == DiHiggsEvent::MU ) makeStats(smpName+"_genmu_mu");
            if( diHiggsEvt.type == DiHiggsEvent::E ) makeStats(smpName+"_gene_e");
        }

        if(!passBase) return false;

        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!passWjjSel) return false;
        if(!hbbCand) return false;
        if(!passWlnuDR) return false;
        if(!passWWDM) return false;




        auto pltSet = [&](const TString& prefix){
            plotter.getOrMake1DPre(prefix,"hh_mass",";HH mass [TeV]; N. events / 5 GeV",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(hBBTM) plotter.getOrMake1DPre(prefix + "_hbbTM","hh_mass",";HH mass [TeV]; N. events / 5 GeV",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);

            plotter.getOrMake1DPre(prefix, "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 700 && hh.mass() < 900)
                plotter.getOrMake1DPre(prefix+"_hh700to900", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 900 && hh.mass() <1100)
                plotter.getOrMake1DPre(prefix+"_hh900to1100", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 1400 && hh.mass() < 1800)
                plotter.getOrMake1DPre(prefix+"_hh1400to1800", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 2500 && hh.mass() < 3500)
                plotter.getOrMake1DPre(prefix+"_hh2500to3000", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
        };

        auto pltRecoSet = [&](const TString& prefix){
            pltSet(prefix);
            if(passHbbSel) pltSet(prefix +"_LI");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MF) pltSet(prefix +"_L");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_ML) pltSet(prefix +"_M");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MM) pltSet(prefix +"_T");
        };

        auto pltSignalBTagSet = [&](const TString& prefix){
            plotter.getOrMake2DPre(prefix, "sigMass_hbbTag" ,";m(X) mass [TeV];H(bb) tagging; arbitrary units",40,0.550,4.550,7,-0.5,6.5)->Fill(float(signal_mass)/1000.,float(hbbCSVCat),weight);
            if(hBBTM)
                plotter.getOrMake2DPre(prefix+"_hbbTM", "sigMass_hbbTag" ,";m(X) mass [TeV];H(bb) tagging; arbitrary units",40,0.550,4.550,7,-0.5,6.5)->Fill(float(signal_mass)/1000.,float(hbbCSVCat),weight);
        };
        auto pltBKGBTagSet = [&](const TString& prefix){
            plotter.getOrMake2DPre(prefix, "hhMass_hbbTag" ,";m(X) mass [TeV];H(bb) tagging; arbitrary units",45,0.500,5.000,7,-0.5,6.5)->Fill(hh.mass()/1000.,float(hbbCSVCat),weight);
            if(hBBTM)
                plotter.getOrMake2DPre(prefix+"_hbbTM", "hhMass_hbbTag" ,";m(X) mass [TeV];H(bb) tagging; arbitrary units",45,0.500,5.000,7,-0.5,6.5)->Fill(hh.mass()/1000.,float(hbbCSVCat),weight);
        };

        pltRecoSet(smpName + "_emu");
        if(selectedLepton->isMuon()) pltRecoSet(smpName+"_mu");
        else pltRecoSet(smpName+"_e");

        pltBKGBTagSet(smpName + "_emu");
        if(selectedLepton->isMuon()) pltBKGBTagSet(smpName+"_mu");
        else pltBKGBTagSet(smpName+"_e");

        if(isSignal()){
            pltSignalBTagSet("signal_emu");
            if(selectedLepton->isMuon()) pltSignalBTagSet("signal_mu");
            else pltSignalBTagSet("signal_e");
        }


        if(isSignal() && diHiggsEvt.type >= DiHiggsEvent::MU){
            pltRecoSet(smpName+"_genemu_emu");
            if( diHiggsEvt.type == DiHiggsEvent::MU)if(selectedLepton->isMuon()) pltRecoSet(smpName+"_genmu_mu");
            if( diHiggsEvt.type == DiHiggsEvent::E)if(!selectedLepton->isMuon()) pltRecoSet(smpName+"_gene_e");
        }




        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

    std::unique_ptr<LeptonProcessor>     leptonNoISOProc ;




};

#endif

void getCutflow(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getCutflow(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
