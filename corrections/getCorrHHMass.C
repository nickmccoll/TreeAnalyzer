
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "Processors/Variables/interface/FatJetSelection.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/HiggsSolver.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }

    MomentumF getInvisible2(const MomentumF& met, const MomentumF& vis, const double hMass, const MomentumF& lep){
        const double a = hMass*hMass - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
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

            ASTypes::CartLorentzVector neutP(met.x(),met.y(),pos,std::sqrt(met.x()*met.x()+met.y()*met.y()+pos*pos));
            ASTypes::CartLorentzVector neutN(met.x(),met.y(),neg,std::sqrt(met.x()*met.x()+met.y()*met.y()+neg*neg));

            if((neutP + lep.p4()).mass() < (neutN + lep.p4()).mass())metZ = pos;
            else metZ = neg;
        }
        ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
        return MomentumF(neutrino);
    }

    void SRPlots(TString prefix, const MomentumF& neutrino,const MomentumF& wjj,   bool tight) {
        auto wlnu = neutrino.p4() + selectedLepton->p4();
        auto hWW = wlnu + wjj.p4();
        auto hh = selectedLepton->p4() + neutrino.p4()+ wjj.p4() + hbbCand->p4();

        ASTypes::CartLorentzVector vM(reader_event->met.x(),reader_event->met.y(),0,std::sqrt(reader_event->met.x()*reader_event->met.x()+reader_event->met.y()*reader_event->met.y()));
        auto plt = [&](const TString& pre){
            plotter.getOrMake1DPre(pre, "wlnu_mass"   ,";W(l#nu) mass [GeV]; arbitrary units",50,0,500)->Fill(wlnu.mass(),weight);
            plotter.getOrMake1DPre(pre, "wjj_mass"   ,";W(jj) mass [GeV]; arbitrary units",50,0,500)->Fill(wjj.mass(),weight);
            plotter.getOrMake1DPre(pre, "hWW_mass"   ,";HWW mass [GeV]; arbitrary units",50,0,500)->Fill(hWW.mass(),weight);
            plotter.getOrMake1DPre(pre, "hWWv_mass"   ,";HWW mass [GeV]; arbitrary units",50,0,500)->Fill((selectedLepton->p4() + wjj.p4()).mass(),weight);
            plotter.getOrMake1DPre(pre, "hWWvm_mass"   ,";HWW mass [GeV]; arbitrary units",50,0,500)->Fill((selectedLepton->p4() + wjj.p4()+vM ).mass(),weight);
            plotter.getOrMake1DPre(pre, "hh_mass"   ,";HH mass [GeV]; arbitrary units",1000,0,5000)->Fill(hh.mass(),weight);
            plotter.getOrMake1DPre(pre, "wjj_wlnu_dR" ,";#DeltaR(Wjj,Wl#nu); arbitrary units",64,0,3.2)->Fill(PhysicsUtilities::deltaR(wlnu,wjj.p4()),weight);
        };

        if(!tight){ plt(prefix + "_hbbL");}
        if(tight){ plt(prefix + "_hbbT");}
        plt(prefix + "_hbbI");
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!isSignal()) return false;
        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;

        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!passWjjSel || !passHbbSel) return false;

        MomentumF sdW = wjjCand->sdMom();
        MomentumF neutrinoFromWlnu =  HiggsSolver::getInvisible(reader_event->met, *selectedLepton, 80);
        MomentumF neutrinoFromWlnu2 =  HiggsSolver::getInvisible(reader_event->met, *selectedLepton, std::max(125. - sdW.mass(),0. ));
        MomentumF neutrinoFromH    =  HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + wjjCand->p4()) ) ;
        MomentumF neutrinoFromHSD    =  HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + sdW.p4()) ) ;

        MomentumF neutrinoFromH2    =  getInvisible2(reader_event->met,(selectedLepton->p4() + wjjCand->p4()),125,selectedLepton->p4() ) ;



        bool passTight = passHbbTSel;

        SRPlots(smpName + "_nuH2" ,neutrinoFromH2  ,wjjCand->p4(), passTight);
        SRPlots(smpName + "_nuH" ,neutrinoFromH  ,wjjCand->p4(), passTight);
        SRPlots(smpName + "_nuW",neutrinoFromWlnu,wjjCand->p4(), passTight);
        SRPlots(smpName + "_nuW2",neutrinoFromWlnu2,wjjCand->p4(), passTight);
        SRPlots(smpName + "_nuHsdW" ,neutrinoFromHSD  ,sdW, passTight);
        SRPlots(smpName + "_nuWsdW",neutrinoFromWlnu,sdW, passTight);



        return true;


    }

    void write(TString fileName){ plotter.write(fileName);}


    HistGetter plotter;

};

#endif

void getCorrHHMass(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getCorrHHMass(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
