
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
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
#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }


    void makeCRPlots(const TString& pre, const bool passHBB, const bool passHBBT,  const float hbbMass, const MomentumF& wjjMom,const MomentumF& hbbMom
            ){

        const MomentumF neutrino = HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + wjjMom.p4()) ) ;
        const MomentumF hh = (selectedLepton->p4() + neutrino.p4()+ wjjMom.p4() +hbbMom.p4()) ;
        const float hhMass = hh.mass();
        const float DR = PhysicsUtilities::deltaR(selectedLepton->p4() + neutrino.p4(), wjjMom);

        auto pltExtras =[&](const TString& pre){
            plotter.getOrMake1DPre(pre, "W_W_dR"  ,";#DeltaR(W,W); arbitrary units",50,0,5.0)->Fill(DR,weight);
            plotter.getOrMake1DPre(pre, "hbb_pt"  ,";H(bb) p_{T} [GeV]; arbitrary units",300,0,3000)->Fill(hbbMom.pt(),weight);
        };
        auto plt = [&](const TString& pre){
            plotter.getOrMake1DPre(pre, "hh_mass" ,";HH mass [GeV]; arbitrary units",500,0,5000)->Fill(hhMass,weight);
            pltExtras(pre);

            if(hhMass >= 900 && hhMass < 1100)  pltExtras(pre +"_hh900to1100");
            if(hhMass >= 1400 && hhMass < 1800) pltExtras(pre +"_hh1400to1800");
            if(hhMass >= 2500 && hhMass < 3500) pltExtras(pre +"_hh2500to3500");


            if(hbbMass >= 90 && hbbMass < 140){
                plotter.getOrMake1DPre(pre + "_hbb90to140", "hh_mass" ,";HH mass [GeV]; arbitrary units",500,0,5000)->Fill(hhMass,weight);
                pltExtras(pre +"_hbb90to140");

                if(hhMass >= 900 && hhMass < 1100) pltExtras(pre +"_hbb90to140_hh900to1100");
                if(hhMass >= 1400 && hhMass < 1800)  pltExtras(pre +"_hbb90to140_hh1400to1800");
                if(hhMass >= 2500 && hhMass < 3500)  pltExtras(pre +"_hbb90to140_hh2500to3500");
            }
        };

        if(passHBB && !passHBBT){ plt(pre + "_hbbL");}
        if(passHBBT){ plt(pre + "_hbbT");}
        plt(pre);

    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size() != 1) return false;
        if(hbbCand == 0 || wjjCand== 0) return false;
        if(!fjProc->passWjjSel()) return false;

        const bool passHBB = fjProc->passHbbSel();
        const bool passHBBT = fjProc->passHbbSelTightBTag();

        const auto wjjSDRawMom  = wjjCand->rawSdMom();
        const auto wjjSDCorrMom = wjjCand->sdMom();
        const auto hbbSDRawMom  = hbbCand->rawSdMom();
        const auto hbbSDCorrMom = hbbCand->sdMom();

        makeCRPlots(smpName+"_sdCorr", passHBB,passHBBT,hbbSDCorrMom.mass(),wjjSDCorrMom.p4(),hbbSDCorrMom.p4());
        makeCRPlots(smpName+"_sdRaw", passHBB,passHBBT,hbbSDCorrMom.mass(),wjjSDRawMom.p4(),hbbSDRawMom.p4());
        makeCRPlots(smpName+"_fj", passHBB,passHBBT,hbbSDCorrMom.mass(),wjjCand->p4(),hbbCand->p4());
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void checkBosonMom(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void checkBosonMom(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
