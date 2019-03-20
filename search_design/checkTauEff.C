
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
        DefaultFatJetSelections::setDefaultFatJetProcessor(noTauW);
        noTauW.wjj_maxT2oT1 = -1;
        DefaultFatJetSelections::setDefaultFatJetProcessor(noTauH);
        noTauH.hbb_maxT2oT1 = -1;
        DefaultFatJetSelections::setDefaultFatJetProcessor(noTauWH);
        noTauWH.wjj_maxT2oT1 = -1;
        noTauWH.hbb_maxT2oT1 = -1;
    }


    void makeCRPlots(const TString& pre, const FatJetParameters& param){
        if(!fjProc->passWjjSel(param)) return;
        bool passHBB = fjProc->passHbbSel(param);
        bool passHBBT = fjProc->passHbbSelTightBTag(param);
        const float hbbMass = hbbCand->sdMom().mass();

        auto plt = [&](const TString& pre){
            plotter.getOrMake1DPre(pre, "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 700 && hh.mass() < 900)
                plotter.getOrMake1DPre(pre, "hh700to900_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 900 && hh.mass() < 1100)
                plotter.getOrMake1DPre(pre, "hh900to1100_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 1400 && hh.mass() < 1800)
                plotter.getOrMake1DPre(pre, "hh1400to1800_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 2500 && hh.mass() < 3000)
                plotter.getOrMake1DPre(pre, "hh2500to3000_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
        };

        if(passHBB && !passHBBT){ plt(pre + "_hbbL");}
        if(passHBBT){ plt(pre + "_hbbT");}


    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size() != 1) return false;
        if(hbbCand == 0 || wjjCand== 0) return false;
        MomentumF Wlnu = selectedLepton->p4() + neutrino.p4();
        if(PhysicsUtilities::deltaR(Wlnu, wjjCand->sdMom().p4()) > 0.5) return false;
        makeCRPlots(smpName+"_tHtW", fjProc->params);
        makeCRPlots(smpName+"_tH"  , noTauW);
        makeCRPlots(smpName+"_tW"  , noTauH);
        makeCRPlots(smpName        , noTauWH);

        if(reader_event->process != FillerConstants::SIGNAL){
            makeCRPlots("bkg_tHtW", fjProc->params);
            makeCRPlots("bkg_tH"  , noTauW);
            makeCRPlots("bkg_tW"  , noTauH);
            makeCRPlots("bkg"     , noTauWH);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

    FatJetParameters noTauW;
    FatJetParameters noTauH;
    FatJetParameters noTauWH;


};

#endif

void checkTauEff(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void checkTauEff(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
