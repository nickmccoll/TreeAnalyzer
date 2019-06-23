#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"

using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){

//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    bool runEvent() override {

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
//        if(!passEventFilters) return false;
        float mtt = float(reader_genpart->mTTBar.val());

        plotter.getOrMake1D("mtt",";M_{tt}",3000,0,3000)->Fill(mtt);
        plotter.getOrMake1D("nGenTops",";N_{gen}",5,-0.5,4.5)->Fill(smDecayEvt.topDecays.size());

        float genMtt = 0;
        if (smDecayEvt.topDecays.size() == 2) {
        	genMtt = (smDecayEvt.topDecays[0].top->p4() + smDecayEvt.topDecays[1].top->p4()).mass();

        	int nLepDecays = 0;
        	if (smDecayEvt.topDecays[0].type > TopDecay::HAD) nLepDecays++;
        	if (smDecayEvt.topDecays[1].type > TopDecay::HAD) nLepDecays++;

        	double sgn = reader_event->weight > 0.0 ? 1.0 : -1.0;
        	double defWt = sgn * xsec() * parameters.event.lumi * 1000 / nSampEvt();

            plotter.getOrMake1D("nLepTops",";N_{lep}",3,-0.5,2.5)->Fill(nLepDecays);

            TString nlep = TString::Format("%dlep",nLepDecays);
            plotter.getOrMake1D("mtt_hand_"+nlep,";M_{tt}",5000,0,5000)->Fill(genMtt,defWt);
            plotter.getOrMake1D("mtt_hand",";M_{tt}",5000,0,5000)->Fill(genMtt,defWt);

            plotter.getOrMake1D("mhh_defwt",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),defWt);
            plotter.getOrMake1D("mhh_norwt",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),weight);
            if (genMtt < 710) {
                plotter.getOrMake1D("mhh_defwt_0to710",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),defWt);
                plotter.getOrMake1D("mhh_norwt_0to710",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),weight);
            } else if (genMtt < 960) {
                plotter.getOrMake1D("mhh_defwt_710to960",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),defWt);
                plotter.getOrMake1D("mhh_norwt_710to960",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),weight);
            } else if (genMtt < 1010) {
                plotter.getOrMake1D("mhh_defwt_960to1010",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),defWt);
                plotter.getOrMake1D("mhh_norwt_960to1010",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),weight);
            } else {
                plotter.getOrMake1D("mhh_defwt_1010toInf",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),defWt);
                plotter.getOrMake1D("mhh_norwt_1010toInf",";M_{HH}",100,700,4000)->Fill(hh_2l.mass(),weight);
            }
        }


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    size64 triggerAccepts=0;

};

#endif

void getMttXSec(std::string fileName, int treeInt,  int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
