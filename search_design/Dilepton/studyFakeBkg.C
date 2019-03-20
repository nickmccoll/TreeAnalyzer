#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;
using namespace std;

class Analyzer : public DileptonSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event);
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiJet",isRealData());       load(reader_fatjet);
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false); load(reader_jet);
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon);

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");          load(reader_genpart   );
        }
    }
    void debugEvt(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
    	cout<<endl<<sn<<endl;
    	printf("hbbjet (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",hbbCSVCat,hbb->E(),hbb->pt(),hbb->eta(),hbb->phi());
    	printf("lep1 (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",lep1->q(),lep1->E(),lep1->pt(),lep1->eta(),lep1->phi());
    	printf("lep2 (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",lep2->q(),lep2->E(),lep2->pt(),lep2->eta(),lep2->phi());
    	cout<<endl;
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if(hbbCSVCat < BTagging::CSVSJ_MF) return false;
        if (hbbMass < 30 || hbbMass > 210) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName+getDilepChan(selectedDileptons[0],selectedDileptons[1]);

//      std::cout<<std::endl<<sn<<std::endl;
//    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        if (sn.Contains("_ee_")) plotter.getOrMake1DPre(sn,"fakeE","",5,0,5)->Fill(isSignal()?signal_mass:ht_puppi,weight);
        if (sn.Contains("_mumu_")) plotter.getOrMake1DPre(sn,"fakeMu","",5,0,5)->Fill(isSignal()?signal_mass:ht_puppi,weight);
        if (sn.Contains("_emu_")) {
        	bool isFakeMu = false;
        	if (smDecayEvt.bosonDecays.size() == 1) {
        		if (smDecayEvt.bosonDecays.front().type == BosonDecay::W_MU) plotter.getOrMake1DPre(sn,"fakeE","",100,0,3000)->Fill(isSignal()?signal_mass:ht_puppi,weight);
        		else if (smDecayEvt.bosonDecays.front().type == BosonDecay::W_E) plotter.getOrMake1DPre(sn,"fakeMu","",100,0,3000)->Fill(isSignal()?signal_mass:ht_puppi,weight);
        	}
        	if (smDecayEvt.bosonDecays.size() > 1) printf("SM event has > 1 bosons\n");
        }

        return true;
    }
    HistGetter plotter;
    void write(TString fileName){
    	plotter.write(fileName);
    }
};

#endif

void studyFakeBkg(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void studyFakeBkg(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
