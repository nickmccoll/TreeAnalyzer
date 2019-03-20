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


	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

	void studyLeps(TString sn) {
		bool W1isHeavier = (diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).mass() > (diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).mass();
		const GenParticle *genlepON = W1isHeavier ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
		const Lepton *lepON;
		bool highestPtIsON = false;
		if (PhysicsUtilities::deltaR2(*genlepON,*selectedDileptons[0]) < PhysicsUtilities::deltaR2(*genlepON,*selectedDileptons[1])) {
			lepON = selectedDileptons[0];
			highestPtIsON = true;
		} else lepON = selectedDileptons[1];

        plotter.getOrMake1DPre(sn, "all", "",100,0,300)->Fill(selectedDileptons[0]->pt(), weight);

		cout<<"highestPtIsON: "<<highestPtIsON<<endl;
		cout<<"lepON index: "<<lepON->index()<<endl;
		cout<<"lep0 index: "<<selectedDileptons[0]->index()<<endl;
		cout<<"lepON isMuon: "<<lepON->isMuon()<<endl;
		cout<<"lep0 isMuon: "<<selectedDileptons[0]->isMuon()<<endl;
		cout<<endl;
		printf("genlepON: pt = %.3f, eta = %.3f, phi = %.3f\n",genlepON->pt(),genlepON->eta(),genlepON->phi());
		printf("lepON: pt = %.3f, eta = %.3f, phi = %.3f\n",lepON->pt(),lepON->eta(),lepON->phi());
		printf("lep_0: pt = %.3f, eta = %.3f, phi = %.3f\n",selectedDileptons[0]->pt(),selectedDileptons[0]->eta(),selectedDileptons[0]->phi());
		printf("lep_1: pt = %.3f, eta = %.3f, phi = %.3f\n",selectedDileptons[1]->pt(),selectedDileptons[1]->eta(),selectedDileptons[1]->phi());
		cout<<endl;

		if (highestPtIsON) plotter.getOrMake1DPre(sn, "higherPt", "",100,0,300)->Fill(selectedDileptons[0]->pt(), weight);
	}
    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        // current selection
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if(hbbCSVCat < BTagging::CSVSJ_MF) return false;
        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (PhysicsUtilities::deltaR2(*selectedDileptons[0],*selectedDileptons[1]) > 1.6*1.6) return false;
        double mll = (selectedDileptons[0]->p4() + selectedDileptons[1]->p4()).mass();
        if (mll > 75 || mll < 12) return false;
        if (nMedBTags_HbbV != 0) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName+getDilepChan(selectedDileptons[0],selectedDileptons[1]);
        if (sn.Contains("_ee_")) {
        	if (!((const Electron*)selectedDileptons[0])->passMedID_noISO() && ((const Electron*)selectedDileptons[1])->passMedID_noISO()) return false;
        }

        studyLeps(sn);
        return true;
    }
    HistGetter plotter;
    void write(TString fileName) {plotter.write(fileName);}
    void outputResults() {
    	vector<TString> chans = {"_ee_","_mumu_","_emu_"};
    	TString samp = TString::Format("m%d",signal_mass);
    	cout<<endl;
    	for (const auto& ch : chans) {
            double num = plotter.getOrMake1DPre(samp+ch, "higherPt", "",100,0,300)->GetEntries();
            double den = plotter.getOrMake1DPre(samp+ch, "all", "",100,0,300)->GetEntries();

            cout<<samp+ch<<":  "<<num/den<<endl;
    	}
    }

};

#endif

void studyLeptons(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
    a.outputResults();
}
void studyLeptons(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
    a.outputResults();
}
