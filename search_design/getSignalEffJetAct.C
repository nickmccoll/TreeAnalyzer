
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
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
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/EventSelection/interface/EventSelection.h"

using namespace TAna;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){

        leptonProcNoISO .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProcNoISO);
        leptonProcNoISO->lepSelParams.el_maxISO =-1;
        leptonProcNoISO->lepSelParams.mu_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.el_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.mu_maxISO =-1;
    }

    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
    }

    const Lepton * getMatchedLepton(const GenParticle& genLepton,const std::vector<const Muon *> muons, const std::vector<const Electron*> electrons){
        if(genLepton.absPdgId() == ParticleInfo::p_muminus){
            double nearestDR =10;
            int idx = PhysicsUtilities::findNearestDRDeref(genLepton,muons,nearestDR,0.2);
            if(idx < 0) return 0;
            else return muons[idx];
        } else {
            double nearestDR =10;
            int idx = PhysicsUtilities::findNearestDRDeref(genLepton,electrons,nearestDR,0.2);
            if(idx < 0) return 0;
            else return electrons[idx];
        }
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        TString sN;
        if(reader_event->process == FillerConstants::SIGNAL){
            if(diHiggsEvt.type != DiHiggsEvent::TAU_MU && diHiggsEvt.type != DiHiggsEvent::MU) return false;

            if(diHiggsEvt.type == DiHiggsEvent::TAU_MU || diHiggsEvt.type == DiHiggsEvent::MU) sN = smpName + "_mu_";
            else if (diHiggsEvt.type == DiHiggsEvent::TAU_E || diHiggsEvt.type == DiHiggsEvent::E) sN = smpName + "_e_";
            else {printf("WARNING1\n"); return false;}
        }
        const std::vector<const Muon    *> muons     = PhysicsUtilities::selObjsMom(reader_muon->muons,26,2.4);
        const std::vector<const Electron*> electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,30,2.4);

        // get reco lepton matched to diHiggsEvt genLepton
        const auto* recoL = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);
        if (recoL == 0) return false;

        // grab leptons that pass ID (does this include IP cuts??? check, Brent)
    	const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
    	if (!passID) return false;

    	Double_t dR_binedges[] = {0, 0.2, 0.4, 0.55, 0.7, 0.8, 1.0, 1.25, 1.5, 3.0, 3.5};
    	Double_t ptratio_binedges[] = {0.0, 0.2, 0.3, 0.4, 0.5, 2.0, 2.5};
    	Double_t pt_binedges[] = {20, 25, 30, 40, 50, 60, 120, 200, 250};
    	Double_t abseta_binedges[] = {0, 0.9, 1.2, 2.1, 2.4, 2.7};

    	float dRnorm = recoL->dRnorm();
    	float ptratio = recoL->PtRatioLepAct();
    	float ptlep = recoL->pt();
    	float absetalep = abs(recoL->eta());

    	// these lines create a visible overflow set of bins
    	if (dRnorm > 3) dRnorm = 3.1;
    	if (ptratio > 2) ptratio = 2.1;
    	if (ptlep > 200) ptlep = 210;
    	if (absetalep > 2.4) absetalep = 2.5;

    	// fill histogram of leptons passing ID
    	plotter.getOrMake2DPre(sN+"jetact", "passID","; #DeltaR(JetAct, lep) / R_{miniIso};P_{T} (JetAct) / P_{T} (lep)", 10, dR_binedges, 6, ptratio_binedges)->Fill(dRnorm, ptratio, weight);
    	plotter.getOrMake2DPre(sN+"pteta", "passID","; P_{T}; |#eta|", 8, pt_binedges, 5, abseta_binedges)->Fill(ptlep, absetalep, weight);

		// fill histogram of leptons passing Isolation
		const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);
		if (passISO) {
			plotter.getOrMake2DPre(sN+"jetact", "PASSiso","; #DeltaR(JetAct, lep) / R_{miniIso};P_{T} (JetAct) / P_{T} (lep)", 10, dR_binedges, 6, ptratio_binedges)->Fill(dRnorm, ptratio, weight);
	    	plotter.getOrMake2DPre(sN+"pteta", "PASSiso","; P_{T}; |#eta|", 8, pt_binedges, 5, abseta_binedges)->Fill(ptlep, absetalep, weight);
		} else {
			plotter.getOrMake2DPre(sN+"jetact", "FAILiso","; #DeltaR(JetAct, lep) / R_{miniIso};P_{T} (JetAct) / P_{T} (lep)", 10, dR_binedges, 6, ptratio_binedges)->Fill(dRnorm, ptratio, weight);
	    	plotter.getOrMake2DPre(sN+"pteta", "FAILiso","; P_{T}; |#eta|", 8, pt_binedges, 5, abseta_binedges)->Fill(ptlep, absetalep, weight);
		}

		//some printouts for study
		if (ptratio > 0.5 && (dRnorm<0.5 || dRnorm>1.5)) {
			if (passISO) {
				std::cout<<std::endl;
				printf("Event PASSES isolation ----> dRnorm = %.2f \n", dRnorm);
				std::cout<<"run:lumi:event = "<< reader_event->run<<":"<<reader_event->lumi<<":"<<reader_event->event<<std::endl;
				std::cout<<std::endl;
			} else {
				std::cout<<std::endl;
				printf("Event FAILS isolation ----> dRnorm = %.2f \n", dRnorm);
				std::cout<<"run:lumi:event = "<< reader_event->run<<":"<<reader_event->lumi<<":"<<reader_event->event<<std::endl;
				std::cout<<std::endl;
			}
		}

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;

};

#endif

void getSignalEffJetAct(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getSignalEffJetAct(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
