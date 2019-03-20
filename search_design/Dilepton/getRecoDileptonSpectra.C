#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
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

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2) {
    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);

    	plotter.getOrMake1DPre(sn,"SIP1",";SIP",100,0,20)->Fill(recolep1->sip3D(),weight);
    	plotter.getOrMake1DPre(sn,"SIP2",";SIP",100,0,20)->Fill(recolep2->sip3D(),weight);
    	plotter.getOrMake1DPre(sn,"dz1",";|#Deltaz|",100,0,2)->Fill(abs(recolep1->dz()),weight);
    	plotter.getOrMake1DPre(sn,"dz2",";|#Deltaz|",100,0,2)->Fill(abs(recolep2->dz()),weight);
    	plotter.getOrMake1DPre(sn,"d01",";|#Delta0|",100,0,2)->Fill(abs(recolep1->d0()),weight);
    	plotter.getOrMake1DPre(sn,"d02",";|#Delta0|",100,0,2)->Fill(abs(recolep2->d0()),weight);
    	plotter.getOrMake1DPre(sn,"miniIso1",";miniIso",100,0,1)->Fill(abs(recolep1->miniIso()),weight);
    	plotter.getOrMake1DPre(sn,"miniIso2",";miniIso",100,0,1)->Fill(abs(recolep2->miniIso()),weight);

    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);
    }
    void printDebugInfo(TString sn, const GenParticle* genlep1, const GenParticle* genlep2, int idx1, int idx2, std::vector<const Lepton*> leps) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
		std::cout << sn << std::endl;
		printf("gen1 = %i --> idx1 = %i; gen2 = %i --> idx2 = %i\n",genlep1->pdgId(),idx1,genlep2->pdgId(),idx2);
    	for (const auto& lep : leps) {
    		printf("lepton %d (%i): (E= %f pT=   %f eta= %f phi=  %f)\n",lep->index(),lep->isMuon() ? (-1)*lep->q()*13:(-1)*lep->q()*11,lep->E(),lep->pt(),lep->eta(),lep->phi());
    	}
		printf("\n");
    }

    const Lepton* getMatchedLepton(const GenParticle* genLep,const std::vector<const Muon*> muons, const std::vector<const Electron*> electrons, double maxDR, bool chargeMatch) {
       if(genLep->absPdgId() == ParticleInfo::p_muminus){
           double nearestDR =10;
           int idx = -1;
       	   for (const auto& mu : muons) {
       	       if (chargeMatch) {
       		       if ((mu->q() > 0) == (genLep->pdgId() > 0)) continue;
       	       }
       	       double dr = PhysicsUtilities::deltaR(*genLep,*mu);
       	       if (dr < nearestDR) {
       	    	   nearestDR = dr;
           	       idx = mu->index();
       	       }
       	   }
       	   if (nearestDR > maxDR) return 0;
       	   if (idx < 0) return 0;
       	   return &reader_muon->muons[idx];
       } else {
           double nearestDR =10;
           int idx = -1;
       	   for (const auto& el : electrons) {
       	       if (chargeMatch) {
       		       if ((el->q() > 0) == (genLep->pdgId() > 0)) continue;
       	       }
       	       double dr = PhysicsUtilities::deltaR(*genLep,*el);
       	       if (dr < nearestDR) {
       	    	   nearestDR = dr;
           	       idx = el->index();
       	       }
       	   }
       	   if (nearestDR > maxDR) return 0;
       	   if (idx < 0) return 0;
       	   return &reader_electron->electrons[idx];
       }
    }
    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 400) return false;
        TString sn = smpName;

        // SIGNAL
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type == DiHiggsEvent::DILEP) {

			plotter.getOrMake1DPre(sn+"_full_signal","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
            const GenParticle* lep1 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
            const GenParticle* lep2 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;
			int lep1id = lep1->pdgId();
			int lep2id = lep2->pdgId();
			if (lep1id<0 == lep2id<0) return false;
			if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
			else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "_ee";
			else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "_mumu";
			else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "_emu";
			else std::cout<<"Error: d1 not a charged lepton"<<std::endl;

			plotter.getOrMake1DPre(sn+"_gen_dilep_notau","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			// GEN lepton pt cuts
            bool passPt1 = lep1->absPdgId() == 13 ? (lep1->pt() > 26) : (lep1->pt() > 30);
            bool passPt2 = (lep2->pt() > 10);
            if (!(passPt1 && passPt2)) return false;

			plotter.getOrMake1DPre(sn+"_gen_dilep_passPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            // get the matched RECO Dileptons
        	// WARNING: have not implemented any check to ensure that these are different
            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
            const Lepton* matchLep1 = getMatchedLepton(lep1,muons,electrons,0.1,true);
            const Lepton* matchLep2 = getMatchedLepton(lep2,muons,electrons,0.1,true);

            if (!(matchLep1 && matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"_foundLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			if ((matchLep1->isMuon() == matchLep2->isMuon()) && (matchLep1->index() == matchLep2->index())) return false; // discard if these are the same RECO lep
			sn += "_pt2_10_";
			plotter.getOrMake1DPre(sn+"goodLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			plotSpectra(sn,matchLep1,matchLep2);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
        	plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
    }
    HistGetter plotter;
};

#endif

void getRecoDileptonSpectra(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getRecoDileptonSpectra(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
