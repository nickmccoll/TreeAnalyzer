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
#include "Processors/Corrections/interface/FatJetScaleFactors.h"

#include "TSystem.h"
using namespace TAna;
using namespace std;


class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
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
    bool passIPcuts(const Lepton* lep1, const Lepton* lep2) {
    	bool pass1 = (fabs(lep1->d0()) < 0.05) && (fabs(lep1->dz()) < 0.1) && (lep1->sip3D() < 4);
    	bool pass2 = (fabs(lep2->d0()) < 0.05) && (fabs(lep2->dz()) < 0.1) && (lep2->sip3D() < 4);

    	if (pass1 && pass2) return true;
    	else {
//    		printf("lep1: d0 = %.3f, dz = %.3f, SIP = %.3f\n",lep1->d0(),lep1->dz(),lep1->sip3D());
//    		printf("lep2: d0 = %.3f, dz = %.3f, SIP = %.3f\n\n",lep2->d0(),lep2->dz(),lep2->sip3D());
    		return false;
    	}
    }
    bool passID(const Lepton *lep1, const Lepton *lep2) {
    	bool pass1 = lep1->isMuon() ? ((const Muon*)lep1)->passMed16ID() : ((const Electron*)lep1)->passMedID_noISO();
    	bool pass2 = lep2->isMuon() ? ((const Muon*)lep2)->passMed16ID() : ((const Electron*)lep2)->passMedID_noISO();

    	if (pass1 && pass2) return true;
    	else return false;
    }

    bool passISO(const Lepton *lep1, const Lepton *lep2) {
    	bool pass1 = (lep1->miniIso() < 0.2);
    	bool pass2 = (lep2->miniIso() < 0.2);

    	if (pass1 && pass2) return true;
    	else return false;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if (ht_chs < 400) return false;
        TString sn = "";

        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type == DiHiggsEvent::DILEP) {

			plotter.getOrMake1DPre(sn+"_full_signal","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
            const GenParticle* lep1 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
            const GenParticle* lep2 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;
			int lep1id = lep1->pdgId();
			int lep2id = lep2->pdgId();
			if (lep1id<0 == lep2id<0) return false;
			if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
			else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "ee_";
			else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "mumu_";
			else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "emu_";
			else std::cout<<"Error: d1 not a charged lepton"<<std::endl;

			plotter.getOrMake1DPre(sn+"gen_dilep_notau","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			// GEN lepton pt cuts
            bool passPt1 = lep1->absPdgId() == 13 ? (lep1->pt() > 26) : (lep1->pt() > 30);
            bool passPt2 = (lep2->pt() > 10);
            if (!(passPt1 && passPt2)) return false;

			plotter.getOrMake1DPre(sn+"gen_dilep_passPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            // get the matched RECO Dileptons
        	// WARNING: have not implemented any check to ensure that these are different
            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
            const Lepton* matchLep1 = getMatchedLepton(lep1,muons,electrons,0.1,true);
            const Lepton* matchLep2 = getMatchedLepton(lep2,muons,electrons,0.1,true);

            if (!(matchLep1 && matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"foundLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			if ((matchLep1->isMuon() == matchLep2->isMuon()) && (matchLep1->index() == matchLep2->index())) return false; // discard if these are the same RECO lep
			plotter.getOrMake1DPre(sn+"goodLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			// RECO lepton pt cuts
			bool passRecoPt1 = matchLep1->isMuon() ? (matchLep1->pt() > 26) : (matchLep1->pt() > 30);
			bool passRecoPt2 = (matchLep2->pt() > 10);
			if (!(passRecoPt1 && passRecoPt2)) return false;
			plotter.getOrMake1DPre(sn+"passRecoPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			if (!passIPcuts(matchLep1,matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"passIP","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			if (!passID(matchLep1,matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"passID","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        	if (!passISO(matchLep1,matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"passISO","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        }

        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getDilepSels(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getDilepSels(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
