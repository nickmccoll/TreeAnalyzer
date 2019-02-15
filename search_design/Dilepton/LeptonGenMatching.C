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

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    void plot(TString sn) {
    	plotter.getOrMake1DPre(sn,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    }
    void debug(TString sn) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
    	std::cout<<sn<<std::endl;
    	for (const auto& e : reader_electron->electrons) {
    		printf("electron %d (%i): (E= %f pT=   %f eta= %f phi=  %f)\n",e.index(),e.isMuon() ? (-1)*e.q()*13:(-1)*e.q()*11,e.E(),e.pt(),e.eta(),e.phi());
    	}
    	for (const auto& mu : reader_muon->muons) {
    		printf("muon %d (%i): (E= %f pT=   %f eta= %f phi=  %f)\n",mu.index(),mu.isMuon() ? (-1)*mu.q()*13:(-1)*mu.q()*11,mu.E(),mu.pt(),mu.eta(),mu.phi());
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
    void testGenMatch(TString sn, const GenParticle* lep1, const GenParticle* lep2) {
        const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
        const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);

        if (muons.size() + electrons.size() < 2) return;
        plot(sn+"_recoLepCheck");

        const Lepton* match1_02 = getMatchedLepton(lep1,muons,electrons,0.2,false);
        const Lepton* match2_02 = getMatchedLepton(lep2,muons,electrons,0.2,false);

        const Lepton* match1_01 = getMatchedLepton(lep1,muons,electrons,0.1,false);
        const Lepton* match2_01 = getMatchedLepton(lep2,muons,electrons,0.1,false);

        const Lepton* match1_02_plusQ = getMatchedLepton(lep1,muons,electrons,0.2,true);
        const Lepton* match2_02_plusQ = getMatchedLepton(lep2,muons,electrons,0.2,true);

        const Lepton* match1_01_plusQ = getMatchedLepton(lep1,muons,electrons,0.1,true);
        const Lepton* match2_01_plusQ = getMatchedLepton(lep2,muons,electrons,0.1,true);

        if (match1_02 && match2_02) plot(sn+"_02");
//        else debug(sn);
        if (match1_01 && match2_01) plot(sn+"_01");
//        else debug(sn);
        if (match1_02_plusQ && match2_02_plusQ) plot(sn+"_02_Q");
//        else debug(sn);
        if (match1_01_plusQ && match2_01_plusQ) plot(sn+"_01_Q");
//        else debug(sn);

        if (match1_02 && match2_02 && !(match1_02_plusQ && match2_02_plusQ)) debug(sn);
        if (match1_01 && match2_01 && !(match1_01_plusQ && match2_01_plusQ)) debug(sn);
    }
    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        TString sn = smpName;
        const GenParticle* lep1 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
        const GenParticle* lep2 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;
        int lep1id = lep1->pdgId();
        int lep2id = lep2->pdgId();

        // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
        if (lep1id<0 == lep2id<0) return false;

        if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
        else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "_ee";
        else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "_mumu";
        else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "_emu";

        // lepton pt cuts
        bool passPt1 = lep1->absPdgId() == 13 ? (lep1->pt() > 26) : (lep1->pt() > 30);
        bool passPt2 = (lep2->pt() > 10);
        if (!(passPt1 && passPt2)) return false;

        sn += "_ept10_mupt10";
        plot(sn);
        testGenMatch(sn,lep1,lep2);
        return true;
    }
    void write(TString fileName){
    	plotter.write(fileName);
    	double ee_den = plotter.getOrMake1DPre(TString::Format("m%i_ee_ept10_mupt10",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double emu_den = plotter.getOrMake1DPre(TString::Format("m%i_emu_ept10_mupt10",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double mumu_den = plotter.getOrMake1DPre(TString::Format("m%i_mumu_ept10_mupt10",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double ee_num01 = plotter.getOrMake1DPre(TString::Format("m%i_ee_ept10_mupt10_01",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double emu_num01 = plotter.getOrMake1DPre(TString::Format("m%i_emu_ept10_mupt10_01",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double mumu_num01 = plotter.getOrMake1DPre(TString::Format("m%i_mumu_ept10_mupt10_01",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double ee_num02 = plotter.getOrMake1DPre(TString::Format("m%i_ee_ept10_mupt10_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double emu_num02 = plotter.getOrMake1DPre(TString::Format("m%i_emu_ept10_mupt10_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double mumu_num02 = plotter.getOrMake1DPre(TString::Format("m%i_mumu_ept10_mupt10_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double ee_num02_Q = plotter.getOrMake1DPre(TString::Format("m%i_ee_ept10_mupt10_02_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double emu_num02_Q = plotter.getOrMake1DPre(TString::Format("m%i_emu_ept10_mupt10_02_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double mumu_num02_Q = plotter.getOrMake1DPre(TString::Format("m%i_mumu_ept10_mupt10_02_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double ee_num01_Q = plotter.getOrMake1DPre(TString::Format("m%i_ee_ept10_mupt10_01_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double emu_num01_Q = plotter.getOrMake1DPre(TString::Format("m%i_emu_ept10_mupt10_01_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
    	double mumu_num01_Q = plotter.getOrMake1DPre(TString::Format("m%i_mumu_ept10_mupt10_01_Q",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();

    	printf("\n");
    	std::cout<<"ee efficiency:"<<std::endl;
    	std::cout<<"01:   "<<ee_num01/ee_den<<std::endl;
    	std::cout<<"02:   "<<ee_num02/ee_den<<std::endl;
    	std::cout<<"01+Q: "<<ee_num01_Q/ee_den<<std::endl;
    	std::cout<<"02+Q: "<<ee_num02_Q/ee_den<<std::endl;
    	printf("\n");
    	std::cout<<"emu efficiency:"<<std::endl;
    	std::cout<<"01:   "<<emu_num01/emu_den<<std::endl;
    	std::cout<<"02:   "<<emu_num02/emu_den<<std::endl;
    	std::cout<<"01+Q: "<<emu_num01_Q/emu_den<<std::endl;
    	std::cout<<"02+Q: "<<emu_num02_Q/emu_den<<std::endl;
    	printf("\n");
    	std::cout<<"mumu efficiency:"<<std::endl;
    	std::cout<<"01:   "<<mumu_num01/mumu_den<<std::endl;
    	std::cout<<"02:   "<<mumu_num02/mumu_den<<std::endl;
    	std::cout<<"01+Q: "<<mumu_num01_Q/mumu_den<<std::endl;
    	std::cout<<"02+Q: "<<mumu_num02_Q/mumu_den<<std::endl;
    	printf("\n");
    }
    HistGetter plotter;
};

#endif

void LeptonGenMatching(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void LeptonGenMatching(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
