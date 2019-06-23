
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/DileptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Processors/Variables/interface/FatJetSelection.h"


using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    bool passGEN_requirements(const GenParticle* lep, const GenParticle* q1, const GenParticle* q2) {
    	const GenParticle* b1 = diHiggsEvt.b1;
    	const GenParticle* b2 = diHiggsEvt.b2;

    	if(b1 && b2 && q1 && q2 && lep) {
        	if(!ParticleInfo::isQuark(q1->pdgId()) || !ParticleInfo::isQuark(q2->pdgId())) {
        		cout << "one of the quarks is not actually a quark" << endl;
        		cout << "q1 and q2 = " << q1->pdgId() << ", " << q2->pdgId() << endl;
        		cout << "lep and nu = " << lep->pdgId() << ", " << diHiggsEvt.w1_d2->pdgId() << endl << endl;
        		return false;
        	}
        	if(!ParticleInfo::isLepton(lep->pdgId())) {
        		cout << "lep is not a lepton" << endl;
        		return false;
        	}

        	if (PhysicsUtilities::deltaR2(*b1,*b2) > 0.8*0.8) return false;
        	if (PhysicsUtilities::deltaR2(*q1,*q2) > 0.8*0.8) return false;
        	if (PhysicsUtilities::deltaR2((q1->p4()+q2->p4()),*lep) > 1.2*1.2) return false;
        	return true;

    	} else return false;
    }

    void plot(TString sn, const GenParticle* lep1, const GenParticle* lep2) {

    	const GenParticle *l1, *l2;
    	if (lep1->pt() > lep2->pt()) {
    		l1 = lep1;
    		l2 = lep2;
    	} else {
    		l1 = lep2;
    		l2 = lep1;
    	}

        plotter.getOrMake1DPre(sn,"pt",";pt",500,0,3000)->Fill(lep1->pt(),weight);
        plotter.getOrMake1DPre(sn,"pt",";pt",500,0,3000)->Fill(lep2->pt(),weight);
        plotter.getOrMake1DPre(sn,"pt1",";pt",500,0,3000)->Fill(l1->pt(),weight);
        plotter.getOrMake1DPre(sn,"pt2",";pt",500,0,3000)->Fill(l2->pt(),weight);

        double drll = PhysicsUtilities::deltaR(*l1,*l2);
        double dphill = PhysicsUtilities::absDeltaPhi(*l1,*l2);
        double mll  = (l1->p4()+l2->p4()).mass();

        plotter.getOrMake1DPre(sn,"drll",";#DeltaR",100,0,5)->Fill(drll,weight);
        plotter.getOrMake1DPre(sn,"dphill",";#Delta#phi",100,0,3.14)->Fill(dphill,weight);
        plotter.getOrMake1DPre(sn,"mll",";M_{ll}",50,0,200)->Fill(mll,weight);

    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;

        if (!isSignal()) return false;
        if (diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        plot(prefix,diHiggsEvt.w1_d1,diHiggsEvt.w2_d1);


        return true;
    }

    std::unique_ptr<FatJetProcessor>        proc     ;
    size64 triggerAccepts = 0;
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getGenVars(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
