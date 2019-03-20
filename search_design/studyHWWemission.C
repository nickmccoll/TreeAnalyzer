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

#include "TLorentzVector.h"
#include "TVector3.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }

    void plotSpectra(TString sn) {

    	//gen
    	TLorentzVector genhww(diHiggsEvt.hww->px(),diHiggsEvt.hww->py(),diHiggsEvt.hww->pz(),diHiggsEvt.hww->E());
    	TLorentzVector genwqq;
    	TLorentzVector genwlnu;
    	const GenParticle *gennu;

    	if (ParticleInfo::isHadron(diHiggsEvt.w1_d1->pdgId())) {
    		if (diHiggsEvt.w1) {genwqq.SetPxPyPzE(diHiggsEvt.w1->px(), diHiggsEvt.w1->py(), diHiggsEvt.w1->pz(), diHiggsEvt.w1->E());}
    		else {
    			TLorentzVector v1(diHiggsEvt.w1_d1->px(), diHiggsEvt.w1_d1->py(), diHiggsEvt.w1_d1->pz(), diHiggsEvt.w1_d1->E());
    			TLorentzVector v2(diHiggsEvt.w1_d2->px(), diHiggsEvt.w1_d2->py(), diHiggsEvt.w1_d2->pz(), diHiggsEvt.w1_d2->E());
    			genwqq = v1+v2;
    		}
    		if (diHiggsEvt.w2) {genwlnu.SetPxPyPzE(diHiggsEvt.w2->px(), diHiggsEvt.w2->py(), diHiggsEvt.w2->pz(), diHiggsEvt.w2->E());}
    		else {
    			TLorentzVector v1(diHiggsEvt.w2_d1->px(), diHiggsEvt.w2_d1->py(), diHiggsEvt.w2_d1->pz(), diHiggsEvt.w2_d1->E());
    			TLorentzVector v2(diHiggsEvt.w2_d2->px(), diHiggsEvt.w2_d2->py(), diHiggsEvt.w2_d2->pz(), diHiggsEvt.w2_d2->E());
    			genwlnu = v1+v2;
    		}
    		gennu = diHiggsEvt.w2_d2;
    	} else {
        	if (diHiggsEvt.w1) {genwlnu.SetPxPyPzE(diHiggsEvt.w1->px(), diHiggsEvt.w1->py(), diHiggsEvt.w1->pz(), diHiggsEvt.w1->E());}
        	else {
        		TLorentzVector v1(diHiggsEvt.w1_d1->px(), diHiggsEvt.w1_d1->py(), diHiggsEvt.w1_d1->pz(), diHiggsEvt.w1_d1->E());
        		TLorentzVector v2(diHiggsEvt.w1_d2->px(), diHiggsEvt.w1_d2->py(), diHiggsEvt.w1_d2->pz(), diHiggsEvt.w1_d2->E());
        		genwlnu = v1+v2;
        	}
        	if (diHiggsEvt.w2) {genwqq.SetPxPyPzE(diHiggsEvt.w2->px(), diHiggsEvt.w2->py(), diHiggsEvt.w2->pz(), diHiggsEvt.w2->E());}
        	else {
        		TLorentzVector v1(diHiggsEvt.w2_d1->px(), diHiggsEvt.w2_d1->py(), diHiggsEvt.w2_d1->pz(), diHiggsEvt.w2_d1->E());
        		TLorentzVector v2(diHiggsEvt.w2_d2->px(), diHiggsEvt.w2_d2->py(), diHiggsEvt.w2_d2->pz(), diHiggsEvt.w2_d2->E());
        		genwqq = v1+v2;
        	}
        	gennu = diHiggsEvt.w1_d2;
    	}
    	TVector3 genboost(genhww.Px()/genhww.E(),genhww.Py()/genhww.E(),genhww.Pz()/genhww.E());

    	TLorentzVector genwqq_rest(genwqq);
    	TLorentzVector genwlep_rest(genwlnu);
    	genwqq_rest.Boost((-1)*genboost);
    	genwlep_rest.Boost((-1)*genboost);

    	double phi = genwqq_rest.Vect().DeltaPhi(genboost);
    	plotter.getOrMake1DPre(sn+"_gen","deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);
    	plotter.getOrMake1DPre(sn+"_gen","had_deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);

    	phi = genwlep_rest.Vect().DeltaPhi(genboost);
    	plotter.getOrMake1DPre(sn+"_gen","deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);
    	plotter.getOrMake1DPre(sn+"_gen","lep_deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);

    	//reco
    	TLorentzVector wqq(wjjCand->px(),wjjCand->py(),wjjCand->pz(),wjjCand->E());
    	TLorentzVector wlep(wlnu.px(),wlnu.py(),wlnu.pz(),wlnu.E());

    	TVector3 boost(hWW.px()/hWW.E(), hWW.py()/hWW.E(), hWW.pz()/hWW.E());

    	TLorentzVector wqq_rest(wqq);
    	TLorentzVector wlep_rest(wlep);
    	wqq_rest.Boost((-1)*boost);
    	wlep_rest.Boost((-1)*boost);

    	phi = wqq_rest.Vect().DeltaPhi(boost);
    	plotter.getOrMake1DPre(sn+"_reco","deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);
    	plotter.getOrMake1DPre(sn+"_reco","had_deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);

    	phi = wlep_rest.Vect().DeltaPhi(boost);
    	plotter.getOrMake1DPre(sn+"_reco","deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);
    	plotter.getOrMake1DPre(sn+"_reco","lep_deltaphi",";#Delta#phi",40,-3.14,3.14)->Fill(phi,weight);

    	// difference in phi between gen and reco
    	double ddphi_had = genwqq_rest.DeltaPhi(wqq_rest);
    	double ddphi_lep = genwlep_rest.DeltaPhi(wlep_rest);
    	plotter.getOrMake1DPre(sn,"delta_deltaphi_lep",";#Delta#phi",40,-3.14,3.14)->Fill(ddphi_lep,weight);
    	plotter.getOrMake1DPre(sn,"delta_deltaphi_had",";#Delta#phi",40,-3.14,3.14)->Fill(ddphi_had,weight);

    	float wwdr = wqq.DeltaR(wlep);
    	plotter.getOrMake2DPre(sn+"_reco","hwwpt_wwdr",";#DeltaR(W,W);H#rightarrowWW P_{T}",40,0,1.5,100,100,1500)->Fill(wwdr,hWW.pt(),weight);
    	plotter.getOrMake2DPre(sn+"_genreco","MET",";recoMET;genMET",50,0,1000,50,0,1000)->Fill(gennu->Et(),neutrino.Et(),weight);

    	plotter.getOrMake1DPre(sn+"_gen","MET",";MET",50,0,1000)->Fill(gennu->Et(),weight);
    	plotter.getOrMake1DPre(sn+"_reco","MET",";MET",50,0,1000)->Fill(neutrino.Et(),weight);
    	float dphi_met = PhysicsUtilities::deltaPhi(*gennu, neutrino);
    	plotter.getOrMake1DPre(sn,"deltaPhi_MET","#Delta#Phi",40,-4,4)->Fill(dphi_met,weight);

    	//debugging plots
    	double dphi_genww = genwqq_rest.DeltaPhi(genwlep_rest);
    	double dphi_recoww = wqq_rest.DeltaPhi(wlep_rest);
    	double ddphi_ww = abs(dphi_recoww) - abs(dphi_genww);
    	plotter.getOrMake1DPre(sn+"_gen","dphi_genww",";",50,-4,4)->Fill(dphi_genww,weight);
    	plotter.getOrMake1DPre(sn+"_gen","dphi_recoww",";",50,-4,4)->Fill(dphi_recoww,weight);
    	plotter.getOrMake1DPre(sn+"_gen","ddphi_ww",";",50,0,4)->Fill(ddphi_ww,weight);
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";
        TString sn = smpName;

        // take only single-lepton events without tau contributions
        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;
        if (selectedLeptons.size() == 0) return false;

        // require lepton to pass pt WP and to exist a WjjCand
//        if (selectedLepton->pt() < (selectedLepton->isMuon() ? 26:30)) return false;
        if (!wjjCand) return false;

        // impose an HT cut
//        if (ht_chs < 400) return false;

        plotSpectra(sn);
        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void studyHWWemission(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void studyHWWemission(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
