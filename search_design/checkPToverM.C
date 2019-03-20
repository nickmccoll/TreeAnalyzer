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

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "Processors/Variables/interface/FatJetSelection.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/HiggsSolver.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"

#include "TSystem.h"
using namespace TAna;

// Macro Description:
// Plotting Higgs (Pt / M_HH) for DiHiggs event in the Single Lepton + Jets channel {X -> HH -> bbWW (lnu jj)}

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    void SRPlots(TString prefix) {

    		// variables for plotting
    		float Xmass = hh.mass();

    		float PToM_hWW = hWW.pt() / Xmass;
    		float PToM_hbb = hbbCand->pt() / Xmass;
    		float PToM_Wlnu = wlnu.pt() / Xmass;
    		float PToM_Wjj = wjjCand->pt() / Xmass;

    		float hWW_rap = hWW.rap();
    		float hbb_rap = hbbCand->rap();
    		float Wjj_rap = wjjCand->rap();
    		float Wlnu_rap = wlnu.rap();

    		// lambda function that will make all desired histograms when called
        auto plt = [&](const TString& pre){
        		plotter.getOrMake1DPre(pre, "hWW_pt_o_hh_mass", ";H#rightarrowWW p_{T} / M_{HH}; a.u.",100,0,1.1)->Fill(PToM_hWW,weight);
        		plotter.getOrMake1DPre(pre, "hbb_pt_o_hh_mass", ";H#rightarrowbb p_{T} / M_{HH}; a.u.",100,0,1.1)->Fill(PToM_hbb,weight);
        		plotter.getOrMake1DPre(pre, "Wlnu_pt_o_hh_mass", ";Wl#nu p_{T} / M_{HH}; a.u.",100,0,1.1)->Fill(PToM_Wlnu,weight);
        		plotter.getOrMake1DPre(pre, "Wjj_pt_o_hh_mass", ";Wqq p_{T} / M_{HH}; a.u.",100,0,1.1)->Fill(PToM_Wjj,weight);

        		plotter.getOrMake1DPre(pre, "hWW_AbsRapidity", ";H#rightarrowWW Abs(y); a.u.",100,0,4.)->Fill(fabs(hWW_rap),weight);
        		plotter.getOrMake1DPre(pre, "hbb_AbsRapidity", ";H#rightarrowbb Abs(y); a.u.",100,0,4.)->Fill(fabs(hbb_rap),weight);
        		plotter.getOrMake1DPre(pre, "Wlnu_AbsRapidity", ";Wl#nu Abs(y); a.u.",100,0,4.)->Fill(fabs(Wlnu_rap),weight);
        		plotter.getOrMake1DPre(pre, "Wjj_AbsRapidity", ";Wqq Abs(y); a.u.",100,0,4.)->Fill(fabs(Wjj_rap),weight);

        		plotter.getOrMake1DPre(pre, "X_Rapidity", ";y_{HH}; a.u.",100,-4.,4.)->Fill(hh.rap(), weight);

        		plotter.getOrMake2DPre(pre, "hWW_Cos_theta_PToM", ";H#rightarrowWW p_{T} / M_{HH}; Cos #theta;",100,0,1.1,40,-1.1,1.1)->Fill(PToM_hWW, TMath::Cos(hWW.theta()), weight);
			plotter.getOrMake2DPre(pre, "hbb_Cos_theta_PToM", ";H#rightarrowbb p_{T} / M_{HH}; Cos #theta;",100,0,1.1,40,-1.1,1.1)->Fill(PToM_hbb, TMath::Cos(hbbCand->theta()), weight);
			plotter.getOrMake2DPre(pre, "Wlnu_Cos_theta_PToM", ";Wl#nu p_{T} / M_{HH}; Cos #theta;",100,0,1.1,40,-1.1,1.1)->Fill(PToM_Wlnu, TMath::Cos(wlnu.theta()), weight);
			plotter.getOrMake2DPre(pre, "Wjj_Cos_theta_PToM", ";Wqq p_{T} / M_{HH}; Cos #theta;",100,0,1.1,40,-1.1,1.1)->Fill(PToM_Wjj, TMath::Cos(wjjCand->theta()), weight);

			plotter.getOrMake2DPre(pre, "hWW_Eta_PToM", ";H#rightarrowWW p_{T} / M_{HH}; #eta;",100,0,1.1,40,-2.5,2.5)->Fill(PToM_hWW, hWW.eta(), weight);
			plotter.getOrMake2DPre(pre, "hbb_Eta_PToM", ";H#rightarrowbb p_{T} / M_{HH}; #eta;",100,0,1.1,40,-2.5,2.5)->Fill(PToM_hbb, hbbCand->eta(), weight);
			plotter.getOrMake2DPre(pre, "Wlnu_Eta_PToM", ";Wl#nu p_{T} / M_{HH}; #eta;",100,0,1.1,40,-2.5,2.5)->Fill(PToM_Wlnu, wlnu.eta(), weight);
			plotter.getOrMake2DPre(pre, "Wjj_Eta_PToM", ";Wqq p_{T} / M_{HH}; #eta;",100,0,1.1,40,-2.5,2.5)->Fill(PToM_Wjj, wjjCand->eta(), weight);

        };

        // plots for different hh invariant mass regions (inclusive and partitioned)
        plt(prefix+"_hhInc");
        if(hh.mass() >= 800) plt(prefix+"_hhgt800");
        if(hh.mass() >= 800 && hh.mass() < 1200) plt(prefix+"_hh800to1200");
        if(hh.mass() >= 2000) plt(prefix+"_hhgt2000");

    }

    bool runEvent() override {

    	// Make plots for the cases where: hbb mass is ON (within) the mass peak and OFF (outside), INClusive
    	   auto make_plots = [&](const TString& pref, float hbbMass) {
    		   SRPlots(pref+"_hbbINC");
    		   if(hbbMass >= 105 && hbbMass < 145){SRPlots(pref+"_hbbON");}
    		   else {SRPlots(pref+"_hbbOFF");}
       };

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";

        // cuts for appropriate lepton+jets channel
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size() != 1) return false;

        if(!wjjCand || !hbbCand) return false;
        if(nMedBTags_HbbV != 0) return false;
//        if(!passHbbSel) return false;

//        if(!passWlnuDR) return false;
        if(!passWWDM) return false;

        // Some code to analyze aspects of W+jets MC sample
/*        if (hh.mass() > 2000 && hbbCSVCat >= BTagging::CSVSJ_MF) {
        		ParticleInfo::printGenInfo(reader_genpart->genParticles, -1);
        		std::cout << "Num fatjets: " << fatjetCands.size() << std::endl;

        		for (const auto& fatjet : fatjetCands) {
        			float fj_ptom = fatjet->pt()/hh.mass();
				std::cout << "FatJet " << fatjet->index() << ": " << *fatjet << std::endl;
				std::cout << "    Rapidity: " << fatjet->rap() << std::endl;
				std::cout << "    Pt / M_hh : " << fj_ptom << std::endl;
        		}

        		std::cout << "Hbb jet: " << hbbCand->index() << " " << *hbbCand << std::endl;
        		std::cout << "Wjj jet: " << wjjCand->index() << " " << *wjjCand << std::endl;
        		std::cout << "hh mass is " << hh.mass() << std::endl;
        		std::cout << " " << std::endl;
        }
*/
        // B tagging cuts

        if(hbbCSVCat < BTagging::CSVSJ_MF) return false;
        make_plots(smpName+"_gtoeLooseB", hbbMass);

        if(hbbCSVCat < BTagging::CSVSJ_ML) return false;
        	make_plots(smpName+"_gtoeMediumB", hbbMass);

        	if(hbbCSVCat < BTagging::CSVSJ_MM) return false;
        	make_plots(smpName+"_gtoeTightB", hbbMass);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}

    HistGetter plotter;

};

#endif

void checkPToverM(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void checkPToverM(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
