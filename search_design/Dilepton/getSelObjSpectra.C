#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
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
    void plotBySelection(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	plot(sn,lep1,lep2,hbb);
    	if (PhysicsUtilities::deltaR(*lep1,*lep2) < 1.6 && nMedBTags_HbbV == 0) plot(sn+"dr_1p6_bVetoAK4_",lep1,lep2,hbb);
    }

    void plotWithID(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	bool passM_e1, passM_e2, passH_e1, passH_e2;
    	if (sn.Contains("_ee_")) {
    		passM_e1 = ((const Electron*)lep1)->passMedID_noISO();
    		passM_e2 = ((const Electron*)lep2)->passMedID_noISO();
    		passH_e1 = ((const Electron*)lep1)->passHEEPID_noISO();
    		passH_e2 = ((const Electron*)lep2)->passHEEPID_noISO();

    		if (passM_e1 && passM_e2) {
    			plotBySelection(sn+"ID_MM_",lep1,lep2,hbb);
    			if (lep1->Et() > 35) plot(sn+"ID_MM_etgt35",lep1,lep2,hbb);
    		}
    		if ((passM_e1 && passH_e2) || (passM_e2 && passH_e1)) plotBySelection(sn+"ID_HM_",lep1,lep2,hbb);
    	} else plotBySelection(sn,lep1,lep2,hbb);
    }
    void plotBjets(TString sn) {
    	int numB_ak8 = 0;
    	for (const auto& fj : reader_fatjet->jets) {
    		bool btag = false;
    		for (const auto& sj : fj.subJets()) {
    			if (sj.csv() > 0.8484) btag = true;
    		}
    		if (btag) numB_ak8++;
    	}
    	plotter.getOrMake1DPre(sn,"numB_ak4_sepHbb",";numB",5,0,5)->Fill(nMedBTags_HbbV,weight);
    	plotter.getOrMake1DPre(sn,"numB_ak8",";numB",5,0,5)->Fill(numB_ak8,weight);
    }
    void plot(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	const MomentumF dilepmom = lep1->p4() + lep2->p4();
    	if (dilepmom.mass() < 12 || dilepmom.mass() > 75) return; // mass cut on 12 < Mll < 75

    	const MomentumF bbllmom = dilepmom.p4() + hbb->p4();
    	const MomentumF metLL = dilepmom.p4() + reader_event->met.p4();
    	double dR_ll = PhysicsUtilities::deltaR(*lep1,*lep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*lep1,*lep2);
    	double dR_bbll = PhysicsUtilities::deltaR(dilepmom,*hbb);
    	double dPhi_bbll = PhysicsUtilities::deltaPhi(dilepmom,*hbb);
    	double dPhi_metLL = PhysicsUtilities::deltaPhi(reader_event->met,dilepmom);

    	plotter.getOrMake1DPre(sn,"Mll",";M_{ll}",100,0,200)->Fill(dilepmom.mass(),weight);
    	plotter.getOrMake1DPre(sn,"Mbb",";M_{bb}",100,30,210)->Fill(hbbMass,weight);
    	plotter.getOrMake1DPre(sn,"Mbbll",";M_{bbll}",200,0,3000)->Fill(bbllmom.mass(),weight);
    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",100,0,4)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";|#Delta#Phi_{ll}|",50,0,3.14)->Fill(abs(dPhi_ll),weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(lep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(lep2->pt(),weight);
    	plotter.getOrMake1DPre(sn,"met",";MET",50,0,2000)->Fill(reader_event->met.Et(),weight);
    	plotter.getOrMake1DPre(sn,"metLL_Et",";GeV",50,0,2000)->Fill(metLL.Et(),weight);
    	plotter.getOrMake1DPre(sn,"metLL_Pt",";GeV",50,0,2000)->Fill(metLL.pt(),weight);
    	plotter.getOrMake1DPre(sn,"ptbb",";p_{T} bb",150,0,1500)->Fill(hbb->pt(),weight);
    	plotter.getOrMake1DPre(sn,"dR_bbll",";#DeltaR_{bb,ll}",100,0,6)->Fill(dR_bbll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_bbll",";#DeltaPhi_{bb,ll}",100,-3.14,3.14)->Fill(dPhi_bbll,weight);
    	plotter.getOrMake1DPre(sn,"ht",";H_{T}",100,0,3000)->Fill(ht_puppi,weight);
    	plotter.getOrMake2DPre(sn,"metLL_x_ht",";MET+ll;H_{T}",100,0,2000,100,0,3000)->Fill(metLL.pt(),ht_puppi,weight);
    	plotter.getOrMake1DPre(sn,"metLL_o_ht",";",100,0,5)->Fill(metLL.pt()/ht_puppi,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_metLL",";|#Delta#Phi_{ll,miss}|",50,0,3.14)->Fill(abs(dPhi_metLL),weight);
    	plotter.getOrMake1DPre(sn,"met_frac_scal",";met / (met pt + ptll)",50,0,1)->Fill(reader_event->met.pt()/(reader_event->met.pt()+dilepmom.pt()),weight);
    	plotter.getOrMake1DPre(sn,"met_frac_vect",";met / (met+ll) pt",100,0,4)->Fill(reader_event->met.pt()/(metLL.pt()),weight);

    	plotBjets(sn);
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
        if (ht_puppi < 400) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName+getDilepChan(selectedDileptons[0],selectedDileptons[1]);

        plotWithID(sn,selectedDileptons[0],selectedDileptons[1],hbbCand);
        return true;
    }
    HistGetter plotter;
    void write(TString fileName){
    	plotter.write(fileName);
    }
};

#endif

void getSelObjSpectra(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getSelObjSpectra(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
