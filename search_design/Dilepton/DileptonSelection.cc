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
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2, const FatJet* hbb) {
    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	MomentumF bbllMOM = dilepMOM.p4() + hbb->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);
    	double dR_bbll = PhysicsUtilities::deltaR(dilepMOM,*hbb);

    	plotter.getOrMake1DPre(sn,"Mll",";m_{ll}",100,0,200)->Fill(dilepMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"Mbbll",";m_{bbll}",200,0,3000)->Fill(bbllMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);
    	plotter.getOrMake1DPre(sn,"ptbb",";p_{T} bb",150,0,1500)->Fill(hbb->pt(),weight);
    	plotter.getOrMake1DPre(sn,"dR_bbll",";#DeltaR_{bb,ll}",100,0,6)->Fill(dR_bbll,weight);
    }
    void printDebugInfo(TString sn, const GenParticle* genlep1, const GenParticle* genlep2, int idx1, int idx2) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
		std::cout << sn << std::endl;
		printf("gen1 = %i --> idx1 = %i; gen2 = %i --> idx2 = %i\n",genlep1->pdgId(),idx1,genlep2->pdgId(),idx2);
    	for (const auto& mu : reader_muon->muons) {
    		printf("muon %d (%i): (E= %f pT= %f eta= %f phi=  %f)\n",mu.index(),mu.q()*(-13),mu.E(),mu.pt(),mu.eta(),mu.phi());
    	}
    	for (const auto& el : reader_electron->electrons) {
    		printf("electron %d (%i): (E= %f pT= %f eta= %f phi=  %f)\n",el.index(),el.q()*(-11),el.E(),el.pt(),el.eta(),el.phi());
    	}
		printf("\n");
    }
    const FatJet* findHbbCand(const Lepton* lep1, const Lepton* lep2) {
    	// lambda function to determine if a FJ is LMT b-tagged (has at least one subjet that passes medium CSV WP)
    	auto isBtag = [&] (const FatJet* fj) {
    		bool hasBtag = false;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.csv() > 0.8484) hasBtag = true;
    		}
    		return hasBtag;
    	};
    	// lambda function to determine if a FJ has two SJs each with pt > 20 and eta < 2.4
    	auto hasGoodSJs = [&] (const FatJet* fj) {
    		bool goodSJs = false;
    		int nGoodSJ = 0;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.pt() > 20 && sj.absEta() < 2.4) nGoodSJ++;
    		}
    		if (nGoodSJ > 1) goodSJs = true;
    		return goodSJs;
    	};
        // only consider the top two fatjets in pt, provided they are above 200 GeV, then take furthest
    	double minDPhi = 2.0;
    	double minPt = 200;
    	int idx = -1;

        const MomentumF recodilepton = lep1->p4() + lep2->p4();
        std::vector<const FatJet*> fatjets;
        for (const auto& fj : reader_fatjet->jets) {
        	if (fj.pt() > minPt) fatjets.push_back(&fj);
        }
        std::sort(fatjets.begin(),fatjets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());

        if (fatjets.size() == 1) {
        	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[0],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[0],*lep1) > 0.8
        			&& PhysicsUtilities::deltaR(*fatjets[0],*lep2) > 0.8;
        	if (separatedFJ && isBtag(fatjets[0]) && hasGoodSJs(fatjets[0])) idx = fatjets[0]->index();

        } else if (fatjets.size() > 1) {
            double fj_dr = 0;
            for (int k=0; k<2; k++) {
            	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[k],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[k],*lep1) > 0.8
                	    && PhysicsUtilities::deltaR(*fatjets[k],*lep2) > 0.8;

                if (separatedFJ && isBtag(fatjets[k]) && hasGoodSJs(fatjets[k])) {
                    if (PhysicsUtilities::deltaR(*fatjets[k],recodilepton) > fj_dr) {
                	    fj_dr = PhysicsUtilities::deltaR(*fatjets[k],recodilepton);
                	    idx = fatjets[k]->index();
                	}
                }
            }
        }
        if (idx < 0) return 0;
        return &reader_fatjet->jets[idx];
    }
    bool matchDileptons(const Lepton* matchReco1, const Lepton* matchReco2, const Lepton* recoCand1, const Lepton* recoCand2) {
    	// check if the two selected reco leptons are matches to the DR-matched leptons (degeneracy = 2)
    	if (matchReco1 && matchReco2 && recoCand1 && recoCand2) {
			if (matchReco1->isMuon() == recoCand1->isMuon() && matchReco2->isMuon() == recoCand2->isMuon()) {
				if (matchReco1->index() == recoCand1->index() && matchReco2->index() == recoCand2->index()) return true;
			}
			if (matchReco1->isMuon() == recoCand2->isMuon() && matchReco2->isMuon() == recoCand1->isMuon()) {
				if (matchReco1->index() == recoCand2->index() && matchReco2->index() == recoCand1->index()) return true;
			}
    	}
    	return false;
    }
	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}
    void testISO_Bkg(TString sn, std::vector<const Lepton*> bkgleps) {
    	std::vector<const Lepton*> leps;
        for (const auto& lep : bkgleps) {
        	if (passSel(lep,0.1,0.05,4.0,"M","M",9999,9999)) leps.push_back(lep);
        }
        if (leps.size() < 2) return;
        if (leps.size() > 2) {printf("More than 2 reco leps pass Sel minus ISO\n"); return;}
        if (leps.front()->isMuon() ? leps.front()->pt() < 26 : leps.front()->pt() < 30) return;

        const FatJet* hbbjet = findHbbCand(leps.front(),leps[1]);
        if (!hbbjet) return;

    	static const std::vector<double> isoWPs = {0.1,0.2,0.3};
        static const std::vector<TString> isoTypeName = {"miniIso","relIso"};

        TString constID = "ID_eM_muM_";
        TString pref = getDilepChan(leps.front(),leps[1]) + "pt2_10_";
        for (unsigned long type=0; type<isoTypeName.size(); type++) {
            if (type == 0) {
                for (const auto& iso : isoWPs) {
                    TString name = TString::Format("miniIso_%.1f",iso);
                    name.ReplaceAll(".","");
                    bool passIso1 = (leps.front()->miniIso() < iso);
                    bool passIso2 = (leps[1]->miniIso() < iso);
                    if (passIso1 && passIso2) {
                    plotter.getOrMake1DPre(sn+pref+constID+name,"evts",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                    plotSpectra(sn+pref+constID+name,leps.front(),leps[1],hbbjet);
                    }
                    if (passIso1) plotter.getOrMake1DPre(sn+pref+constID+name,"evts_lep1",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                    if (passIso2) plotter.getOrMake1DPre(sn+pref+constID+name,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                }
            } else if (type == 1) {
                for (const auto& iso : isoWPs) {
                    TString name = TString::Format("relIso_%.1f",iso);
                    name.ReplaceAll(".","");
                    bool passIso1 = leps.front()->isMuon() ? (((const Muon*)leps.front())->dbRelISO() < iso) : (((const Electron*)leps.front())->eaRelISO() < iso);
                    bool passIso2 = leps[1]->isMuon() ? (((const Muon*)leps[1])->dbRelISO() < iso) : (((const Electron*)leps[1])->eaRelISO() < iso);
                    if (passIso1 && passIso2) {
                        plotter.getOrMake1DPre(sn+pref+constID+name,"evts",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                        plotSpectra(sn+pref+constID+name,leps.front(),leps[1],hbbjet);
                    }
                    if (passIso1) plotter.getOrMake1DPre(sn+pref+constID+name,"evts_lep1",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                    if (passIso2) plotter.getOrMake1DPre(sn+pref+constID+name,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
                }
            }
        }
    }
    void testISO_Sig(TString sn, const Lepton* recoLep1, const Lepton* recoLep2, const FatJet* hbbjet) {

    	static const std::vector<double> isoWPs = {0.1,0.2,0.3};
        static const std::vector<TString> isoTypeName = {"miniIso","relIso"};
    	// require the RECO leps to pass the current (single-lep channel) ID WPs
    	bool passID1 = recoLep1->isMuon() ? ((const Muon*)recoLep1)->passMed16ID() : ((const Electron*)recoLep1)->passTightID_noISO();
    	bool passID2 = recoLep2->isMuon() ? ((const Muon*)recoLep2)->passMed16ID() : ((const Electron*)recoLep2)->passTightID_noISO();
    	if (!(passID1 && passID2)) return;

    	TString constID = "ID_eT_muM_";
        for (unsigned long type=0; type<isoTypeName.size(); type++) {
            if (type == 0) {
            	for (const auto& iso : isoWPs) {
            		TString name = TString::Format("miniIso_%.1f",iso);
            		name.ReplaceAll(".","");
            		bool passIso1 = (recoLep1->miniIso() < iso);
            		bool passIso2 = (recoLep2->miniIso() < iso);
            		if (passIso1 && passIso2) {
                        plotter.getOrMake1DPre(sn+constID+name,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                        plotSpectra(sn+constID+name,recoLep1,recoLep2,hbbjet);
            		} else {
//            			printf("Max miniIso: %.3f: lep1 miniIso = %.3f, lep2 miniIso = %.3f\n\n",iso,recoLep1->miniIso(),recoLep2->miniIso());
            		}
            		if (passIso1) plotter.getOrMake1DPre(sn+constID+name,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            		if (passIso2) plotter.getOrMake1DPre(sn+constID+name,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            	}
            } else if (type == 1) {
            	for (const auto& iso : isoWPs) {
            		TString name = TString::Format("relIso_%.1f",iso);
            		name.ReplaceAll(".","");
            		bool passIso1 = recoLep1->isMuon() ? (((const Muon*)recoLep1)->dbRelISO() < iso) : (((const Electron*)recoLep1)->eaRelISO() < iso);
            		bool passIso2 = recoLep2->isMuon() ? (((const Muon*)recoLep2)->dbRelISO() < iso) : (((const Electron*)recoLep2)->eaRelISO() < iso);
            		if (passIso1 && passIso2) {
                        plotter.getOrMake1DPre(sn+constID+name,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                        plotSpectra(sn+constID+name,recoLep1,recoLep2,hbbjet);
            		} else {
/*            				printf("Max relIso: %.3f: lep1 relIso = %.3f, lep2 relIso = %.3f\n\n",iso,
            						recoLep1->isMuon()?((const Muon*)recoLep1)->dbRelISO():((const Electron*)recoLep1)->eaRelISO(),
            								recoLep2->isMuon()?((const Muon*)recoLep2)->dbRelISO():((const Electron*)recoLep2)->eaRelISO());
*/            		}
            		if (passIso1) plotter.getOrMake1DPre(sn+constID+name,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            		if (passIso2) plotter.getOrMake1DPre(sn+constID+name,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
            	}
            }
        }
    }
    void testID_Bkg(TString sn, std::vector<const Lepton*> bkgleps) {
        static const std::vector<TString> ids = {"L","M","T","H"};
        std::vector<const Lepton*> leps;
        for (const auto& lep : bkgleps) {
        	if (passSel(lep,0.1,0.05,4.0,"I","I",0.1,0.2)) leps.push_back(lep);
        }
        if (leps.size() < 2) return;
        if (leps.size() > 2) {printf("More than 2 reco leps pass Sel minus ID\n"); return;}
        if (leps.front()->isMuon() ? leps.front()->pt() < 26 : leps.front()->pt() < 30) return;
    	const FatJet* hbbjet = findHbbCand(leps[0],leps[1]);
        if (!hbbjet) return;

        TString constISO = "miniIso_e01_mu02";
        for (unsigned long id=0; id<ids.size(); id++) {
        	bool passId1 = false;
        	bool passId2 = false;;
        	if (id==0) {
        		passId1 = leps.front()->isMuon() ? ((const Muon*)leps.front())->passLooseID() : ((const Electron*)leps.front())->passLooseID_noISO();
        		passId2 = leps[1]->isMuon() ? ((const Muon*)leps[1])->passLooseID() : ((const Electron*)leps[1])->passLooseID_noISO();
        	} else if (id==1) {
        		passId1 = leps.front()->isMuon() ? ((const Muon*)leps.front())->passMed16ID() : ((const Electron*)leps.front())->passMedID_noISO();
        		passId2 = leps[1]->isMuon() ? ((const Muon*)leps[1])->passMed16ID() : ((const Electron*)leps[1])->passMedID_noISO();
        	} else if (id==2) {
        		passId1 = leps.front()->isMuon() ? ((const Muon*)leps.front())->passTightID() : ((const Electron*)leps.front())->passTightID_noISO();
        		passId2 = leps[1]->isMuon() ? ((const Muon*)leps[1])->passTightID() : ((const Electron*)leps[1])->passTightID_noISO();
        	} else if (id==3) {
        		passId1 = leps.front()->isMuon() ? ((const Muon*)leps.front())->passHighPT() : ((const Electron*)leps.front())->passHEEPID_noISO();
        		passId2 = leps[1]->isMuon() ? ((const Muon*)leps[1])->passHighPT() : ((const Electron*)leps[1])->passHEEPID_noISO();
        	} else printf("check your IDs vector\n");

        	TString name = "ID_"+ids[id]+"_";
        	TString pref = getDilepChan(leps[0],leps[1]) + "pt2_10_";
        	if (passId1 && passId2) {
        	    plotter.getOrMake1DPre(sn+pref+name+constISO,"evts",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
        	    plotSpectra(sn+pref+name+constISO,leps.front(),leps[1],hbbjet);
        	}
        	if (passId1) plotter.getOrMake1DPre(sn+pref+name+constISO,"evts_lep1",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
        	if (passId2) plotter.getOrMake1DPre(sn+pref+name+constISO,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
        }
    }
    void testID_Sig(TString sn, const Lepton* recoLep1, const Lepton* recoLep2, const FatJet* hbbjet) {
        static const std::vector<TString> ids = {"L","M","T","H"};
        // require the RECO leps to pass the current (single-lep) ISO WPs
        bool passISO1 = recoLep1->isMuon() ? recoLep1->miniIso() < 0.2 : recoLep1->miniIso() < 0.1;
        bool passISO2 = recoLep2->isMuon() ? recoLep2->miniIso() < 0.2 : recoLep2->miniIso() < 0.1;
        if (!(passISO1 && passISO2)) return;

        TString constISO = "miniIso_e01_mu02";
        for (unsigned long id=0; id<ids.size(); id++) {
        	bool passId1 = false;
        	bool passId2 = false;;
        	if (id==0) {
        		passId1 = recoLep1->isMuon() ? ((const Muon*)recoLep1)->passLooseID() : ((const Electron*)recoLep1)->passLooseID_noISO();
        		passId2 = recoLep2->isMuon() ? ((const Muon*)recoLep2)->passLooseID() : ((const Electron*)recoLep2)->passLooseID_noISO();
        	} else if (id==1) {
        		passId1 = recoLep1->isMuon() ? ((const Muon*)recoLep1)->passMed16ID() : ((const Electron*)recoLep1)->passMedID_noISO();
        		passId2 = recoLep2->isMuon() ? ((const Muon*)recoLep2)->passMed16ID() : ((const Electron*)recoLep2)->passMedID_noISO();
        	} else if (id==2) {
        		passId1 = recoLep1->isMuon() ? ((const Muon*)recoLep1)->passTightID() : ((const Electron*)recoLep1)->passTightID_noISO();
        		passId2 = recoLep2->isMuon() ? ((const Muon*)recoLep2)->passTightID() : ((const Electron*)recoLep2)->passTightID_noISO();
        	} else if (id==3) {
        		passId1 = recoLep1->isMuon() ? ((const Muon*)recoLep1)->passHighPT() : ((const Electron*)recoLep1)->passHEEPID_noISO();
        		passId2 = recoLep2->isMuon() ? ((const Muon*)recoLep2)->passHighPT() : ((const Electron*)recoLep2)->passHEEPID_noISO();
        	} else printf("check your IDs vector\n");

        	TString name = "ID_"+ids[id]+"_";
        	if (passId1 && passId2) {
           	    plotter.getOrMake1DPre(sn+name+constISO,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
           	    plotSpectra(sn+name+constISO,recoLep1,recoLep2,hbbjet);
        	}
        	if (passId1) plotter.getOrMake1DPre(sn+name+constISO,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        	if (passId2) plotter.getOrMake1DPre(sn+name+constISO,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        }
    }
    void testIP_Bkg(TString sn, std::vector<const Lepton*> bkgleps) {
    	std::vector<const Lepton*> leps;
    	for (const auto& lep : bkgleps) {
    		if (passSel(lep,999,999,999,"M","M",0.1,0.2)) leps.push_back(lep);
    	}
    	if (leps.size() < 2) return;
    	TString constIDISO = "ID_eM_muM_miniIS0_e01_mu02";

    	// N-1 plots for each parameter in IP selection
        static const std::vector<double> vec_dz = {0.01,0.05,0.1,0.2,0.3};
        static const std::vector<double> vec_d0 = {0.01,0.05,0.1,0.2,0.3};
        static const std::vector<double> vec_SIP = {2,3,4,5,6,7,8};

        std::vector<const Lepton*> lepCandsDZ;
        std::vector<const Lepton*> lepCandsD0;
        std::vector<const Lepton*> lepCandsSIP;
        for (const auto& lep : leps) {
        	if (lep->d0() < 0.05 && lep->sip3D() < 4.0) lepCandsDZ.push_back(lep);
        	if (lep->dz() < 0.1 && lep->sip3D() < 4.0) lepCandsD0.push_back(lep);
        	if (lep->dz() < 0.1 && lep->d0() < 0.05) lepCandsSIP.push_back(lep);
        }
        if (lepCandsDZ.size() == 2) {
        	if (lepCandsDZ.front()->isMuon() ? lepCandsDZ.front()->pt() > 26 : lepCandsDZ.front()->pt() > 30) {
				const FatJet* hbbjet = findHbbCand(lepCandsDZ.front(),lepCandsDZ[1]);
				if (!hbbjet) return;

				TString pref = getDilepChan(lepCandsDZ.front(),lepCandsDZ[1]) + "pt2_10_";
				for (const auto& dz : vec_dz) {
					bool passDZ1 = lepCandsDZ.front()->dz() < dz;
					bool passDZ2 = lepCandsDZ[1]->dz() < dz;

					TString name = TString::Format("IP_dz_%.2f_",dz);
					name.ReplaceAll(".","");
					if (passDZ1 && passDZ2) {
						plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
						plotSpectra(sn+pref+name+constIDISO,lepCandsDZ.front(),lepCandsDZ[1],hbbjet);
					}
					if (passDZ1) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep1",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
					if (passDZ2) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
				}
			}
        } else if (lepCandsDZ.size() > 2) printf("More than 2 RECO leps pass Sel minus DZ\n");

        if (lepCandsD0.size() == 2) {
        	if (lepCandsD0.front()->isMuon() ? lepCandsD0.front()->pt() > 26 : lepCandsD0.front()->pt() > 30) {
				const FatJet* hbbjet = findHbbCand(lepCandsD0.front(),lepCandsD0[1]);
				if (!hbbjet) return;

				TString pref = getDilepChan(lepCandsD0.front(),lepCandsD0[1]) + "pt2_10_";
				for (const auto& d0 : vec_d0) {
					bool passD01 = lepCandsDZ.front()->d0() < d0;
					bool passD02 = lepCandsDZ[1]->d0() < d0;

					TString name = TString::Format("IP_d0_%.2f_",d0);
					name.ReplaceAll(".","");
					if (passD01 && passD02) {
						plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
						plotSpectra(sn+pref+name+constIDISO,lepCandsD0.front(),lepCandsD0[1],hbbjet);
					}
					if (passD01) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep1",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
					if (passD02) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
				}
        	}
        } else if (lepCandsD0.size() > 2) printf("More than 2 RECO leps pass Sel minus D0\n");

        if (lepCandsSIP.size() == 2) {
        	if (lepCandsSIP.front()->isMuon() ? lepCandsSIP.front()->pt() > 26 : lepCandsSIP.front()->pt() > 30) {
				const FatJet* hbbjet = findHbbCand(lepCandsSIP.front(),lepCandsSIP[1]);
				if (!hbbjet) return;

				TString pref = getDilepChan(lepCandsSIP.front(),lepCandsSIP[1]) + "pt2_10_";
				for (const auto& sip : vec_SIP) {
					bool passSIP1 = lepCandsSIP.front()->sip3D() < sip;
					bool passSIP2 = lepCandsSIP[1]->sip3D() < sip;

					TString name = TString::Format("IP_SIP_%.0f_",sip);
					name.ReplaceAll(".","");
					if (passSIP1 && passSIP2) {
						plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts",";T_{T}",50,600,4600)->Fill(ht_chs,weight);
						plotSpectra(sn+pref+name+constIDISO,lepCandsSIP.front(),lepCandsSIP[1],hbbjet);
					}
					if (passSIP1) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep1",";T_{T}",50,600,4600)->Fill(ht_chs,weight);
					if (passSIP2) plotter.getOrMake1DPre(sn+pref+name+constIDISO,"evts_lep2",";H_{T}",50,600,4600)->Fill(ht_chs,weight);
				}
        	}
        } else if (lepCandsSIP.size() > 2) printf("More than 2 RECO leps pass Sel minus SIP\n");
    }

    void testIP_Sig(TString sn, const Lepton* recoLep1, const Lepton* recoLep2, const FatJet* hbbjet) {
        static const std::vector<double> vec_dz = {0.01,0.05,0.1,0.2,0.3};
        static const std::vector<double> vec_d0 = {0.01,0.05,0.1,0.2,0.3};
        static const std::vector<double> vec_SIP = {2,3,4,5,6,7,8};

        // require the RECO leps to pass the current (single-lep) ID+ISO WPs
        bool passIDISO1 = recoLep1->isMuon() ? (recoLep1->miniIso() < 0.2) && (((const Muon*)recoLep1)->passMed16ID()) :
        			                           (recoLep1->miniIso() < 0.1) && (((const Electron*)recoLep1)->passMedID_noISO());
        bool passIDISO2 = recoLep2->isMuon() ? (recoLep2->miniIso() < 0.2) && (((const Muon*)recoLep2)->passMed16ID()) :
        			                           (recoLep2->miniIso() < 0.1) && (((const Electron*)recoLep2)->passMedID_noISO());
        if (!(passIDISO1 && passIDISO2)) return;

        TString constIDISO = "ID_eM_muM_miniIso_e01_mu02";
        bool passNomDZ = recoLep1->dz() < 0.1 && recoLep2->dz() < 0.1;
        bool passNomD0 = recoLep1->d0() < 0.05 && recoLep2->d0() < 0.05;
        bool passNomSIP = recoLep1->sip3D() < 4.0 && recoLep2->sip3D() < 4.0;

        // N-1 plots for each parameter in the IP selection
        if (passNomDZ && passNomD0) {
        	for (const auto& sip : vec_SIP) {
        		bool passSIP1 = recoLep1->sip3D() < sip;
        		bool passSIP2 = recoLep2->sip3D() < sip;

        		TString name = TString::Format("IP_SIP_%.0f_",sip);
        		name.ReplaceAll(".","");
        		if (passSIP1 && passSIP2) {
                   	plotter.getOrMake1DPre(sn+name+constIDISO,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                   	plotSpectra(sn+name+constIDISO,recoLep1,recoLep2,hbbjet);
        		}
        		if (passSIP1) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		if (passSIP2) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        	}
        }
        if (passNomDZ && passNomSIP) {
        	for (const auto& d0 : vec_d0) {
        		bool passD01 = recoLep1->d0() < d0;
        		bool passD02 = recoLep2->d0() < d0;

       			TString name = TString::Format("IP_d0_%.2f_",d0);
       			name.ReplaceAll(".","");
       			if (passD01 && passD02) {
                   	plotter.getOrMake1DPre(sn+name+constIDISO,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                   	plotSpectra(sn+name+constIDISO,recoLep1,recoLep2,hbbjet);
       			}
        		if (passD01) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		if (passD02) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
       		}
       	}
       	if (passNomD0 && passNomSIP) {
       		for (const auto& dz : vec_dz) {
       			bool passDZ1 = recoLep1->dz() < dz;
       			bool passDZ2 = recoLep2->dz() < dz;

       			TString name = TString::Format("IP_dz_%.2f_",dz);
       			name.ReplaceAll(".","");
       			if (passDZ1 && passDZ2) {
                   	plotter.getOrMake1DPre(sn+name+constIDISO,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                   	plotSpectra(sn+name+constIDISO,recoLep1,recoLep2,hbbjet);
       			}
       			if (passDZ1) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
       			if (passDZ2) plotter.getOrMake1DPre(sn+name+constIDISO,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
       		}
       	}
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
    	bool pass1 = (abs(lep1->d0()) < 0.05) && (abs(lep1->dz()) < 0.1) && (lep1->sip3D() < 4);
    	bool pass2 = (abs(lep2->d0()) < 0.05) && (abs(lep2->dz()) < 0.1) && (lep2->sip3D() < 4);

    	if (pass1 && pass2) return true;
    	else {
//    		printf("lep1: d0 = %.3f, dz = %.3f, SIP = %.3f\n",lep1->d0(),lep1->dz(),lep1->sip3D());
//    		printf("lep2: d0 = %.3f, dz = %.3f, SIP = %.3f\n\n",lep2->d0(),lep2->dz(),lep2->sip3D());
    		return false;
    	}
    }
    bool passSel(const Lepton* lep, double dz, double d0, double sip, TString eID, TString muID, double eISO, double muISO) {
    	bool passIP = lep->dz() < dz && lep->d0() < d0 && lep->sip3D() < sip;
    	if (!passIP) return false;

    	bool passID = false;
    	bool passISO = false;
    	if (lep->isMuon()) {
    		passISO = lep->miniIso() < muISO;
    		if (muID == "I") passID = true;
    		else if (muID == "L") passID = ((const Muon*)lep)->passLooseID();
    		else if (muID == "M") passID = ((const Muon*)lep)->passMed16ID();
    		else if (muID == "T") passID = ((const Muon*)lep)->passTightID();
    		else if (muID == "H") passID = ((const Muon*)lep)->passHighPT();
    	} else {
    		passISO = lep->miniIso() < eISO;
    		if (eID == "I") passID = true;
    		else if (eID == "L") passID = ((const Electron*)lep)->passLooseID_noISO();
    		else if (eID == "M") passID = ((const Electron*)lep)->passMedID_noISO();
    		else if (eID == "T") passID = ((const Electron*)lep)->passTightID_noISO();
    		else if (eID == "H") passID = ((const Electron*)lep)->passHEEPID_noISO();
    	}
    	if (passID && passISO) return true;
    	else return false;
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

			// RECO lepton pt cuts
			bool passRecoPt1 = matchLep1->isMuon() ? (matchLep1->pt() > 26) : (matchLep1->pt() > 30);
			bool passRecoPt2 = (matchLep2->pt() > 10);
			if (!(passRecoPt1 && passRecoPt2)) {
//				printDebugInfo(sn,lep1,lep2,matchLep1->index(),matchLep2->index());
				return false;
			}
			plotter.getOrMake1DPre(sn+"passRecoPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			const FatJet* hbbjet = findHbbCand(matchLep1,matchLep2);
			if (!hbbjet) return false;
			plotter.getOrMake1DPre(sn+"validHbb","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			testIP_Sig(sn,matchLep1,matchLep2,hbbjet);
			if (!passIPcuts(matchLep1,matchLep2)) return false;
			plotter.getOrMake1DPre(sn+"passIP","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			testISO_Sig(sn,matchLep1,matchLep2,hbbjet);
        	testID_Sig(sn,matchLep1,matchLep2,hbbjet);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
        	plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(ht_chs,weight);
            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
//printf("debug0\n");
        	// collect the muons and electrons together and then sort the leptons by pt
        	std::vector<const Lepton*> leps;
        	for (const auto* mu : muons) leps.push_back(mu);
        	for (const auto* e : electrons) leps.push_back(e);
            std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());
//printf("debug1\n");
            testIP_Bkg(sn,leps);
//            printf("debug2\n");
            testID_Bkg(sn,leps);
//            printf("debug3\n");
            testISO_Bkg(sn,leps);
//            printf("debug4\n");
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
    }
    HistGetter plotter;
};

#endif

void DileptonSelection(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void DileptonSelection(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
