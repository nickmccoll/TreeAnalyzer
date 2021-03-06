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
#include "Processors/Variables/interface/JetKinematics.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "Processors/Variables/interface/HiggsSolver.h"
#include "/Users/brentstone/Dropbox/Physics/GitRepos/HistoPlotting/include/Plotter.h"

#include "TSystem.h"
#include "TFitter.h"
using namespace TAna;
using namespace std;

class HwwInfo {
public:
	HwwInfo(){};
	double testStat = -1;
	ASTypes::CylLorentzVector nu1;
	ASTypes::CylLorentzVector nu2;
	ASTypes::CylLorentzVector W;
	ASTypes::CylLorentzVector Wstar;
	ASTypes::CylLorentzVector hww;
};

class HwwSolver {
public:

	TFitter *minimizer;

	HwwSolver();
	~HwwSolver();

	static double HwwFunction(const double lep1X, const double lep1Y, const double lep1Z, const double lep2X, const double lep2Y, const double lep2Z,
					const double nu1X, const double nu1Y, const double nu1Z, const double nu2X, const double nu2Y, const double nu2Z,
					const double metX, const double metY, const double WmassON, const double WmassOFF, HwwInfo* info = 0
					);

	static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);
	double HwwMinimization(const ASTypes::CylLorentzVectorF& lep1, const ASTypes::CylLorentzVectorF& lep2, const ASTypes::CylLorentzVectorF& met, HwwInfo * info, DiHiggsEvent diHiggsEvt);

};  // end of class HwwSolver

HwwSolver::HwwSolver() : minimizer(new TFitter(6)) {
	cout.precision(11);
	double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
}

HwwSolver::~HwwSolver(){
    delete minimizer;
}

void HwwSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par, int flg) {
	if (nDim) {}
	if (gout) {}
	if (flg) {}

	result = HwwFunction(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],
            		par[8],par[9],par[10],par[11],par[12],par[13],par[14],par[15]);

}

double HwwSolver::HwwMinimization(const ASTypes::CylLorentzVectorF& lep1, const ASTypes::CylLorentzVectorF& lep2, const ASTypes::CylLorentzVectorF& met, HwwInfo * info, DiHiggsEvent diHiggsEvt) {
	minimizer->SetParameter(0,"lep1x"     ,lep1.x(),0,lep1.x()-0.001,lep1.x()+0.001);
	minimizer->SetParameter(1,"lep1y"     ,lep1.y(),0,lep1.y()-0.001,lep1.y()+0.001);
	minimizer->SetParameter(2,"lep1z"     ,lep1.z(),0,lep1.z()-0.001,lep1.z()+0.001);
	minimizer->SetParameter(3,"lep2x"     ,lep2.x(),0,lep2.x()-0.001,lep2.x()+0.001);
	minimizer->SetParameter(4,"lep2y"     ,lep2.y(),0,lep2.y()-0.001,lep2.y()+0.001);
	minimizer->SetParameter(5,"lep2z"     ,lep2.z(),0,lep2.z()-0.001,lep2.z()+0.001);

	minimizer->SetParameter(6,"nu1x"      ,met.x()/2,500,-3000,3000);
	minimizer->SetParameter(7,"nu1y"      ,met.y()/2,500,-3000,3000);
//	minimizer->SetParameter(8,"nu1z"      ,0,500,-3000,3000);
	minimizer->SetParameter(9,"nu2x"      ,met.x()/2,500,-3000,3000);
	minimizer->SetParameter(10,"nu2y"     ,met.y()/2,500,-3000,3000);
//	minimizer->SetParameter(11,"nu2z"     ,0,500,-3000,3000);

	double massON, massOFF, nuON_pz(0), nuOFF_pz(0), nuON_px, nuON_py, nuOFF_px, nuOFF_py;
	if (diHiggsEvt.hww) {
		double mass1 = diHiggsEvt.w1 ? diHiggsEvt.w1->mass() : (diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).mass();
		double mass2 = diHiggsEvt.w2 ? diHiggsEvt.w2->mass() : (diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).mass();
		if (mass1 > mass2) {
			massON = mass1;
			massOFF = mass2;
			nuON_pz = diHiggsEvt.w1_d2->pz();
			nuOFF_pz= diHiggsEvt.w2_d2->pz();
			nuON_px = diHiggsEvt.w1_d2->px();
			nuOFF_px= diHiggsEvt.w2_d2->px();
			nuON_py = diHiggsEvt.w1_d2->py();
			nuOFF_py= diHiggsEvt.w2_d2->py();
		} else {
			massON = mass2;
			massOFF = mass1;
			nuON_pz = diHiggsEvt.w2_d2->pz();
			nuOFF_pz= diHiggsEvt.w1_d2->pz();
			nuON_px = diHiggsEvt.w2_d2->px();
			nuOFF_px= diHiggsEvt.w1_d2->px();
			nuON_py = diHiggsEvt.w2_d2->py();
			nuOFF_py= diHiggsEvt.w1_d2->py();
		}
	}
	MomentumF dilepton; dilepton.setP4(lep1+lep2);
	double pz_init = 0.5*met.Pt() / TMath::Tan(dilepton.theta());

	ASTypes::CylLorentzVectorF v1(met.Pt()*lep1.P()/(lep1.P()+lep2.P()), lep1.Eta(), lep1.Phi(), 0);
	ASTypes::CylLorentzVectorF v2(met.Pt()*lep2.P()/(lep1.P()+lep2.P()), lep2.Eta(), lep2.Phi(), 0);

	minimizer->SetParameter(8,"nu1z"      ,v1.Pz(),1000,-3000,3000);
	minimizer->SetParameter(11,"nu2z"     ,v2.Pz(),1000,-3000,3000);
	minimizer->SetParameter(14,"massON"   ,80,0,80-0.001,80+0.001);
	minimizer->SetParameter(15,"massOFF"  ,40,0,40-0.001,40+0.001);
//	minimizer->SetParameter(6,"nu1x"      ,nuON_px,0,nuON_px-100,nuON_px+100);
//	minimizer->SetParameter(7,"nu1y"      ,nuON_py,0,nuON_py-100,nuON_py+100);
//	minimizer->SetParameter(9,"nu2x"      ,nuOFF_px,0,nuOFF_px-100,nuOFF_px+100);
//	minimizer->SetParameter(10,"nu2y"      ,nuOFF_py,0,nuOFF_py-100,nuOFF_py+100);


	minimizer->SetParameter(12,"metx"     ,met.x(),0,met.x()-0.001,met.x()+0.001);
	minimizer->SetParameter(13,"mety"     ,met.y(),0,met.y()-0.001,met.y()+0.001);

	minimizer->FixParameter(0);
	minimizer->FixParameter(1);
	minimizer->FixParameter(2);
	minimizer->FixParameter(3);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(12);
	minimizer->FixParameter(13);

//	minimizer->FixParameter(8);
//	minimizer->FixParameter(11);
	minimizer->FixParameter(14);
	minimizer->FixParameter(15);

//	minimizer->FixParameter(6);
//	minimizer->FixParameter(7);
//	minimizer->FixParameter(9);
//	minimizer->FixParameter(10);

	// Run the simplex minimizer to get close to the minimum [no good precision]
	minimizer->ExecuteCommand("SIMPLEX",0,0);

	// Run the migrad minimizer to precisely estimate the minimum
	//    minimizer->ExecuteCommand("MIGRAD",0,0);

	return  HwwFunction(
	  minimizer->GetParameter(0 ),
	  minimizer->GetParameter(1 ),
	  minimizer->GetParameter(2 ),
	  minimizer->GetParameter(3 ),
	  minimizer->GetParameter(4 ),
	  minimizer->GetParameter(5 ),
	  minimizer->GetParameter(6 ),
	  minimizer->GetParameter(7 ),
	  minimizer->GetParameter(8 ),
	  minimizer->GetParameter(9 ),
	  minimizer->GetParameter(10),
	  minimizer->GetParameter(11),
	  minimizer->GetParameter(12),
	  minimizer->GetParameter(13),
	  minimizer->GetParameter(14),
	  minimizer->GetParameter(15),
	  info);
}


double HwwSolver::HwwFunction(const double lep1X, const double lep1Y, const double lep1Z, const double lep2X, const double lep2Y, const double lep2Z,
		const double nu1X, const double nu1Y, const double nu1Z, const double nu2X, const double nu2Y, const double nu2Z,
		const double metX, const double metY, const double WmassON, const double WmassOFF, HwwInfo* info) {

	const double hwwPar_x = lep1X + lep2X + metX;
	const double hwwPar_y = lep1Y + lep2Y + metY;
	const double hwwMag = std::sqrt(hwwPar_x*hwwPar_x + hwwPar_y*hwwPar_y);

	const double PARA_x = hwwPar_x / hwwMag;
	const double PARA_y = hwwPar_y / hwwMag;

	const double PERP_x = -1*PARA_y;
	const double PERP_y = PARA_x;

	// lambda functions for getting parallel and perp components relative to the rudimentary hww unit vector axis
	auto getPerp =[&](const double momx, const double momy)->double{ return momx*PERP_x+momy*PERP_y;};
	auto getPar =[&](const double momx, const double momy)->double{ return momx*PARA_x+momy*PARA_y;};

	const double metPerp = getPerp(metX,metY);
	const double metPara = getPar(metX,metY);
	const double nusPerp = getPerp(nu1X+nu2X, nu1Y+nu2Y);
	const double nusPara = getPar(nu1X+nu2X, nu1Y+nu2Y);

	const ASTypes::CartLorentzVector lep1(lep1X,lep1Y,lep1Z,sqrt(lep1X*lep1X+lep1Y*lep1Y+lep1Z*lep1Z));
	const ASTypes::CartLorentzVector lep2(lep2X,lep2Y,lep2Z,sqrt(lep2X*lep2X+lep2Y*lep2Y+lep2Z*lep2Z));
	const ASTypes::CartLorentzVector nu1(nu1X,nu1Y,nu1Z,sqrt(nu1X*nu1X+nu1Y*nu1Y+nu1Z*nu1Z));
	const ASTypes::CartLorentzVector nu2(nu2X,nu2Y,nu2Z,sqrt(nu2X*nu2X+nu2Y*nu2Y+nu2Z*nu2Z));

	const ASTypes::CartLorentzVector w1 = lep1 + nu1;
	const ASTypes::CartLorentzVector w2 = lep2 + nu2;
	const ASTypes::CartLorentzVector hWW = w1 + w2;

	// extra met PERP term
	const double extraMetPerp = (metPerp - nusPerp) / hwwMag;
	const double metPerpError = 0.1;
	const double metPerpTestStat = pow( (extraMetPerp/metPerpError), 2);

	// extra met PARA term
	const double extraMetPara = (metPara - nusPara) / hwwMag;
	const double metParError = 0.15;
	const double metParaTestStat = pow( (extraMetPara/metParError), 2);

	// on-shell W term
	const double extraWmass = w1.mass() - WmassON;
	const double W_on_error = 3;
	const double wOnTestStat = pow( (extraWmass/W_on_error), 2);

	// off-shell W term
	const double WstarMass = WmassOFF;
	const double extraWStarMass = w2.mass() - WstarMass;
	const double WStarError = w2.mass() > WstarMass ? 3 : 19;
	const double wStarTestStat = pow( (extraWStarMass/WStarError), 2);

	// Higgs mass term
	const double HmassDiff = hWW.mass() - 125;
	const double HmassError = HmassDiff > 0 ? 4 : 7;
	const double HiggsTestStat = pow( (HmassDiff/HmassError), 2);

	MomentumF w1mom; w1mom.setP4(w1);
	MomentumF w2mom; w2mom.setP4(w2);
	const double dTheta_ww = w1.theta() - w2.theta();
	const double dThetaError = 0.15;
	const double newTerm = pow( (dTheta_ww/dThetaError), 2 );
	MomentumF nunu; nunu.setP4(nu1+nu2);
	MomentumF ll; ll.setP4(lep1+lep2);
//	const double newTerm = pow((PhysicsUtilities::deltaR(nu1,nu2) - PhysicsUtilities::deltaR(lep1,lep2)),2);
//	const double newTerm = 10*PhysicsUtilities::deltaR2(nunu,ll);
//	const double newTerm = pow( (nunu.theta()-ll.theta())/0.5,2);


	const double testStat = metPerpTestStat + metParaTestStat + wOnTestStat + wStarTestStat + HiggsTestStat + newTerm;
//	const double testStat = wOnTestStat + wStarTestStat + HiggsTestStat;
	if (info) {
		info->testStat = testStat;
		info->nu1 = nu1;
		info->nu2 = nu2;
		info->W = w1;
		info->Wstar = w2;
		info->hww = hWW;
	}

	return testStat;
}

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void makePlots(TString sn, const MomentumF hww) {
        plotter.getOrMake1DPre(sn, "hww_mass", ";M_{Hww} (GeV)",100,0,400)->Fill(hww.mass(), weight);
        plotter.getOrMake1DPre(sn, "hh_mass", ";M_{hh}",1000,0,5000)->Fill((hbbCand_2l->p4()+hww.p4()).mass(),weight);
    }

    const GenParticle *getGenLepFromTau(const GenParticle* tau) {
    	if (tau->absPdgId() != ParticleInfo::p_tauminus) return 0;
    	if (!ParticleInfo::isLastInChain(tau)) {
    		cout<<"tau not last in chain"<<endl;
    		return 0;
    	}

    	const GenParticle* gp=0;
    	int nlepsfromtau = 0;
    	for (unsigned int k = 0; k < tau->numberOfDaughters(); k++) {
    		const auto* dau = tau->daughter(k);
    		if (dau->absPdgId() == ParticleInfo::p_muminus || dau->absPdgId() == ParticleInfo::p_eminus) {

    			if (nlepsfromtau >= 1) {
    				cout<<"already a lep in the tau decay!!"<<endl;
    				continue;
    			}
    			gp = dau;
    			nlepsfromtau++;
    		}
    	}
    	if (!gp) return 0;
    	if (nlepsfromtau > 1) {
    		cout<<"Found "<<nlepsfromtau<<" leps in tau decay!!!!"<<endl;
    		return 0;
    	}

    	const GenParticle *ptcl = gp;
    	return ptcl;
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

	void testSolution0(TString sn, const Lepton* lep1, const Lepton* lep2) {
		// assume the MET perfectly encapsulates neutrino momenta and no pz component
		const MomentumF hww = lep1->p4() + lep2->p4() + reader_event->met.p4();
		makePlots(sn+"Sol0_",hww);
	}

	void testSolution1(TString sn, const Lepton* lep1, const Lepton* lep2) {
		if (!isSignal()) return;

		const GenParticle *glep1=0;
		const GenParticle *glep2=0;

		if (diHiggsEvt.w1_d1->absPdgId() == 15) glep1 = getGenLepFromTau(diHiggsEvt.w1_d1);
		else glep1 = diHiggsEvt.w1_d1;

		if (diHiggsEvt.w2_d1->absPdgId() == 15) glep2 = getGenLepFromTau(diHiggsEvt.w2_d1);
		else glep2 = diHiggsEvt.w2_d1;

		if (!glep1 || !glep2) return;

		// use GEN-level invisible momentum
		const MomentumF hww = lep1->p4() + lep2->p4() +
				diHiggsEvt.hww->p4() - glep1->p4() - glep2->p4();
		makePlots(sn+"SolGEN_",hww);
	}

	void testSolution2(TString sn, const Lepton* lep1, const Lepton* lep2) {

		// assume the missing momentum has same theta as the dileptons and retain the pt from MET
		// mvv = 0
		double pz = reader_event->met.pt() / TMath::Tan((lep1->p4()+lep2->p4()).theta());
		pz = (pz < 0 == (lep1->p4()+lep2->p4()).pz() < 0) ? pz : (-1)*pz;
		ASTypes::CartLorentzVector pinv(reader_event->met.px(),reader_event->met.py(),pz,
				sqrt(pow(reader_event->met.px(),2)+pow(reader_event->met.py(),2)+pz*pz));
		const MomentumF inv(pinv);
		const MomentumF hww = lep1->p4() + lep2->p4() + inv.p4();
		makePlots(sn+"Sol2_",hww);
	}

	void testSolution3(TString sn, const Lepton* lep1, const Lepton* lep2) {

		// assume the missing momentum has same theta as the dileptons and retain the pt from MET
		// mvv = 55
		double pz = reader_event->met.pt() / TMath::Tan((lep1->p4()+lep2->p4()).theta());
		pz = (pz < 0 == (lep1->p4()+lep2->p4()).pz() < 0) ? pz : (-1)*pz;
		ASTypes::CartLorentzVector pinv(reader_event->met.px(),reader_event->met.py(),pz,
				sqrt(pow(reader_event->met.px(),2)+pow(reader_event->met.py(),2)+pz*pz+55*55));
		const MomentumF inv(pinv);
		const MomentumF hww = lep1->p4() + lep2->p4() + inv.p4();
		makePlots(sn+"Sol3_",hww);
	}

	void testSolution4(TString sn, const Lepton* lep1, const Lepton* lep2) {

		MomentumF pll = lep1->p4() + lep2->p4();

		// get pz of missing momentum with the delta Theta, then solve for minv by imposing higgs mass
		double pxinv = reader_event->met.px();
		double pyinv = reader_event->met.py();

		double pzinv = reader_event->met.pt() / TMath::Tan(pll.theta());
//		pzinv = (pzinv < 0 == pll.pz() < 0) ? pzinv : (-1)*pzinv;

		double A = 125*125 - pll.mass()*pll.mass() + 2*(pll.px()*pxinv + pll.py()*pyinv + pll.pz()*pzinv);
		double B = 2 * pll.E();
		double C = pxinv*pxinv + pyinv*pyinv + pzinv*pzinv;

		double b = (2*A + B*B);
		double c = A*A - B*B*C;

		double delta = b*b - 4*c;

		double minv = 0;
		if (delta < 0) {
			cout<<"delta is less than 0"<<endl;
			minv = TMath::Sqrt(b/2);
		} else {

			if (sqrt(delta) > b) {
				minv = sqrt( (b + sqrt(delta)) / 2);
				if (minv > 120) {
					cout<<"minv gt 120!"<<endl;
					cout<<"pos m^2 is "<<((b + sqrt(delta)) / 2)<<endl;
					cout<<"neg m^2 is "<<((b - sqrt(delta)) / 2)<<endl;
					printf("ll  mom: px = %f, py = %f, pz = %f, E = %f\n",pll.px(),pll.py(),pll.pz(),pll.E());
					printf("inv mom: px = %f, py = %f, pz = %f, m = %f\n",pxinv,pyinv,pzinv,minv);
					printf("hbb mom: px = %f, py = %f, pz = %f, E = %f\n",hbbCand_2l->px(),hbbCand_2l->py(),
							hbbCand_2l->pz(),hbbCand_2l->E());

					cout<<endl;

					double sol1 = sqrt( (b + sqrt(delta)) / 2);
					double sol2 = sqrt( fabs((b - sqrt(delta)) / 2));
					minv = sol2;
				}
			}
			else {
				cout << "two viable solutions here" << endl;
				double sol1 = sqrt( (b + sqrt(delta)) / 2);
				double sol2 = sqrt( (b - sqrt(delta)) / 2);
				cout<<sol1<<" or "<<sol2<<endl<<endl;;

				minv = sol2;
			}
		}

		cout<<endl<<"MINV = "<<minv<<endl;
		double Einv = sqrt(pxinv*pxinv + pyinv*pyinv + pzinv*pzinv + minv*minv);
		ASTypes::CartLorentzVector pinv(pxinv,pyinv,pzinv,Einv);

		const MomentumF hww = lep1->p4() + lep2->p4() + pinv;
		makePlots(sn+"Sol4_",hww);
	}

	void testMinimization(TString sn, const Lepton* lep1, const Lepton* lep2) {
		HwwInfo info, infoHigh, infoLow;
		double testStatTruth,testStatHighPt,testStatLowPt;
		if (isSignal()) {
			double w1mass = diHiggsEvt.w1 ? diHiggsEvt.w1->mass() : (diHiggsEvt.w1_d1->p4() + diHiggsEvt.w1_d2->p4()).mass();
			double w2mass = diHiggsEvt.w2 ? diHiggsEvt.w2->mass() : (diHiggsEvt.w2_d1->p4() + diHiggsEvt.w2_d2->p4()).mass();

			bool lep1ToW1 = PhysicsUtilities::deltaR2(*lep1,*diHiggsEvt.w1_d1) < PhysicsUtilities::deltaR2(*lep1,*diHiggsEvt.w2_d1) ? true : false;
			double massON, massOFF;
			const Lepton *lepON, *lepOFF;
			if (w1mass > w2mass) {
				massON = w1mass;
				massOFF = w2mass;
				if (lep1ToW1) {
					lepON = lep1;
					lepOFF = lep2;
				} else {
					lepON = lep2;
					lepOFF = lep1;
				}
			} else {
				massON = w2mass;
				massOFF = w1mass;
				if (lep1ToW1) {
					lepON = lep2;
					lepOFF = lep1;
				} else {
					lepON = lep1;
					lepOFF = lep2;
				}
			}
			testStatHighPt = solver.HwwMinimization(lep1->p4(),lep2->p4(),reader_event->met.p4(),&infoHigh,diHiggsEvt);
			testStatLowPt = solver.HwwMinimization(lep2->p4(),lep1->p4(),reader_event->met.p4(),&infoLow,diHiggsEvt);
			testStatTruth = solver.HwwMinimization(lepON->p4(),lepOFF->p4(),reader_event->met.p4(),&info,diHiggsEvt);
		} else {
			testStatHighPt = solver.HwwMinimization(lep1->p4(),lep2->p4(),reader_event->met.p4(),&infoHigh,diHiggsEvt);
			testStatLowPt = solver.HwwMinimization(lep2->p4(),lep1->p4(),reader_event->met.p4(),&infoLow,diHiggsEvt);
		}
//		makePlots(sn+"SolMIN_",info.hww);
//        plotter.getOrMake1DPre(sn, "testStat", ";test Statistic",300,0,30)->Fill(testStatTruth, weight);

		makePlots(sn+"SolMIN_highPt_",infoHigh.hww);
        plotter.getOrMake1DPre(sn, "testStat_highPt", ";test Statistic",300,0,30)->Fill(testStatHighPt, weight);


        if (infoHigh.testStat < infoLow.testStat) {
    		makePlots(sn+"SolMIN_minTestStat_",infoHigh.hww);
            plotter.getOrMake1DPre(sn, "testStat_minTestStat", ";test Statistic",300,0,30)->Fill(testStatHighPt, weight);
        } else {
    		makePlots(sn+"SolMIN_minTestStat_",infoLow.hww);
            plotter.getOrMake1DPre(sn, "testStat_minTestStat", ";test Statistic",300,0,30)->Fill(testStatLowPt, weight);
        }
	}

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(isSignal() && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        TString sn = smpName + "_";

        // current selection
        if(!hbbCand_2l) return false;
        if(selectedDileptons.size() != 2) return false;
        if(hbbCSVCat_2l < BTagging::CSVSJ_MF) return false;
        if (hbbMass_2l < 30 || hbbMass_2l > 210) return false;
        if (ht_puppi < 400) return false;

        testSolution0(sn,dilep1,dilep2);
        testSolution1(sn,dilep1,dilep2);
        testSolution2(sn,dilep1,dilep2);
        testSolution3(sn,dilep1,dilep2);
        testSolution4(sn,dilep1,dilep2);

        if (PhysicsUtilities::deltaR2(*dilep1,*dilep2) > 1.0) return false;
        double mll = (dilep1->p4() + dilep2->p4()).mass();
        if (mll > 75 || mll < 6) return false;
        if (nMedBTags_HbbV != 0) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection

        testSolution0(sn+"fullSel_",dilep1,dilep2);
        testSolution1(sn+"fullSel_",dilep1,dilep2);
        testSolution2(sn+"fullSel_",dilep1,dilep2);
        testSolution3(sn+"fullSel_",dilep1,dilep2);
        testSolution4(sn+"fullSel_",dilep1,dilep2);
//        testMinimization(sn,dilep1,dilep2);

        return true;
    }
    HistGetter plotter;
    void write(TString fileName) {plotter.write(fileName);}
    void outputResults() {
    	TString samp = TString::Format("m%d_",signal_mass);

        Plotter *p = new Plotter();
        p->addHistLine(plotter.getOrMake1DPre(samp+"SolGEN_","hh_mass", ";M_{hh}",100,0,5000),"Perfect");
        p->addHistLine(plotter.getOrMake1DPre(samp+"Sol2_","hh_mass", ";M_{hh}",100,0,5000),"best simple");
        p->addHistLine(plotter.getOrMake1DPre(samp+"SolMIN_minTestStat_","hh_mass", ";M_{hh}",100,0,5000),"Min chi2");
        p->addHistLine(plotter.getOrMake1DPre(samp+"SolMIN_highPt_","hh_mass", ";M_{hh}",100,0,5000),"high Pt");

        p->draw(false,samp);

    }
    HwwSolver solver;
};

#endif

void testHww2lReco(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
    a.outputResults();
}
void testHww2lReco(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
    a.outputResults();
}
