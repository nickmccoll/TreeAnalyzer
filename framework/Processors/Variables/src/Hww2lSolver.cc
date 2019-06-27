#include "Processors/Variables/interface/Hww2lSolver.h"
#include "DataFormats/interface/Momentum.h"
#include "TMath.h"
#include "TFitter.h"

using namespace std;
namespace TAna {

Hww2lSolver::Hww2lSolver() : minimizer(new TFitter(6)) {
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
}

Hww2lSolver::~Hww2lSolver(){
  delete minimizer;
}

void Hww2lSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par, int flg) {
	if (nDim) {}
	if (gout) {}
	if (flg) {}

	result = HwwFunction(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],
            		par[8],par[9],par[10],par[11],par[12],par[13],par[14],par[15]);

}

double Hww2lSolver::HwwFunction(const double lep1X, const double lep1Y, const double lep1Z, const double lep2X, const double lep2Y, const double lep2Z,
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
		info->w_on = w1;
		info->w_star = w2;
		info->hWW = hWW;
	}

	return testStat;
}

double Hww2lSolver::HwwMinimization(const ASTypes::CylLorentzVectorF& lep1, const ASTypes::CylLorentzVectorF& lep2, const ASTypes::CylLorentzVectorF& met, HwwInfo * info) {
	minimizer->SetParameter(0,"lep1x"     ,lep1.x(),0,lep1.x()-0.001,lep1.x()+0.001);
	minimizer->SetParameter(1,"lep1y"     ,lep1.y(),0,lep1.y()-0.001,lep1.y()+0.001);
	minimizer->SetParameter(2,"lep1z"     ,lep1.z(),0,lep1.z()-0.001,lep1.z()+0.001);
	minimizer->SetParameter(3,"lep2x"     ,lep2.x(),0,lep2.x()-0.001,lep2.x()+0.001);
	minimizer->SetParameter(4,"lep2y"     ,lep2.y(),0,lep2.y()-0.001,lep2.y()+0.001);
	minimizer->SetParameter(5,"lep2z"     ,lep2.z(),0,lep2.z()-0.001,lep2.z()+0.001);

	minimizer->SetParameter(6,"nu1x"      ,met.x()/2,500,-3000,3000);
	minimizer->SetParameter(7,"nu1y"      ,met.y()/2,500,-3000,3000);
	minimizer->SetParameter(9,"nu2x"      ,met.x()/2,500,-3000,3000);
	minimizer->SetParameter(10,"nu2y"     ,met.y()/2,500,-3000,3000);

//	MomentumF dilepton; dilepton.setP4(lep1+lep2);
//	double pz_init = 0.5*met.Pt() / TMath::Tan(dilepton.theta());

	ASTypes::CylLorentzVectorF v1(met.Pt()*lep1.P()/(lep1.P()+lep2.P()), lep1.Eta(), lep1.Phi(), 0);
	ASTypes::CylLorentzVectorF v2(met.Pt()*lep2.P()/(lep1.P()+lep2.P()), lep2.Eta(), lep2.Phi(), 0);

	minimizer->SetParameter(8,"nu1z"      ,v1.Pz(),1000,-3000,3000);
	minimizer->SetParameter(11,"nu2z"     ,v2.Pz(),1000,-3000,3000);
	minimizer->SetParameter(14,"massON"   ,80,0,80-0.001,80+0.001);
	minimizer->SetParameter(15,"massOFF"  ,40,0,40-0.001,40+0.001);

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
	minimizer->FixParameter(14);
	minimizer->FixParameter(15);

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

MomentumF Hww2lSolver::getSimpleHiggsMom(MomentumF llMom, MomentumF met, float massinv) {

	double pz = met.pt() / TMath::Tan(llMom.theta());
	pz = ( (pz < 0) == (llMom.pz() < 0) ) ? pz : (-1)*pz;
	double E = sqrt( pow(met.px(),2) + pow(met.py(),2) + pz*pz + massinv*massinv);
	ASTypes::CartLorentzVector pnunu(met.px(),met.py(),pz,E);

	return (llMom.p4() + pnunu);

}

}
