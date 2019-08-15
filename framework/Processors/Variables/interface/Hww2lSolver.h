#ifndef PROCESSORS_VARIABLES_HWW2LSOLVER_H
#define PROCESSORS_VARIABLES_HWW2LSOLVER_H

#include "DataFormats/interface/Momentum.h"
#include "TFitter.h"

namespace TAna {

class HwwInfo {
	// this class holds reconstructed objects from minimization and value of the test statistic
public:
	HwwInfo() {};
	float testStat = -1;
	ASTypes::CylLorentzVectorF nu1;
	ASTypes::CylLorentzVectorF nu2;
	ASTypes::CylLorentzVectorF w_on;
	ASTypes::CylLorentzVectorF w_star;
	ASTypes::CylLorentzVectorF hWW;
};

class Hww2lSolver {
public:
	Hww2lSolver ();
	~Hww2lSolver();
	TFitter *minimizer;

	// MINUIT function wrapper, std form for minimizations
	static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);

	// call function to do minimization which fills HwwInfo and outputs test statistic val
	double HwwMinimization(const ASTypes::CylLorentzVectorF& lep1, const ASTypes::CylLorentzVectorF& lep2, const ASTypes::CylLorentzVectorF& met, HwwInfo *info);

	// define the loss function to be minimized
	static double HwwFunction(const double lep1X, const double lep1Y, const double lep1Z, const double lep2X, const double lep2Y, const double lep2Z,
						const double nu1X, const double nu1Y, const double nu1Z, const double nu2X, const double nu2Y, const double nu2Z,
						const double metX, const double metY, const double Wmass, const double WstarMass, HwwInfo* info = 0
						);

	static MomentumF getSimpleHiggsMom(const MomentumF llMom, const MomentumF met, const float massinv);

};

}
#endif
