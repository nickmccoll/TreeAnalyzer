#ifndef PROCESSORS_VARIABLES_HWW2LSOLVER_H
#define PROCESSORS_VARIABLES_HWW2LSOLVER_H

#include "DataFormats/interface/Momentum.h"
#include "TFitter.h"

using namespace std;
namespace TAna {

class HwwInfo {
public:
	HwwInfo();
	float testStat;
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

	static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);
	double HwwMinimization(const ASTypes::CylLorentzVectorF& lep1, const ASTypes::CylLorentzVectorF& lep2, const ASTypes::CylLorentzVectorF& met, HwwInfo *info);

	static double HwwFunction(const double lep1X, const double lep1Y, const double lep1Z, const double lep2X, const double lep2Y, const double lep2Z,
						const double nu1X, const double nu1Y, const double nu1Z, const double nu2X, const double nu2Y, const double nu2Z,
						const double metX, const double metY, const double Wmass, const double WstarMass, HwwInfo* info = 0
						);
	TFitter *minimizer;

};

}
#endif
