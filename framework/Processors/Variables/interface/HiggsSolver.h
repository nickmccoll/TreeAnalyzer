
#ifndef PROCESSORS_VARIABLES_HIGGSSOLVER_H
#define PROCESSORS_VARIABLES_HIGGSSOLVER_H


#include "DataFormats/interface/Momentum.h"
#include "TFitter.h"


namespace TAna {

class HiggsSolverInfo {
public:
    HiggsSolverInfo()
{};
    float chiSq = -1;
    float SF = -1;
    ASTypes::CylLorentzVectorF neutrino;
    ASTypes::CylLorentzVectorF wlnu;
    ASTypes::CylLorentzVectorF hWW;
    ASTypes::CylLorentzVectorF wqqjet;
};





class HiggsSolver {

public:

    HiggsSolver();
    ~HiggsSolver();


    static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);
    double hSolverMinimization(const ASTypes::CylLorentzVectorF& lep, const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met, bool jetIsVirtual, HiggsSolverInfo * info);


    static double hSolverFunction( const double leptonX, const double leptonY, const double leptonZ, const double neutrinoX, const double neutrinoY,const double neutrinoZ,
            const double jetX,    const double jetY,    const double jetZ,    const double jetM, const double jetSF,
            const double metX,    const double metY, HiggsSolverInfo * info = 0
    );


    static MomentumF getInvisible(const MomentumF& met, const MomentumF& vis, const double hMass = 125);

    TFitter *minimizer;
};



}


#endif

