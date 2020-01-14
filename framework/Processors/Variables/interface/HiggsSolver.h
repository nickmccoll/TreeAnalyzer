
#ifndef PROCESSORS_VARIABLES_HIGGSSOLVER_H
#define PROCESSORS_VARIABLES_HIGGSSOLVER_H

#include "DataFormats/interface/Momentum.h"
#include "TFitter.h"
#include "TH1F.h"
#include "TFile.h"
#include "Math/Minimizer.h"
#include "Fit/Fitter.h"

#include "Minuit2/Minuit2Minimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"

namespace TAna {
struct HWWParameters;

//--------------------------------------------------------------------------------------------------
// Solve a second order equation
//--------------------------------------------------------------------------------------------------
namespace HSolverBasic {
MomentumF getInvisible(const MomentumF& met, const MomentumF& vis,
            const double hMass = 125);
};


//--------------------------------------------------------------------------------------------------
// Higgs solver based on a chi-square fit
//--------------------------------------------------------------------------------------------------
struct HSolverChiInfo {
    HSolverChiInfo(){};
    float chiSq = -1;
    float SF = -1;
    ASTypes::CylLorentzVectorF neutrino;
    ASTypes::CylLorentzVectorF wlnu;
    ASTypes::CylLorentzVectorF hWW;
    ASTypes::CylLorentzVectorF wqqjet;
};
//--------------------------------------------------------------------------------------------------
class HSolverChi {
public:
    HSolverChi();
    ~HSolverChi();

    static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par,int flg);
    double hSolverMinimization(const ASTypes::CylLorentzVectorF& lep,
            const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met,
            bool jetIsVirtual,const HWWParameters& params, HSolverChiInfo * info);


    static double hSolverFunction( const double leptonX, const double leptonY, const double leptonZ,
            const double neutrinoX, const double neutrinoY,const double neutrinoZ,
            const double jetX,    const double jetY,    const double jetZ,    const double jetM,
            const double jetSF, const double metX, const double metY,
            HSolverChiInfo * info = 0
    );

    TFitter *minimizer;

    //conststants for the hSolver
    static double posMETParErr     ;
    static double negMETParErr     ;
    static double metPerpErr       ;
    static double jetErr           ;
    static double onWlnuMeanJet    ;
    static double offWlnuMeanJet   ;
    static double onWlnuMeanWlnu   ;
    static double offWlnuMeanWlnu  ;
    static double offWlnuPosWlnuErr;
    static double offWnluNegWlnuErr;
    static double onWlnuWlnuErr    ;
    static double onWlnuHWWErr     ;
    static double offWlnuHWWErr    ;
};
//--------------------------------------------------------------------------------------------------
// PDF support for the likelihood HSolvers
//--------------------------------------------------------------------------------------------------
class BASEPDF {
public:
    virtual ~BASEPDF(){}
    virtual double getProbability(const double x,const double y) const {return x+y;}
    virtual double getProbability(const double x)  const {return x;}


    //Helper functions
    static double interpolate(const double x,const double x1,const double x2,
            const double y1,const double y2);
    static double extapDown(const double x, const double bound, const double probAtBound,
            const double integBound);
    static double extapUp(const double x, const double bound, const double probAtBound,
            const double integBound);
};
//--------------------------------------------------------------------------------------------------
// 1D PDF that supports interpolation between two PDFs
// We assume that the underflow and overflow make sense.....
// Bin contents are the probability density
//--------------------------------------------------------------------------------------------------
class OneDimPDFWInterp : public BASEPDF {
public:
    virtual ~OneDimPDFWInterp(){}
    virtual void setup(TFile * inFile, const std::string& hLowName, const std::string& hHighName,
            const double lowValue, const double highValue, bool verbose= false );
    void setInterp(double i);
    using BASEPDF::getProbability;
    virtual double getProbability(const double x) const override;


    double iLow = 0;
    double iHigh = 0;
    double bW = 0;
    int nB = 0;

    std::unique_ptr<TH1> hLow  = 0;
    std::unique_ptr<TH1> hHigh = 0;
    std::unique_ptr<TH1> hCur  = 0;
    std::unique_ptr<TH1> hM    = 0;
    std::unique_ptr<TH1> hB    = 0;
};

//--------------------------------------------------------------------------------------------------
// 1D PDF that supports extrapolation above and below the histogram edges
// Histogram assumptions: Equal bin widths
// bin contents are the probability density (divided by bin widths)
// That way we can interpolate between bin centers and not mess up the integral
// the overflow and underflow contain the INTEGRAL of probability from the last bins bin center to
// infinity. We do an extrapolation assuming an exponential form.
//--------------------------------------------------------------------------------------------------
class OneDimPDFWExtrap : public BASEPDF {
public:
    virtual void setup(TFile * inFile, const std::string& hName, bool verbose= false );
    using BASEPDF::getProbability;
    double getProbability(const double x) const override;
    double bW = 0;
    int nB = 0;
    double minX = 0;
    double maxX = 0;
    std::unique_ptr<TH1> hCur  = 0;
};

//--------------------------------------------------------------------------------------------------
// 1D PDF that supports extrapolation and interpolation (follows Extrap assumptions and req.)
//--------------------------------------------------------------------------------------------------
class OneDimPDFWInterpAndExtrap : public OneDimPDFWInterp {
public:
    virtual void setup(TFile * inFile, const std::string& hLowName, const std::string& hHighName,
            const double lowValue, const double highValue, bool verbose= false ) override;
    using BASEPDF::getProbability;
    double getProbability(const double x) const override;
    double minX = 0;
    double maxX = 0;
};

//--------------------------------------------------------------------------------------------------
// 2D PDF that supports extrapolation...same requirements as OneDimPDFWExtrap
//--------------------------------------------------------------------------------------------------
class TwoDimPDF : public BASEPDF {
public:
    void setup(TFile * inFile, const std::string& hName,  bool verbose= false );
    using BASEPDF::getProbability;
    double getProbability(const double x,const double y) const override;
    double bW = 0;
    int nBX = 0;
    double minX = 0;
    double maxX = 0;
    int nBY = 0;
    double minY = 0;
    double maxY = 0;
    std::unique_ptr<TH1> h    = 0;
};
//--------------------------------------------------------------------------------------------------
// Solve with a likelihood based method
//--------------------------------------------------------------------------------------------------
struct HSolverLiInfo {
    HSolverLiInfo() {};

    //minimizer output
    int minOut = -2;
    float likeli = -1;
    float rawLikeli = -1; //Not normalized
    float noSDLikli = -1; // No SD mass (raw)

    //fitted values
    std::vector<float> params = {0,0,0,0,0};

    //For convenience
    ASTypes::CylLorentzVectorF neutrino;
    ASTypes::CylLorentzVectorF wlnu;
    ASTypes::CylLorentzVectorF hWW;
    ASTypes::CylLorentzVectorF wqqjet;
};
//--------------------------------------------------------------------------------------------------
// H->WW assumption function
//--------------------------------------------------------------------------------------------------
class HSolverFunction {
public:
    virtual ~HSolverFunction(){};
    virtual double operator()(const double * p) =0;

    //constants for this class
    std::vector<std::shared_ptr<BASEPDF>> pdfs;

    //output
    double LL = 0;
};

class HSolverLiFunction : public HSolverFunction {
public:
    enum PARAMList {EMET_PERP,EMET_PAR, NEUT_Z,WQQ_RES,DUMMY, NUMPARAMS};

    HSolverLiFunction() {}
    virtual ~HSolverLiFunction(){};

    void setObservables(const MomentumF& inL, const MomentumF& inM,
            const MomentumF& inJ);

    void setIterationStorage(const double * p);

    double operator()(const double * p);

    //constants for this class
    double jetM =0;

    //constants for every run (set in setObservables)
    ASTypes::CartLorentzVector lepton;
    ASTypes::CartLorentzVector met;
    ASTypes::CartLorentzVector qqJet;
    double origQQPT = 0;

    double hwwParX =0;
    double hwwParY =0;
    double hwwParNormX =0;
    double hwwParNormY =0;
    double hwwPerpNormX =0;
    double hwwPerpNormY =0;
    double hwwMag  =0;
    double metPerp =0;
    double metPar  =0;

    //per iteration storage (set in set iteration vars)
    double neutE = 0;
    double neutPerp = 0;
    double neutPar = 0;
    double neutX = 0;
    double neutY = 0;
    double extraMetPerp = 0;
    double extraMetPar = 0;
    double wqqSF = 0;
    double jetE = 0;
    ASTypes::CartLorentzVector neutrino;
    ASTypes::CartLorentzVector scaledQQJet;
    ASTypes::CartLorentzVector scaledLep;
    ASTypes::CartLorentzVector wlnu;
    ASTypes::CartLorentzVector hww;

};
//--------------------------------------------------------------------------------------------------
// H->WW assumption function normalization
//--------------------------------------------------------------------------------------------------
class HSolverLiAltFunction : public HSolverFunction {
public:
public:
    enum PARAMList {EMET_PERP,EMET_PAR,WQQ_RES,WLNU_MASS,HWW_MASS,NUMPARAMS};

    HSolverLiAltFunction() {}
    virtual ~HSolverLiAltFunction(){};
    double operator()(const double * p);

};

//--------------------------------------------------------------------------------------------------
// H->WW assumption solver
//--------------------------------------------------------------------------------------------------
class HSolverLi {
public:
    typedef TMinuitMinimizer Minimizer;
        enum PDFList {EMET_PERP,EMET_PAR,WQQ_RES, WQQ_SDMASS, WLNU_MASS, HWW_MASS, NPDFS};

    //Setup
    HSolverLi(const std::string& dataDir);
    ~HSolverLi() {}
    void setParamters(const HWWParameters& hwwParam, bool verbose=false);

    // Run
    double minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
            double qqSDMass, HSolverLiInfo& out,
            HSolverLiInfo * osqq_sol = 0, HSolverLiInfo * vqq_sol = 0,
            HSolverLiInfo * nom_sol=0);



    //Support
    void interpolatePDFs(const MomentumF& lepton, const MomentumF& met,const MomentumF& qqJet);
    double fit(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
            double qqSDMass, HSolverLiInfo& out,
            HSolverLiInfo * osqq_sol = 0, HSolverLiInfo * vqq_sol = 0);
    double fitNormalization(double qqSDMass,HSolverLiInfo * out=0);
    void resetParameters(Minimizer& min, const double neutZ=0);
    void resetAltParameters(Minimizer& min);
    double getCorrHWWPT(const double recoPT) const;

    const std::string dataDir;
    double ptCorB=0;
    double ptCorM=0;
    std::vector<std::shared_ptr<BASEPDF>> osqq_pdfs;
    std::vector<std::shared_ptr<BASEPDF>> vqq_pdfs;
    std::unique_ptr<TH1> hwwPT;

    HSolverLiFunction osqq_function;
    HSolverLiFunction vqq_function;
    HSolverLiAltFunction osqq_alt_function;
    HSolverLiAltFunction vqq_alt_function;

    ROOT::Math::Functor osqq_functor;
    ROOT::Math::Functor vqq_functor;
    ROOT::Math::Functor osqq_alt_functor;
    ROOT::Math::Functor vqq_alt_functor;

    Minimizer minimizer;
};


}







#endif

