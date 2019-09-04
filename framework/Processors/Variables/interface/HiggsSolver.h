
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


    static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par,int flg);
    double hSolverMinimization(const ASTypes::CylLorentzVectorF& lep,
            const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met,
            bool jetIsVirtual,const HWWParameters& params, HiggsSolverInfo * info);


    static double hSolverFunction( const double leptonX, const double leptonY, const double leptonZ,
            const double neutrinoX, const double neutrinoY,const double neutrinoZ,
            const double jetX,    const double jetY,    const double jetZ,    const double jetM,
            const double jetSF, const double metX, const double metY,
            HiggsSolverInfo * info = 0
    );


    static MomentumF getInvisible(const MomentumF& met, const MomentumF& vis,
            const double hMass = 125);

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

class OneDimPDFWInterp : public BASEPDF {
public:
    virtual ~OneDimPDFWInterp(){}
    virtual void setup(TFile * inFile, const std::string& hLowName, const std::string& hHighName,
            const double lowValue, const double highValue, bool verbose= false );
    void setInterp(double i);
    virtual double getProbability(const double x) const;


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
//Histogram assumptions: Equal bin widths
//bin contents are the probability density (divided by bin widths)
//That way we can interpolate between bin centers and not mess up the integral
//the overflow and underflow contain the INTEGRAL of probability from the last bins bin center to
//infinity. We do an extrapolation assuming an exponential form.
class OneDimPDFWInterpAndExtrap : public OneDimPDFWInterp {
public:
    virtual void setup(TFile * inFile, const std::string& hLowName, const std::string& hHighName,
            const double lowValue, const double highValue, bool verbose= false );
    double getProbability(const double x) const;
    double minX = 0;
    double maxX = 0;
};

//Histogram assumptions: Equal bin widths
//bin contents are the probability density (divided by bin widths)
//That way we can interpolate between bin centers and not mess up the integral
//the overflow and underflow contain the INTEGRAL of probability from the last bins bin center to
//infinity. We do an extrapolation assuming an exponential form.
class OneDimPDFWExtrap : public BASEPDF {
public:
    virtual void setup(TFile * inFile, const std::string& hName, bool verbose= false );
    double getProbability(const double x) const;
    double bW = 0;
    int nB = 0;
    double minX = 0;
    double maxX = 0;
    std::unique_ptr<TH1> hCur  = 0;
};

class TwoDimPDF : public BASEPDF {
public:
    void setup(TFile * inFile, const std::string& hName,  bool verbose= false );
    double getProbability(const double x,const double y) const;
    double bW = 0;
    int nBX = 0;
    double minX = 0;
    double maxX = 0;
    int nBY = 0;
    double minY = 0;
    double maxY = 0;
    std::unique_ptr<TH1> h    = 0;
};


class HiggsLiFunction {
public:
//    enum PARAMList {NEUT_X,NEUT_Y,NEUT_Z, QQJET_RES};

    enum PARAMList {EMET_PERP,EMET_PAR, NEUT_Z,WQQ_RES};

    HiggsLiFunction() {}

    void setObservables(const MomentumF& inL, const MomentumF& inM,
            const MomentumF& inJ);

    void setIterationStorage(const double * p);

    double operator()(const double * p);
    double test(const double * p);

    //constants for this class
    std::vector<std::shared_ptr<BASEPDF>> pdfs;
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
    ASTypes::CartLorentzVector wlnu;
    ASTypes::CartLorentzVector hww;

    //output
    double LL = 0;
};
class HiggsSolInfoDebug {
public:
    HiggsSolInfoDebug()
{};
    int minOut = -2;
    float likeli = -1;
    float noSDLikli = -1;
    float ptRes = -1;
    ASTypes::CylLorentzVectorF neutrino;
    ASTypes::CylLorentzVectorF wlnu;
    ASTypes::CylLorentzVectorF hWW;
    ASTypes::CylLorentzVectorF wqqjet;

    float emetperp=0;
    float emetpar=0;
};


class HiggsSolverInfoDebug {
public:
    HiggsSolverInfoDebug()
{};
    HiggsSolInfoDebug osqq_sol;
    HiggsSolInfoDebug vqq_sol;
    HiggsSolInfoDebug min_sol;
};
class HiggsLi {
public:
    typedef TMinuitMinimizer Minimizer;

//    typedef ROOT::Minuit2::Minuit2Minimizer Minimizer;

    HiggsLi(const std::string& dataDir);
    ~HiggsLi() {

    }
//    enum PDFList {EMET_PERP,EMET_PAR,WQQ_RES, WQQ_SDMASS, HWW_WLNU_MASS, NPDFS};
    enum PDFList {EMET_PERP,EMET_PAR,WQQ_RES, WQQ_SDMASS, WLNU_MASS, HWW_MASS, NPDFS};



    void setup(std::string fileName, const double ptCorB_,const double ptCorM_,
            bool verbose= false );

//    void resetParameters(ROOT::Minuit2::Minuit2Minimizer& min,
//            const double metX=0, const double metY=0, const double metZ=0);

    void resetParameters(Minimizer& min,
            const double neutZ=0,const double bounds=1000);

    double getCorrHWWPT(const double recoPT) const;

    void minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
            double qqSDMass, HiggsSolverInfoDebug& out);

    const std::string dataDir;
    double ptCorB=0;
    double ptCorM=0;
    std::vector<std::shared_ptr<BASEPDF>> osqq_pdfs;
    std::vector<std::shared_ptr<BASEPDF>> vqq_pdfs;
    std::unique_ptr<TH1> hwwPT;

    HiggsLiFunction osqq_function;
    HiggsLiFunction vqq_function;

    ROOT::Math::Functor osqq_functor;
    ROOT::Math::Functor vqq_functor;


    Minimizer osqq_minimizer;
    Minimizer vqq_minimizer;
};



class BkgLiFunction {
public:

    enum PARAMList {EMET_PERP,EMET_PAR, NEUT_Z};

    BkgLiFunction() {}

    void setObservables(const MomentumF& inL, const MomentumF& inM,
            const MomentumF& inJ);

    void setIterationStorage(const double * p);

    double operator()(const double * p);

    //constants for this class
    std::vector<std::shared_ptr<BASEPDF>> pdfs;

    //constants for every run (set in setObservables)
    ASTypes::CartLorentzVector lepton;
    ASTypes::CartLorentzVector met;
    ASTypes::CartLorentzVector qqJet;

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
    ASTypes::CartLorentzVector neutrino;
    ASTypes::CartLorentzVector wlnu;
    ASTypes::CartLorentzVector hww;

    //output
    double LL = 0;
};


class BkgLi {
public:
    typedef TMinuitMinimizer Minimizer;

    BkgLi(const std::string& dataDir);
    ~BkgLi() {

    }
    enum PDFList {EMET_PERP,EMET_PAR, WLNU_MASS, HWW_MASS, NPDFS};

    void setup(std::string fileName,bool verbose= false );
    void resetParameters(Minimizer& min);
    void minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
            HiggsSolverInfoDebug& out);

    const std::string dataDir;
    double ptCorB=0;
    double ptCorM=0;
    std::vector<std::shared_ptr<BASEPDF>> pdfs;
    std::unique_ptr<TH1> hwwPT;
    BkgLiFunction function;
    ROOT::Math::Functor functor;
    Minimizer minimizer;
};



}







#endif

