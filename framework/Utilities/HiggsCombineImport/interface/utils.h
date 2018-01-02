#ifndef HiggsAnalysis_CombinedLimit_utils_h
#define HiggsAnalysis_CombinedLimit_utils_h

#include <vector>
#include <string>
#include <unordered_map>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <RooHistError.h>
#include <RooFitResult.h>
#include <TH1.h>
struct RooDataHist;
struct RooAbsData;
struct RooAbsPdf;
struct RooAbsReal;
struct RooAbsArg;
struct RooArgSet;
struct RooArgList;
struct RooSimultaneous;
struct RooAbsCollection;
struct RooWorkspace;
struct RooPlot;
struct RooRealVar;
struct RooProduct;
namespace RooStats { class ModelConfig; }
namespace utils {
    void printRDH(RooAbsData *data) ;

}
#endif
