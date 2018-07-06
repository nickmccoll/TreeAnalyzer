#include "TFile.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooAbsPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooAbsRealLValue.h"
#include "RooRealVar.h"
#include "TH1D.h"
#include "TH2.h"
#include <vector>
#include <string>

struct Contribution {
    Contribution(const std::string& name, const bool isSignal, const std::string postFix)
    : name(name), isSignal(isSignal), postFix(postFix)
    {}
    std::string name = "";
    bool isSignal = false;
    std::string postFix = "";
};

class PostFitter {
public:
    PostFitter(const std::string& inputFileName);
    ~PostFitter();

    void addSignal(const std::string& name, const std::string postFix="") {contributions.emplace_back(name,true,postFix);}
    void addBkg(const std::string& name,const std::string postFix=""){contributions.emplace_back(name,false,postFix);}
    void doFit(const std::string& model = "s",const bool minos=false,const bool weighted = false, const bool verbose = false);
    TH2 * get2DHistogram(const std::string& var1,const std::string& var2,const std::string& cat, const std::string& contribution);
    const Contribution* getContribution(const std::string& name);
    void fix(const std::string& parameter, const double val);
    const std::vector<Contribution>& getContributions() const {return contributions;}
private:
    TFile * rootFile=0;
    RooWorkspace* w=0;
    RooFitResult* fitResult = 0;
    std::vector<Contribution> contributions;
};


PostFitter::PostFitter(const std::string& inputFileName){
    rootFile = new TFile(inputFileName.c_str(),"READ");
    rootFile->GetObject("w",w);
}
PostFitter::~PostFitter(){
    rootFile->Close();
    delete rootFile;
}
void PostFitter::doFit(const std::string& model,const bool minos,const bool weighted, const bool verbose){
    fitResult = w->pdf((std::string("model_")+model).c_str())->fitTo(*w->data("data_obs"),
            RooFit::NumCPU(8),RooFit::SumW2Error(weighted),RooFit::Minos(minos),RooFit::Verbose(verbose),RooFit::Save(true) );
}
const Contribution* PostFitter::getContribution(const std::string& name){
    for(const auto& c: contributions){
        if(c.name == name) return &c;
    }
    return 0;
}

TH2 * PostFitter::get2DHistogram(const std::string& var1,const std::string& var2,const std::string& cat, const std::string& contribution){
    auto * c = getContribution(contribution);
    if(!c) return 0;

    const std::string prefix =  c->isSignal ? "shapeSig" : "shapeBkg";
    const std::string pdfName = prefix + "_" + c->name + "_"+  cat +c->postFix ;
    auto histogram = w->pdf( pdfName.c_str())->createHistogram((c->name+"_"+cat).c_str(),*w->var(var1.c_str()),RooFit::YVar(*w->var(var2.c_str())),RooFit::IntrinsicBinning());
    auto frame = w->var(var1.c_str())->frame();
    const std::string cutStr = std::string("CMS_channel==CMS_channel::") + cat;
    auto reducedDataset=w->data("data_obs")->reduce(cutStr.c_str());
    reducedDataset->plotOn(frame,RooFit::Name("datapoints"),RooFit::Invisible());
    ((RooSimultaneous*)w->pdf("model_s"))->getPdf(cat.c_str())->plotOn(frame,RooFit::Components(pdfName.c_str()),RooFit::Name("tmp"),RooFit::Invisible(),RooFit::Normalization(1.0,RooAbsReal::RelativeExpected));
    auto curve=frame->getCurve("tmp");
    auto binArray = w->var(var1.c_str())->getBinning().array();
    auto nBins = w->var(var1.c_str())->getBinning().numBins();
    auto histo = new TH1D("tmp","histo",nBins,binArray);
    histo->SetDirectory(0);
    for(int iJ = 1; iJ <= histo->GetNbinsX(); ++iJ ){
        auto x=histo->GetXaxis()->GetBinCenter(iJ);
        histo->SetBinContent(iJ,curve->Eval(x));
    }
    histogram->Scale(histo->Integral()/histogram->Integral());
    delete histo;
    return (TH2*)histogram;
}

void PostFitter::fix(const std::string& parameter, const double val){
    w->var(parameter.c_str())->setVal(val);
    w->var(parameter.c_str())->setConstant(1);
}


