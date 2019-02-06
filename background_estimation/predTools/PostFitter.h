#include "TFile.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooAbsPdf.h"
#include <RooStats/ModelConfig.h>
#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/ToyMCSamplerOpt.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooAbsRealLValue.h"
#include "RooRealVar.h"
#include "RooUniform.h"
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

typedef std::pair<std::string,std::string> Category; //plotting name, workspace name
typedef std::pair<std::string,RooRealVar*> Variable; //plotting name, workspace variable

class PostFitter {
public:
    PostFitter(const std::string& inputFileName, double rvalue = 0, double mhvalue = 2000, const std::string& modConfig = "ModelConfig", const std::string& modConfigBonly = "ModelConfig_bonly");
    ~PostFitter();

    void addSignal(const std::string& name, const std::string postFix="") {contributions.emplace_back(name,true,postFix);}
    void addBkg(const std::string& name,const std::string postFix=""){contributions.emplace_back(name,false,postFix);}
    void addCategory(const std::string& name,const std::string& wsname){categories.emplace_back(name,wsname);}
    void addVariable(const std::string& name){variables.emplace_back(name,w->var(name.c_str()));}


    void doDataFit(bool doBonly = false);
    void doToys(int nToys);
    void write(const std::string& outname){plotter.write(outname.c_str());}

private:
    void plot2DPDFComponents(RooAbsPdf* func, const std::string& name);
    void plotTotal2DPDF(RooAbsPdf* func,RooAbsData* data,  const std::string& name);
    void plot2DData(RooAbsData* data, const std::string& name);
    RooFitResult* fit(RooStats::ModelConfig *mc, RooAbsData &data, bool saveFitResult = false  );
    double getNormalization(RooRealVar* var, RooAbsPdf* pdf, RooAbsData * data, const std::string& compString = "");
    RooDataSet * copyDataWOZeros(RooDataSet* indata);

    TFile * rootFile=0;
    RooWorkspace* w=0;
    RooFitResult* fitResult = 0;
    std::vector<Contribution> contributions;
    std::vector<Category> categories;
    std::vector<Variable> variables;


    RooStats::ModelConfig * mc = 0;
    RooStats::ModelConfig * mc_bonly = 0;
    const RooArgSet *POI=0;
    const RooArgSet * observables  = 0;
    const RooArgSet * nuisances = 0;
    bool isExtended = false;
    RooRealVar *MH = 0;
    RooAbsData *dobs =  0;
    RooRealVar *weightVar_ = 0;
    std::auto_ptr<RooAbsPdf> nuisancePdf;

    HistGetter plotter;
};


PostFitter::PostFitter(const std::string& inputFileName, double rvalue, double mhvalue, const std::string& modConfig, const std::string& modConfigBonly){
    rootFile = new TFile(inputFileName.c_str(),"READ");
    rootFile->GetObject("w",w);
    mc       = dynamic_cast<RooStats::ModelConfig *>(w->genobj(modConfig.c_str()));
    mc_bonly = dynamic_cast<RooStats::ModelConfig *>(w->genobj(modConfigBonly.c_str()));


    if(dynamic_cast<RooSimultaneous *>(mc->GetPdf()) && !dynamic_cast<RooSimultaneousOpt *>(mc->GetPdf())  ) {
        RooSimultaneousOpt *optpdf = new RooSimultaneousOpt(static_cast<RooSimultaneous&>(*mc->GetPdf()), TString(mc->GetPdf()->GetName())+"_opt");
        w->import(*optpdf);
        mc->SetPdf(*optpdf);
    }

    observables = mc->GetObservables();
    POI = mc->GetParametersOfInterest();
    nuisances = mc->GetNuisanceParameters();

    // Always reset the POIs to floating (post-fit workspaces can actually have them frozen in some cases, in any case they can be re-frozen in the next step
    TIterator *pois = POI->createIterator();
    while (RooRealVar *arg = (RooRealVar*)pois->Next()) {
        arg->setConstant(0);
    }

    //If we want a constant
    //    std::cout << "Will set nuisance parameters to constants: " ;
    //    utils::setAllConstant(*nuisances, true);

    if(mc->GetGlobalObservables()) utils::setAllConstant(*mc->GetGlobalObservables(), true);
    if(mc_bonly && mc_bonly->GetGlobalObservables()) utils::setAllConstant(*mc_bonly->GetGlobalObservables(), true);


    isExtended = mc->GetPdf()->canBeExtended();
    MH = w->var("MH");
    dobs = w->data("data_obs");
    //Remove zeros....this will fix a warning
    //but has no actual effect
    if(dobs->isNonPoissonWeighted()){
        auto newData = copyDataWOZeros((RooDataSet*)dobs);
        dobs =newData;
    }

    

    // Generate with signal model if r or other physics model parameters are defined
    //signal settings...may want to change
    if(rvalue >= 0){
        ((RooRealVar*)POI->find("r"))->setVal(rvalue);
        ((RooRealVar*)POI->find("r"))->setConstant(true);
    }

    MH->setVal(mhvalue);

    // Ok now we're ready to go lets save a "clean snapshot" for the current parameters state
    // w->allVars() misses the RooCategories, useful for some things - so need to include them. Set up a utils function for that
    w->saveSnapshot("prefit", utils::returnAllVars(w));

}


RooFitResult* PostFitter::fit(RooStats::ModelConfig *mc, RooAbsData &data, bool saveFitResult ){
    RooAbsPdf *pdf = mc->GetPdf();
    const RooCmdArg &constrainCmdArg = RooFit::Constrain(*mc->GetNuisanceParameters());
    return pdf->fitTo(data, constrainCmdArg, RooFit::Extended(pdf->canBeExtended()),RooFit::Save(saveFitResult),RooFit::Verbose(0) );
}

void PostFitter::doDataFit(bool doBonly){
    w->loadSnapshot("prefit");
    auto model = doBonly ? mc_bonly : mc;
    plotTotal2DPDF(model->GetPdf(),dobs,"prefit");
    fit(model,*dobs,false);
    plotTotal2DPDF(model->GetPdf(),dobs,"postfit");
    plot2DPDFComponents(model->GetPdf(),"postfit");
    plot2DData(dobs,"data");
}


void PostFitter::plot2DPDFComponents(RooAbsPdf* func, const std::string& name){

    auto components = func->getComponents();
    auto simPdf = dynamic_cast<RooSimultaneous *>(func);
       for(const auto& c :contributions){
           const std::string prefix =  c.isSignal ? "shapeSig" : "shapeBkg";
           if(simPdf){
               for(const auto& cat : categories){
                   const std::string pdfName = prefix + "_" + c.name+"_"+cat.second+c.postFix;
                   const std::string normName = "n_exp_bin" +cat.second +"_proc_"+c.name;
                   auto pdf = (RooAbsPdf*)components->find(pdfName.c_str());
                   auto norm = (RooAbsReal*)components->find(normName.c_str());
                   auto histogram = pdf->createHistogram((name+"_"+c.name+"_"+cat.first).c_str(),*variables[0].second,RooFit::YVar(*variables[1].second),RooFit::IntrinsicBinning());
                   histogram->Scale(norm->getVal()/histogram->Integral());
                   histogram->SetDirectory(0);
                   plotter.add1D(histogram);
               }
           } else {
               const std::string pdfName = prefix + "_" + c.name+c.postFix;
               const std::string normName = "n_exp_bin_proc_"+c.name;
               auto pdf = (RooAbsPdf*)components->find(pdfName.c_str());
               auto norm = (RooAbsReal*)components->find(normName.c_str());
               auto histogram = pdf->createHistogram((name+"_"+c.name).c_str(),*variables[0].second,RooFit::YVar(*variables[1].second),RooFit::IntrinsicBinning());
               histogram->Scale(norm->getVal()/histogram->Integral());
               histogram->SetDirectory(0);
               plotter.add1D(histogram);
           }
       }
       delete components;
}

void PostFitter::plotTotal2DPDF(RooAbsPdf* func,RooAbsData* data,  const std::string& name){
    auto simPdf = dynamic_cast<RooSimultaneous *>(func);
    if(simPdf){
        for(const auto& cat : categories){
            auto pdf = simPdf->getPdf(cat.second.c_str());
            auto obs = pdf->getObservables(*data);
            double norm = pdf->expectedEvents(*obs);
            auto histogram = pdf->createHistogram((name+"_"+cat.first).c_str(),*variables[0].second,RooFit::YVar(*variables[1].second),RooFit::IntrinsicBinning());
            histogram->Scale(norm/histogram->Integral());
            histogram->SetDirectory(0);
            plotter.add1D(histogram);
            delete obs;
        }
    } else {
        auto obs = func->getObservables(*data);
        double norm = func->expectedEvents(*obs);
        auto histogram = func->createHistogram(name.c_str(),*variables[0].second,RooFit::YVar(*variables[1].second),RooFit::IntrinsicBinning());
        histogram->Scale(norm/histogram->Integral());
        histogram->SetDirectory(0);
        plotter.add1D(histogram);
        delete obs;
    }
}

double PostFitter::getNormalization(RooRealVar* var, RooAbsPdf* pdf, RooAbsData * data, const std::string& compString){
    auto compArg = compString.size() ? RooFit::Components(compString.c_str()) : RooCmdArg();
    auto frame = var->frame();
    data->plotOn(frame,RooFit::Name("datapoints"),RooFit::Invisible());
    pdf->plotOn(frame,RooFit::Name("tmp"),compArg,RooFit::Invisible(),RooFit::Normalization(1.0,RooAbsReal::RelativeExpected));
    auto curve=frame->getCurve("tmp");
    auto binArray = var->getBinning().array();
    auto nBins =var->getBinning().numBins();
    auto histo = new TH1D("tmp","histo",nBins,binArray);
    histo->SetDirectory(0);
    for(int iJ = 1; iJ <= histo->GetNbinsX(); ++iJ ){
        auto x=histo->GetXaxis()->GetBinCenter(iJ);
        histo->SetBinContent(iJ,curve->Eval(x));
    }
    double norm = histo->Integral();
    delete histo;
    delete frame;
    return norm;
}


void PostFitter::plot2DData(RooAbsData* data, const std::string& name){
    if(categories.size()){
        for(const auto& cat : categories){
            const std::string cutStr = std::string("CMS_channel==CMS_channel::") + cat.second;
            auto reducedDataset=data->reduce(cutStr.c_str());
            auto histogram = reducedDataset->createHistogram((name+"_"+cat.first).c_str(),*variables[0].second,RooFit::YVar(*variables[1].second));
            histogram->SetDirectory(0);
            plotter.add1D(histogram);
            delete reducedDataset;
        }
    } else {
        auto histogram = data->createHistogram(name.c_str(),*variables[0].second,RooFit::YVar(*variables[1].second));
        histogram->SetDirectory(0);
        plotter.add1D(histogram);
    }
}


void PostFitter::doToys(int nToys){
    w->loadSnapshot("prefit");

    auto* toyModel = mc;
    std::unique_ptr<toymcoptutils::SimPdfGenInfo> newToyMC;

    RooArgSet allFloatingParameters = w->allVars();
    allFloatingParameters.remove(*POI);
    int nFloatingNonPoiParameters = utils::countFloating(allFloatingParameters);

    //if there are floating parameters, fit to the data
    RooFitResult* fitres = 0;
    RooAbsReal* cloneFunc = 0;
    RooAbsPdf* paramPdf = 0;
    RooDataSet* toyParams=0;
    RooArgSet* cloneParams = 0;
    RooArgSet* cloneObservables = 0;
    if (nFloatingNonPoiParameters) {
        utils::setAllConstant(*POI, true);
        fitres = fit(toyModel,*dobs,true);
        utils::setAllConstant(*POI, false);
        w->saveSnapshot("postfit", utils::returnAllVars(w));
        //Now generate toy models about the fit errors
        cloneFunc =  (RooAbsReal*)toyModel->GetPdf()->cloneTree() ;
        cloneParams = cloneFunc->getObservables(fitres->floatParsFinal()) ;
        cloneObservables = cloneFunc->getObservables(dobs);
        paramPdf = fitres->createHessePdf(*cloneParams) ;
        toyParams = paramPdf->generate(*cloneParams,nToys) ;
        newToyMC.reset(new toymcoptutils::SimPdfGenInfo((*(RooAbsPdf*)cloneFunc),*cloneObservables,true) );
    } else {
        newToyMC.reset(new toymcoptutils::SimPdfGenInfo(*toyModel->GetPdf(),*observables,true) );
    }

    for(int iToy = 1; iToy <= nToys; ++iToy) {
        std::cout <<" STARTING TOY: "<< iToy <<std::endl;
        RooAbsData *absdata_toy = 0;
        if(fitres){
            *cloneParams = (*toyParams->get(iToy-1)) ;
            plotTotal2DPDF((RooAbsPdf*)cloneFunc,dobs,std::string("toyModel_")+ASTypes::int2Str(iToy));
            if (isExtended) {
                absdata_toy = newToyMC->generate(weightVar_);
            } else {
                absdata_toy = ((RooAbsPdf*)cloneFunc)->generate(*observables,1);
            }
        } else {
            w->loadSnapshot("prefit");
            if (isExtended) {
                absdata_toy = newToyMC->generate(weightVar_);
            } else {
                absdata_toy = toyModel->GetPdf()->generate(*observables,1);
            }
        }
        auto newData = copyDataWOZeros((RooDataSet*)absdata_toy);
        plot2DData(newData,std::string("toyData_")+ASTypes::int2Str(iToy));
        w->loadSnapshot("prefit");
        fit(toyModel,*newData,false);
        plotTotal2DPDF(toyModel->GetPdf(),newData,std::string("toyDataFit_")+ASTypes::int2Str(iToy));
        delete absdata_toy;
        delete newData;
    }

    if (nFloatingNonPoiParameters) {
        delete fitres;
        delete paramPdf;
        delete cloneFunc;
        delete cloneParams;
        delete cloneObservables;

    }


}

RooDataSet * PostFitter::copyDataWOZeros(RooDataSet* indata){
    RooArgSet vars(*observables), varsPlusWeight(*observables);

    if (indata->isWeighted()) {
        if (weightVar_ == 0) weightVar_ = new RooRealVar("_weight_","",1.0);
        varsPlusWeight.add(*weightVar_);
    }

    auto ret = new RooDataSet(TString::Format("%sDataCopy", indata->GetName()), "", varsPlusWeight, (weightVar_ ? weightVar_->GetName() : 0));
    RooAbsArg::setDirtyInhibit(true); // don't propagate dirty flags while filling histograms
    for (unsigned int i = 0, n = indata->numEntries(); i < n; ++i) {
        vars = *indata->get(i);
        if(indata->weight() ==0) continue;
        ret->add(vars,indata->weight());
    }
    RooAbsArg::setDirtyInhibit(false); // restore proper propagation of dirty flags
    return ret;
}


PostFitter::~PostFitter(){
    rootFile->Close();
    delete rootFile;
}
