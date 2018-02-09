
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "Utilities/HiggsCombineImport/interface/VerticalInterpHistPdf.h"
#include "Utilities/HiggsCombineImport/interface/HZZ2L2QRooPdfs.h"
#include <string.h>
#include <regex>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "TCanvas.h"

class FunctionFitter {
public:
    std::string postFix;
    std::unique_ptr<RooWorkspace> w;
    std::vector<std::string> vars;
    std::vector<std::string> params;
    std::vector<std::unique_ptr<TF1>> varDep;
    //for making histograms out of a json
    int    plot_bins = -1;
    double plot_min = -1;
    double plot_max = -1;

    FunctionFitter(const TH1* iH, const std::string postFix = "T", const std::vector<std::string>&  plotVars = {"M"}):
        postFix(postFix)
    {
        w.reset(new RooWorkspace("w",false));
        for(const auto& v: plotVars){
            vars.push_back(v);
            w->factory((v+"[0,10000]").c_str());
        }
        if(iH) setupHist(iH);
    }
    void setVar(std::string var, double val) {
        w->var(var.c_str())->setVal(val);
    }
    void setVar(std::string var, double val, double min, double max) {
        w->var(var.c_str())->setVal(val);
        w->var(var.c_str())->setRange(min,max);
    }
    void setConst(std::string var, double val) {
        w->var(var.c_str())->setConstant(val);
    }
    double getVal(std::string var) {
        return w->var(var.c_str())->getVal();
    }

    void fit(const std::vector<RooCmdArg>& options){
        std::string thisModel = std::string("model")+postFix;
        std::string thisData  = "data";
        if(options.size()==0) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()));
        else if(options.size()==1) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0]);
        else if(options.size()==2) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1]);
        else if(options.size()==3) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1],options[2]);
        else if(options.size()==4) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1],options[2],options[3]);
        else if(options.size()==5) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1],options[2],options[3],options[4]);
    }

    TCanvas* projection(const std::string& name, double& chi2){
        std::string model = std::string("model")+postFix;
        std::string data  = std::string("data");
        auto frame = w->var(vars[0].c_str())->frame();
        w->data(data.c_str())->plotOn(frame);
        w->pdf(model.c_str())->fixAddCoefRange("fit");
        w->pdf(model.c_str())->plotOn(frame, RooFit::NormRange("fit"));
        TCanvas* can = new TCanvas(name.c_str());
        can->cd();
        frame->Draw();
        frame->GetYaxis()->SetTitle(name.c_str());
        frame->GetXaxis()->SetTitle(vars[0].c_str());
        frame->SetTitle("");
        chi2 = frame->chiSquare();
        return can;
    }
    TCanvas* projection2D(const std::string& name, double& chi2){
        std::string model = std::string("model")+postFix;
        auto * pdf = w->pdf(model.c_str())->createHistogram(model.c_str(), *w->var(vars[0].c_str()),RooFit::YVar(*w->var(vars[1].c_str()))   );
        TCanvas* can = new TCanvas(name.c_str());
        can->cd();
        pdf->Draw("COLZ");
        pdf->GetYaxis()->SetTitle(vars[1].c_str());
        pdf->GetXaxis()->SetTitle(vars[0].c_str());
        pdf->SetTitle("");
        return can;
    }
    TH1* data2D(const std::string& name){
        std::string data  = std::string("data");
        auto * hist = w->data(data.c_str())->createHistogram(name.c_str(), *w->var(vars[0].c_str()),RooFit::YVar(*w->var(vars[1].c_str()))   );
        hist->GetYaxis()->SetTitle(vars[1].c_str());
        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }
    TH1* pdf2D(const std::string& name){
        std::string model = std::string("model")+postFix;
        auto * hist = w->pdf(model.c_str())->createHistogram(name.c_str(), *w->var(vars[0].c_str()),RooFit::YVar(*w->var(vars[1].c_str()))   );
        hist->GetYaxis()->SetTitle(vars[1].c_str());
        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }
    TH1* pdf1D(const std::string& name){
        std::string model = std::string("model")+postFix;
        auto * hist = w->pdf(model.c_str())->createHistogram(name.c_str(), *w->var(vars[0].c_str()));
        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }


protected:
    void addParam(std::string name, std::string bounds){
        params.push_back(name);
        w->factory((name + bounds).c_str());
    }

    void setupHist(const TH1* iH,const std::string& name="data"){
        const unsigned int nD = vars.size();
        RooArgList args;
        auto doBinning =[&](const std::string& var, const TAxis* ax) {
            args.add(*w->var(var.c_str()));
            w->var(var.c_str())->setMin(ax->GetXmin());
            w->var(var.c_str())->setMax(ax->GetXmax());
            w->var(var.c_str())->setBins(ax->GetNbins());
            w->var(var.c_str())->setVal((ax->GetXmin()+ax->GetXmax())/2.);
        };
        doBinning(vars[0],iH->GetXaxis());
        if(nD > 1)doBinning(vars[1],iH->GetYaxis());
        if(nD > 2)doBinning(vars[2],iH->GetZaxis());
        RooDataHist dataHist(name.c_str(),name.c_str(),args,iH);
        w->import(dataHist);
    }


    void addCB(const std::string& postFix, const std::string& modelName, const std::string& varName, bool doExpo){
        auto pn =[&] (const std::string& v) ->std::string {return v+postFix;};
        auto pv =[&] (const std::string& v) ->RooRealVar* {return w->var((v+postFix).c_str());};

        std::string model = std::string("model")+postFix;
        addParam(pn("mean"  ),"[90,0,300]");
        addParam(pn("sigma" ),"[8,0,100]");
        addParam(pn("alpha" ),"[1,0.001,100]");
        addParam(pn("alpha2"),"[1,0.001,100]");

        addParam(pn("n")  ,"[5,1,100]");
        addParam(pn("n2"),"[5,1,100]");

        std::string cbName = doExpo ? modelName +"P" : modelName;
        RooDoubleCB modelP(cbName .c_str(),cbName.c_str(),*w->var(varName.c_str())  ,
                *pv("mean") ,*pv("sigma"),*pv("alpha"),*pv("n") ,*pv("alpha2"),*pv("n2"));
        w->import(modelP);

        if(doExpo){
            addParam(pn("slope"),"[-1,-10,0]");
            std::string expName =  modelName +"E";
            RooExponential modelE(expName.c_str(),expName.c_str(),*w->var(varName.c_str()),*pv("slope"));
            w->import(modelE);
            addParam(pn("fE"),"[0.1,0,1]");
            RooAddPdf modelC(modelName.c_str(),modelName.c_str(),*w->pdf(expName.c_str()),*w->pdf(cbName.c_str()),*pv("fE"));
            w->import(modelC);
        }
    }

};
class CBFunctionFitter : public FunctionFitter{
public:
    CBFunctionFitter(const TH1* iH, bool doExpo, const std::string postFix, const std::vector<std::string>&  plotVars = {"M"}) :
    FunctionFitter(iH,postFix,plotVars){
        addCB(postFix,std::string("model")+postFix,vars[0],doExpo);
    }
};


class CBFunctionFitter2D : public FunctionFitter{
public:
    CBFunctionFitter2D(const TH2* iH,bool doExpoX, const std::string postFix, const std::vector<std::string>&  plotVars = {"MJJ","MHH"}) :
    FunctionFitter(iH,postFix,plotVars){
        std::string model = std::string("model")+postFix;
        addCB(postFix+plotVars[0],std::string("model")+postFix+plotVars[0],vars[0],doExpoX);

        auto pn0 =[&] (std::string v) ->std::string{return v+postFix+plotVars[0];};
        auto pn1 =[&] (std::string v) ->std::string{return v+postFix+plotVars[1];};
        auto pn0V =[&] (std::string v) ->RooRealVar*{return w->var((v+postFix+plotVars[0]).c_str());};
        auto pn1V =[&] (std::string v) ->RooRealVar*{return w->var((v+postFix+plotVars[1]).c_str());};
        auto pn1F =[&] (std::string v) ->RooAbsReal*{return w->function((v+postFix+plotVars[1]).c_str());};
        auto pn1P =[&] (std::string v) ->RooAbsPdf*{return w->pdf((v+postFix+plotVars[1]).c_str());};
        auto pn0P =[&] (std::string v) ->RooAbsPdf*{return w->pdf((v+postFix+plotVars[0]).c_str());};

        addParam(pn1("mean"  ),"[90,0,5000]");
        addParam(pn1("sigma" ),"[8,0,1000]");
        addParam(pn1("alpha" ),"[1,0.001,100]");
        addParam(pn1("alpha2"),"[1,0.001,100]");
        addParam(pn1("n")        ,"[5,1,100]");
        addParam(pn1("n2")       ,"[5,1,100]");
        addParam(pn1("maxS")     ,"[2.5,0,5]");
        addParam(pn1("mean_p1"  ),"[0,-5000,5000]");
        addParam(pn1("sigma_p1" ),"[0,-1000,1000]");

        RooFormulaVar xSig  (pn1("xSig"  ).c_str(),"(@0-@1)/@2",RooArgList(*w->var((plotVars[0]).c_str()),*pn0V("mean"),*pn0V("sigma")));
        w->import(xSig  );
        RooFormulaVar xSigC  (pn1("xSigC"  ).c_str(),"max(-1*@0,min(@1,@0))",RooArgList(*pn1V("maxS"  ),*pn1F("xSig"  )));
        w->import(xSigC  );
        RooFormulaVar mean  (pn1("meanF"  ).c_str(),"@0*(1+@1*@2)"     ,RooArgList(*pn1V("mean"  ),*pn1V("mean_p1"  ),*pn1F("xSigC"  )));
        RooFormulaVar sigma (pn1("sigmaF" ).c_str(),"@0*(1+@1*abs(@2))",RooArgList(*pn1V("sigma" ),*pn1V("sigma_p1" ),*pn1F("xSigC"  )));
        w->import(mean  );
        w->import(sigma );
        std::cout <<"IMPORT!2"<<std::endl;
        RooDoubleCB modelMVV(pn1("model").c_str(),pn1("model").c_str(),*w->var(plotVars[1].c_str())  ,
                *pn1F("meanF"),*pn1F("sigmaF") ,*pn1V("alpha"),*pn1V("n"),*pn1V("alpha2"),*pn1V("n2"));
        std::cout <<"IMPORT!3"<<std::endl;
        w->import(modelMVV);
        std::cout <<"IMPORT!4"<<std::endl;
        RooProdPdf condProd((std::string("model")+postFix).c_str(), (pn0("model")+"*"+pn1("model")).c_str(),*pn0P("model"),RooFit::Conditional(*pn1P("model"), *w->var(plotVars[1].c_str())) );
        std::cout <<"IMPORT!5"<<std::endl;
        w->import(condProd);
        std::cout <<"IMPORT!6"<<std::endl;
    }
};

class ThreeGausExpoFunctionFitter : public FunctionFitter{
public:
    ThreeGausExpoFunctionFitter(const TH1* iH, const std::string postFix = "", const std::vector<std::string>&  plotVars = {"M"}) :
    FunctionFitter(iH,postFix,plotVars){
        std::string model = std::string("model")+postFix;
        addParam("meanW","[90,80,110]");
        addParam("sigmaW","[8,2,15]");
        addParam("fW","[0.2,0,1]");
        addParam("meanTop","[180,170,195]");
        addParam("sigmaTop","[15,5,20]");
        addParam("meanMix","[100,95,150]");
        addParam("sigmaMix","[35,30,45]");
        addParam("fMix","[0.5,0,1]");
        addParam("slopeD","[-0.05,-0.2,-0.001]");
        addParam("fD","[.03,0,1]");

        w->factory((std::string("RooExponential:") + model +"D("+vars[0]+",slopeD)").c_str());
        w->factory((std::string("RooGaussian:") + model +"W("+vars[0]+",meanW,sigmaW)").c_str());
        w->factory((std::string("RooGaussian:") + model +"Top("+vars[0]+",meanTop,sigmaTop)").c_str());
        w->factory((std::string("RooGaussian:") + model +"Mix("+vars[0]+",meanMix,sigmaMix)").c_str());
        RooAddPdf addModel(model.c_str(),model.c_str(),RooArgList(*w->pdf((model+"D").c_str()),*w->pdf((model+"W").c_str()),*w->pdf((model+"Mix").c_str()),*w->pdf((model+"Top").c_str())),
                RooArgList(*w->var("fD"),*w->var("fW"),*w->var("fMix")),true);
        w->import(addModel);
    }
};


class FunctionParameterPlotter{
public:
    std::vector<std::unique_ptr<TCanvas>> cans;
    std::vector<std::unique_ptr<TGraphErrors>> graphs;
    std::vector<std::unique_ptr<TH1>> hists;
    int curPt = 0;

    FunctionParameterPlotter(){}

void addFit(FunctionFitter* fitter, double pt, std::string name){
    if(!curPt){
        for(const auto& p: fitter->params){
            graphs.emplace_back(new TGraphErrors);
            graphs.back()->SetTitle(p.c_str());
            graphs.back()->SetName(p.c_str());
        }
        graphs.emplace_back(new TGraphErrors);
        graphs.back()->SetTitle("chi2");
        graphs.back()->SetName("chi2");
    }


    double chi2 = 0;
    cans.emplace_back(fitter->projection(std::string("can_")+name,chi2));

    for(unsigned int iV = 0; iV< fitter->params.size(); ++iV ){
        graphs[iV]->SetPoint(curPt,pt,fitter->w->var(fitter->params[iV].c_str())->getVal());
        graphs[iV]->SetPointError(curPt,0.0,fitter->w->var(fitter->params[iV].c_str())->getError());
    }
    graphs.back()->SetPoint(curPt,pt,chi2);
    graphs.back()->SetPointError(curPt,0.0,0.0);
    curPt++;
}
void addFit2D(FunctionFitter* fitter, double pt, std::string name){
    if(!curPt){
        for(const auto& p: fitter->params){
            graphs.emplace_back(new TGraphErrors);
            graphs.back()->SetTitle(p.c_str());
            graphs.back()->SetName(p.c_str());
        }
    }
    hists.emplace_back(fitter->data2D(std::string("data_")+name));
    hists.emplace_back(fitter->pdf2D(std::string("pdf_")+name));

    for(unsigned int iV = 0; iV< fitter->params.size(); ++iV ){
        graphs[iV]->SetPoint(curPt,pt,fitter->w->var(fitter->params[iV].c_str())->getVal());
        graphs[iV]->SetPointError(curPt,0.0,fitter->w->var(fitter->params[iV].c_str())->getError());
    }
    curPt++;
}

    void write(const std::string& outFileName){
        TFile * oF = new TFile(outFileName.c_str(),"recreate");
        oF->cd();
        for(auto& c : cans) c->Write();
        for(auto& c : graphs) c->Write();
        for(auto& c : hists) c->Write();
        oF->Close();
    }
};


#endif
