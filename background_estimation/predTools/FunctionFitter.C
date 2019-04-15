
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "PDFAdder.h"
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
    std::unique_ptr<TAxis> xAxis;
    std::unique_ptr<TAxis> yAxis;
    std::unique_ptr<TAxis> zAxis;
    //for making histograms out of a json
    int    plot_bins = -1;
    double plot_min = -1;
    double plot_max = -1;

    RooFitResult *fitResultCache=0;
    std::vector<std::string> secondaryModelNames;

    //If the input histogram is not defined, get the axes from the TAxis vector
    FunctionFitter(const TH1* iH, const std::vector<TAxis*>& axs,
            const std::string postFix = "T",const std::vector<std::string>&  plotVars = {"M"}):
        postFix(postFix)
    {
        w.reset(new RooWorkspace("w",false));
        for(const auto& v: plotVars){
            vars.push_back(v);
            w->factory((v+"[0,10000]").c_str());
        }
        if(iH){
            PDFAdder::addBinnedData(w.get(),iH,"data",vars);
            if(vars.size()>0) xAxis.reset((TAxis*)iH->GetXaxis()->Clone());
            if(vars.size()>1) yAxis.reset((TAxis*)iH->GetYaxis()->Clone());
            if(vars.size()>2) zAxis.reset((TAxis*)iH->GetZaxis()->Clone());
        } else {
            if(vars.size()>0){
                PDFAdder::setBinning(w.get(),vars[0],axs[0]);
                xAxis.reset((TAxis*)axs[0]->Clone());
            }
            if(vars.size()>1){
                PDFAdder::setBinning(w.get(),vars[1],axs[2]);
                yAxis.reset((TAxis*)axs[1]->Clone());
            }
            if(vars.size()>2){
                PDFAdder::setBinning(w.get(),vars[2],axs[2]);
                zAxis.reset((TAxis*)axs[2]->Clone());
            }
        }
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
    double getError(std::string var) {
        return w->var(var.c_str())->getError();
    }
    double getFVal(std::string var) {
        return w->function(var.c_str())->getVal();
    }
    double getFError(std::string var) {
        if(fitResultCache)
            return w->function(var.c_str())->getPropagatedError(*fitResultCache);
        else return -1;
    }

    void fit(const std::vector<RooCmdArg>& options){
        std::string thisModel = std::string("model")+postFix;
        std::string thisData  = "data";
        if(options.size()==0) fitResultCache =  w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()));
        else if(options.size()==1) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0]);
        else if(options.size()==2) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0],options[1]);
        else if(options.size()==3) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0],options[1],options[2]);
        else if(options.size()==4) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0],options[1],options[2],options[3]);
        else if(options.size()==5) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0],options[1],options[2],options[3],options[4]);
        else if(options.size()==6) fitResultCache = w->pdf(thisModel.c_str())->fitTo(
                *w->data(thisData.c_str()),options[0],options[1],options[2],options[3],options[4],
                options[5]);
    }

    TCanvas* projection(const std::string& name, double& chi2){
        std::string model = std::string("model")+postFix;
        std::string data  = std::string("data");

        auto* dataH =  w->data(data.c_str())->
                createHistogram((name+"_data").c_str(),*w->var(vars[0].c_str()));
        auto * pdf = createTH1FromPDF( w->pdf(model.c_str()),
                w->var(vars[0].c_str()),model.c_str(),"",xAxis.get());

        pdf->Scale(dataH->Integral()/pdf->Integral());
//        if(secondaryModelNames.size())
//            w->pdf(model.c_str())->plotOn(frame, RooFit::NormRange("fit"),RooFit::Binning(xBins),
//                 RooFit::Components(secondaryModelNames[0].c_str()), RooFit::LineStyle(kDashed));
        TCanvas* can = new TCanvas(name.c_str());
        can->cd();
        dataH->Draw();
        dataH->GetYaxis()->SetTitle(name.c_str());
        dataH->GetXaxis()->SetTitle(vars[0].c_str());
        dataH->SetTitle("");
        pdf->SetLineColor(kBlue);
        pdf->Draw("SAME HIST");
//        chi2 = frame->chiSquare();
        return can;
    }
    TCanvas* projection2D(const std::string& name, double& chi2){
        std::string model = std::string("model")+postFix;
//        w->pdf(model.c_str())->fixAddCoefRange("coef");
        auto * pdf = createTH2FromPDF( w->pdf(model.c_str()),
                w->var(vars[0].c_str()),w->var(vars[1].c_str()),model.c_str(),"",
                xAxis.get(),yAxis.get());

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
        auto * hist = w->data(data.c_str())->createHistogram(name.c_str(),
                *w->var(vars[0].c_str()),RooFit::YVar(*w->var(vars[1].c_str()))   );
        hist->GetYaxis()->SetTitle(vars[1].c_str());
        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }
    TH1* pdf2D(const std::string& name){
        std::string model = std::string("model")+postFix;
//        w->pdf(model.c_str())->fixAddCoefRange("coef");

        auto * hist = createTH2FromPDF( w->pdf(model.c_str()),
                w->var(vars[0].c_str()),w->var(vars[1].c_str()),name.c_str(),"",
                xAxis.get(),yAxis.get());
        hist->GetYaxis()->SetTitle(vars[1].c_str());
        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }
    TH1* pdf1D(const std::string& name){
        std::string model = std::string("model")+postFix;
//        w->pdf(model.c_str())->fixAddCoefRange("coef");
        auto * hist = createTH1FromPDF( w->pdf(model.c_str()), w->var(vars[0].c_str()),
                name.c_str(),"",xAxis.get());

        hist->GetXaxis()->SetTitle(vars[0].c_str());
        return hist;
    }


protected:
    void addParam(std::string name, std::string bounds){
        params.push_back(name);
        w->factory((name + bounds).c_str());
    }


    void addCB(const std::string& modelName, const std::string& pPFix,  const std::string& varName,
            bool addParams = true, RooAbsReal* meanFunc =0, RooAbsReal* sigmaFunc =0){
        auto pn =[&] (const std::string& v) ->std::string {return v+pPFix;};
        auto pv =[&] (const std::string& v) ->RooRealVar* {return w->var((v+pPFix).c_str());};

        if(addParams){
            if(meanFunc==0) addParam(pn("mean"  ),"[90,0,5000]");
            if(sigmaFunc==0) addParam(pn("sigma" ),"[8,0,5000]");
            addParam(pn("alpha" ),"[1,0.001,100]");
            addParam(pn("alpha2"),"[1,0.001,100]");
            addParam(pn("n")  ,"[5,1,100]");
            addParam(pn("n2"),"[5,1,100]");
        }



        RooDoubleCB modelP(modelName .c_str(),modelName.c_str(),*w->var(varName.c_str())  ,
                meanFunc  ? *meanFunc : *pv("mean"), sigmaFunc  ? *sigmaFunc : *pv("sigma"),
                *pv("alpha"),*pv("n") ,*pv("alpha2"),*pv("n2"));
        w->import(modelP);
    }

    void addExpo(const std::string& modelName, const std::string& pPFix, const std::string& varName,
            bool addParams = true){
        auto pn =[&] (const std::string& v) ->std::string {return v+pPFix;};
        if(addParams){
            addParam(pn("slope"),"[-1,-10,0]");
        }
        RooExponential modelE(modelName.c_str(),modelName.c_str(),*w->var(varName.c_str()),
                *w->var(pn("slope").c_str()));
        w->import(modelE);
    }

    void addGaus(const std::string& modelName, const std::string& pPFix, const std::string& varName,
            bool addParams = true, RooAbsReal* meanFunc =0, RooAbsReal* sigmaFunc =0){
        auto pn =[&] (const std::string& v) ->std::string {return v+pPFix;};
        auto pv =[&] (const std::string& v) ->RooRealVar* {return w->var((v+pPFix).c_str());};
        if(addParams){
            if(meanFunc==0) addParam(pn("mean"  ),"[90,0,5000]");
            if(sigmaFunc==0) addParam(pn("sigma" ),"[8,0,1000]");;
        }
        RooGaussian model(modelName.c_str(),modelName.c_str(),*w->var(varName.c_str()),
                meanFunc  ? *meanFunc : *pv("mean"),sigmaFunc  ? *sigmaFunc : *pv("sigma"));
        w->import(model);
    }



};
class CBFunctionFitter : public FunctionFitter{
public:
    CBFunctionFitter(const TH1* iH, const std::vector<TAxis*>& axs,
            bool doExpo,
            const std::string postFix,const std::vector<std::string>&  plotVars = {"M"}) :
        FunctionFitter(iH,axs,postFix,plotVars){
        //        addCB(postFix,std::string("model")+postFix,vars[0],doExpo);
        const std::string modelName = std::string("model")+postFix;
        const std::string cbName = doExpo ?modelName +"P" : modelName;
        addCB(cbName,postFix,vars[0]);
        if(doExpo){
            auto pen =[&] (const std::string& v) ->std::string {return v+postFix;};
            addParam(pen("fE"),"[0.1,0,1]");
            const std::string expName =  modelName +"E";
            addExpo(expName,postFix,vars[0]);
            RooAddPdf modelC(modelName.c_str(),modelName.c_str(),
                    *w->pdf(expName.c_str()),*w->pdf(cbName.c_str()),*w->var(pen("fE").c_str()));
            secondaryModelNames.push_back(expName);
            w->import(modelC);
        }
    }
};


class CBFunctionFitter2DWithSepExpo : public FunctionFitter{
public:
    CBFunctionFitter2DWithSepExpo(const TH2* iH, const std::vector<TAxis*>& axs,
            bool doExpoX,
            const std::string postFix,const std::vector<std::string>&  plotVars = {"MJJ","MHH"}) :
        FunctionFitter(iH,axs,postFix,plotVars){

        const std::string xPF = postFix+plotVars[0];
        const std::string yPF = postFix+plotVars[1];
        const std::string mN = std::string("model")+postFix;
        const std::string mNP  = mN+ (doExpoX ? "P" : "");
        const std::string mNE  = mN+ "E";
        const std::string mXNP = std::string("model") + xPF + (doExpoX ? "P" : "");
        const std::string mYNP = std::string("model") + yPF + (doExpoX ? "P" : "");
        const std::string mXNE = std::string("model") + xPF + "E" ;
        const std::string mYNE = std::string("model") + yPF + "E" ;

        auto pXN =[&] (std::string v) ->std::string{return v+xPF;};
        auto pXV =[&] (std::string v) ->RooRealVar*{return w->var(pXN(v).c_str());};
        auto pYN =[&] (std::string v) ->std::string{return v+yPF;};
        auto pYV =[&] (std::string v) ->RooRealVar*{return w->var(pYN(v).c_str());};
        auto pYF =[&] (std::string v) ->RooAbsReal*{return w->function(pYN(v).c_str());};

        addCB(mXNP,xPF,vars[0]);

        addParam(pYN("mean"  ),"[90,0,5000]");
        addParam(pYN("sigma" ),"[8,0,1000]");
        addParam(pYN("alpha" ),"[1,0.001,100]");
        addParam(pYN("alpha2"),"[1,0.001,100]");
        addParam(pYN("n")        ,"[5,1,100]");
        addParam(pYN("n2")       ,"[5,1,100]");
        addParam(pYN("maxS")     ,"[2.5,0,5]");
        addParam(pYN("mean_p1"  ),"[0,-5000,5000]");
        addParam(pYN("sigma_p1" ),"[0,-1000,1000]");
        RooFormulaVar xSig  (pYN("xSig"  ).c_str(),"(@0-@1)/@2",
                RooArgList(*w->var((plotVars[0]).c_str()),*pXV("mean"),*pXV("sigma")));
        w->import(xSig  );
        RooFormulaVar xSigC  (pYN("xSigC"  ).c_str(),"max(-1*@0,min(@1,@0))",
                RooArgList(*pYV("maxS"  ),*pYF("xSig"  )));
        w->import(xSigC  );
        RooFormulaVar meanYDepX  (pYN("meanYDepX"  ).c_str(),"@0*(1+@1*@2)"     ,
                RooArgList(*pYV("mean"  ),*pYV("mean_p1"  ),*pYF("xSigC"  )));
        RooFormulaVar sigmaYDepX (pYN("sigmaYDepX" ).c_str(),"@0*(1+(@2>=0?0:@1*abs(@2)))",
                RooArgList(*pYV("sigma" ),*pYV("sigma_p1" ),*pYF("xSigC"  )));
        w->import(meanYDepX  );
        w->import(sigmaYDepX );

        addCB(mYNP,yPF,vars[1],false,pYF("meanYDepX"),pYF("sigmaYDepX"));
        RooProdPdf condProdP(mNP.c_str(), (mXNP+"*"+mYNP).c_str(),*w->pdf(mXNP.c_str()),
                RooFit::Conditional(*w->pdf(mYNP.c_str()), *w->var(plotVars[1].c_str())) );
        w->import(condProdP);

        if(doExpoX){
            addParam(pXN("fE"),"[0.1,0,1]");
            addExpo(mXNE,xPF,vars[0]);
            addGaus(mYNE,std::string("E")+yPF,vars[1],true);

            RooProdPdf condProdE(mNE.c_str(), (mXNE+"*"+mYNE).c_str(),
                    *w->pdf(mXNE.c_str()),*w->pdf(mYNE.c_str()));
            w->import(condProdE);

            RooAddPdf model(mN.c_str(),mN.c_str(),*w->pdf(mNE.c_str()),
                    *w->pdf(mNP.c_str()),*pXV("fE"));
            w->import(model);
        }
    }
};

class CBFunctionFitter2D : public FunctionFitter{
public:
    CBFunctionFitter2D(const TH2* iH, const std::vector<TAxis*>& axs,
            bool doExpoX,
            const std::string postFix, const std::vector<std::string>&  plotVars = {"MJJ","MHH"}) :
        FunctionFitter(iH,axs,postFix,plotVars){

        const std::string xPF = postFix+plotVars[0];
        const std::string yPF = postFix+plotVars[1];
        const std::string mN = std::string("model")+postFix;
        const std::string mNX  = std::string("model") + xPF;
        const std::string mNY  = std::string("model") + yPF;
        const std::string mXNP = std::string("model") + xPF + (doExpoX ? "P" : "");
        const std::string mXNE = std::string("model") + xPF + "E" ;

        auto pXN =[&] (std::string v) ->std::string{return v+xPF;};
        auto pXV =[&] (std::string v) ->RooRealVar*{return w->var(pXN(v).c_str());};
        auto pYN =[&] (std::string v) ->std::string{return v+yPF;};
        auto pYV =[&] (std::string v) ->RooRealVar*{return w->var(pYN(v).c_str());};
        auto pYF =[&] (std::string v) ->RooAbsReal*{return w->function(pYN(v).c_str());};

        addCB(mXNP,xPF,vars[0]);
        if(doExpoX){
            addParam(pXN("fE"),"[0.1,0,1]");
            addExpo(mXNE,xPF,vars[0]);;
            RooAddPdf modelC(mNX.c_str(),mNX.c_str(),
                    *w->pdf(mXNE.c_str()),*w->pdf(mXNP.c_str()),*pXV("fE"));
            w->import(modelC);
        }

        addParam(pYN("mean"  ),"[90,0,5000]");
        addParam(pYN("sigma" ),"[8,0,1000]");
        addParam(pYN("alpha" ),"[1,0.001,100]");
        addParam(pYN("alpha2"),"[1,0.001,100]");
        addParam(pYN("n")        ,"[5,1,100]");
        addParam(pYN("n2")       ,"[5,1,100]");
        addParam(pYN("maxS")     ,"[2.5,0,5]");
        addParam(pYN("mean_p1"  ),"[0,-5000,5000]");
        addParam(pYN("sigma_p1" ),"[0,-1000,1000]");
        addParam(pYN("sigma_p2" ),"[0,-1,1]");
        RooFormulaVar xSig  (pYN("xSig"  ).c_str(),"(@0-@1)/@2",
                RooArgList(*w->var((plotVars[0]).c_str()),*pXV("mean"),*pXV("sigma")));
        w->import(xSig  );
        RooFormulaVar xSigC  (pYN("xSigC"  ).c_str(),"max(-1*@0,min(@1,@0))",
                RooArgList(*pYV("maxS"  ),*pYF("xSig"  )));
        w->import(xSigC  );


        RooFormulaVar meanYDepX  (pYN("meanYDepX"  ).c_str(),"@0*(1+@1*@2)",
                RooArgList(*pYV("mean"  ),*pYV("mean_p1"  ),*pYF("xSigC"  )));
        RooFormulaVar sigmaYDepX (pYN("sigmaYDepX" ).c_str(),"@0*(1+(@2>=0?0:@1*abs(@2)))",
                RooArgList(*pYV("sigma" ),*pYV("sigma_p1" ),*pYF("xSigC"  )));
        w->import(meanYDepX  );
        w->import(sigmaYDepX );
        addCB(mNY,yPF,vars[1],false,pYF("meanYDepX"),pYF("sigmaYDepX"));
        RooProdPdf condProdP(mN.c_str(), (mNX+"*"+mNY).c_str(),*w->pdf(mNX.c_str()),
                RooFit::Conditional(*w->pdf(mNY.c_str()), *w->var(plotVars[1].c_str())) );
        w->import(condProdP);
    }

};


class CBFunctionFitter2DNoCond : public FunctionFitter{
public:
    CBFunctionFitter2DNoCond(const TH2* iH, const std::vector<TAxis*>& axs,
            bool doExpoX,
            const std::string postFix,const std::vector<std::string>&  plotVars = {"MJJ","MHH"}) :
        FunctionFitter(iH,axs,postFix,plotVars){

        const std::string xPF = postFix+plotVars[0];
        const std::string yPF = postFix+plotVars[1];
        const std::string mN = std::string("model")+postFix;
        const std::string mNX  = std::string("model") + xPF;
        const std::string mNY  = std::string("model") + yPF;
        const std::string mXNP = std::string("model") + xPF + (doExpoX ? "P" : "");
        const std::string mXNE = std::string("model") + xPF + "E" ;

        auto pXN =[&] (std::string v) ->std::string{return v+xPF;};
        auto pXV =[&] (std::string v) ->RooRealVar*{return w->var(pXN(v).c_str());};
        auto pYN =[&] (std::string v) ->std::string{return v+yPF;};
        auto pYV =[&] (std::string v) ->RooRealVar*{return w->var(pYN(v).c_str());};
        auto pYF =[&] (std::string v) ->RooAbsReal*{return w->function(pYN(v).c_str());};

        addCB(mXNP,xPF,vars[0]);
        if(doExpoX){
            addParam(pXN("fE"),"[0.1,0,1]");
            addExpo(mXNE,xPF,vars[0]);;
            RooAddPdf modelC(mNX.c_str(),mNX.c_str(),*w->pdf(mXNE.c_str()),
                    *w->pdf(mXNP.c_str()),*pXV("fE"));
            w->import(modelC);
        }

        addCB(mNY,yPF,vars[1]);
        RooProdPdf condProdP(mN.c_str(), (mNX+"*"+mNY).c_str(),
                RooArgSet(*w->pdf(mNX.c_str()), *w->pdf(mNY.c_str())) );
        w->import(condProdP);
    }

};

class CBFunctionFitter2DCondMVV : public FunctionFitter{
public:
    CBFunctionFitter2DCondMVV(const TH2* iH, const std::vector<TAxis*>& axs,
            bool doExpoX,
            const std::string postFix,const std::vector<std::string>&  plotVars = {"MJJ","MHH"}) :
        FunctionFitter(iH,axs,postFix,plotVars){

        const std::string xPF = postFix+plotVars[0];
        const std::string yPF = postFix+plotVars[1];
        const std::string mN = std::string("model")+postFix;
        const std::string mNX  = std::string("model") + xPF;
        const std::string mNY  = std::string("model") + yPF;
        const std::string mXNP = std::string("model") + xPF + (doExpoX ? "P" : "");
        const std::string mXNE = std::string("model") + xPF + "E" ;

        auto pXN =[&] (std::string v) ->std::string{return v+xPF;};
        auto pXV =[&] (std::string v) ->RooRealVar*{return w->var(pXN(v).c_str());};
        auto pYN =[&] (std::string v) ->std::string{return v+yPF;};
        auto pYV =[&] (std::string v) ->RooRealVar*{return w->var(pYN(v).c_str());};
        auto pXF =[&] (std::string v) ->RooAbsReal*{return w->function(pXN(v).c_str());};

        addCB(mNY,yPF,vars[1]);

        addParam(pXN("mean"  ),"[90,0,5000]");
        addParam(pXN("sigma" ),"[8,0,1000]");
        addParam(pXN("alpha" ),"[1,0.001,100]");
        addParam(pXN("alpha2"),"[1,0.001,100]");
        addParam(pXN("n")        ,"[5,1,100]");
        addParam(pXN("n2")       ,"[5,1,100]");
        addParam(pXN("maxS")     ,"[2.5,0,5]");
        addParam(pXN("mean_p1"  ),"[0,-5000,5000]");

        RooFormulaVar ySig  (pXN("ySig"  ).c_str(),"(@0-@1)/@2",
                RooArgList(*w->var((plotVars[1]).c_str()),*pYV("mean"),*pYV("sigma")));
        w->import(ySig  );
        RooFormulaVar ySigC  (pXN("ySigC"  ).c_str(),"max(-1*@0,min(@1,@0))",
                RooArgList(*pXV("maxS"  ),*pXF("ySig"  )));
        w->import(ySigC  );

        RooFormulaVar meanXDepY  (pXN("meanXDepY"  ).c_str(),"@0*(1+@1*@2)",
                RooArgList(*pXV("mean"  ),*pXV("mean_p1"  ),*pXF("ySigC"  )));
        w->import(meanXDepY  );


        addCB(mNX,xPF,vars[0],false,pXF("meanXDepY"));

        if(doExpoX){
            addParam(pXN("fE"),"[0.1,0,1]");
            addExpo(mXNE,xPF,vars[0]);;
            RooAddPdf modelC(mNX.c_str(),mNX.c_str(),*w->pdf(mXNE.c_str()),*w->pdf(mXNP.c_str()),
                    *pXV("fE"));
            w->import(modelC);
        }
        RooProdPdf condProdP(mN.c_str(), (mNY+"*"+mNX).c_str(),*w->pdf(mNY.c_str()),
                RooFit::Conditional(*w->pdf(mNX.c_str()), *w->var(plotVars[0].c_str())) );
        w->import(condProdP);
    }

};

class ThreeGausExpoFunctionFitter : public FunctionFitter{
public:
    ThreeGausExpoFunctionFitter(const TH1* iH,  const std::vector<TAxis*>& axs,
            const std::string postFix = "",const std::vector<std::string>&  plotVars = {"M"}) :
        FunctionFitter(iH,axs,postFix,plotVars){
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

        w->factory((
                std::string("RooExponential:") + model +"D("+vars[0]+",slopeD)"
                ).c_str());
        w->factory((
                std::string("RooGaussian:") + model +"W("+vars[0]+",meanW,sigmaW)"
                ).c_str());
        w->factory((
                std::string("RooGaussian:") + model +"Top("+vars[0]+",meanTop,sigmaTop)"
                ).c_str());
        w->factory((
                std::string("RooGaussian:") + model +"Mix("+vars[0]+",meanMix,sigmaMix)"
                ).c_str());
        RooAddPdf addModel(model.c_str(),model.c_str(),
                RooArgList(*w->pdf((model+"D").c_str()),*w->pdf((model+"W").c_str()),
                        *w->pdf((model+"Mix").c_str()),*w->pdf((model+"Top").c_str())),
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
            double val = fitter->w->var(fitter->params[iV].c_str())->getVal();
            double error = fitter->w->var(fitter->params[iV].c_str())->getError();
            double max = fitter->w->var(fitter->params[iV].c_str())->getMax();
            double min = fitter->w->var(fitter->params[iV].c_str())->getMin();
            if( std::fabs(max - val)/(max-min) <= 0.0001) error = (max-min);
            if( std::fabs(val-min)/(max-min) <= 0.0001) error = (max-min);
            graphs[iV]->SetPoint(curPt,pt,val);
            graphs[iV]->SetPointError(curPt,0.0,error);
        }
        graphs.back()->SetPoint(curPt,pt,chi2);
        graphs.back()->SetPointError(curPt,0.0,0.0);
        curPt++;
    }
    void addFit(const TGraphErrors* graph, const std::string& name){
        graphs.emplace_back((TGraphErrors*)graph->Clone(name.c_str()));
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
        hists.emplace_back(fitter->pdf2D(std::string("pdf_")+name+"__"+MOD_MJ+"_"+MOD_MR));

        for(unsigned int iV = 0; iV< fitter->params.size(); ++iV ){
            double val = fitter->w->var(fitter->params[iV].c_str())->getVal();
            double error = fitter->w->var(fitter->params[iV].c_str())->getError();
            double max = fitter->w->var(fitter->params[iV].c_str())->getMax();
            double min = fitter->w->var(fitter->params[iV].c_str())->getMin();
            if( std::fabs(max - val)/(max-min) <= 0.0001) error = (max-min);
            if( std::fabs(val-min)/(max-min) <= 0.0001) error = (max-min);
            graphs[iV]->SetPoint(curPt,pt,val);
            graphs[iV]->SetPointError(curPt,0.0,error);
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
