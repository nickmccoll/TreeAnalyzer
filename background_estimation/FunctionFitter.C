
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

    void loadJSON(const std::string& jsonFN,const std::string& var, const int nBinsX, const double minX, const double maxX){
        std::ifstream file (jsonFN);
        if (!file.is_open())  throw std::invalid_argument("FunctionFitter::loadJSON() -> Bad file");
        varDep.resize(params.size());
        std::stringstream strStream;
        strStream << file.rdbuf();
        std::string str = strStream.str();
        str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
        str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
        std::vector<std::string> paramFits(std::sregex_token_iterator(str.begin(), str.end(), std::regex("\", \""), -1), std::sregex_token_iterator());
        if(paramFits.size() != params.size()) throw std::invalid_argument("FunctionFitter::loadJSON() -> not the correct number of parameters");
        for(const auto& s :paramFits){
            std::vector<std::string> ps(std::sregex_token_iterator(s.begin(), s.end(), std::regex("\": \""), -1), std::sregex_token_iterator());
            if(ps.size() != 2) {
                for(auto& p :ps) std::cout << p <<" ";
                std::cout <<std::endl;
                throw std::invalid_argument("FunctionFitter::loadJSON() -> Bad parsing");
            }
            ps[0].erase(std::remove(ps[0].begin(), ps[0].end(), '"'), ps[0].end());
            ps[1].erase(std::remove(ps[1].begin(), ps[1].end(), '"'), ps[1].end());
            int pInd=-1;
            for(unsigned int iP = 0; iP<params.size();iP++){
                if(params[iP] != ps[0]) continue;
                pInd =iP; break;
            }
            if(pInd < 0) {
                for(auto& p :ps) std::cout << p <<" ";
                std::cout <<std::endl;
                throw std::invalid_argument("FunctionFitter::loadJSON() -> Bad parameter");
            }
            std:size_t index = 0;
            while (true) {
                 index = ps[1].find(var, index);
                 if (index == std::string::npos) break;
                 ps[1].replace(index, 2, "x");
                 index += 1;
            }
            varDep[pInd].reset(new TF1(ps[0].c_str(),ps[1].c_str(),1,13000));
        }

        plot_bins = nBinsX;
        plot_min  = minX;
        plot_max  = maxX;
        w->var(vars[0].c_str())->setMin(plot_min);
        w->var(vars[0].c_str())->setMax(plot_max);
        w->var(vars[0].c_str())->setBins(plot_bins);
        w->var(vars[0].c_str())->setVal((plot_min+plot_max)/2.);

    }
    double evalFitFromJSON(std::string paramName, double varV){
        for(unsigned int iP = 0; iP < params.size(); ++iP){
            if(params[iP] == paramName) return varDep[iP]->Eval(varV);
        }
        throw std::invalid_argument("FunctionFitter::evalFitFromJSON() -> Bad name");

    }
    TH1* getHistFromJSON(std::string histName,const double varV){
        for(unsigned int iP = 0; iP < params.size(); ++iP){
            setVar(params[iP],varDep[iP]->Eval(varV));
        }
        std::string thisModel = std::string("model")+postFix;
        return w->pdf(thisModel.c_str())->createHistogram(histName.c_str(),*w->var(vars[0].c_str()),RooFit::Binning(RooBinning (plot_bins,plot_min,plot_max))) ;
    }


    void fit(const std::vector<RooCmdArg>& options){
        std::string thisModel = std::string("model")+postFix;
        std::string thisData  = "data";
        if(options.size()==0) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()));
        else if(options.size()==1) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0]);
        else if(options.size()==2) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1]);
        else if(options.size()==3) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1],options[2]);
        else if(options.size()==4) w->pdf(thisModel.c_str())->fitTo(*w->data(thisData.c_str()),options[0],options[1],options[2],options[3]);
    }

    TCanvas* projection(const std::string& name, double& chi2){
        std::string model = std::string("model")+postFix;
        std::string data  = std::string("data");
        auto frame = w->var(vars[0].c_str())->frame();
        w->data(data.c_str())->plotOn(frame);
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

};

class DoubleCBFunctionFitter : public FunctionFitter{
public:
    DoubleCBFunctionFitter(const TH1* iH, const std::string postFix = "", const std::vector<std::string>&  plotVars = {"M"}) :
    FunctionFitter(iH,postFix,plotVars){
        std::string model = std::string("model")+postFix;
        addParam("meanW"  ,"[90,80,110]");
        addParam("sigmaW" ,"[8,5,20]");
        addParam("alphaW" ,"[1,0.1,10]");
        addParam("alphaW2","[1,0.1,10]");

        addParam("meanTop"  ,"[180,140,200]");
        addParam("sigmaTop" ,"[15,5,30]");
        addParam("alphaTop" ,"[1,0.1,10]");
        addParam("alphaTop2","[1,0.1,10]");

        addParam("n","[5]");
        addParam("fW","[0.5,0,1]");

        RooDoubleCB modelW((model+"W").c_str(),(model+"W").c_str(),*w->var(vars[0].c_str()),*w->var("meanW"),*w->var("sigmaW"),*w->var("alphaW"),*w->var("n"),*w->var("alphaW2"),*w->var("n"));
        RooDoubleCB modelTop((model+"Top").c_str(),(model+"Top").c_str(),*w->var(vars[0].c_str()),*w->var("meanTop"),*w->var("sigmaTop"),*w->var("alphaTop"),*w->var("n"),*w->var("alphaTop2"),*w->var("n"));
        w->import(modelW);
        w->import(modelTop);
        w->factory((std::string("SUM::")+model+"(fW*" + model+"W"+","+model+"Top)").c_str());
    }
};

class CBFunctionFitter : public FunctionFitter{
public:
    CBFunctionFitter(const TH1* iH, const std::string postFix = "", const std::vector<std::string>&  plotVars = {"M"}) :
    FunctionFitter(iH,postFix,plotVars){
        std::string model = std::string("model")+postFix;
        addParam(std::string("mean"  )+postFix,"[90,0,300]");
        addParam(std::string("sigma" )+postFix,"[8,0,100]");
        addParam(std::string("alpha" )+postFix,"[1,0.001,100]");
        addParam(std::string("alpha2")+postFix,"[1,0.001,100]");

        addParam(std::string("n")+postFix ,"[5,1,100]");
        addParam(std::string("n2")+postFix,"[5,1,100]");
        RooDoubleCB modelW(model.c_str(),model.c_str(),*w->var(vars[0].c_str())  ,  *w->var((std::string("mean")  +postFix).c_str()) ,*w->var((std::string("sigma") +postFix).c_str()) ,
                                                                                *w->var((std::string("alpha") +postFix).c_str()) ,*w->var((std::string("n")     +postFix).c_str()) ,
                                                                                *w->var((std::string("alpha2")+postFix).c_str()) ,*w->var((std::string("n2")    +postFix).c_str())
        );
        w->import(modelW);
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

    void write(const std::string& outFileName){
        TFile * oF = new TFile(outFileName.c_str(),"recreate");
        oF->cd();
        for(auto& c : cans) c->Write();
        for(auto& c : graphs) c->Write();
        oF->Close();
    }
};


#endif
