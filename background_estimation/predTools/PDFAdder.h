#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_PARAMSHAPEADDER_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_PARAMSHAPEADDER_H
#include "makeJSON.C"
#include "CutConstants.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "RooExponential.h"

#include <utility>


namespace PDFAdder{

typedef std::pair<std::string,double> StrFlt;
typedef std::vector<StrFlt> StrFlts;

typedef std::pair<std::string,std::string> StrStr;
typedef std::vector<StrStr> StrStrs;

struct InterpSyst {
    InterpSyst(const std::string& hName, const StrStrs& systPs) :
        hName(hName), systPs(systPs){}
    std::string hName; //Histogram name for interpolation
    StrStrs systPs;    //list of systematic parameters
};

struct InterpSysts : public std::vector<InterpSyst> {
    void addSyst(const std::string& hName, const StrStrs& systPs) {
        emplace_back(hName,systPs);
    }
};

class HistRebin{
public:
    virtual ~HistRebin() {}
    virtual TH1 * rebin1D(TH1* inH,const std::string& newName,bool isYAxis1D = false) = 0;
    virtual TH1 * rebin2D(TH1* inH,const std::string& newName) = 0;
    TH1 * process(TH1* inH, const std::string& newName,bool isYAxis1D = false) {
        TH1 * oH =0;
        if (  dynamic_cast<const TH2*>(inH) ){
            oH =rebin2D(inH,newName);
        } else
            oH= rebin1D(inH,newName,isYAxis1D);
        oH->SetDirectory(0);
        return oH;
    }

    TString getTitle(const TH1* inH, bool is2D){
        return is2D ? TString("; ")+  inH->GetXaxis()->GetTitle() + " ; "+inH->GetYaxis()->GetTitle() :
                TString("; ")+  inH->GetXaxis()->GetTitle();
    }
    void add1DData(const TH1* inH, TH1* outH){
        outH->Sumw2(true);
        for(int inX = 1; inX <= inH->GetNbinsX(); ++inX){
            int outX = outH->FindFixBin(inH->GetBinCenter(inX));
            if(outX == 0 || outX > outH->GetNbinsX()) continue;
            outH->SetBinContent(outX, outH->GetBinContent(outX)+inH->GetBinContent(inX)  );
            (*outH->GetSumw2())[outX] += (*inH->GetSumw2())[inX];
        }
    }
    void add2DData(const TH1* inH, TH1* outH){
        outH->Sumw2(true);
        for(int inX = 1; inX <= inH->GetNbinsX(); ++inX){
            int outX = outH->GetXaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(inX));
            if(outX == 0 || outX > outH->GetNbinsX()) continue;
            for(int inY = 1; inY <= inH->GetNbinsY(); ++inY){
                int outY = outH->GetYaxis()->FindFixBin(inH->GetYaxis()->GetBinCenter(inY));
                if(outY == 0 || outY > outH->GetNbinsY()) continue;
                outH->SetBinContent(outX,outY, outH->GetBinContent(outX,outY)+inH->GetBinContent(inX,inY)  );
                int outBin = outH->GetBin(outX,outY);
                int inBin = inH->GetBin(inX,inY);
                (*outH->GetSumw2())[outBin] += (*inH->GetSumw2())[inBin];
            }
        }
    }

};

class HistRebinX : public HistRebin {
public:

    HistRebinX(int nBins, double min, double max) : nBins(nBins), min(min), max(max) {}
    HistRebinX(const  std::vector<double>& bins): bins(bins),nBins(bins.size()-1){}
    virtual ~HistRebinX() {}

    virtual TH1 * rebin1D(TH1* inH,const std::string& newName,bool isYAxis1D = false){
        if(isYAxis1D) return (TH1*)inH->Clone(newName.c_str());
        TH1 * outH = 0;
        if(bins.size())
            outH = new TH1F(newName.c_str(),getTitle(inH,false),nBins,&bins[0]);
        else
            outH = new TH1F(newName.c_str(),getTitle(inH,false),nBins,min,max);
        add1DData(inH,outH);
        return outH;
    }
    virtual TH1 * rebin2D(TH1* inH,const std::string& newName){
        TH1 * outH = 0;
        if(bins.size()){
            if(inH->GetYaxis()->GetXbins()->GetSize())
                outH = new TH2F(newName.c_str(),getTitle(inH,true),nBins,&bins[0],inH->GetNbinsY(),inH->GetYaxis()->GetXbins()->GetArray());
            else
                outH = new TH2F(newName.c_str(),getTitle(inH,true),nBins,&bins[0],inH->GetNbinsY(),inH->GetYaxis()->GetXmin(),inH->GetYaxis()->GetXmax());
        } else {
            if(inH->GetYaxis()->GetXbins()->GetSize())
                outH = new TH2F(newName.c_str(),getTitle(inH,true),nBins,min,max,inH->GetNbinsY(),inH->GetYaxis()->GetXbins()->GetArray());
            else
                outH = new TH2F(newName.c_str(),getTitle(inH,true),nBins,min,max,inH->GetNbinsY(),inH->GetYaxis()->GetXmin(),inH->GetYaxis()->GetXmax());
        }
        add2DData(inH,outH);
        return outH;
    }

    std::vector<double> bins; //set to size 0 if using fixed width bins
    int nBins = -1;
    double min =-1;
    double max =-1;
};

class HistRebinY : public HistRebin {
public:

    HistRebinY(int nBins, double min, double max) : nBins(nBins), min(min), max(max) {}
    HistRebinY(const  std::vector<double>& bins):  bins(bins),nBins(bins.size()-1){}
    virtual ~HistRebinY() {}


    virtual TH1 * rebin1D(TH1* inH,const std::string& newName,bool isYAxis1D = false){
        if(!isYAxis1D) return (TH1*)inH->Clone(newName.c_str());
        TH1 * outH = 0;
        if(bins.size())
            outH = new TH1F(newName.c_str(),getTitle(inH,false),nBins,&bins[0]);
        else
            outH = new TH1F(newName.c_str(),getTitle(inH,false),nBins,min,max);
        add1DData(inH,outH);
        return outH;
    }
    virtual TH1 * rebin2D(TH1* inH,const std::string& newName){
        TH1 * outH = 0;
        if(bins.size()){
            if(inH->GetXaxis()->GetXbins()->GetSize())
                outH = new TH2F(newName.c_str(),getTitle(inH,true),inH->GetNbinsX(),inH->GetXaxis()->GetXbins()->GetArray(),nBins,&bins[0]);
            else
                outH = new TH2F(newName.c_str(),getTitle(inH,true),inH->GetNbinsX(),inH->GetXaxis()->GetXmin(),inH->GetXaxis()->GetXmax(),nBins,&bins[0]);
        } else {
            if(inH->GetXaxis()->GetXbins()->GetSize())
                outH = new TH2F(newName.c_str(),getTitle(inH,true),inH->GetNbinsX(),inH->GetXaxis()->GetXbins()->GetArray(),nBins,min,max);
            else
                outH = new TH2F(newName.c_str(),getTitle(inH,true),inH->GetNbinsX(),inH->GetXaxis()->GetXmin(),inH->GetXaxis()->GetXmax(),nBins,min,max);
        }
        add2DData(inH,outH);
        return outH;
    }
    std::vector<double> bins; //set to size 0 if using fixed width bins
    int nBins = -1;
    double min =-1;
    double max =-1;
};

class HistRebinXY : public HistRebin {
public:

    HistRebinXY(int nBinsX, double minX, double maxX,int nBinsY, double minY, double maxY) :
        nBinsX(nBinsX), minX(minX), maxX(maxX),nBinsY(nBinsY), minY(minY), maxY(maxY) {}
    HistRebinXY(const std::vector<double>& binsX,int nBinsY, double minY, double maxY) :
         binsX(binsX),nBinsX(binsX.size()-1),nBinsY(nBinsY), minY(minY), maxY(maxY) {}
    HistRebinXY(int nBinsX, double minX, double maxX,std::vector<double>& binsY) :
        nBinsX(nBinsX), minX(minX), maxX(maxX), binsY(binsY),nBinsY(binsY.size() -1) {}
    HistRebinXY(const std::vector<double>& binsX,const std::vector<double>& binsY) :
        binsX(binsX),nBinsX(binsX.size()-1),binsY(binsY),nBinsY(binsY.size()-1) {}

    virtual TH1 * rebin1D(TH1* inH,const std::string& newName,bool isYAxis1D = false){
        TH1 * outH = 0;
        if(isYAxis1D){
            if(binsY.size())
                outH = new TH1F(newName.c_str(),getTitle(inH,false),nBinsY,&binsY[0]);
            else
                outH = new TH1F(newName.c_str(),getTitle(inH,false),nBinsY,minY,maxY);
        } else {
            if(binsX.size())
                outH = new TH1F(newName.c_str(),getTitle(inH,false),nBinsX,&binsX[0]);
            else
                outH = new TH1F(newName.c_str(),getTitle(inH,false),nBinsX,minX,maxX);
        }
        add1DData(inH,outH);
        return outH;
    }
    virtual TH1 * rebin2D(TH1* inH,const std::string& newName){
        TH1 * outH = 0;
        if(binsX.size() && binsY.size())
            outH = new TH2F(newName.c_str(),getTitle(inH,true),nBinsX,&binsX[0],nBinsY,&binsY[0]);
        else if(binsX.size() && !binsY.size())
            outH = new TH2F(newName.c_str(),getTitle(inH,true),nBinsX,&binsX[0],nBinsY,minY,maxY);
        else if(!binsX.size() && binsY.size())
            outH = new TH2F(newName.c_str(),getTitle(inH,true),nBinsX,minX,maxX,nBinsY,&binsY[0]);
        else
            outH = new TH2F(newName.c_str(),getTitle(inH,true),nBinsX,minX,maxX,nBinsY,minY,maxY);
        add2DData(inH,outH);
        return outH;
    }


    std::vector<double> binsX; //set to size 0 if using fixed width bins
    int nBinsX = -1;
    double minX =-1;
    double maxX =-1;


    std::vector<double> binsY; //set to size 0 if using fixed width bins
    int nBinsY = -1;
    double minY =-1;
    double maxY =-1;
};


//Add a systematic to the worspace, it makes a new formula variable based on the initial expression and the systematics
//systematics are a list of the systematic variables and their multiplicative scales
//extra params are the parameters that need to be added to the parameter list (usually what is included in the variable expression)
void makeParamFormula (RooWorkspace* w, const std::string& newFormulaName,const std::string& varExp, const StrFlts& systematics, const std::vector<std::string>& extraParams ) {
    //First make the systematic string and the paramlist
    RooArgList paramList;
    for(const auto& p : extraParams) paramList.add(*w->var(p.c_str()));
    std::string systMathStr = "0";
    for(const auto& s : systematics) {
        systMathStr = systMathStr +"+"+ASTypes::flt2Str(s.second)+"*"+s.first;
        paramList.add(*w->var(s.first.c_str()));
    }
    if(systematics.size()) systMathStr = "*(1+"+systMathStr+")";
    else systMathStr = "";

    RooFormulaVar varF(newFormulaName.c_str(),newFormulaName.c_str(),(varExp+systMathStr).c_str(),paramList);
    w->import(varF);
}



void conditionalProduct(RooWorkspace* w, const std::string& pdfName, const std::string& pdfName_cxoy, const std::string& varName_x, const std::string& pdfName_y){
    std::string title = pdfName_cxoy+"("+varName_x+"|y)*"+pdfName_y+"(y)";
    RooProdPdf condProd(pdfName.c_str(),title.c_str(),*w->pdf(pdfName_y.c_str()),RooFit::Conditional(*w->pdf(pdfName_cxoy.c_str()), *w->var(varName_x.c_str()))  );
    w->import(condProd);
}


void addExpo(RooWorkspace* w, const std::string& pdfName,const std::string& pVar,const std::string& varName, const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json){
    auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
    auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
    auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};

    makeParamFormula(w,pN("slope")    ,pJS("slope"),StrFlts(),{pVar});
    RooExponential modelE(pdfName.c_str(),pdfName.c_str(),*w->var(varName.c_str()), *pF("slope"));
    w->import(modelE);
}

void addCB(RooWorkspace* w,const std::string& pdfName,const std::string& pVar, const std::string& varName, const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json,const StrFlts& scale,const StrFlts& resolution,const StrFlts& alpha1 = StrFlts()){
    auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
    auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
    auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};

    makeParamFormula(w,pN("mean"),pJS("mean"),scale,{pVar});
    makeParamFormula(w,pN("sigma"),pJS("sigma"),resolution,{pVar});
    makeParamFormula(w,pN("alpha") ,pJS("alpha"),alpha1,{pVar});
    makeParamFormula(w,pN("n")     ,pJS("n"),StrFlts(),{pVar});
    makeParamFormula(w,pN("alpha2"),pJS("alpha2"),StrFlts(),{pVar});
    makeParamFormula(w,pN("n2")    ,pJS("n2"),StrFlts(),{pVar});

    RooDoubleCB modelP(pdfName.c_str(),pdfName.c_str(),*w->var(varName.c_str()),
            *pF("mean"),*pF("sigma"),*pF("alpha"),*pF("n") ,*pF("alpha2"),*pF("n2"));
    w->import(modelP);
}


void addCondCB(RooWorkspace* w, const std::string& pdfName,const std::string& pVar, const std::string& varName,  const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json,const StrFlts& scale,const StrFlts& resolution,
        const std::string& condVarName, const std::string& condParPostFix
){
    auto cPN = [&](const std::string& v)->std::string{return v +"_"+condParPostFix;};
    auto cPF = [&](const std::string& v)->RooAbsReal*{return w->function(cPN(v).c_str());};

    auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
    auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
    auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};

    makeParamFormula(w,pN("mean"),pJS("mean"),scale,{pVar});
    makeParamFormula(w,pN("sigma"),pJS("sigma"),resolution,{pVar});
    makeParamFormula(w,pN("alpha") ,pJS("alpha"),StrFlts(),{pVar});
    makeParamFormula(w,pN("n")     ,pJS("n"),StrFlts(),{pVar});
    makeParamFormula(w,pN("alpha2"),pJS("alpha2"),StrFlts(),{pVar});
    makeParamFormula(w,pN("n2")    ,pJS("n2"),StrFlts(),{pVar});

    makeParamFormula(w,pN("maxS")     ,pJS("maxS"),StrFlts(),{pVar});
    makeParamFormula(w,pN("mean_p1")  ,pJS("mean_p1"),StrFlts(),{pVar});
    makeParamFormula(w,pN("sigma_p1") ,pJS("sigma_p1"),StrFlts(),{pVar});

    RooFormulaVar xSig  (pN("xSig"  ).c_str(),"(@0-@1)/@2",RooArgList(*w->var(condVarName.c_str()),*cPF("mean"),*cPF("sigma")));
    w->import(xSig  );
    RooFormulaVar xSigC  (pN("xSigC"  ).c_str(),"max(-1*@0,min(@1,@0))",RooArgList(*pF("maxS"  ),*pF("xSig"  )));
    w->import(xSigC  );

    RooFormulaVar meanYDepX  (pN("meanYDepX"  ).c_str(),"@0*(1+@1*@2)"     ,RooArgList(*pF("mean"  ),*pF("mean_p1"  ),*pF("xSigC"  )));
    RooFormulaVar sigmaYDepX (pN("sigmaYDepX" ).c_str(),"@0*(1+(@2>=0?0:@1*abs(@2)))",RooArgList(*pF("sigma" ),*pF("sigma_p1" ),*pF("xSigC"  )));
    w->import(meanYDepX  );
    w->import(sigmaYDepX );

    RooDoubleCB modelP(pdfName.c_str(),pdfName.c_str(),*w->var(varName.c_str()),
            *pF("meanYDepX"),*pF("sigmaYDepX"),*pF("alpha"),*pF("n") ,*pF("alpha2"),*pF("n2"));
    w->import(modelP);
}


void add2DCBNoCond(RooWorkspace* w, const std::string& name,const std::string& PF,const std::string& variableX,const std::string& variableY,const std::string& jsonFile,const StrFlts& scale_X,const StrFlts& resolution_X
        ,const StrFlts& scale_Y,const StrFlts& resolution_Y, bool exponential_X, const std::string& pVar){
    CJSON json( jsonFile);
    const std::string varXCBPDFName = exponential_X ? name+"_"+variableX+"_"+"CB"+"_"+PF : name+"_"+variableX+"_"+PF;
    const std::string varXEPDFName = name+"_"+variableX+"_"+"E"+"_"+PF;
    const std::string varXPF = name+"_"+variableX+"_"+PF;
    const std::string varYPF = name+"_"+variableY+"_"+PF;
    const std::string varXPDFName = name+"_"+variableX+"_"+PF;
    const std::string varYPDFName = name+"_"+variableY+"_"+PF;
    const std::string PDFName = name+"_"+PF;

    addCB(w,varXCBPDFName,pVar,variableX,varXPF,variableX,json,scale_X,resolution_X);
    if(exponential_X){
        std::string vN = std::string("fE") +"_"+varXPF;
        std::string expr = std::string("min(")+ pVar+"*0+"+json.getP(std::string("fE")+variableX)+",1)";
        RooFormulaVar varF(vN.c_str(),vN.c_str(),expr.c_str(),RooArgList(*w->var(pVar.c_str())));
        w->import(varF);
        PDFAdder::addExpo(w,varXEPDFName,pVar,variableX,varXPF,variableX,json);

        RooAddPdf modelC(varXPDFName.c_str(),varXPDFName.c_str(),*w->pdf(varXEPDFName.c_str()),*w->pdf(varXCBPDFName.c_str()),*w->function(vN.c_str()));
        w->import(modelC);
    }

    addCB(w,varYPDFName,pVar,variableY,varYPF,variableY,json,scale_Y,resolution_Y);
    RooProdPdf prodP(PDFName.c_str(), (varXPDFName+"*"+varYPDFName).c_str(),RooArgSet(*w->pdf(varXPDFName.c_str()), *w->pdf(varYPDFName.c_str())) );
    w->import(prodP);


}

void add2DCB(RooWorkspace* w, const std::string& name,const std::string& PF, const std::string& variableX,const std::string& variableY,const std::string& jsonFile,const StrFlts& scale_X,const StrFlts& resolution_X
        ,const StrFlts& scale_Y,const StrFlts& resolution_Y, bool exponential_X, const std::string& pVar){
    CJSON json( jsonFile);
    const std::string varXCBPDFName = exponential_X ? name+"_"+variableX+"_"+"CB"+"_"+PF : name+"_"+variableX+"_"+PF;
    const std::string varXEPDFName = name+"_"+variableX+"_"+"E"+"_"+PF;
    const std::string varXPF = name+"_"+variableX+"_"+PF;
    const std::string varYPF = name+"_"+variableY+"_"+PF;
    const std::string varXPDFName = name+"_"+variableX+"_"+PF;
    const std::string varYPDFName = name+"_"+variableY+"_"+PF;
    const std::string PDFName = name+"_"+PF;

    PDFAdder::addCB(w,varXCBPDFName,pVar,variableX,varXPF,variableX,json,scale_X,resolution_X);
    if(exponential_X){
        std::string vN = std::string("fE") +"_"+varXPF;
        std::string expr = std::string("min(")+ pVar+"*0+"+json.getP(std::string("fE")+variableX)+",1)";
        RooFormulaVar varF(vN.c_str(),vN.c_str(),expr.c_str(),RooArgList(*w->var(pVar.c_str())));
        w->import(varF);
        PDFAdder::addExpo(w,varXEPDFName,pVar,variableX,varXPF,variableX,json);

        RooAddPdf modelC(varXPDFName.c_str(),varXPDFName.c_str(),*w->pdf(varXEPDFName.c_str()),*w->pdf(varXCBPDFName.c_str()),*w->function(vN.c_str()));
        w->import(modelC);
    }
    addCondCB(w,varYPDFName,pVar,variableY,varYPF,variableY,json,scale_Y,resolution_Y,variableX,varXPF);
    conditionalProduct(w,PDFName,varYPDFName,variableY,varXPDFName);
}

//Add a FastVerticleInterp
//Name is the prefix of everything
//PF is the post fix to everything...the final PDF name will be Name_PF
//obs is the number of observables
//filename is the filename to find input histograms...hName is the nominal shape
//Systs are the shape systematics...the first entry is the name in the file....hName_FIRST(Up|Down)...the second is the vector of of the systematics and the relative scales
//Conditional and order are parameters in the roofit classes
//All systematics are assumed to have already been added to the worksapce (def for shapes: addVar(s.second,0,-1,1,false);)
void addHistoShapeFromFile(RooWorkspace* w, std::string name,std::string PF, const std::vector<std::string>& obs, const std::string& fileName, const std::string& hName, const InterpSysts& systs,
        const bool conditional = false, const int order = 0, HistRebin* rebinner= 0, bool isY = false){
    //Save a few steps
    if(PF.size()) PF = "_"+PF;
    const std::string PDFName = name+PF;

    //Argsets and such for the observables
    RooArgSet varset;
    RooArgList varlist;
    std::vector<RooRealVar*> varPointers;
    for(const auto& v : obs){
        varPointers.push_back(w->var(v.c_str()));
        varset.add(*varPointers.back());
        varlist.add(*varPointers.back());
    }
    auto* iF =  TObjectHelper::getFile(fileName);
    auto inh = TObjectHelper::getObject<TH1>(iF,hName);
    TH1 * h = &*inh;
    if(rebinner) h = rebinner->process(&*h,std::string(h->GetName()) +"_rebin",isY);

    const std::string nominalHistName = name+"_Nominal_HIST"+PF;
    const std::string nominalPDFName  = systs.size() ? name+"_Nominal"+PF : PDFName;

    RooDataHist nHist(nominalHistName.c_str(),nominalHistName.c_str(),varlist,h);
    w->import(nHist);
    RooHistPdf nPDF(nominalPDFName.c_str(),nominalPDFName.c_str(),varset,nHist,order);
    w->import(nPDF);

    RooArgList coefList;
    RooArgList pdfList(*w->pdf(nominalPDFName.c_str()));
    for(const auto& s : systs ){

        std::string coefName = name+"_"+s.hName+PF;
        RooArgList paramList;
        std::string systMathStr = "0";
        for(const auto& es : s.systPs) {
            systMathStr = systMathStr +"+("+es.second+")*"+es.first;
            paramList.add(*w->var(es.first.c_str()));
        }
        RooFormulaVar varF(coefName.c_str(),coefName.c_str(),systMathStr.c_str(),paramList);
        w->import(varF);
        coefList.add(*w->arg(coefName.c_str()));

        auto addVar = [&](const std::string& var){
            auto inh =  TObjectHelper::getObject<TH1>(iF,hName+"_"+s.hName+var);
            const std::string vHistName = name+"_" +s.hName+var+"_HIST"+PF;
            const std::string vPDFName = name +"_"+s.hName+var+PF;
            TH1 * sh = &*inh;
            if(rebinner) sh =  rebinner->process(&*sh,std::string(sh->GetName()) +"_rebin",isY);
            RooDataHist vRooHist(vHistName.c_str(),vHistName.c_str(),varlist,sh);
            w->import(vRooHist);
            RooHistPdf pdf(vPDFName.c_str(),vPDFName.c_str(),varset,vRooHist,order);
            w->import(pdf);
            pdfList.add(*w->pdf(vPDFName.c_str()));
            if(rebinner) delete sh;
        };
        addVar("Up");
        addVar("Down");

    }


    if(systs.size()){
        if(obs.size() == 1){
            FastVerticalInterpHistPdf total(PDFName.c_str(),PDFName.c_str(),*w->var(obs[0].c_str()),pdfList,coefList);
            w->import(total);
        } else if(obs.size() == 2){
            FastVerticalInterpHistPdf2D total(PDFName.c_str(),PDFName.c_str(),*w->var(obs[0].c_str()),*w->var(obs[1].c_str()),conditional,pdfList,coefList);
            w->import(total);
        }
    }

    iF->Close();
    delete iF;
    if(rebinner) delete h;
}


void addBinnedData(RooWorkspace* w, const TH1* iH, const std::string& name, const std::vector<std::string>& variables ){
    const unsigned int nD = variables.size();
    RooArgList args;
    auto doBinning =[&](const std::string& var, const TAxis* ax) {
        args.add(*w->var(var.c_str()));
        w->var(var.c_str())->setMin(ax->GetXmin());
        w->var(var.c_str())->setMax(ax->GetXmax());
        w->var(var.c_str())->setBins(ax->GetNbins());
        w->var(var.c_str())->setVal((ax->GetXmin()+ax->GetXmax())/2.);
    };
    doBinning(variables[0],iH->GetXaxis());
    if(nD > 1)doBinning(variables[1],iH->GetYaxis());
    if(nD > 2)doBinning(variables[2],iH->GetZaxis());
    if(nD > 3)throw std::invalid_argument("makeCard::importBinnedData() -> Too many observables!");
    RooDataHist dataHist(name.c_str(),name.c_str(),args,iH);
    w->import(dataHist);
}




};





#endif
