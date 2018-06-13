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

typedef std::pair<std::string,StrFlts> InterpSyst; // [histogram name, vector<systematics>]
typedef std::vector<InterpSyst> InterpSysts;

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

//Add a FastVerticleInterp
//Name is the prefix of everything
//PF is the post fix to everything...the final PDF name will be Name_PF
//obs is the number of observables
//filename is the filename to find input histograms...hName is the nominal shape
//Systs are the shape systematics...the first entry is the name in the file....hName_FIRST(Up|Down)...the second is the vector of of the systematics and the relative scales
//Conditional and order are parameters in the roofit classes
//All systematics are assumed to have already been added to the worksapce (def for shapes: addVar(s.second,0,-1,1,false);)
void addHistoShapeFromFile(RooWorkspace* w, std::string name,std::string PF, const std::vector<std::string>& obs, const std::string& fileName, const std::string& hName, const InterpSysts& systs, const bool conditional = false, const int order = 0){
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
    auto h = TObjectHelper::getObject<TH1>(iF,hName);

    const std::string nominalHistName = name+"_Nominal_HIST"+PF;
    const std::string nominalPDFName  = systs.size() ? name+"_Nominal"+PF : PDFName;

    RooDataHist nHist(nominalHistName.c_str(),nominalHistName.c_str(),varlist,&*h);
    w->import(nHist);
    RooHistPdf nPDF(nominalPDFName.c_str(),nominalPDFName.c_str(),varset,nHist,order);
    w->import(nPDF);

    RooArgList coefList;
    RooArgList pdfList(*w->pdf(nominalPDFName.c_str()));
    for(const auto& s : systs ){

        std::string coefName = name+"_"+s.first+PF;
        RooArgList paramList;
        std::string systMathStr = "0";
        for(const auto& es : s.second) {
            systMathStr = systMathStr +"+"+ASTypes::flt2Str(es.second)+"*"+es.first;
            paramList.add(*w->var(es.first.c_str()));
        }
        RooFormulaVar varF(coefName.c_str(),coefName.c_str(),systMathStr.c_str(),paramList);
        w->import(varF);
        coefList.add(*w->arg(coefName.c_str()));

        auto addVar = [&](const std::string& var){
            auto sh =  TObjectHelper::getObject<TH1>(iF,hName+"_"+s.first+var);
            const std::string vHistName = name+"_" +s.first+var+"_HIST"+PF;
            const std::string vPDFName = name +"_"+s.first+var+PF;
            RooDataHist vRooHist(vHistName.c_str(),vHistName.c_str(),varlist,&*sh);
            w->import(vRooHist);
            RooHistPdf pdf(vPDFName.c_str(),vPDFName.c_str(),varset,vRooHist,order);
            w->import(pdf);
            pdfList.add(*w->pdf(vPDFName.c_str()));
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
}


void conditionalProduct(RooWorkspace* w, const std::string& pdfName, const std::string& pdfName1, const std::string& varName, const std::string& pdfName2){
    w->factory((std::string("PROD::")+pdfName+"("+pdfName1+"|"+varName+","+pdfName2+")").c_str());
}




};





#endif
