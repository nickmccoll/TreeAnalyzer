#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_DATACARDMAKER_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_DATACARDMAKER_H
#include "makeJSON.C"
#include "CutConstants.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "Utilities/HiggsCombineImport/interface/VerticalInterpHistPdf.h"
#include "Utilities/HiggsCombineImport/interface/HZZ2L2QRooPdfs.h"
#include "RooExponential.h"

#include <utility>



typedef std::pair<std::string,std::string> StrStr;
typedef std::vector<StrStr> StrStrs;

typedef std::pair<std::string,double> StrFlt;
typedef std::vector<StrFlt> StrFlts;

struct Systematic{
    //For parameterizaed uncertainties
    Systematic(const std::string& name,const double pV0,const double pV1):
        name(name),kind("param"),isParam(true),pV0(pV0),pV1(pV1){}
    //for others
    Systematic(const std::string& name,const std::string& kind,const StrFlts& values):
        name(name),kind(kind),isParam(false),pV0(0),pV1(0),values(values){}
    std::string name="";
    std::string kind="";
    bool isParam = false;
    double pV0 = 0;
    double pV1 = 0;
    StrFlts values;
};
struct Contribution{
    Contribution(const std::string& name,const std::string& pdf,const int id, const double yield):
        name(name),pdf(pdf),id(id),yield(yield){}
    std::string name="";
    std::string pdf="";
    int id=-1;
    double yield=1.0;
};

class DataCardMaker{
public:



    DataCardMaker(const std::string& finalstate,const std::string& category,const std::string& period, const double luminosity, const std::string& physics)
: finalstate(finalstate),category(category),period(period),luminosity(luminosity),physics(physics),lumiV(physics+"_"+period+"_lumi"),tag(physics +"_"+ finalstate+"_"+category+"_"+period)
{
        rootFile = new TFile((std::string("datacardInputs_")+tag+".root").c_str(),"RECREATE");
        rootFile->cd();
        w = new RooWorkspace("w","w");
        w->factory((lumiV+"["+ASTypes::flt2Str(luminosity)+"]").c_str());
}
    void addSystematic(const std::string& name,const std::string& kind, const StrFlts& values){
        systematics.emplace_back(name,kind,values);
    }
    void addParamSystematic(const std::string& name,const double pV0,const double pV1){
        systematics.emplace_back(name,pV0,pV1);
    }

    void addVar(const std::string& var, const bool constant){
        w->factory((var).c_str());
        w->var(var.c_str())->setConstant(constant);
    }
    void addVar(const std::string& var, const double val, const bool constant){
        w->factory((var + "["+ASTypes::flt2Str(val)+"]").c_str());
        w->var(var.c_str())->setConstant(constant);
    }
    void addVar(const std::string& var, const double val,const double min, const double max, const bool constant){
        w->factory((var + "["+ASTypes::flt2Str(val)+","+ASTypes::flt2Str(min)+","+ASTypes::flt2Str(max)+"]").c_str());
        w->var(var.c_str())->setConstant(constant);
    }

    void add1DBKGParametricShape(const std::string& name,const std::string& variableX, const std::string& jsonFile,  const StrFlts& scale_X,const StrFlts& resolution_X
            , const std::string& pVar="MR",const std::string& PDFPF="",const std::string& modelPF=""){
        CJSON json( jsonFile);
        const std::string varXPF =  name+"_"+PDFPF+(PDFPF.size() ? "_":"")+tag;
        const std::string PDFName = name+"_"+PDFPF+(PDFPF.size() ? "_":"")+tag;
        addCB(PDFName,pVar,variableX,varXPF,modelPF,json,scale_X,resolution_X);
    }

    void add2DSignalParametricShape(const std::string& name,const std::string& variableX,const std::string& variableY,const std::string& jsonFile,const StrFlts& scale_X,const StrFlts& resolution_X
            ,const StrFlts& scale_Y,const StrFlts& resolution_Y, bool exponential_X, const std::string& pVar="MH"){
        CJSON json( jsonFile);
        const std::string varXCBPDFName = exponential_X ? name+"_"+variableX+"_"+"CB"+"_"+tag : name+"_"+variableX+"_"+tag;
        const std::string varXEPDFName = name+"_"+variableX+"_"+"E"+"_"+tag;
        const std::string varXPF = name+"_"+variableX+"_"+tag;
        const std::string varYPF = name+"_"+variableY+"_"+tag;
        const std::string varXPDFName = name+"_"+variableX+"_"+tag;
        const std::string varYPDFName = name+"_"+variableY+"_"+tag;
        const std::string PDFName = name+"_"+tag;

        addCB(varXCBPDFName,pVar,variableX,varXPF,variableX,json,scale_X,resolution_X);
        if(exponential_X){
            std::string vN = std::string("fE") +"_"+varXPF;
            std::string expr = std::string("min(")+ pVar+"*0+"+json.getP(std::string("fE")+variableX)+",1)";
            RooFormulaVar varF(vN.c_str(),vN.c_str(),expr.c_str(),RooArgList(*w->var(pVar.c_str())));
            w->import(varF);
            addExpo(varXEPDFName,pVar,variableX,varXPF,variableX,json);

            RooAddPdf modelC(varXPDFName.c_str(),varXPDFName.c_str(),*w->pdf(varXEPDFName.c_str()),*w->pdf(varXCBPDFName.c_str()),*w->function(vN.c_str()));
            w->import(modelC);
            w->pdf(varXPDFName.c_str())->fixAddCoefRange("coef");
        }

        addCondCB(varYPDFName,pVar,variableY,varYPF,variableY,json,scale_Y,resolution_Y,variableX,varXPF);
        RooProdPdf condProdP(PDFName.c_str(), (varXPDFName+"*"+varYPDFName).c_str(),*w->pdf(varXPDFName.c_str()),RooFit::Conditional(*w->pdf(varYPDFName.c_str()), *w->var(variableY.c_str())) );
        w->import(condProdP);
        w->pdf(PDFName.c_str())->fixAddCoefRange("coef");
    }
    void addParametricYieldWithUncertainty(const std::string& name,const unsigned int ID,const std::string& jsonFile, const double constant, const std::string& uncName,const std::string& uncFormula,const double& uncValue,  const std::string& pVar="MS"){
        CJSON json( jsonFile);
        const std::string PDFName = name+"_"+tag;
        const std::string PDFNorm = PDFName+"_norm";
        std::string expr = std::string("(")+ json.getP("yield") +")*"+lumiV+"*("+ASTypes::flt2Str(constant)+"+"+uncName+"*"+uncFormula+")";
        addVar(uncName,0,-1,1,false);

        RooFormulaVar varF(PDFNorm.c_str(),PDFNorm.c_str(),expr.c_str(),RooArgList(*w->var(pVar.c_str()),*w->var(lumiV.c_str()),*w->var(uncName.c_str())));
        w->import(varF);
        addParamSystematic(uncName,0.0,uncValue);
        contributions.emplace_back(name,PDFName,ID,1.0);
    }

    void addHistoShapeFromFile(const std::string& name, const std::vector<std::string>& obs, const std::string& fileName, const std::string& hName, const StrStrs& systs, const bool conditional = false, const int order = 0,
            const std::string& PDFPF = "", const std::string& newTag = ""  ){
        RooArgSet varset;
        RooArgList varlist;
        std::vector<RooRealVar*> varPointers;

        for(const auto& v : obs){
            varPointers.push_back(w->var(v.c_str()));
            varset.add(*varPointers.back());
            varlist.add(*varPointers.back());
        }

        const std::string modelTag = newTag.size() ? newTag : name + "_"+tag;
        auto* iF =  TObjectHelper::getFile(fileName);
        auto h = TObjectHelper::getObject<TH1>(iF,hName);

        const std::string nominalHistName = name + (systs.size()?"NominalHIST_":"HIST_") +modelTag;
        const std::string nominalPDFName  = name + (systs.size()?"Nominal_":"_") +tag;

        RooDataHist nHist(nominalHistName.c_str(),nominalHistName.c_str(),varlist,&*h);
        w->import(nHist);
        RooHistPdf nPDF(nominalPDFName.c_str(),nominalPDFName.c_str(),varset,nHist,order);
        w->import(nPDF);

        RooArgList coefList;
        RooArgList pdfList(*w->pdf(nominalPDFName.c_str()));
        for(const auto& s : systs ){
            addVar(s.second,0,-1,1,false);
            coefList.add(*w->var(s.second.c_str()));
            auto sh = TObjectHelper::getObject<TH1>(iF,hName);

            auto addVar = [&](const std::string& var){
                auto sh =  TObjectHelper::getObject<TH1>(iF,hName+"_"+s.first+var);
                const std::string vHistName = name +"_"+s.first+var+"HIST_"+modelTag;
                const std::string vPDFName = name +"_"+s.first+var+"_"+tag;
                RooDataHist vRooHist(vHistName.c_str(),vHistName.c_str(),varlist,&*sh);
                w->import(vRooHist);
                RooHistPdf pdf(vPDFName.c_str(),vPDFName.c_str(),varset,vRooHist,order);
                w->import(pdf);
                pdfList.add(*w->pdf(vPDFName.c_str()));
            };
            addVar("Up");
            addVar("Down");

        }

        const std::string& PDFName = name+"_"+PDFPF+(PDFPF.size() ? "_":"")+tag;
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

    void addFixedYieldFromFile(const std::string& name, const int id, const std::string& fileName, const std::string& hName, const double constant = 1.0){
        auto* iF =  TObjectHelper::getFile(fileName);
        auto h = TObjectHelper::getObject<TH1>(iF,hName);
        const double nEvts = h->Integral()*luminosity*constant;
        contributions.emplace_back(name,name+"_"+tag,id,nEvts);
        iF->Close();
        delete iF;

    }
    void conditionalProduct(const std::string& name, const std::string& pdf1, const std::string& varName, const std::string& pdf2, const std::string& tag1="", const std::string& tag2=""){
        const std::string& pdfName = name +"_"+tag;
        const std::string& pdfName1 = pdf1 +"_"+ (tag1.size() ? tag1 : tag);
        const std::string& pdfName2 = pdf2 +"_"+ (tag2.size() ? tag2 : tag);
        w->factory((std::string("PROD::")+pdfName+"("+pdfName1+"|"+varName+","+pdfName2+")").c_str());
    }

    void importBinnedData(const std::string& fileName,const std::string& hName, const std::vector<std::string>& variables, const std::string& name="data_obs", const double scale = 1){
        auto* iF =  TObjectHelper::getFile(fileName);
        auto h = TObjectHelper::getObject<TH1>(iF,hName);
        if(scale != 1) h->Scale(scale);

        const unsigned int nD = variables.size();
        RooArgList args;
        auto doBinning =[&](const std::string& var, const TAxis* ax) {
            args.add(*w->var(var.c_str()));
            w->var(var.c_str())->setMin(ax->GetXmin());
            w->var(var.c_str())->setMax(ax->GetXmax());
            w->var(var.c_str())->setBins(ax->GetNbins());
            w->var(var.c_str())->setVal((ax->GetXmin()+ax->GetXmax())/2.);
        };
        doBinning(variables[0],h->GetXaxis());
        if(nD > 1)doBinning(variables[1],h->GetYaxis());
        if(nD > 2)doBinning(variables[2],h->GetZaxis());
        if(nD > 3)throw std::invalid_argument("makeCard::importBinnedData() -> Too many observables!");
        RooDataHist dataHist(name.c_str(),name.c_str(),args,&*h);
        w->import(dataHist);

        iF->Close();
        delete iF;
    }

    void makeCard(){
        std::ofstream f (std::string("datacard_")+tag+".txt",std::ios::out|std::ios::trunc);
        const std::string rfn = std::string("datacardInputs_")+tag+".root";
        const std::string dn = "data_obs";
        f<< "imax 1\n";
        f<<"jmax "<< (int(contributions.size()) -1)<<"\n";
        f<<"kmax *\n";
        f<<"-------------------------\n";
        for(const auto& c : contributions)
            f<<"shapes "<<c.name<<" "<<tag<<" "<<rfn<<" w:"<<c.pdf<<"\n";
        f<<"shapes "<<dn<<" "<<tag<<" "<<rfn<<" w:"<<dn<<"\n";
        f<<"-------------------------\n";
        f<<"bin "<<tag<<"\n";
        f<<"observation  -1\n";
        f<<"-------------------------\n";
        f<<"bin\t";

        for(const auto& c : contributions)
            f<<tag<<"\t";
        f<<"\n";

        //Sort the shapes by ID
        std::sort(contributions.begin(),contributions.end(), [](const Contribution& a,const Contribution& b){return a.id < b.id;});

        //print names
        f<<"process\t";
        for(const auto& c : contributions)
            f<<c.name<<"\t";
        f<<"\n";

        //Print ID
        f<<"process\t";
        for(const auto& c : contributions)
            f<< c.id<<"\t";
        f<<"\n";

        //print rates
        f<<"rate\t";
        for(const auto& c : contributions)
            f<< c.yield<<"\t";
        f<<"\n";

        //Now systematics
        auto getSContrib=[&](const Systematic& syst){
            for(const auto& c : contributions){
                bool has = false;
                for(const auto& sV : syst.values){
                    if(c.name == sV.first){
                        f<<sV.second<<"\t";
                        has = true;
                        break;
                    }
                }
                if(!has) f <<"-\t";
            }
        };


        for(const auto& s : systematics){
            f <<s.name<<"\t"<<s.kind;
            if(s.isParam)
                f<<"\t"<<s.pV0<<"\t"<<s.pV1;
            else if(s.kind == "discrete"){}
            else if(s.kind == "lnN"||s.kind == "lnU"){
                f<<"\t";
                getSContrib(s);
            }
            else throw std::invalid_argument( std::string("makeCard::makeCard() -> Bad Systematic kind: ")+ s.kind  );
            f<<"\n";
        }
        f.close();
        rootFile->cd();
        w->Write();
        rootFile->Close();
        delete rootFile;
    }


private:


    void addSystToWksp (const StrFlts& inSysts, const std::vector<std::string>& extraParams, std::string& systMathStr, RooArgList& paramList) {
        for(const auto& p : extraParams) paramList.add(*w->var(p.c_str()));
        systMathStr = "0";
        for(const auto& s : inSysts) {
            systMathStr = systMathStr +"+"+ASTypes::flt2Str(s.second)+"*"+s.first;
            paramList.add(*w->var(s.first.c_str()));
        }
    }

    void addParam (const std::string& vN,const std::string& expr, const RooArgList& pList, const std::string& systStr = ""){
        std::string sS ="";
        if(systStr.size()){sS = "*(1+"+systStr+")";}
        RooFormulaVar varF(vN.c_str(),vN.c_str(),(expr+sS).c_str(),pList);
        w->import(varF);
    }

    void addExpo(const std::string& pdfName,const std::string& pVar,const std::string& varName, const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json){
        auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
        auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
        auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};

        RooArgList paramList_std; paramList_std.add(*w->var(pVar.c_str()));
        addParam(pN("slope")    ,pJS("slope"),paramList_std);
        RooExponential modelE(pdfName.c_str(),pdfName.c_str(),*w->var(varName.c_str()), *pF("slope"));
        w->import(modelE);
    }

    void addCB(const std::string& pdfName,const std::string& pVar, const std::string& varName, const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json,const StrFlts& scale,const StrFlts& resolution){
        auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
        auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
        auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};


        RooArgList paramList_std; paramList_std.add(*w->var(pVar.c_str()));
        std::string systStr_scale;
        std::string systStr_res;
        RooArgList paramList_scale;
        RooArgList paramList_res;
        addSystToWksp(scale,{pVar}, systStr_scale,paramList_scale);
        addSystToWksp(resolution,{pVar},systStr_res,paramList_res);

        addParam(pN("mean")  ,pJS("mean"),paramList_scale,systStr_scale);
        addParam(pN("sigma") ,pJS("sigma"),paramList_res,systStr_res);
        addParam(pN("alpha") ,pJS("alpha"),paramList_std);
        addParam(pN("n")     ,pJS("n"),paramList_std);
        addParam(pN("alpha2"),pJS("alpha2"),paramList_std);
        addParam(pN("n2")    ,pJS("n2"),paramList_std);

        RooDoubleCB modelP(pdfName.c_str(),pdfName.c_str(),*w->var(varName.c_str()),
                *pF("mean"),*pF("sigma"),*pF("alpha"),*pF("n") ,*pF("alpha2"),*pF("n2"));
        w->import(modelP);
    }


    void addCondCB(const std::string& pdfName,const std::string& pVar, const std::string& varName,  const std::string& parPostFix,const std::string& JSONPostFix, const CJSON& json,const StrFlts& scale,const StrFlts& resolution,
            const std::string& condVarName, const std::string& condParPostFix
    ){
        auto cPN = [&](const std::string& v)->std::string{return v +"_"+condParPostFix;};
        auto cPF = [&](const std::string& v)->RooAbsReal*{return w->function(cPN(v).c_str());};

        auto pN = [&](const std::string& v)->std::string{return v +"_"+parPostFix;};
        auto pJS = [&](const std::string& v)->std::string{return std::string("(0*")+pVar+"+"  +json.getP(v+JSONPostFix)+")";};
        auto pF = [&](const std::string& v)->RooAbsReal*{return w->function(pN(v).c_str());};


        RooArgList paramList_std; paramList_std.add(*w->var(pVar.c_str()));
        std::string systStr_scale;
        std::string systStr_res;
        RooArgList paramList_scale;
        RooArgList paramList_res;
        addSystToWksp(scale,{pVar}, systStr_scale,paramList_scale);
        addSystToWksp(resolution,{pVar},systStr_res,paramList_res);

        addParam(pN("mean")  ,pJS("mean"),paramList_scale,systStr_scale);
        addParam(pN("sigma") ,pJS("sigma"),paramList_res,systStr_res);
        addParam(pN("alpha") ,pJS("alpha"),paramList_std);
        addParam(pN("n")     ,pJS("n"),paramList_std);
        addParam(pN("alpha2"),pJS("alpha2"),paramList_std);
        addParam(pN("n2")    ,pJS("n2"),paramList_std);

        addParam(pN("maxS")     ,pJS("maxS"),paramList_std);
        addParam(pN("mean_p1")  ,pJS("mean_p1"),paramList_std);
        addParam(pN("sigma_p1") ,pJS("sigma_p1"),paramList_std);

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

    const std::string finalstate;
    const std::string category;
    const std::string period;
    const double luminosity;
    const std::string physics;
    const std::string lumiV;
    const std::string tag;
    const double sqrt_s= 13000;


    TFile * rootFile=0;
public:
    RooWorkspace* w=0;
    std::vector<Contribution> contributions;
    std::vector<Systematic> systematics;
};
#endif
