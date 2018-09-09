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
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "RooExponential.h"
#include "PDFAdder.h"

#include <utility>


using PDFAdder::StrFlt;
using PDFAdder::StrFlts;


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
    void rebinX(int nBins, double min, double max) {rebinner.reset(new PDFAdder::HistRebinX(nBins,min,max)); }
    void rebinX(const  std::vector<double>& bins){rebinner.reset(new PDFAdder::HistRebinX(bins)); }
    void rebinY(int nBins, double min, double max) {rebinner.reset(new PDFAdder::HistRebinY(nBins,min,max)); }
    void rebinY(const  std::vector<double>& bins){rebinner.reset(new PDFAdder::HistRebinY(bins)); }
    void rebinXY(int nBinsX, double minX, double maxX,int nBinsY, double minY, double maxY)
    {rebinner.reset(new PDFAdder::HistRebinXY(nBinsX,minX,maxX,nBinsY,minY,maxY)); }
    void rebinXY(const  std::vector<double>& binsX,int nBinsY, double minY, double maxY)
    {rebinner.reset(new PDFAdder::HistRebinXY(binsX,nBinsY,minY,maxY)); }
    void rebinXY(int nBinsX, double minX, double maxX, std::vector<double>& binsY)
    {rebinner.reset(new PDFAdder::HistRebinXY(nBinsX,minX,maxX,binsY)); }
    void rebinXY(std::vector<double>& binsX, std::vector<double> binsY)
    {rebinner.reset(new PDFAdder::HistRebinXY(binsX,binsY)); }

    void addSystematic(const std::string& name,const std::string& kind, const StrFlts& values){
        systematics.emplace_back(name,kind,values);
    }
    void addParamSystematic(const std::string& name,const double pV0,const double pV1){
        addVar(name,pV0,pV0-pV1*4,pV0+pV1*4,false);
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
        PDFAdder::addCB(w,PDFName,pVar,variableX,varXPF,modelPF,json,scale_X,resolution_X);
    }


    void add2DSignalParametricShape(const std::string& name,const std::string& variableX,const std::string& variableY,const std::string& jsonFile,const StrFlts& scale_X,const StrFlts& resolution_X
            ,const StrFlts& scale_Y,const StrFlts& resolution_Y, bool exponential_X, const std::string& pVar="MH"){
        PDFAdder::add2DCB(w,name,tag,variableX,variableY,jsonFile,
                scale_X,resolution_X,scale_Y,resolution_Y,exponential_X,pVar);
    }

    void add2DSignalParametricShapeNoCond(const std::string& name,const std::string& variableX,const std::string& variableY,const std::string& jsonFile,const StrFlts& scale_X,const StrFlts& resolution_X
            ,const StrFlts& scale_Y,const StrFlts& resolution_Y, bool exponential_X, const std::string& pVar="MH"){
        PDFAdder::add2DCBNoCond(w,name,tag,variableX,variableY,jsonFile,
                scale_X,resolution_X,scale_Y,resolution_Y,exponential_X,pVar);
    }

    void addParametricYieldWithUncertainty(const std::string& name,const unsigned int ID,const std::string& jsonFile, const double constant, const std::string& uncFormula, const  std::vector<std::string>& uncPList, const std::string& pVar="MS"){
        CJSON json( jsonFile);
        const std::string PDFName = name+"_"+tag;
        const std::string PDFNorm = PDFName+"_norm";
        std::string expr = std::string("(")+ json.getP("yield") +")*"+lumiV+"*"+ASTypes::flt2Str(constant);
        RooArgList args(*w->var(pVar.c_str()),*w->var(lumiV.c_str()));
        if(uncFormula.size()){
            expr += "*(" + uncFormula+")";
            for(const auto& p: uncPList) args.add(*w->var(p.c_str()));
        }
        RooFormulaVar varF(PDFNorm.c_str(),PDFNorm.c_str(),expr.c_str(),args);
        w->import(varF);
        contributions.emplace_back(name,PDFName,ID,1.0);
    }


    void addHistoShapeFromFile(const std::string& name,const std::vector<std::string>& obs, const std::string& fileName,const std::string& hName,
            const PDFAdder::InterpSysts& systs, const bool conditional = false, const int order = 0,const std::string& PDFPF = "",bool isY = false){
        std::string PF = PDFPF+(PDFPF.size() ? "_":"")+tag;
        PDFAdder::addHistoShapeFromFile(w,name,PF,obs,fileName,hName,systs,conditional,order,&*rebinner,isY);
    }

    double addFixedYieldFromFile(const std::string& name, const int id, const std::string& fileName, const std::string& hName, const double constant = 1.0, bool isY = false){
        auto* iF =  TObjectHelper::getFile(fileName);
        auto h = TObjectHelper::getObject<TH1>(iF,hName);
        double nEvents =luminosity*constant;
        if(rebinner){
            auto * nh = rebinner->process(&*h,std::string(h->GetName())+"_rebin",isY);
            nEvents *= nh->Integral();
            delete nh;
        }
        else nEvents *= h->Integral();
        contributions.emplace_back(name,name+"_"+tag,id,nEvents);
        iF->Close();
        delete iF;
        return nEvents;

    }
    void conditionalProduct(const std::string& name, const std::string& pdfName_cxoy, const std::string& varName_x, const std::string& pdfName_y, const std::string& tag_pdf_cxoy="", const std::string& tag_pdf_y=""){
        const std::string& pdfName = name +"_"+tag;
        const std::string& pdfName1 = pdfName_cxoy +"_"+ (tag_pdf_cxoy.size() ? tag_pdf_cxoy : tag);
        const std::string& pdfName2 = pdfName_y +"_"+ (tag_pdf_y.size() ? tag_pdf_y : tag);
        PDFAdder::conditionalProduct(w,pdfName,pdfName1,varName_x,pdfName2);
    }

    void importBinnedData(const std::string& fileName,const std::string& hName, const std::vector<std::string>& variables, const std::string& name="data_obs", const double scale = 1, bool isY=false){
        auto* iF =  TObjectHelper::getFile(fileName);
        auto inh = TObjectHelper::getObject<TH1>(iF,hName);
        TH1 * h = &*inh;
        if(rebinner){
            h = rebinner->process(&*h,std::string(h->GetName())+"_rebin",isY);
        }
        if(scale != 1) h->Scale(scale);
        PDFAdder::addBinnedData(w,h,name,variables);
        iF->Close();
        delete iF;
        if(rebinner) delete h;
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
    const std::string finalstate;
    const std::string category;
    const std::string period;
    const double luminosity;
    const std::string physics;
    const std::string lumiV;
    const std::string tag;
    const double sqrt_s= 13000;
    TFile * rootFile=0;
    std::unique_ptr<PDFAdder::HistRebin> rebinner;
public:
    RooWorkspace* w=0;
    std::vector<Contribution> contributions;
    std::vector<Systematic> systematics;
};
#endif
