
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/KDEProducer.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "thread"
#include <string.h>


using namespace TAna;

class make1DTemplateWithAdaKernAnalyzer : public BaseTreeAnalyzer {
public:

    make1DTemplateWithAdaKernAnalyzer(std::string fileName, std::string treeName,std::string arguments ) : BaseTreeAnalyzer(fileName,treeName,2)
{


        ParParser p;
        name     = p.addString("n","histogram name",true);
        auto v   = p.addString("x","variable",true);
        auto g   = p.addString("g","gen variable",true);
        auto b   = p.addVFloat("xb","x-variable binning",true);
        auto s   = p.addString("s","selection",false,"1.0");
        auto w   = p.addString("w","weight",false,"1.0");
        auto sw  = p.addBool("sw","add systematic weight");
        auto swUp= p.addString("swUp","systematic weight up",false,"1.0");
        auto swDown= p.addString("swDown","systematic weight down",false,"1.0");
        kt       = p.addString("t","Types of templates to calculate",false, "nrRsSo");
        ks       = p.addFloat("ks","KDE scaling: scale");
        kr       = p.addFloat("kr","KDE scaling: resolution");
        khs      = p.addFloat("khs","KDE h-scale factor",false,1.);
        khc      = p.addFloat("khc","KDE adaptive bandwidth cutoff",false,5);
        kss      = p.addBool ("kss","KDE sigma scaling");


        hs       = p.addFloat("hs","Histogram scaling: x proportional scale",true);
        hr       = p.addFloat("hr","Histogram scaling: 1/x proportional scale",true);

        doS      = p.addBool("doS","Apply exponential smoothing");
        emin     = p.addFloat ("emin","Exponential fit min");
        emax     = p.addFloat ("emax","Exponential fit max");

        auto asf       = p.addString("vsf","Average scale file",true);
        auto ash       = p.addString("vsh","Average scale hist",true);
        auto asv       = p.addString("vsv","Average scale variable",true);
        p.parse(arguments);

        if(b->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        vAxis.reset(new TAxis((*b)[0],(*b)[1],(*b)[2]));

        tree.getTree()->SetBranchStatus("*",1);
        sForm.reset(new TTreeFormula("sForm", TString::Format("%s*(%s)",w->c_str(),s->c_str()),tree.getTree()));
        if(*sw){
            sFormUp.reset(new TTreeFormula("sFormUp", TString::Format("%s*(%s)",swUp->c_str(),s->c_str()),tree.getTree()));
            sFormDown.reset(new TTreeFormula("sFormDown", TString::Format("%s*(%s)",swDown->c_str(),s->c_str()),tree.getTree()));
        }
        vForm.reset(new TTreeFormula("vForm", v->c_str(),tree.getTree()));

        const int nEntries =  tree.getTree()->GetEntries(TString::Format("%s*(%s)",w->c_str(),s->c_str()));

        nominalX .reset(new std::vector<double>);
        upSX     .reset(new std::vector<double>);
        downSX   .reset(new std::vector<double>);
        weight   .reset(new std::vector<double>);
        if(*sw){
        systWeightUp   .reset(new std::vector<double>);
        systWeightDown   .reset(new std::vector<double>);
        }

        nominalX ->reserve(nEntries);
        upSX     ->reserve(nEntries);
        downSX   ->reserve(nEntries);
        weight   ->reserve(nEntries);
        if(*sw){
        systWeightUp   ->reserve(nEntries);
        systWeightDown   ->reserve(nEntries);
        }

        needScaleFile = (kt->find('R') != std::string::npos)||(kt->find('r') != std::string::npos) ;
        if( needScaleFile ){
            auto file =  TObjectHelper::getFile(*asf);
            avgScale.reset(new TObjectHelper::Hist1DContainer(file,*ash));
            delete file;

            upRX     .reset(new std::vector<double>);
            downRX   .reset(new std::vector<double>);
            upRX     ->reserve(nEntries);
            downRX   ->reserve(nEntries);

            gForm.reset(new TTreeFormula("gForm", g->c_str(),tree.getTree()));
            vsForm.reset(new TTreeFormula("vsForm", asv->c_str(),tree.getTree()));
        }
}


    void loadVariables() override{};
    bool runEvent() override {
        double s   = sForm->EvalInstance();

        if(s==0) return false;

        double v   = vForm->EvalInstance();

        double upScale = *ks*v;
        double downScale = (*ks-1.)*v;

        nominalX ->push_back(v);
        upSX     ->push_back(upScale);
        downSX   ->push_back(downScale);
        weight   ->push_back(s);
        if(systWeightUp){
        double sUp   = sFormUp->EvalInstance();
        double sDown   = sFormDown->EvalInstance();
        systWeightUp   ->push_back(sUp);
        systWeightDown   ->push_back(sDown);
        }

        if(needScaleFile){
            double g   = gForm->EvalInstance();
            double sv   = vsForm->EvalInstance();

            double avgS = avgScale->getBinContentByValue(sv).val();
            double upRes= std::max(0.0,avgS*g + *kr* (v -avgS*g));
            double downRes = std::max(0.0,avgS*g + (*kr-1.)*(v -avgS*g));

            upRX     ->push_back(upRes);
            downRX   ->push_back(downRes);
        }

        return true;
    }

    const TH1* makeKDE(std::string name, const std::vector<double>& xvals, const std::vector<double>& weights){
        const int   nBsX   = vAxis->GetNbins();
        const float minX   = vAxis->GetXmin();
        const float maxX   = vAxis->GetXmax();

        KDEProducer pdfProd(&xvals,&weights,*khs,nBsX,minX,maxX,*khc,*kss);

        plotter.add1D(pdfProd.getPDF(name+"_debug_KDE0","",nBsX,minX,maxX));

        TH1 * dataH = new TH1F((name+"_data").c_str(),"",nBsX,minX,maxX);
        for(unsigned int iP = 0; iP < xvals.size(); ++iP){
            dataH->Fill((xvals)[iP],(weights)[iP]);
        }
        plotter.add1D(dataH);

        auto * pilot = pdfProd.getPilotPDF();
        pilot->SetName((name + "_debug_pilotKDE").c_str());
        plotter.add1D(pilot);

        plotter.add1D(pdfProd.getABandwidths(name+"_debug_bandwidths","",nBsX,minX,maxX));
        plotter.add1D(pdfProd.getLocalVariance(name   +"_debug_var","",nBsX,minX,maxX)              ) ;

        TH1 * kde = pdfProd.getAPDF(name+"_KDE","",nBsX,minX,maxX);
        kde->Scale(dataH->Integral()/kde->Integral());

        plotter.add1D(kde);
        return kde;
    }

    const TH1 * cloneAndWrite(const std::string& name, const TH1* iHist){
        TH1 * oHist = (TH2*)iHist->Clone(name.c_str());
        oHist->Scale(1.0/oHist->Integral());
        plotter.add1D(oHist);
        return oHist;
    }

    const TH1 * smoothTail(const std::string& name, const TH1* iHist){
        TH1 * oHistD = (TH1*)iHist->Clone((name +"_debugFitKDE") .c_str());
        TH1 * oHist = (TH1*)iHist->Clone((name) .c_str());
        TF1 expo("expo","expo",*emin,*emax);
        oHistD->Fit(&expo,"","",*emin,*emax);
        for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
            const double x = oHist->GetXaxis()->GetBinCenter(iX);
            if(x > *emin+300){
                double fv = expo.Eval(x);
                if(x < *emin +700){
                    const double kFr = ((*emin +700) - x)/400;
                    fv = (1-kFr)*fv + kFr* oHist->GetBinContent(iX);
                }
                oHist->SetBinContent(iX,fv);
            }

        }
        oHist->Scale(1.0/oHist->Integral());
        plotter.add1D(oHist);
        return oHist;
    }


    const TH1* transform ( std::string name,const TH1 * iHist, std::function<double(double)> f) {
        TH1 * outH = (TH1*)iHist->Clone(name.c_str());
        for(int iB = 1; iB <= outH->GetNbinsX(); ++iB){
            outH->SetBinContent(iB,outH->GetBinContent(iB)*f(outH->GetBinCenter(iB)) );
        }
        outH->Scale(1.0/outH->Integral());
        plotter.add1D(outH);
        return outH;
    };


    void process(std::string outFileName) {

        if (kt->find('n') != std::string::npos) {
            auto kde = makeKDE(*name,*nominalX,*weight);
            kde = *doS ? smoothTail(*name,kde) : cloneAndWrite(*name,kde);
            transform(*name+"_PTDown" ,kde, [&](double x){return  1./(1. + *hs*x);})  ;
            transform(*name+"_PTUp"   ,kde, [&](double x){return  (1. + *hs*x);})  ;

            transform(*name+"_PT2Down" ,kde, [&](double x){return  1./(1. + *hs*x*x);})  ;
            transform(*name+"_PT2Up"   ,kde, [&](double x){return  (1. + *hs*x*x);})  ;

            transform(*name+"_OPTDown",kde, [&](double x){return  1./(1. + *hr/x);})  ;
            transform(*name+"_OPTUp"  ,kde, [&](double x){return  (1. + *hr/x);})  ;

            if(systWeightUp){
            auto kdeUp = makeKDE(*name+"_WeightUp",*nominalX,*systWeightUp);
            kde = *doS ? smoothTail(*name+"_WeightUp",kdeUp) : cloneAndWrite(*name+"_WeightUp",kdeUp);
            auto kdeDown = makeKDE(*name+"_WeightDown",*nominalX,*systWeightDown);
            kde = *doS ? smoothTail(*name+"_WeightDown",kdeDown) : cloneAndWrite(*name+"_WeightDown",kdeDown);
            }
        }

        auto justDoNominal = [&](std::string name, const std::vector<double>& xvals ){
            auto kde = makeKDE(name,xvals,*weight);
            *doS ? smoothTail(name,kde) : cloneAndWrite(name,kde);
        };

        if (kt->find('S') != std::string::npos) justDoNominal(*name+"_ScaleUp"  ,*upSX  );
        if (kt->find('s') != std::string::npos) justDoNominal(*name+"_ScaleDown",*downSX);
        if (kt->find('R') != std::string::npos) justDoNominal(*name+"_ResUp"    ,*upRX  );
        if (kt->find('r') != std::string::npos) justDoNominal(*name+"_ResDown"  ,*downRX);

        plotter.write(outFileName);
    }

    std::unique_ptr<TAxis> vAxis;
    std::unique_ptr<TTreeFormula> sForm;
    std::unique_ptr<TTreeFormula> sFormUp;
    std::unique_ptr<TTreeFormula> sFormDown;
    std::unique_ptr<TTreeFormula> vForm;
    std::unique_ptr<TTreeFormula> gForm;
    std::unique_ptr<TTreeFormula> vsForm;
    std::unique_ptr<TObjectHelper::Hist1DContainer>               avgScale;

    std::shared_ptr<std::string> name;
    std::shared_ptr<std::string>  kt;
    std::shared_ptr<double>       kr;
    std::shared_ptr<double>       ks;
    std::shared_ptr<double>       khs;
    std::shared_ptr<double>       khc;
    std::shared_ptr<bool>         kss;
    std::shared_ptr<double>       hs;
    std::shared_ptr<double>       hr;
    std::shared_ptr<bool>         doS;
    std::shared_ptr<double>       emin;
    std::shared_ptr<double>       emax;


    std::unique_ptr<std::vector<double>> nominalX  ;
    std::unique_ptr<std::vector<double>> upSX      ;
    std::unique_ptr<std::vector<double>> downSX    ;
    std::unique_ptr<std::vector<double>> upRX      ;
    std::unique_ptr<std::vector<double>> downRX    ;
    std::unique_ptr<std::vector<double>> weight    ;
    std::unique_ptr<std::vector<double>> systWeightUp;
    std::unique_ptr<std::vector<double>> systWeightDown;

    bool needScaleFile = false;
    HistGetter plotter;



};

#endif

void make1DTemplateWithAdaKern(std::string fileName, std::string outFileName,std::string arguments){
    make1DTemplateWithAdaKernAnalyzer a(fileName,"treeMaker/Events",arguments);
    a.analyze(100000);
    a.process(outFileName);
}
