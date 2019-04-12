\
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/KDEProducer2D.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "thread"
#include <string.h>


using namespace TAna;

class make2DTemplateWithAdaKernAnalyzer : public BaseTreeAnalyzer {
public:

    make2DTemplateWithAdaKernAnalyzer(std::string fileName, std::string treeName,
            std::string arguments ) : BaseTreeAnalyzer(fileName,treeName,2)
    {


        ParParser p;
        name      = p.addString("n","histogram name",true);
        xIsCond   = p.addBool("xIsCond",
                "If true, set x as the conditional variable P(X|Y)*P(Y), if false: P(Y|X)*PY(Y)");
        auto vx   = p.addString("vx","x variable",true);
        auto vy   = p.addString("vy","y variable",true);
        auto xb   = p.addVFloat("xb","x-variable binning",true);
        auto yb   = p.addVFloat("yb","y-variable binning",true);

        auto s   = p.addString("s","selection",false,"1.0");
        auto w   = p.addString("w","weight",false,"1.0");
        khxs     = p.addFloat("khxs","KDE h-x scale factor",false,1.);
        khxc     = p.addFloat ("khxc","KDE adaptive x-bandwidth cutoff",false,5);
        khys     = p.addFloat("khys","KDE h-y scale factor",false,1.);
        khyc     = p.addFloat ("khyc","KDE adaptive y-bandwidth cutoff",false,5);
//        kss      = p.addBool ("kss","KDE sigma scaling");

        eOpt     = p.addString("eopt",
                "Add an exponential fit, set to x if you want to fit x, y otherwise",false,"");
        emin     = p.addFloat ("emin","Exponential fit min");
        emax     = p.addFloat ("emax","Exponential fit max");

        hs       = p.addFloat("hs","Conditional variable hist scaling: x proportional",true);
        hr       = p.addFloat("hr","Conditional variable hist scaling: 1/x proportional",true);

        p.parse(arguments);

        if(xb->size() != 3) throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");
        xAxis.reset(new TAxis((*xb)[0],(*xb)[1],(*xb)[2]));
        if(yb->size() != 3){
            yAxis.reset(new TAxis((*yb).size() -1, &(*yb)[0]));
        } else {
            yAxis.reset(new TAxis((*yb)[0],(*yb)[1],(*yb)[2]));
        }

        tree.getTree()->SetBranchStatus("*",1);
        sForm.reset(new TTreeFormula("sForm",
                TString::Format("%s*(%s)",w->c_str(),s->c_str()),tree.getTree()));
        vxForm.reset(new TTreeFormula("vxForm", vx->c_str(),tree.getTree()));
        vyForm.reset(new TTreeFormula("vyForm", vy->c_str(),tree.getTree()));

        const int nEntries =  tree.getTree()->GetEntries(
                TString::Format("%s*(%s)",w->c_str(),s->c_str()));

        nominalX .reset(new std::vector<double>);
        nominalY .reset(new std::vector<double>);
        weight   .reset(new std::vector<double>);

        nominalX ->reserve(nEntries);
        nominalY ->reserve(nEntries);
        weight   ->reserve(nEntries);
    }


    void loadVariables() override{};
    bool runEvent() override {
        double s   = sForm->EvalInstance();
        if(s==0) return false;

        double vx   = vxForm->EvalInstance();
        double vy   = vyForm->EvalInstance();
        nominalX ->push_back(vx);
        nominalY ->push_back(vy);
        weight   ->push_back(s);
        return true;
    }

    //output the coarse KDE used for later steps
    const TH2* makeKDE(std::string name){
        const int   nBsX   = xAxis->GetNbins();
        const float minX   = xAxis->GetXmin();
        const float maxX   = xAxis->GetXmax();
        const int   nBsY   = yAxis->GetNbins();
        const float minY   = yAxis->GetXmin();
        const float maxY   = yAxis->GetXmax();
        const bool vBinsY   = yAxis->GetXbins()->GetSize();
        const double* fBinY = vBinsY ? yAxis->GetXbins()->GetArray() : 0;

        KDEProducer2D * pdfProd = vBinsY ?
                new KDEProducer2D(nominalX.get(),nominalY.get(),weight.get(),
                        *khxs,nBsX,minX,maxX,*khxc,
                        *khys,nBsY,fBinY,*khyc) :
                new KDEProducer2D(nominalX.get(),nominalY.get(),weight.get(),
                        *khxs,nBsX,minX,maxX,*khxc,
                        *khys,nBsY,minY,maxY,*khyc) ;

//        plotter.add2D(pdfProd.getPDF(name+"_debug_KDE0","",nBsX,minX,maxX,nBsY,minY,maxY));

        TH2 * data = vBinsY ?
                new TH2F((name+"_fine_data").c_str(),";data",nBsX,minX,maxX,nBsY,fBinY)
                : new TH2F((name+"_fine_data").c_str(),";data",nBsX,minX,maxX,nBsY,minY,maxY);
        for(unsigned int iP = 0; iP < nominalX->size(); ++iP){
            data->Fill((*nominalX)[iP],(*nominalY)[iP],(*weight)[iP]);
        }

        plotter.add2D(data);
//        plotter.add2D(pdfProd.getAPDF(name+"_debug_fineKDE","",nBsX,minX,maxX,nBsY,minY,maxY));

//        auto * pilot = pdfProd.getPilotPDF();
//        pilot->SetName((name + "_debug_pilotKDE").c_str());
//        plotter.add2D(pilot);

        if(vBinsY){
            plotter.add2D(pdfProd->getABandwidths(name+"_debug_bandwidthsX","",
                    nBsX,minX,maxX,nBsY,fBinY,true));
            plotter.add2D(pdfProd->getABandwidths(name+"_debug_bandwidthsY","",
                    nBsX,minX,maxX,nBsY,fBinY,false));
        } else {
            plotter.add2D(pdfProd->getABandwidths(name+"_debug_bandwidthsX","",
                    nBsX,minX,maxX,nBsY,minY,maxY,true));
            plotter.add2D(pdfProd->getABandwidths(name+"_debug_bandwidthsY","",
                    nBsX,minX,maxX,nBsY,minY,maxY,false));
        }
//        plotter.add2D(pdfProd.getLocalVarX(name   +"_debug_varX","",nBsX,minX,maxX,nBsY,minY,maxY)              ) ;
//        plotter.add2D(pdfProd.getLocalVarY(name   +"_debug_varY","",nBsX,minX,maxX,nBsY,minY,maxY)              ) ;

        TH2 * kde = vBinsY ? pdfProd->getAPDF(name+"_KDE","",nBsX,minX,maxX,nBsY,fBinY) :
                pdfProd->getAPDF(name+"_KDE","",nBsX,minX,maxX,nBsY,minY,maxY);
        kde->Scale(1.0/kde->Integral());
        plotter.add2D(kde);
        delete pdfProd;
        return kde;
    }

    const TH2 * smoothXTail(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        oHist->Scale(1.0/oHist->Integral());

        for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
            auto * proj = oHist->ProjectionX("q",iY,iY);
            TF1 expo("expo","expo",*emin,*emax);
            proj->Fit(&expo,"","",*emin,*emax);
            for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
                const double x = oHist->GetXaxis()->GetBinCenter(iX);
                if(x > *emin+300){
                    double fv = expo.Eval(x);
                    if(x < *emin +700){
                        const double kFr = ((*emin +700) - x)/400;
                        fv = (1-kFr)*fv + kFr* oHist->GetBinContent(iX,iY);
                    }
                    oHist->SetBinContent(iX,iY,fv);
                }

            }
        }
        oHist->Scale(1.0/oHist->Integral());
        plotter.add2D(oHist);
        return oHist;

    }
    const TH2 * smoothYTail(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
         oHist->Scale(1.0/oHist->Integral());

        for(int iX = 1; iX <= oHist->GetNbinsX(); ++iX){
            auto * proj = oHist->ProjectionY("q",iX,iX);
            TF1 expo("expo","expo",*emin,*emax);
            proj->Fit(&expo,"","",*emin,*emax);
            for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                const double y = oHist->GetYaxis()->GetBinCenter(iY);
                if(y > *emin+300){
                    double fv = expo.Eval(y);
                    if(y < *emin +700){
                        const double kFr = ((*emin +700) - y)/400;
                        fv = (1-kFr)*fv + kFr* oHist->GetBinContent(iX,iY);
                    }
                    oHist->SetBinContent(iX,iY,fv);
                }

            }
        }
        oHist->Scale(1.0/oHist->Integral());
        plotter.add2D(oHist);
        return oHist;
    }
    const TH2 * cloneAndWrite(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        for(int iX = 1; iX <= iHist->GetNbinsX(); ++iX)
            for(int iY = 1; iY <= iHist->GetNbinsY(); ++iY){
                const double binW = iHist->GetXaxis()->GetBinWidth(iX)
                        *iHist->GetYaxis()->GetBinWidth(iY);
                oHist->SetBinContent(iX,iY,binW*oHist->GetBinContent(iX,iY));
                oHist->SetBinError(iX,iY,binW*iHist->GetBinError(iX,iY));
            }
        oHist->Scale(1.0/oHist->Integral());
        plotter.add2D(oHist);
        return oHist;
    }
    const TH2 * XConditionalOnY(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
            auto * proj = oHist->ProjectionX("q",iY,iY);
            double xIntegral = oHist->Integral(1,oHist->GetNbinsX(),iY,iY);
            if(xIntegral)
            for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
                oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)/xIntegral);
            }
        }
        plotter.add2D(oHist);
        return oHist;
    }

    const TH2 * YConditionalOnX(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        for(int iX = 1; iX <= oHist->GetNbinsX(); ++iX){
            auto * proj = oHist->ProjectionY("q",iX,iX);
            double yIntegral = oHist->Integral(iX,iX,1,oHist->GetNbinsY());
            if(yIntegral)
            for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)/yIntegral);
            }
        }
        plotter.add2D(oHist);
        return oHist;
    }

    const TH2 * transformX ( std::string name,const TH1 * iHist, std::function<double(double)> f ) {
      TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
      for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
          const double s = f(oHist->GetXaxis()->GetBinCenter(iX));
          for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
              oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*s);
          }
      }
      plotter.add2D(oHist);
       return oHist;
    };

    const TH2 * transformY ( std::string name,const TH1 * iHist, std::function<double(double)> f ) {
      TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
      for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
          const double s = f(oHist->GetYaxis()->GetBinCenter(iY));
          for(int iX = 1; iX <= oHist->GetNbinsX(); ++iX){
              oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*s);
          }
      }
      plotter.add2D(oHist);
       return oHist;
    };

    void process(std::string outFileName) {
        auto * kde = makeKDE(*name);
        auto * nomH = cloneAndWrite(*name+"_noCond",kde);

        if(ASTypes::strFind(*eOpt,"x")){
            kde = smoothXTail(*name+ "_smoothKDE", kde);
            nomH = cloneAndWrite(*name+"_smoothNoCond",kde);
        }
        else if (ASTypes::strFind(*eOpt,"y")){
            kde = smoothYTail(*name+ "_smoothKDE", kde);
            nomH = cloneAndWrite(*name+"_smoothNoCond",kde);
        }


        if(*xIsCond){
            XConditionalOnY(*name,nomH);

            auto * upScale = transformX(*name+"_PTUp_noCond",nomH,
                    [&](double x){return  (1. + *hs*x);});
            XConditionalOnY(*name+"_PTUp",upScale);
            auto * downScale = transformX(*name+"_PTDown_noCond",nomH,
                    [&](double x){return  1./(1. + *hs*x);});
            XConditionalOnY(*name+"_PTDown",downScale);

            auto * upRes = transformX(*name+"_OPTUp_noCond",nomH,
                    [&](double x){return  (1. + *hr/x);});
            XConditionalOnY(*name+"_OPTUp",upRes);
            auto * downRes = transformX(*name+"_OPTDown_noCond",nomH,
                    [&](double x){return 1./(1. + *hr/x);});
            XConditionalOnY(*name+"_OPTDown",downRes);


        } else {
            YConditionalOnX(*name,nomH);

            auto * upScale = transformY(*name+"_PTUp_debug_noCond",nomH,
                    [&](double x){return  (1. + *hs*x);});
            YConditionalOnX(*name+"_PTUp",upScale);
            auto * downScale = transformY(*name+"_PTDown_debug_noCond",nomH,
                    [&](double x){return  1./(1. + *hs*x);});
            YConditionalOnX(*name+"_PTDown",downScale);

            auto * upRes = transformY(*name+"_OPTUp_debug_noCond",nomH,
                    [&](double x){return  (1. + *hr/x);});
            YConditionalOnX(*name+"_OPTUp",upRes);
            auto * downRes = transformY(*name+"_OPTDown_debug_noCond",nomH,
                    [&](double x){return 1./(1. + *hr/x);});
            YConditionalOnX(*name+"_OPTDown",downRes);
        }

        plotter.write(outFileName);
    }
    std::shared_ptr<std::string> name;
    std::shared_ptr<bool>        xIsCond;

    std::unique_ptr<TAxis> xAxis;
    std::unique_ptr<TAxis> yAxis;
    std::unique_ptr<TTreeFormula> sForm;
    std::unique_ptr<TTreeFormula> vxForm;
    std::unique_ptr<TTreeFormula> vyForm;

    std::shared_ptr<double>       khxs;
    std::shared_ptr<double>       khxc;
    std::shared_ptr<double>       khys;
    std::shared_ptr<double>       khyc;
//    std::shared_ptr<bool>         kss;
    std::shared_ptr<double>       hs;
    std::shared_ptr<double>       hr;

    std::shared_ptr<std::string>  eOpt;
    std::shared_ptr<double>       emin;
    std::shared_ptr<double>       emax;

    std::unique_ptr<std::vector<double>> nominalX  ;
    std::unique_ptr<std::vector<double>> nominalY  ;
    std::unique_ptr<std::vector<double>> weight    ;

    HistGetter plotter;



};

#endif

void make2DTemplateWithAdaKern(std::string fileName, std::string outFileName,std::string arguments){
    std::cout << outFileName <<std::endl;
    make2DTemplateWithAdaKernAnalyzer a(fileName,"treeMaker/Events",arguments);
    a.analyze(100000);
    a.process(outFileName);
}
