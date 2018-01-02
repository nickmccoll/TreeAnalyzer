
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

    make2DTemplateWithAdaKernAnalyzer(std::string fileName, std::string treeName,std::string arguments ) : BaseTreeAnalyzer(fileName,treeName,2)
    {


        ParParser p;
        name      = p.addString("n","histogram name",true);
        auto vx   = p.addString("vx","x variable, the probability variable P(X|Y)",true);
        auto vy   = p.addString("vy","y variable, the conditional variable P(X|Y)",true);
        auto xb   = p.addVFloat("xb","x-variable binning",true);
        auto yb   = p.addVFloat("yb","y-variable binning",true);
        auto ycb  = p.addVFloat("ycb","y-variable coarse binning (explicit bin edges)",true);
        auto s   = p.addString("s","selection",false,"1.0");
        auto w   = p.addString("w","weight",false,"1.0");
        khxs     = p.addFloat("khxs","KDE h-x scale factor",false,1.);
        khxc     = p.addFloat ("khxc","KDE adaptive x-bandwidth cutoff",false,5);
        khys     = p.addFloat("khys","KDE h-y scale factor",false,1.);
        khyc     = p.addFloat ("khyc","KDE adaptive y-bandwidth cutoff",false,5);
        kss      = p.addBool ("kss","KDE sigma scaling");

        emin     = p.addFloat ("emin","Exponential fit min");
        emax     = p.addFloat ("emax","Exponential fit max");

        hs       = p.addFloat("hs","Histogram scaling: x proportional scale",true);
        hr       = p.addFloat("hr","Histogram scaling: 1/x proportional scale",true);

        p.parse(arguments);

        if(xb->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");
        if(yb->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");
        if(ycb->size() < 2)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        xAxis.reset(new TAxis((*xb)[0],(*xb)[1],(*xb)[2]));
        yAxis.reset(new TAxis((*yb)[0],(*yb)[1],(*yb)[2]));
        ycAxis.reset(new TAxis((*ycb).size() -1,&ycb->at(0)));


        tree.getTree()->SetBranchStatus("*",1);
        sForm.reset(new TTreeFormula("sForm", TString::Format("%s*(%s)",w->c_str(),s->c_str()),tree.getTree()));
        vxForm.reset(new TTreeFormula("vxForm", vx->c_str(),tree.getTree()));
        vyForm.reset(new TTreeFormula("vyForm", vy->c_str(),tree.getTree()));

        const int nEntries =  tree.getTree()->GetEntries(TString::Format("%s*(%s)",w->c_str(),s->c_str()));

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
    const TH2* makeKDE(std::string name,bool doCoarse){
        const int   nBsX   = xAxis->GetNbins();
        const float minX   = xAxis->GetXmin();
        const float maxX   = xAxis->GetXmax();
        const int   nBsY   = yAxis->GetNbins();
        const float minY   = yAxis->GetXmin();
        const float maxY   = yAxis->GetXmax();
        const int   nBsYc  = !doCoarse ? nBsY : ycAxis->GetNbins();
        const float minYc  = !doCoarse ? minY : ycAxis->GetXmin();
        const float maxYc  = !doCoarse ? maxY : ycAxis->GetXmax();

        KDEProducer2D pdfProd(nominalX.get(),nominalY.get(),weight.get(),
                *khxs,nBsX,minX,maxX,*khxc,
                *khys,nBsY,minY,maxY,*khyc,
                *kss);


        plotter.add2D(pdfProd.getPDF(name+"_debug_KDE0","",nBsX,minX,maxX,nBsY,minY,maxY));

        TH2 * data = new TH2F((name+"_fine_data").c_str(),";data",nBsX,minX,maxX,nBsY,minY,maxY);
        for(unsigned int iP = 0; iP < nominalX->size(); ++iP){
            data->Fill((*nominalX)[iP],(*nominalY)[iP],(*weight)[iP]);
        }

        plotter.add2D(data);
        plotter.add2D(pdfProd.getAPDF(name+"_debug_fineKDE","",nBsX,minX,maxX,nBsY,minY,maxY));

        auto * pilot = pdfProd.getPilotPDF();
        pilot->SetName((name + "_debug_pilotKDE").c_str());
        plotter.add2D(pilot);

        plotter.add2D(pdfProd.getABandwidthsX(name+"_debug_bandwidthsX","",nBsX,minX,maxX,nBsY,minY,maxY));
        plotter.add2D(pdfProd.getABandwidthsY(name+"_debug_bandwidthsY","",nBsX,minX,maxX,nBsY,minY,maxY));
        plotter.add2D(pdfProd.getLocalVarX(name   +"_debug_varX","",nBsX,minX,maxX,nBsY,minY,maxY)              ) ;
        plotter.add2D(pdfProd.getLocalVarY(name   +"_debug_varY","",nBsX,minX,maxX,nBsY,minY,maxY)              ) ;

        TH2 * kde = pdfProd.getAPDF(name+"_KDE","",nBsX,minX,maxX,nBsYc,minYc,maxYc);
        kde->Scale(data->Integral()/kde->Integral());

        plotter.add2D(kde);
        return kde;
    }

    const TH2 * smoothTail(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
         oHist->Scale(1.0/oHist->Integral());
        TF1 expo("expo","expo",*emin,*emax);

        for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
            auto * proj = oHist->ProjectionX("q",iY,iY);
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
        plotter.add2D(oHist);
        return oHist;

    }
    const TH2 * conditional(const std::string& name, const TH2* iHist){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
            auto * proj = oHist->ProjectionX("q",iY,iY);
            double xIntegral = oHist->Integral(1,oHist->GetNbinsX(),iY,iY);
            if(!xIntegral) { throw std::invalid_argument("Analyzer::Analyzer() -> Making a conditional slice with no events!!!!!");}
            for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
                oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)/xIntegral);
            }
        }
        plotter.add2D(oHist);
        return oHist;

    }

    const TH2 * expandHisto(const std::string& name, const TH2* iHist){
        const int   nBsX   = xAxis->GetNbins();
        const float minX   = xAxis->GetXmin();
        const float maxX   = xAxis->GetXmax();
        const int   nBsY   = yAxis->GetNbins();
        const float minY   = yAxis->GetXmin();
        const float maxY   = yAxis->GetXmax();
        TH2 * oHist = new TH2F(name.c_str(),iHist->GetTitle(),nBsX,minX,maxX,nBsY,minY,maxY);
        for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
            auto * proj = iHist->ProjectionY("q",iX,iX);
            TGraph * graph = new TGraph(proj);
            for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
                const double y = oHist->GetYaxis()->GetBinCenter(iY);
                oHist->SetBinContent(iX,iY,graph->Eval(y,0,"S"));
            }
        }
        plotter.add2D(oHist);
        return oHist;
    }

    const TH2 * transform ( std::string name,const TH1 * iHist, std::function<double(double)> f ) {
      TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
      for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
          const double x = oHist->GetXaxis()->GetBinCenter(iX);
          for(int iY = 1; iY <= oHist->GetNbinsY(); ++iY){
              oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*f(x));
          }
      }
      plotter.add2D(oHist);
       return oHist;
    };

    void process(std::string outFileName) {
//
//        auto * kde = makeKDE(*name,true);
//        kde = smoothTail(*name+ "_coarse_smooth", kde);
//        kde = conditional(*name + "_coarse_conditional", kde);
//        kde = expandHisto(*name + "_expanded",kde);
//        kde = conditional(*name,kde);

        auto * kde = makeKDE(*name,false);
        kde = smoothTail(*name+ "_smooth", kde);
        kde = conditional(*name,kde);

        auto * upScale = transform(*name+"_PTUp_debug_beforeCond",kde,[&](double x){return  (1. + *hs*x);});
        conditional(*name+"_PTUp",upScale);
        auto * downScale = transform(*name+"_PTDown_debug_beforeCond",kde,[&](double x){return  1./(1. + *hs*x);});
        conditional(*name+"_PTDown",downScale);

        auto * upRes = transform(*name+"_OPTUp_debug_beforeCond",kde,[&](double x){return  (1. + *hr/x);});
        conditional(*name+"_OPTUp",upRes);
        auto * downRes = transform(*name+"_OPTDown_debug_beforeCond",kde,[&](double x){return 1./(1. + *hr/x);});
        conditional(*name+"_OPTDown",downRes);


        plotter.write(outFileName);
    }

    std::unique_ptr<TAxis> xAxis;
    std::unique_ptr<TAxis> yAxis;
    std::unique_ptr<TAxis> ycAxis;
    std::unique_ptr<TTreeFormula> sForm;
    std::unique_ptr<TTreeFormula> vxForm;
    std::unique_ptr<TTreeFormula> vyForm;

    std::shared_ptr<std::string> name;
    std::shared_ptr<double>       khxs;
    std::shared_ptr<double>       khxc;
    std::shared_ptr<double>       khys;
    std::shared_ptr<double>       khyc;
    std::shared_ptr<bool>         kss;
    std::shared_ptr<double>       hs;
    std::shared_ptr<double>       hr;
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
