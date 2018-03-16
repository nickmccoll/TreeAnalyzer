
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TTreeFormula.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooNDKeysPdf.h"
#include "RooBinning.h"
#include "TMath.h"


using namespace TAna;

class make1DTemplateWithScaledKernelsAnalyzer : public BaseTreeAnalyzer {
public:

    make1DTemplateWithScaledKernelsAnalyzer(std::string fileName, std::string treeName,std::string arguments ) : BaseTreeAnalyzer(fileName,treeName,2)
, x("x","x",1), w("w","w",1){


        nominalDataset.reset(new RooDataSet("nominal","nominal",RooArgSet(x,w),"w"));
        nominalDataset.reset(new RooDataSet("nominal","nominal",RooArgSet(x,w),"w"));
        nominalDataset.reset(new RooDataSet("nominal","nominal",RooArgSet(x,w),"w"));
        nominalDataset.reset(new RooDataSet("nominal","nominal",RooArgSet(x,w),"w"));

        nominalDataset  .reset(new RooDataSet("nominal","nominal",RooArgSet(x,w),"w"));
        upScaleDataset  .reset(new RooDataSet("upScale","upScale",RooArgSet(x,w),"w"));
        downScaleDataset.reset(new RooDataSet("downScale","nominal",RooArgSet(x,w),"w"));
        upResDataset    .reset(new RooDataSet("upRes","upRes",RooArgSet(x,w),"w"));
        downResDataset  .reset(new RooDataSet("downRes","downRes",RooArgSet(x,w),"w"));

        ParParser p;
        name     = p.addString("n","histogram name",true);
        auto v   = p.addString("x","variable",true);
        auto g   = p.addString("g","gen variable",true);
        auto b   = p.addVFloat("xb","x-variable binning",true);
        auto s   = p.addString("s","selection",false,"1.0");
        auto w   = p.addString("w","weight",false,"1.0");
        kr       = p.addFloat("kr","KDE h-sigma",false,3.0);
        ks       = p.addFloat("ks","KDE h-scale factor",false,1.);
        ken      = p.addBool("ken","Apply effective scale to h");
        ss       = p.addFloat("ss","Scale scale",true);
        rs       = p.addFloat("rs","Resolution scale",true);
        sa       = p.addFloat("sa","Scale alpha",true);
        ra       = p.addFloat("ra","Resolution alpha",true);
        t        = p.addString("t","Types of templates to calculate",false, "nrRsSo");

        auto vsf       = p.addString("vsf","varaible scale file",true);
        auto vsh       = p.addString("vsh","varaible scale hist",true);
        auto vsv       = p.addString("vsv","varaible scale variable",true);
        p.parse(arguments);

        if(b->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        auto file =  TObjectHelper::getFile(*vsf);
        avgScale.reset(new TObjectHelper::Hist1DContainer(file,*vsh));
        delete file;

        vAxis.reset(new TAxis((*b)[0],(*b)[1],(*b)[2]));


        tree.getTree()->SetBranchStatus("*",1);
        sForm.reset(new TTreeFormula("sForm", TString::Format("%s*(%s)",w->c_str(),s->c_str()),tree.getTree()));
        vForm.reset(new TTreeFormula("vForm", v->c_str(),tree.getTree()));
        gForm.reset(new TTreeFormula("gForm", g->c_str(),tree.getTree()));

        vsForm.reset(new TTreeFormula("vsForm", vsv->c_str(),tree.getTree()));



    }


    void loadVariables() override{};
    bool runEvent() override {
        double s   = sForm->EvalInstance();
        if(s==0) return false;

        sumW += s;
        sumW2 += s*s;


        double v   = vForm->EvalInstance();
        double g   = gForm->EvalInstance();
        double sv   = vsForm->EvalInstance();
        double avgS = avgScale->getBinContentByValue(sv).val();

        x.setVal(v);
        nominalDataset->add( RooArgSet(x),s);

        double upScale = (1.+*ss)*v;
        double downScale = (1.-*ss)*v;

        double upRes= std::max(0.0,avgS*g + (1+*rs)* (v -avgS*g));
        double downRes = std::max(0.0,avgS*g + (1-*rs)*(v -avgS*g));

        x.setVal(upScale);
        upScaleDataset    ->add( RooArgSet(x),s);
        x.setVal(downScale);
        downScaleDataset  ->add( RooArgSet(x),s);
        x.setVal(upRes);
        upResDataset      ->add( RooArgSet(x),s);
        x.setVal(downRes);
        downResDataset    ->add( RooArgSet(x),s);
        return true;
    }

    void process(std::string outFileName) {
        TFile * f = new TFile(outFileName.c_str(),"recreate");
        f->cd();

        double hScale = *ks;
        if(*ken && sumW ){
         hScale  *= TMath::Power(sumW/sumW2,-0.2);
        }


        auto mkKernel = [&](const TString name,  RooDataSet& data){
            std::cout <<"Starting: "<< name<<std::endl;
            RooNDKeysPdf pdf(name + "PDF",name+ "PDF",RooArgSet(x)        ,data  ,"m",hScale,*kr) ;
            auto h1 = data.createHistogram(name+"hist",x,RooFit::Binning(RooBinning (vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax()))) ;
            h1->Write();
            auto h2 = pdf.createHistogram(name+"pdf",x,RooFit::Binning(RooBinning (vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax()))) ;
            h2->Write();
        };

        if (t->find('n') != std::string::npos) mkKernel("nom"  ,*nominalDataset  );
        if (t->find('S') != std::string::npos) mkKernel("upS"  ,*upScaleDataset  );
        if (t->find('s') != std::string::npos) mkKernel("downS",*downScaleDataset);
        if (t->find('R') != std::string::npos) mkKernel("upR"  ,*upResDataset    );
        if (t->find('r') != std::string::npos) mkKernel("downR",*downResDataset  );


        f->Close();

    }

    std::unique_ptr<TAxis> vAxis;
    std::unique_ptr<TTreeFormula> sForm;
    std::unique_ptr<TTreeFormula> vForm;
    std::unique_ptr<TTreeFormula> gForm;
    std::unique_ptr<TTreeFormula> vsForm;
    std::unique_ptr<TObjectHelper::Hist1DContainer>               avgScale;

    std::shared_ptr<std::string> name;
    std::shared_ptr<double>       kr;
    std::shared_ptr<double>       ks;
    std::shared_ptr<bool>        ken;
    std::shared_ptr<double>       ss;
    std::shared_ptr<double>       rs;
    std::shared_ptr<double>       sa;
    std::shared_ptr<double>       ra;
    std::shared_ptr<std::string>   t;


    double sumW=0;
    double sumW2=0;

    RooRealVar x;
    RooRealVar w;
    std::unique_ptr<RooDataSet> nominalDataset  ;
    std::unique_ptr<RooDataSet> upScaleDataset  ;
    std::unique_ptr<RooDataSet> downScaleDataset;
    std::unique_ptr<RooDataSet> upResDataset    ;
    std::unique_ptr<RooDataSet> downResDataset  ;


};

#endif

void make1DTemplateWithScaledKernels(std::string fileName, std::string outFileName,std::string arguments){
    make1DTemplateWithScaledKernelsAnalyzer a(fileName,"treeMaker/Events",arguments);
    a.analyze(100000);
    a.process(outFileName);
}
