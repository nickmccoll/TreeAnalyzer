
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
        kt       = p.addString("t","Types of templates to calculate",false, "nrRsSo");
        ks       = p.addFloat("ks","KDE scaling: scale");
        kr       = p.addFloat("kr","KDE scaling: resolution");
        khs      = p.addFloat("khs","KDE h-scale factor",false,1.);
        khc      = p.addFloat ("khc","KDE adaptive bandwidth cutoff",false,5);
        kss      = p.addBool ("kss","KDE sigma scaling");


        hs       = p.addFloat("hs","Histogram scaling: x proportional scale",true);
        hr       = p.addFloat("hr","Histogram scaling: 1/x proportional scale",true);

        auto asf       = p.addString("vsf","Average scale file",true);
        auto ash       = p.addString("vsh","Average scale hist",true);
        auto asv       = p.addString("vsv","Average scale variable",true);
        p.parse(arguments);

        if(b->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        auto file =  TObjectHelper::getFile(*asf);
        avgScale.reset(new TObjectHelper::Hist1DContainer(file,*ash));
        delete file;

        vAxis.reset(new TAxis((*b)[0],(*b)[1],(*b)[2]));


        tree.getTree()->SetBranchStatus("*",1);
        sForm.reset(new TTreeFormula("sForm", TString::Format("%s*(%s)",w->c_str(),s->c_str()),tree.getTree()));
        vForm.reset(new TTreeFormula("vForm", v->c_str(),tree.getTree()));
        gForm.reset(new TTreeFormula("gForm", g->c_str(),tree.getTree()));
        vsForm.reset(new TTreeFormula("vsForm", asv->c_str(),tree.getTree()));

        const int nEntries =  tree.getTree()->GetEntries(TString::Format("%s*(%s)",w->c_str(),s->c_str()));

        nominalX .reset(new std::vector<double>);
        upSX     .reset(new std::vector<double>);
        downSX   .reset(new std::vector<double>);
        upRX     .reset(new std::vector<double>);
        downRX   .reset(new std::vector<double>);
        weight   .reset(new std::vector<double>);

        nominalX ->reserve(nEntries);
        upSX     ->reserve(nEntries);
        downSX   ->reserve(nEntries);
        upRX     ->reserve(nEntries);
        downRX   ->reserve(nEntries);
        weight   ->reserve(nEntries);
    }


    void loadVariables() override{};
    bool runEvent() override {
        double s   = sForm->EvalInstance();
        if(s==0) return false;

        double v   = vForm->EvalInstance();
        double g   = gForm->EvalInstance();
        double sv   = vsForm->EvalInstance();
        double avgS = avgScale->getBinContentByValue(sv).val();

        double upScale = *ks*v;
        double downScale = (*ks-1.)*v;

        double upRes= std::max(0.0,avgS*g + *kr* (v -avgS*g));
        double downRes = std::max(0.0,avgS*g + (*kr-1.)*(v -avgS*g));

        nominalX ->push_back(v);
        upSX     ->push_back(upScale);
        downSX   ->push_back(downScale);
        upRX     ->push_back(upRes);
        downRX   ->push_back(downRes);
        weight   ->push_back(s);
        return true;
    }

    void process(std::string outFileName) {
        TFile * f = new TFile(outFileName.c_str(),"recreate");
        f->cd();


        auto mkKernel = [](const std::string name,  const std::vector<double> * xvals, const std::vector<double> * weights,
                double hSF, int nBinsX, double xMin, double xMax, double trimFactor, bool doSigmaScaling, std::vector<TH1*>& outHists){
            std::cout <<"Starting: "<< name<<std::endl;

            KDEProducer pdfProd(xvals,weights,hSF,nBinsX,xMin,xMax,trimFactor,doSigmaScaling);
            std::string nName = "_pdf0";
            if(nBinsX <= 0) nName = "_pdf"; //in case the non-adaptive is the main pdf

            TH1 * hPDF = pdfProd.getPDF(name+nName,";pdf",nBinsX,xMin,xMax);
            TH1 * hPDFD = pdfProd.convToHist(hPDF);
            hPDFD->SetName( (std::string(hPDF->GetName())+ "D").c_str());

            TH1 * dataH = new TH1F((name+"_data").c_str(),";data",nBinsX,xMin,xMax);
            for(unsigned int iP = 0; iP < xvals->size(); ++iP){
                dataH->Fill((*xvals)[iP],(*weights)[iP]);
            }

            outHists.push_back(hPDF);
            outHists.push_back(hPDFD);
            outHists.push_back(dataH);

            if(nBinsX > 0){
                TH1 * haPDF = pdfProd.getAPDF(name+"_pdf",";pdf",nBinsX,xMin,xMax);
                TH1 * haPDFD = pdfProd.convToHist(haPDF);
                haPDFD->SetName( (std::string(haPDF->GetName())+ "D").c_str());
                TH1 * hpPDF = pdfProd.getPilotPDF();
                hpPDF->SetName( (std::string(haPDF->GetName())+ "P").c_str());
                TH1 * hpPDFD = pdfProd.convToHist(hpPDF);
                hpPDFD->SetName( (std::string(haPDF->GetName())+ "PD").c_str());
                TH1 * hbPDF = pdfProd.getABandwidths(name+"_bandwidths",";bandwidths",nBinsX,xMin,xMax);
                TH1 * hs    = pdfProd.getLocalVariance(name+"_sigmas",";sigmas",nBinsX,xMin,xMax);



                outHists.push_back(haPDF);
                outHists.push_back(haPDFD);
                outHists.push_back(hpPDF);
                outHists.push_back(hpPDFD);
                outHists.push_back(hbPDF);
                outHists.push_back(hs);
            }
        };


        std::vector<std::thread> threads;
        std::vector<TH1*> outHs;
        if (kt->find('n') != std::string::npos) mkKernel("nom"  ,nominalX.get(),weight.get(),*khs,vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax(),*khc,*kss,std::ref(outHs));
        if (kt->find('S') != std::string::npos) mkKernel("upS"  ,upSX    .get(),weight.get(),*khs,vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax(),*khc,*kss,std::ref(outHs));
        if (kt->find('s') != std::string::npos) mkKernel("downS",downSX  .get(),weight.get(),*khs,vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax(),*khc,*kss,std::ref(outHs));
        if (kt->find('R') != std::string::npos) mkKernel("upR"  ,upRX    .get(),weight.get(),*khs,vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax(),*khc,*kss,std::ref(outHs));
        if (kt->find('r') != std::string::npos) mkKernel("downR",downRX  .get(),weight.get(),*khs,vAxis->GetNbins(),vAxis->GetXmin(),vAxis->GetXmax(),*khc,*kss,std::ref(outHs));
        for(auto& t : threads) t.join();

        auto mkTrnformed = [&] (const TH1 * inH, std::string name, std::function<double(double)> f ) -> TH1* {
          TH1 * outH = (TH1*)inH->Clone(name.c_str());
          for(int iB = 1; iB <= inH->GetNbinsX(); ++iB){
              outH->SetBinContent(iB,inH->GetBinContent(iB)*f(inH->GetBinCenter(iB)) );
              outH->SetBinError(iB,0);
          }
           return outH;
        };

        if (kt->find('n') != std::string::npos){
            const TH1 * nomH = 0;
            for(const auto * h : outHs) if(!std::strcmp(h->GetName(), "nom_pdf")) nomH = h;

            outHs.push_back(mkTrnformed(nomH,"downPT_pdf", [&](double x){return  1./(1. + *hs*x/vAxis->GetXmax());})  );
            outHs.push_back(mkTrnformed(nomH,"upPT_pdf", [&](double x){return  (1. + *hs*x/vAxis->GetXmax());})  );
            outHs.push_back(mkTrnformed(nomH,"downOPT_pdf", [&](double x){return  1./(1. + *hr*vAxis->GetXmin()/x);})  );
            outHs.push_back(mkTrnformed(nomH,"upOPT_pdf", [&](double x){return  (1. + *hr*vAxis->GetXmin()/x);})  );
        }

        for(auto* h : outHs) h->Write();
        f->Close();

    }

    std::unique_ptr<TAxis> vAxis;
    std::unique_ptr<TTreeFormula> sForm;
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


    std::unique_ptr<std::vector<double>> nominalX  ;
    std::unique_ptr<std::vector<double>> upSX      ;
    std::unique_ptr<std::vector<double>> downSX    ;
    std::unique_ptr<std::vector<double>> upRX      ;
    std::unique_ptr<std::vector<double>> downRX    ;
    std::unique_ptr<std::vector<double>> weight    ;



};

#endif

void make1DTemplateWithAdaKern(std::string fileName, std::string outFileName,std::string arguments){
    make1DTemplateWithAdaKernAnalyzer a(fileName,"treeMaker/Events",arguments);
    a.analyze(100000);
    a.process(outFileName);
}
