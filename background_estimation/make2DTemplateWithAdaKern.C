
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
        khys     = p.addFloat("khys","KDE h-y scale factor",false,1.);
        khc      = p.addFloat ("khc","KDE adaptive bandwidth cutoff",false,5);
        kss      = p.addBool ("kss","KDE sigma scaling");

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

    void process(std::string outFileName) {
//
        const int   nBsX   = xAxis->GetNbins();
        const float minX   = xAxis->GetXmin();
        const float maxX   = xAxis->GetXmax();
        const int   nBsY   = yAxis->GetNbins();
        const float minY   = yAxis->GetXmin();
        const float maxY   = yAxis->GetXmax();
        std::vector<TH1*> outHists;

        KDEProducer2D pdfProd(nominalX.get(),nominalY.get(),weight.get(),
                *khxs,nBsX,minX,maxX,
                *khys,nBsY,minY,maxY,
                *khc,*kss);
        auto nomName = *name + "_nom";
        TH2 * hPDF = pdfProd.getPDF(nomName+"_pdf0",";pdf",nBsX,minX,maxX,nBsY,minY,maxY);
        TH2 * hPDFD = pdfProd.convToHist(hPDF);
        hPDFD->SetName( (std::string(hPDF->GetName())+ "D").c_str());

        TH2 * dataH = new TH2F((nomName+"_data").c_str(),";data",nBsX,minX,maxX,nBsY,minY,maxY);
        for(unsigned int iP = 0; iP < nominalX->size(); ++iP){
            dataH->Fill((*nominalX)[iP],(*nominalY)[iP],(*weight)[iP]);
        }


        outHists.push_back(hPDF);
        outHists.push_back(hPDFD);
        outHists.push_back(dataH);

        TH2 * haPDF = pdfProd.getAPDF(nomName+"_pdf",";pdf",nBsX,minX,maxX,nBsY,minY,maxY);
        TH2 * haPDFD = pdfProd.convToHist(haPDF);
        haPDFD->SetName( (std::string(haPDF->GetName())+ "D").c_str());
        TH2 * hpPDF = pdfProd.getPilotPDF();
        hpPDF->SetName( (std::string(haPDF->GetName())+ "P").c_str());
        TH2 * hpPDFD = pdfProd.convToHist(hpPDF);
        hpPDFD->SetName( (std::string(haPDF->GetName())+ "PD").c_str());
        TH2 * hbx = pdfProd.getABandwidthsX(nomName+"_bandwidthsX",";bandwidths",nBsX,minX,maxX,nBsY,minY,maxY);
        TH2 * hby = pdfProd.getABandwidthsY(nomName+"_bandwidthsY",";bandwidths",nBsX,minX,maxX,nBsY,minY,maxY);
        TH2 * hsx    = pdfProd.getLocalVarX(nomName+"_varX",";sigmas",nBsX,minX,maxX,nBsY,minY,maxY);
        TH2 * hsy    = pdfProd.getLocalVarY(nomName+"_varY",";sigmas",nBsX,minX,maxX,nBsY,minY,maxY);

        outHists.push_back(haPDF);
        outHists.push_back(haPDFD);
        outHists.push_back(hpPDF);
        outHists.push_back(hpPDFD);
        outHists.push_back(hbx);
        outHists.push_back(hby);
        outHists.push_back(hsx);
        outHists.push_back(hsy);

        TFile * f = new TFile(outFileName.c_str(),"recreate");
        f->cd();
        for(auto* h : outHists){
            h->Write();
        }
        f->Close();

    }

    std::unique_ptr<TAxis> xAxis;
    std::unique_ptr<TAxis> yAxis;
    std::unique_ptr<TAxis> ycAxis;
    std::unique_ptr<TTreeFormula> sForm;
    std::unique_ptr<TTreeFormula> vxForm;
    std::unique_ptr<TTreeFormula> vyForm;

    std::shared_ptr<std::string> name;
    std::shared_ptr<double>       khxs;
    std::shared_ptr<double>       khys;
    std::shared_ptr<double>       khc;
    std::shared_ptr<bool>         kss;
    std::shared_ptr<double>       hs;
    std::shared_ptr<double>       hr;


    std::unique_ptr<std::vector<double>> nominalX  ;
    std::unique_ptr<std::vector<double>> nominalY  ;
    std::unique_ptr<std::vector<double>> weight    ;



};

#endif

void make2DTemplateWithAdaKern(std::string fileName, std::string outFileName,std::string arguments){
    std::cout << outFileName <<std::endl;
    make2DTemplateWithAdaKernAnalyzer a(fileName,"treeMaker/Events",arguments);
    a.analyze(100000);
    a.process(outFileName);
}
