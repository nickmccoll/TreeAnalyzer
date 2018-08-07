
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"

#include "HiggsAnalysis/CombinedLimit/interface/utils.h"

#include <string.h>



class NuisPlotterAnalyzer {
public:

    NuisPlotterAnalyzer(const std::string inputFileName, const std::string& pdfName, const std::vector<std::string>& nuisNames, const std::string outName = "")
{
        auto rootFile = new TFile(inputFileName.c_str(),"READ");
        RooWorkspace * w;
        rootFile->GetObject("w",w);
        w->saveSnapshot("prefit", utils::returnAllVars(w));
        auto pdf = w->pdf(pdfName.c_str()); // shapeBkg_losttw_std_e_L_LP_full_13TeV_opt
        if(pdf==0) return;

        auto hN = pdf->createHistogram((pdfName+"_nom").c_str(),*w->var("MJ"),RooFit::YVar(*w->var("MR")),RooFit::IntrinsicBinning());
        hN->SetName((pdfName+"_nom").c_str());
        plotter.add1D(hN);

        for(const auto& n : nuisNames){
            auto var = w->var(n.c_str()); //CMS_HHlnujj_losttw_PTX_L
            auto nus = w->pdf((n+"_Pdf").c_str());
            double mean = w->var((n+"_In").c_str())->getVal();

            double err=0;
            for(unsigned int iS = 0; iS < 3; ++iS ){
                if(std::strcmp(nus->findServer(iS)->ClassName(), "RooConstVar")!=0) continue;
                err = ((RooConstVar*)nus->findServer(iS))->getVal();
                break;
            }

            var->setVal(mean+err);
            auto hNU= pdf->createHistogram((pdfName+"_"+n+"_up").c_str(),*w->var("MJ"),RooFit::YVar(*w->var("MR")),RooFit::IntrinsicBinning());
            hNU->SetName((pdfName+"_"+n+"_up").c_str());
            plotter.add1D(hNU);
            var->setVal(mean-err);
            auto hND= pdf->createHistogram((pdfName+"_"+n+"_down").c_str(),*w->var("MJ"),RooFit::YVar(*w->var("MR")),RooFit::IntrinsicBinning());
            hND->SetName((pdfName+"_"+n+"_down").c_str());
            plotter.add1D(hND);
            w->loadSnapshot("prefit");
        }

        plotter.write(outName);
        rootFile->Close();

}
    HistGetter plotter;
};

#endif

void NuisPlotter(const std::string inputFileName, const std::string& pdfName, const std::vector<std::string>& nuisNames, const std::string outName = ""){
    NuisPlotterAnalyzer(inputFileName, pdfName, nuisNames, outName);
}
