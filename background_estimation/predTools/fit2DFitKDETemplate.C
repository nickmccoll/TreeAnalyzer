
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include <string.h>
#include <regex>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"

#include "InputsHelper.h"
using ASTypes::flt2Str;



class fit2DFitKDETemplateAnalyzer {
public:

    fit2DFitKDETemplateAnalyzer(std::string outFileName,std::string arguments )
    {
        ParParser p;
        auto fTN = p.addString("fT","KDE template file name",true);
        auto nT  = p.addString("nT","KDE template histogram name",true);
        auto varX= p.addString("varX","x-variable name",true);
        auto varY = p.addString("varY","y-variable name",true);
        auto kdeX = p.addBool("kdeX","True if the KDE is set to be the x axis");
        auto fJS  = p.addString("fJS","Fit JSON filename",true);
        auto sS   = p.addFloat("sS"  ,"Scale systematic size",true);
        auto sR   = p.addFloat("sR"  ,"Res systematic size",true);
        auto sA1  = p.addFloat("sA1" ,
                "Alpha 1 systematic size...negative if you don't want to fit it",false,-1);
        auto mTop = p.addBool("mTop","True if you want a lower bound (100GeV) on the fit");
        auto fHN = p.addString("fH","Fitting histogram file name",true);
        auto nHN  = p.addString("nH","fitting histogram name",true);
        p.parse(arguments);


        auto * fH =  TObjectHelper::getFile(*fHN);
        auto oH =TObjectHelper::getObject<TH2F>(fH,*nHN);
        std::unique_ptr<TAxis> oXAxis((TAxis*)oH->GetXaxis()->Clone());
        std::unique_ptr<TAxis> oYAxis((TAxis*)oH->GetYaxis()->Clone());

        std::unique_ptr<TAxis> nXAxis((TAxis*)oH->GetXaxis()->Clone());
        std::unique_ptr<TAxis> nYAxis((TAxis*)oH->GetYaxis()->Clone());

        RooWorkspace w("w",false);
        TH2 * nH = 0;
        if(*mTop){
            plotter.add1D((TH2*)oH->Clone("orig_data"));

            std::vector<double> newHbbBins;
            for(int iB = 1; iB <= oXAxis->GetNbins(); ++iB){
                if(oXAxis->GetBinLowEdge(iB) < 100) continue;
                newHbbBins.push_back(oXAxis->GetBinLowEdge(iB));
            }
            newHbbBins.push_back(oXAxis->GetXmax());
            nXAxis.reset(new TAxis(newHbbBins.size()-1, &newHbbBins[0]));
            nH = cutHistogram(oH.get(),"cutHist","",nXAxis.get(),nYAxis.get());
            nH->SetDirectory(0);
        } else {
            nH = &*oH;
        }
        plotter.add1D((TH2*)nH->Clone("data"));


        //Setup axes
        w.factory((*varX+"[0,10000]").c_str());
        w.factory((*varY+"[0,10000]").c_str());

        auto getBinning=[](const TAxis* axis)->RooBinning{
            return axis->GetXbins()->GetSize()
                    ? RooBinning(axis->GetNbins(),axis->GetXbins()->GetArray())
                    : RooBinning(axis->GetNbins(),axis->GetXmin(),axis->GetXmax());
        };

        RooBinning oXBins =  getBinning(oXAxis.get());
        RooBinning oYBins =  getBinning(oYAxis.get());
        RooBinning nXBins =  getBinning(nXAxis.get());
        RooBinning nYBins =  getBinning(nYAxis.get());


        w.var(varX->c_str())->setBinning(nXBins);
        w.var(varY->c_str())->setBinning(nYBins);
        RooArgSet varset;
        RooArgList varlist;
        varset.add(*w.var(varX->c_str()));
        varset.add(*w.var(varY->c_str()));
        varlist.add(*w.var(varX->c_str()));
        varlist.add(*w.var(varY->c_str()));

        //Now for binning use the non-cut version:
        const std::string modelName = "model";
        const std::string fitName   = "model_fit";
        const std::string kdeName   = "model_kde";
        const std::string fitVar    = *kdeX ? *varY : *varX;
        const std::string kdeVar    = *kdeX ? *varX : *varY;

        //Data to fit to
        RooDataHist fitDataHist((*nHN+"DH").c_str(),(*nHN+"DH").c_str(),varlist,&*nH);
        w.import(fitDataHist);

        //KDE dimension
        PDFAdder::addHistoShapeFromFile(&w,"model","kde",{kdeVar},*fTN,*nT,PDFAdder::InterpSysts());

        //Fit dimension
        CJSON json( *fJS);
        w.factory((
                std::string("scaleSyst[0,-")+flt2Str(*sS)+","+flt2Str(*sS)+"]"
                ).c_str());
        w.factory((
                std::string("resSyst[0,-")+flt2Str(*sR)+","+flt2Str(*sR)+"]"
                ).c_str());

        if(*sA1 > 0){
            w.factory((
                    std::string("alpha1Syst[0,-")+flt2Str(*sA1)+","+flt2Str(*sA1)+"]"
                    ).c_str());
            PDFAdder::addCB(&w,fitName, kdeVar,fitVar,"","",json,
                    {{"scaleSyst",1}},{{"resSyst",1}},{{"alpha1Syst",1}});
        }else
            PDFAdder::addCB(&w,fitName, kdeVar,fitVar,"","",json,
                    {{"scaleSyst",1}},{{"resSyst",1}});

        //2D PDF
        PDFAdder::conditionalProduct(&w,modelName,fitName,fitVar,kdeName);

        auto original2D = w.pdf(modelName.c_str())->createHistogram("originalPDF",
                *w.var(varX->c_str()),RooFit::Binning(oXBins),
                RooFit::YVar(*w.var(varY->c_str()),RooFit::Binning(oYBins))) ;
        original2D->SetName("originalPDF");
        plotter.add1D(original2D);

        auto original1D = w.pdf(modelName.c_str())->createHistogram("nominalPDF",
                *w.var(fitVar.c_str()),*kdeX ? RooFit::Binning(oYBins):RooFit::Binning(oXBins));
        original1D->SetName(("originalPDF_"+ fitVar).c_str());
        plotter.add1D(original1D);

        w.pdf(modelName.c_str())->fitTo(*w.data((*nHN+"DH").c_str()),RooFit::SumW2Error(kTRUE));

        std::cout <<"Doing "<< outFileName <<" Fit"<<std::endl;
        std::cout << "scaleSyst -> "<< w.var("scaleSyst")->getVal()<<std::endl;
        std::cout << "resSyst -> "<< w.var("resSyst")->getVal()<<std::endl;
        if(*sA1 > 0)std::cout << "alpha1Syst -> "<< w.var("alpha1Syst")->getVal()<<std::endl;


        std::string newMean = "(" +
                json.getP("mean")+")*(1+"+flt2Str(w.var("scaleSyst")->getVal())
                +")";
        json.replaceEntry("mean", newMean );
        std::string newSigma = "(" +
                json.getP("sigma")+")*(1+"+flt2Str(w.var("resSyst")->getVal())
                +")";
        json.replaceEntry("sigma", newSigma );
        if(*sA1 > 0){
            std::string newAlpha1 = "(" +
                    json.getP("alpha")+")*(1+"+flt2Str(w.var("alpha1Syst")->getVal())
                    +")";
            json.replaceEntry("alpha", newAlpha1 );
        }
        json.write(outFileName);

        auto fit2D = w.pdf(modelName.c_str())->createHistogram("histo",
                *w.var(varX->c_str()),RooFit::Binning(oXBins),
                RooFit::YVar(*w.var(varY->c_str()),RooFit::Binning(oYBins))) ;
        fit2D->SetName("histo");
        plotter.add1D(fit2D);

        auto fit1D = w.pdf(modelName.c_str())->createHistogram("histo",
                *w.var(fitVar.c_str()),*kdeX ? RooFit::Binning(oYBins) : RooFit::Binning(oXBins));
        fit1D->SetName(("histo_"+ fitVar).c_str());
        plotter.add1D(fit1D);

        plotter.write(outFileName +".root");
        fH->Close();
        delete fH;
        if(*mTop) delete nH;
    }


    HistGetter plotter;
};

#endif

void fit2DFitKDETemplate(std::string outFileName,std::string arguments){
    fit2DFitKDETemplateAnalyzer a(outFileName, arguments);
}
