
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

#include "PDFAdder.h"




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
        auto sA1  = p.addFloat("sA1"  ,"Alpha 1 systematic size...negative if you don't want to fit it",false,-1);
        auto mTop = p.addBool("mTop","True if you want a lower bound (100GeV) on the fit");
        auto fHN = p.addString("fH","Fitting histogram file name",true);
        auto nH  = p.addString("nH","fitting histogram name",true);
        p.parse(arguments);


        RooWorkspace w("w",false);

        auto * fH =  TObjectHelper::getFile(*fHN);
        auto oH =TObjectHelper::getObject<TH2F>(fH,*nH);
        TH2 * inH = 0;
        if(*mTop){
            plotter.add1D((TH2*)oH->Clone("orig_data"));
            float binW = (maxHbbMass-minHbbMass)/float(nHbbMassBins);
            int nBins = (maxHbbMass-100.0)/binW;
            inH = cutHist(&*oH,nBins,100,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
            inH->SetDirectory(0);
        } else {
            inH = &*oH;
        }
        plotter.add1D((TH2*)inH->Clone("data"));


        //Setup axes
        auto * xAxis = inH->GetXaxis();
        auto * yAxis = inH->GetYaxis();
        w.factory((*varX+"[0,10000]").c_str());
        w.factory((*varY+"[0,10000]").c_str());
        w.var(varX->c_str())->setBinning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax()));
        w.var(varY->c_str())->setBinning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax()));
        RooArgSet varset;
        RooArgList varlist;
        varset.add(*w.var(varX->c_str()));
        varset.add(*w.var(varY->c_str()));
        varlist.add(*w.var(varX->c_str()));
        varlist.add(*w.var(varY->c_str()));

        //Now for binning use the non-cut version:
        xAxis = oH->GetXaxis();
        yAxis = oH->GetYaxis();

        const std::string modelName = "model";
        const std::string fitName   = "model_fit";
        const std::string kdeName   = "model_kde";
        const std::string fitVar = *kdeX ? *varY : *varX;
        const std::string kdeVar = *kdeX ? *varX : *varY;


        //Data to fit to
        RooDataHist fitDataHist((*nH+"DH").c_str(),(*nH+"DH").c_str(),varlist,&*inH);
        w.import(fitDataHist);



        //KDE dimension
        PDFAdder::addHistoShapeFromFile(&w,"model","kde",{kdeVar},*fTN,*nT,PDFAdder::InterpSysts());

        //Fit dimension
        CJSON json( *fJS);
        w.factory((std::string("scaleSyst[0,-")+ASTypes::flt2Str(*sS)+","+ASTypes::flt2Str(*sS)+"]").c_str());
        w.factory((std::string("resSyst[0,-")+ASTypes::flt2Str(*sR)+","+ASTypes::flt2Str(*sR)+"]").c_str());

        if(*sA1 > 0){
            w.factory((std::string("alpha1Syst[0,-")+ASTypes::flt2Str(*sA1)+","+ASTypes::flt2Str(*sA1)+"]").c_str());
            PDFAdder::addCB(&w,fitName, kdeVar,fitVar,"","",json,{{"scaleSyst",1}},{{"resSyst",1}},{{"alpha1Syst",1}});
        }else
            PDFAdder::addCB(&w,fitName, kdeVar,fitVar,"","",json,{{"scaleSyst",1}},{{"resSyst",1}});

        //2D PDF
        PDFAdder::conditionalProduct(&w,modelName,fitName,fitVar,kdeName);


        auto original2D = w.pdf(modelName.c_str())->createHistogram("originalPDF",
                *w.var(varX->c_str()),RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())),
                RooFit::YVar(*w.var(varY->c_str()),RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax())))) ;
        original2D->SetName("originalPDF");
        plotter.add1D(original2D);
        auto original1D = w.pdf(modelName.c_str())->createHistogram("nominalPDF",
                *w.var(fitVar.c_str()),   *kdeX ? RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax())) :
                        RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())));
        original1D->SetName(("originalPDF_"+ fitVar).c_str());
        plotter.add1D(original1D);

        w.pdf(modelName.c_str())->fitTo(*w.data((*nH+"DH").c_str()),RooFit::SumW2Error(kTRUE));

        std::cout <<"Doing "<< outFileName <<" Fit"<<std::endl;
        std::cout << "scaleSyst -> "<< w.var("scaleSyst")->getVal()<<std::endl;
        std::cout << "resSyst -> "<< w.var("resSyst")->getVal()<<std::endl;
        if(*sA1 > 0)std::cout << "alpha1Syst -> "<< w.var("alpha1Syst")->getVal()<<std::endl;


        std::string newMean = "(" + json.getP("mean")+")*(1+"+ASTypes::flt2Str(w.var("scaleSyst")->getVal())+")";
        json.replaceEntry("mean", newMean );
        std::string newSigma = "(" + json.getP("sigma")+")*(1+"+ASTypes::flt2Str(w.var("resSyst")->getVal())+")";
        json.replaceEntry("sigma", newSigma );
        if(*sA1 > 0){
            std::string newAlpha1 = "(" + json.getP("alpha")+")*(1+"+ASTypes::flt2Str(w.var("alpha1Syst")->getVal())+")";
            json.replaceEntry("alpha", newAlpha1 );
        }
        json.write(outFileName);

        auto fit2D = w.pdf(modelName.c_str())->createHistogram("histo",
                *w.var(varX->c_str()),RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())),
                RooFit::YVar(*w.var(varY->c_str()),RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax())))) ;
        fit2D->SetName("histo");
        plotter.add1D(fit2D);
        auto fit1D = w.pdf(modelName.c_str())->createHistogram("histo",
                *w.var(fitVar.c_str()),   *kdeX ? RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax())) :
                        RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())));
        fit1D->SetName(("histo_"+ fitVar).c_str());
        plotter.add1D(fit1D);

        plotter.write(outFileName +".root");
        fH->Close();
        delete fH;
        if(*mTop) delete inH;
    }

    TH2 * cutHist(const TH2* inHist, int nBinX, float minX, float maxX, int nBinY, float minY, float maxY){
        auto outH = new TH2F ("cutHist","",nBinX,minX,maxX,nBinY,minY,maxY);

        //first cut up conditional template
        for(int iX = 1; iX <= inHist->GetNbinsX(); ++iX){
            int oBinX = outH->GetXaxis()->FindFixBin(inHist->GetXaxis()->GetBinCenter(iX));
            if(oBinX < 1 || oBinX > outH->GetNbinsX() ) continue;
            for(int iY = 1; iY <= inHist->GetNbinsY(); ++iY){
                int oBinY = outH->GetYaxis()->FindFixBin(inHist->GetYaxis()->GetBinCenter(iY));
                if(oBinY < 1 || oBinY > outH->GetNbinsY() ) continue;
                outH->SetBinContent(oBinX,oBinY,inHist->GetBinContent(iX,iY));
                outH->SetBinError(oBinX,oBinY,inHist->GetBinError(iX,iY));
            }
        }

        return outH;
    }


    HistGetter plotter;
};

#endif

void fit2DFitKDETemplate(std::string outFileName,std::string arguments){
    fit2DFitKDETemplateAnalyzer a(outFileName, arguments);
}
