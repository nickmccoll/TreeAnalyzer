
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
#include "RooConstVar.h"
#include "RooBinning.h"
#include "RooCmdConfig.h"
#include <RooStats/ModelConfig.h>
#include "TCanvas.h"

#include "HistoPlotting/include/Drawing.h"

    bool isSignal = true;
    const std::string cName = "radHH";
    const std::string l = "mu";
    const std::string b = "M";
    const std::string p = "HP";
    const std::string h = "full";
    const std::string cat  = "std_"+ l +"_"+b+"_"+p +"_"+h+"_13TeV";
    const std::string postFix = true ? "" : "_opt"; //only for losttw and qg
    const std::string prefix =  isSignal ? "shapeSig" : "shapeBkg";
    const std::string pdfName = prefix + "_" + cName+"_"+cat+postFix;
    RooWorkspace * w = 0;

//    int nBinsX = 90;
    int nBinsX = 30;
//    int nBinsY = 72;
//    std::vector<double> binsY = {700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2600,2700,2800,2900,3000,3200,3400,3600,3800,4000};
//    int nBinsY = 132;
//    std::vector<double> binsY =  {700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2025,2050,2075,2100,2125,2150,2175,2200,2225,2250,2275,2300,2325,2350,2375,2400,2425,2450,2475,2500,2525,2550,2575,2600,2625,2650,2675,2700,2725,2750,2775,2800,2825,2850,2875,2900,2925,2950,2975,3000,3025,3050,3075,3100,3125,3150,3175,3200,3225,3250,3275,3300,3325,3350,3375,3400,3425,3450,3475,3500,3525,3550,3575,3600,3625,3650,3675,3700,3725,3750,3775,3800,3825,3850,3875,3900,3925,3950,3975,4000};
//    int nBinsY = 66;
//    std::vector<double> binsY = {700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,3100,3150,3200,3250,3300,3350,3400,3450,3500,3550,3600,3650,3700,3750,3800,3850,3900,3950,4000};
    int nBinsY = 86;
    std::vector<double> binsY = {700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,3100,3175,3250,3325,3400,3475,3550,3625,3700,3775,3850,3925,4000};
    std::string hhResSysName = "jer"; //0.05
    std::string hhScaleSysName = "jes"; //0.01

    std::string getSystName (const std::string& prefix, const std::string& proc,const std::string& name, const std::string& sel){
        std::string sn = prefix;
        if(proc.size() )sn += "_"+proc;
        if(name.size() )sn += "_"+name;
        if(sel.size() )sn += "_"+sel;
        return sn;
    };
    std::string systName(const std::string& name, const std::string& sel = "-1") {
        return getSystName("HHlnujj",cName,name, sel == "-1" ? cat : sel  );
    };


 TH2* createHistogram(const std::string& name,RooAbsPdf* pdf, const TH2* histBinning,
         RooRealVar* xV,RooRealVar* yV,RooArgSet* components ){
     TH2 * oH = (TH2*)histBinning->Clone(name.c_str());
     oH->Reset("M");

     const std::string normName = "n_exp_bin" +cat +"_proc_"+cName;
     auto norm = (RooAbsReal*)components->find(normName.c_str());
     double normV = norm->getVal();
          RooArgSet vars(*xV) ;
          vars.add(*yV) ;

     for( int iInX = 1; iInX <= histBinning->GetNbinsX(); ++iInX ){
       double centerX = histBinning->GetXaxis()->GetBinCenter(iInX);
       double widthX = histBinning->GetXaxis()->GetBinWidth(iInX);
       xV->setVal(centerX);
       for( int iInY = 1; iInY <= histBinning->GetNbinsY(); ++iInY ){
         double centerY = histBinning->GetYaxis()->GetBinCenter(iInY);
         double widthY = histBinning->GetYaxis()->GetBinWidth(iInY);
         yV->setVal(centerY);
         double val = pdf->getVal(vars);
         oH->Fill(centerX,centerY,val*widthX*widthY/normV);
       }
     }
     return oH;
 }

void doOneSet(const std::string& name,RooArgSet* components ){
    std::vector<TObject*> paramPads;
    auto pdf = (RooAbsPdf*)components->find(pdfName.c_str());


    auto fineXBins = RooFit::Binning(360,30,210);
    auto fineYBins = RooFit::Binning(3300,700,4000);

    auto coarseXBins = RooFit::Binning(nBinsX,30,210);
    RooBinning varBinningY(nBinsY,&binsY[0]);
//    auto coarseYBins = RooFit::Binning(132,700,4000);
    auto coarseYBins = RooFit::Binning(varBinningY);

    const std::string fineHistName = "Fine_"+name+"_"+cName+"_"+cat;
    TH2* histF = (TH2*)pdf->createHistogram(fineHistName.c_str(),
            *w->var("MJ"),fineXBins, RooFit::YVar(*w->var("MR"),fineYBins));
    paramPads.push_back(histF);


    //testing out my impementation of create histogram
//    const std::string coarseHist = "Coarse_"+name+"_"+cName+"_"+cat;
//    TH2* histC = (TH2*)pdf->createHistogram(coarseHist.c_str(),
//            *w->var("MJ"),coarseXBins,RooFit::YVar(*w->var("MR"),coarseYBins));
//    paramPads.push_back(histC);

    auto rebin=[](const TH2* hIn, const std::string& outName)->TH2* {
//        TH2 * histRe = new TH2F(outName.c_str(),";MJ;MR",nBinsX,30,210,132,700,4000);
        TH2 * histRe = new TH2F(outName.c_str(),";MJ;MR",nBinsX,30,210,nBinsY,&binsY[0]);
        for( int iInX = 1; iInX <= hIn->GetNbinsX(); ++iInX ){
          double centerX = hIn->GetXaxis()->GetBinCenter(iInX);
          for( int iInY = 1; iInY <= hIn->GetNbinsY(); ++iInY ){
            double centerY = hIn->GetYaxis()->GetBinCenter(iInY);
            histRe->Fill(centerX,centerY,hIn->GetBinContent(iInX,iInY));
          }
        }
        return histRe;
    };



    TH2 * histFRe = rebin(histF,fineHistName+"_Re");
//    paramPads.push_back(histFRe);
    const std::string coarseHist = "Coarse_"+name+"_"+cName+"_"+cat;
    TH2* histC = createHistogram(coarseHist,pdf,histFRe,w->var("MJ"),w->var("MR"),components);
    paramPads.push_back(histC);
//    paramPads.push_back(histCRe);
    TH2* histRat = (TH2*)histC->Clone();
    histRat->Divide(histFRe);
    paramPads.push_back(histRat);


    auto setNSigma=[&](const::std::string& nuisName, float nSigma = 0){
        auto var = w->var(nuisName.c_str()); //CMS_HHlnujj_losttw_PTX_L
        auto nus = w->pdf((nuisName+"_Pdf").c_str());
        double mean = w->var((nuisName+"_In").c_str())->getVal();

        double err=0;
        for(unsigned int iS = 0; iS < 3; ++iS ){
            if(std::strcmp(nus->findServer(iS)->ClassName(), "RooConstVar")!=0) continue;
            err = ((RooConstVar*)nus->findServer(iS))->getVal();
            break;
        }
        var->setVal(mean+err*nSigma);
    };


    setNSigma(hhResSysName,1);
    const std::string resUpName = "ResUp_"+name+"_"+cName+"_"+cat;
    TH2 * histResUp = (TH2*)pdf->createHistogram(resUpName.c_str(),
            *w->var("MJ"),fineXBins, RooFit::YVar(*w->var("MR"),fineYBins));
    TH2 * histResUpRe = rebin(histResUp,resUpName + "_Re");
    setNSigma(hhResSysName,-1);
    const std::string resDownName = "ResDown_"+name+"_"+cName+"_"+cat;
    TH2 * histResDown = (TH2*)pdf->createHistogram(resDownName.c_str(),
            *w->var("MJ"),fineXBins, RooFit::YVar(*w->var("MR"),fineYBins));
    TH2 * histResDownRe = rebin(histResDown,resDownName + "_Re");
    setNSigma(hhResSysName,0);
    setNSigma(hhScaleSysName,1);
    const std::string scaleUpName = "ScaleUp_"+name+"_"+cName+"_"+cat;
    TH2 * histScaleUp = (TH2*)pdf->createHistogram(scaleUpName.c_str(),
            *w->var("MJ"),fineXBins, RooFit::YVar(*w->var("MR"),fineYBins));
    TH2 * histScaleUpRe = rebin(histScaleUp,scaleUpName + "_Re");
    setNSigma(hhScaleSysName,0);
    histResUpRe->Divide(histFRe);
    histResDownRe->Divide(histFRe);
    histScaleUpRe->Divide(histFRe);

    paramPads.push_back(histResDownRe);
    paramPads.push_back(histResUpRe);
    paramPads.push_back(histScaleUpRe);

    Drawing::drawAll(paramPads, name.c_str(),"COLZ");
    std::cout << name<<" -> " << histFRe->Integral() <<" "<<" "<<histC->Integral()<<std::endl;
}


void go(std::string inFile){
    TFile * f = new TFile(inFile.c_str(),"read");
    w = 0;
    f->GetObject("w",w);
    auto model = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
    auto totPDF = model->GetPdf();
    auto components = totPDF->getComponents();

    w->var("MH")->setVal(1000);
    doOneSet("1000",components);

    w->var("MH")->setVal(2000);
    doOneSet("2000",components);

    w->var("MH")->setVal(2200);
    doOneSet("2200",components);

    w->var("MH")->setVal(3100);
    doOneSet("3100",components);

    w->var("MH")->setVal(3500);
    doOneSet("3500",components);

//    const std::string prefix =  isSignal ? "shapeSig" : "shapeBkg";
//    const std::string pdfName = prefix + "_" + cName+"_"+cat+postFix;
//    const std::string normName = "n_exp_bin" +cat +"_proc_"+cName;
//
//






  }

#endif

void testTemplateBinning(std::string inFile= "combined.root"){
    go(inFile);
}
