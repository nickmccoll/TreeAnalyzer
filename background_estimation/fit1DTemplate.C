
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "Utilities/HiggsCombineImport/interface/VerticalInterpHistPdf.h"
#include <string.h>
#include <regex>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"





class fit1DTemplateAnalyzer {
public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    fit1DTemplateAnalyzer(std::string outFileName,std::string arguments )
    {


        ParParser p;
        auto fTN = p.addString("fT","template file name",true);
        auto nT  = p.addString("nT","template histogram base name",true);
        auto s   = p.addString("s" ,"Comma separated list of systematics",true);
        auto fHN = p.addString("fH","Fitting histogram file name",true);
        auto nH  = p.addString("nH","fitting histogram name",true);
        p.parse(arguments);

        std::vector<std::string> systList = getList(*s);

        fT =  TObjectHelper::getFile(*fTN);
        fH =  TObjectHelper::getFile(*fHN);
        auto h = TObjectHelper::getObject<TH1F>(fT,*nT);
        auto * xAxis = h->GetXaxis();

        RooWorkspace w("w",false);
        RooArgSet varset;
        RooArgList varlist;
        w.factory("x[0,10000]");
        w.var("x")->setBinning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax()));
        varset.add(*w.var("x"));
        varlist.add(*w.var("x"));

        RooDataHist rH("nominalHist","nominalHist",varlist,&*h);
        RooHistPdf rP("nominalPDF","nominalPDF",varlist,rH);
        w.import(rH);
        w.import(rP);
        RooArgList coeffList;
        RooArgList pdfList(*w.pdf("nominalPDF"));

        const std::vector<std::string> systVar = {"Up" , "Down"};
        for(const auto& syst : systList){
            w.factory((syst+"[-5,5]").c_str());
            coeffList.add(*w.var(syst.c_str()));
            for(const auto& var : systVar){
                auto sH = TObjectHelper::getObject<TH1F>(fT,*nT+"_"+syst+var);
                RooDataHist rSH((syst+var+"Hist").c_str(),(syst+var+"Hist").c_str(),varlist,&*sH);
                RooHistPdf Spdf((syst+var+"PDF").c_str(),(syst+var+"PDF").c_str(),varset,rSH,0);
                w.import(rSH);
                w.import(Spdf);
                pdfList.add(*w.pdf((syst+var+"PDF").c_str()));
            }

        }

        FastVerticalInterpHistPdf interpPDF("interpPDF","interpPDF",*w.var("x"),pdfList,coeffList);
        w.import(interpPDF);

        auto inH =TObjectHelper::getObject<TH1F>(fH,*nH);
        RooDataHist fitDataHist((*nH+"DH").c_str(),(*nH+"DH").c_str(),varlist,&*inH);
        w.import(fitDataHist);
        w.pdf("interpPDF")->fitTo(*w.data((*nH+"DH").c_str()),RooFit::SumW2Error(kTRUE));

        auto fitHist = w.pdf("interpPDF")->createHistogram(nT->c_str(),
                *w.var("x"),RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax()))) ;
        //                fitHist->SetName(nT->c_str());
        for(const auto& syst : systList){
            std::cout << syst <<" -> "<< w.var(syst.c_str())->getVal()<<std::endl;
        }
        fitHist->SetName(nT->c_str());
        plotter.add1D(fitHist);


        plotter.write(outFileName);
        fT->Close();
        fH->Close();
    }

    std::vector<std::unique_ptr<TH1F> > pdfHistos;
    std::vector<RooRealVar> plottingVars;
    std::vector<RooRealVar> systVars;
    std::vector<RooDataHist> dataHists;
    std::vector<RooHistPdf> pdfs;


    std::vector<std::string> getList(const std::string& inList){
        std::vector<std::string> systList(std::sregex_token_iterator(inList.begin(), inList.end(), std::regex(","), -1), std::sregex_token_iterator());
        return systList;
    }

    HistGetter plotter;
    TFile * fT = 0;
    TFile * fH = 0;
    std::shared_ptr<std::vector<double>>       xb;
    std::shared_ptr<std::vector<double>>       yb;

};

#endif

void fit1DTemplate(std::string outFileName,std::string arguments){
    fit1DTemplateAnalyzer a(outFileName, arguments);
}
