
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

    const TH1* makeFitHist(const std::string& name,const TH1 * iHist, const std::vector<double>& coefList, const std::vector<std::unique_ptr<TH1F>>& upHists
            , const std::vector<std::unique_ptr<TH1F>>& downHists) {

        const double _smoothRegion =1;
        auto smoothStepFunc =[&] (double x) ->double{
            if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
            double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
            return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
        };

        TH1 * hd = (TH1*)iHist->Clone(name.c_str());
        const int nBinsX = hd->GetNbinsX();
        std::vector<double> upHistN;
        std::vector<double> downHistN;
        for(unsigned int iC = 0; iC < coefList.size(); ++iC){
            upHistN.push_back(upHists[iC]->Integral(1,nBinsX));
            downHistN.push_back(downHists[iC]->Integral(1,nBinsX));
        }

        hd->Scale(1./hd->Integral(1,hd->GetNbinsX()) );
        for(int iX = 1; iX <= nBinsX; ++iX){
                const double xCen = hd->GetXaxis()->GetBinCenter(iX);
                const double oldV =  hd->GetBinContent(iX);
                double newV =  oldV;
                for(unsigned int iC = 0; iC < coefList.size(); ++iC){
                    const double x = coefList[iC];
                    const double a = 0.5*x;
                    const double b = smoothStepFunc(x);
                    const double upV = upHists[iC]->GetBinContent(iX)/upHistN[iC];
                    const double downV = downHists[iC]->GetBinContent(iX)/downHistN[iC];
                    newV+= a*((upV-downV) + b*(upV+downV - 2*oldV));
                }
                hd->SetBinContent(iX,newV);
                hd->SetBinError(iX,0);
            }
        hd->Scale(1./hd->Integral(1,hd->GetNbinsX()) );
        plotter.add1D(hd);
        return hd;
    }

    const TH1* systHist(const std::string& name, const TH1* iH1,const TH1* nH1,const TH1* sH1){
        TH1* oH1 = (TH1*)iH1->Clone(name.c_str());
        for(int iX = 1; iX <= oH1->GetNbinsX(); ++iX){
            if(nH1->GetBinContent(iX)){
                const double sf = sH1->GetBinContent(iX)/nH1->GetBinContent(iX);
                oH1->SetBinContent(iX,oH1->GetBinContent(iX)*sf);
                oH1->SetBinError(iX,0);
            }
        }
        oH1->Scale(1./oH1->Integral(1,oH1->GetNbinsX()));
        plotter.add1D(oH1);
        return oH1;
    }

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
        std::vector<std::unique_ptr<TH1F>> downHists;
        std::vector<std::unique_ptr<TH1F>> upHists;
        for(const auto& syst : systList){
            w.factory((syst+"[-5,5]").c_str());
            coeffList.add(*w.var(syst.c_str()));
            for(const auto& var : systVar){
                auto * histV = &(var== "Up" ? upHists : downHists);
                histV->emplace_back(TObjectHelper::getObject<TH1F>(fT,*nT+"_"+syst+var));
                RooDataHist rSH((syst+var+"Hist").c_str(),(syst+var+"Hist").c_str(),varlist,&*histV->back());
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

        auto fitHist = w.pdf("interpPDF")->createHistogram((std::string("debug_") + *nT).c_str(),
                *w.var("x"),RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax()))) ;

        for(const auto& syst : systList){
            std::cout << syst <<" -> "<< w.var(syst.c_str())->getVal()<<std::endl;
        }
//        fitHist->SetName(nT->c_str());
        plotter.add1D(fitHist);

        std::vector<double> coefList;
        for(const auto& syst : systList){ coefList.push_back(w.var(syst.c_str())->getVal());}
        auto * hfit = makeFitHist(*nT,&*h,coefList,upHists,downHists);

        for(unsigned int iS = 0; iS < systList.size(); ++iS){
            systHist(*nT + "_"+systList[iS]+"_Up",hfit,&*h,&*upHists[iS] );
            systHist(*nT + "_"+systList[iS]+"_Down",hfit,&*h,&*downHists[iS] );
        }



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
