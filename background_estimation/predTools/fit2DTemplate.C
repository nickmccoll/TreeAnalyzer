
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "InputsHelper.h"

#include <string.h>
#include <regex>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"





class fit2DTemplateAnalyzer {
public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;



    const TH2* makeFitHist(const std::string& name,const TH2 * iHist
            ,const std::vector<double>& coefList, const std::vector<std::unique_ptr<TH2F>>& upHists
            ,const std::vector<std::unique_ptr<TH2F>>& downHists) {

        const double _smoothRegion =1;
        auto smoothStepFunc =[&] (double x) ->double{
            if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
            double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
            return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
        };

        TH2 * hd = (TH2*)iHist->Clone(name.c_str());
        const int nBinsX = hd->GetNbinsX();
        const int nBinsY = hd->GetNbinsY();
        std::vector<double> upHistN;
        std::vector<double> downHistN;
        for(unsigned int iC = 0; iC < coefList.size(); ++iC){
            upHistN.push_back(upHists[iC]->Integral(1,nBinsX,1,nBinsY));
            downHistN.push_back(downHists[iC]->Integral(1,nBinsX,1,nBinsY));
        }

        hd->Scale(1./hd->Integral(1,hd->GetNbinsX(),1,hd->GetNbinsY()) );
        bool anyNegative = false;
        for(int iX = 1; iX <= nBinsX; ++iX)
            for(int iY = 1; iY <= nBinsY; ++iY){
                const double xCen = hd->GetXaxis()->GetBinCenter(iX);
                const double yCen = hd->GetYaxis()->GetBinCenter(iY);
                const double oldV =  hd->GetBinContent(iX,iY);
                double newV =  oldV;
                for(unsigned int iC = 0; iC < coefList.size(); ++iC){
                    const double x = coefList[iC];
                    const double a = 0.5*x;
                    const double b = smoothStepFunc(x);
                    const double upV = upHists[iC]->GetBinContent(iX,iY)/upHistN[iC];
                    const double downV = downHists[iC]->GetBinContent(iX,iY)/downHistN[iC];
                    newV+= a*((upV-downV) + b*(upV+downV - 2*oldV));
                }
                if(newV < 0){
                    newV *= -1;
                    anyNegative = true;
                }
                hd->SetBinContent(iX,iY,newV);
                hd->SetBinError(iX,iY,0);
            }
        if(anyNegative) std::cout << std::endl <<"ALERT -> Negative values!"<<std::endl;
        hd->Scale(1./hd->Integral(1,hd->GetNbinsX(),1,hd->GetNbinsY()) );
        plotter.add2D(hd);
        return hd;
    }

    const TH2 * conditional(const std::string& name, const TH2* iHist,const TH2* nomHist,
            const TH2* systHist, const bool condOnX = true){
        TH2 * oHist = (TH2*)iHist->Clone(name.c_str());
        if(condOnX){
            for(int iX = 1; iX <= oHist->GetNbinsX(); ++iX){
                for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                    if(nomHist->GetBinContent(iX,iY)){
                        const double sf =
                                systHist->GetBinContent(iX,iY)/nomHist->GetBinContent(iX,iY);
                        oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*sf);
                    }
                }
                double yIntegral = oHist->Integral(iX,iX,1,oHist->GetNbinsY());
                if(!yIntegral) { throw std::invalid_argument(
                        "Analyzer::Analyzer() -> Making a conditional slice with no events!!!!!");}
                for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                    oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)/yIntegral);
                }
            }
        } else {
            for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                for(int iX = 1; iX <= oHist->GetNbinsX(); ++iX){

                    if(nomHist->GetBinContent(iX,iY)){
                        const double sf =
                                systHist->GetBinContent(iX,iY)/nomHist->GetBinContent(iX,iY);
                        oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*sf);
                    }
                }
                double xIntegral = oHist->Integral(1,oHist->GetNbinsX(),iY,iY);
                if(!xIntegral) { throw std::invalid_argument(
                        "Analyzer::Analyzer() -> Making a conditional slice with no events!!!!!");}
                for(int iX =1; iX <= oHist->GetNbinsX(); ++iX ){
                    oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)/xIntegral);
                }
            }
        }
        plotter.add2D(oHist);
        return oHist;
    }

    const TH1* conditionalOn(const std::string& name, const TH2* iHist,const TH2* nomHist,
            const TH2* systHist, const bool condOnX = true){
        TH1* oH1 = condOnX
                ?iHist->ProjectionX(name.c_str(),1,iHist->GetNbinsY())
                :iHist->ProjectionY(name.c_str(),1,iHist->GetNbinsX());
        TH1* sH1 = condOnX
                ?systHist->ProjectionX((name +"_syst").c_str(),1,iHist->GetNbinsY())
                :systHist->ProjectionY((name +"_syst").c_str(),1,iHist->GetNbinsX());
        TH1* nH1 = condOnX
                ?nomHist->ProjectionX((name +"_nom").c_str(),1,iHist->GetNbinsY())
                :nomHist->ProjectionY((name +"_nom").c_str(),1,iHist->GetNbinsX());
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

    const TH2* mergeConditionalHistos(const std::string& name, const TH2* iConHist,
            const TH1* iConOnHist, const bool condOnX = true){
        TH2 * oHist = (TH2*)iConHist->Clone(name.c_str());
        if(condOnX){
            for(int iX = 0; iX <= oHist->GetNbinsX(); ++iX){
                const double condOnNorm = iConOnHist->GetBinContent(iX);
                const double originalNorm = iConHist->Integral(iX,iX,0,-1);
                for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                    oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*condOnNorm/originalNorm);
                }
            }
        } else {
            for(int iY =1; iY <= oHist->GetNbinsY(); ++iY ){
                const double condOnNorm = iConOnHist->GetBinContent(iY);
                const double originalNorm = iConHist->Integral(0,-1,iY,iY);
                for(int iX = 0; iX <= oHist->GetNbinsX(); ++iX){
                    oHist->SetBinContent(iX,iY,oHist->GetBinContent(iX,iY)*condOnNorm/originalNorm);
                }
            }
        }
        oHist->Scale(1./oHist->Integral());
        plotter.add2D(oHist);
        return oHist;
    }



    fit2DTemplateAnalyzer(std::string outFileName,std::string arguments )
    {
        ParParser p;
        auto fTN = p.addString("fT","template file name",true);
        auto nT  = p.addString("nT","template histogram base name",true);
        auto s   = p.addString("s" ,"Comma separated list of systematics",true);
        auto sA  = p.addString("sA",
                "Comma separated list of systematics to be included in the output but not fit to."
                ,false,"");
        auto fHN = p.addString("fH","Fitting histogram file name",true);
        auto nH  = p.addString("nH","fitting histogram name",true);
        auto xCy  = p.addBool("xCy",
                "True if x is conditional on y (P(x|y)), otherwise assume P(y|x)");
        p.parse(arguments);

        std::vector<std::string> systList = getList(*s);
        std::vector<std::string> extraSysts;
        if(sA->size()) extraSysts = getList(*sA);


        fT =  TObjectHelper::getFile(*fTN);
        fH =  TObjectHelper::getFile(*fHN);
        auto h = TObjectHelper::getObject<TH2F>(fT,*nT);
        plotter.add1D((TH2*)h->Clone("originalPDF"));

        auto * xAxis = h->GetXaxis();
        auto * yAxis = h->GetYaxis();
        auto xBins =  xAxis->GetXbins()->GetSize()
                ? RooBinning(xAxis->GetNbins(),xAxis->GetXbins()->GetArray())
                : RooBinning(xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax());
        auto yBins =  yAxis->GetXbins()->GetSize()
                ? RooBinning(yAxis->GetNbins(),yAxis->GetXbins()->GetArray())
                : RooBinning(yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax());


        RooWorkspace w("w",false);
        RooArgSet varset;
        RooArgList varlist;
        w.factory("x[0,10000]");
        w.factory("y[0,10000]");
        w.var("x")->setBinning(xBins);
        w.var("y")->setBinning(yBins);
        varset.add(*w.var("x"));
        varset.add(*w.var("y"));
        varlist.add(*w.var("x"));
        varlist.add(*w.var("y"));


        RooDataHist rH("nominalHist","nominalHist",varlist,&*h);
        RooHistPdf rP("nominalPDF","nominalPDF",varlist,rH);
        w.import(rH);
        w.import(rP);
        RooArgList coeffList;
        RooArgList pdfList(*w.pdf("nominalPDF"));

        const std::vector<std::string> systVar = {"Up" , "Down"};

        std::vector<std::unique_ptr<TH2F>> downHists;
        std::vector<std::unique_ptr<TH2F>> upHists;
        for(const auto& syst : systList){
            w.factory((syst+"[-1,1]").c_str());
            coeffList.add(*w.var(syst.c_str()));
            for(const auto& var : systVar){
                auto * histV = &(var== "Up" ? upHists : downHists);
                histV->emplace_back(TObjectHelper::getObject<TH2F>(fT,*nT+"_"+syst+var));
                RooDataHist rSH((syst+var+"Hist").c_str(),(syst+var+"Hist").c_str()
                        ,varlist,&*histV->back());
                RooHistPdf Spdf((syst+var+"PDF").c_str(),(syst+var+"PDF").c_str()
                        ,varset,rSH,0);
                w.import(rSH);
                w.import(Spdf);
                pdfList.add(*w.pdf((syst+var+"PDF").c_str()));
            }

        }

        FastVerticalInterpHistPdf2D interpPDF("interpPDF","interpPDF"
                ,*w.var("x"),*w.var("y"),false,pdfList,coeffList);
        w.import(interpPDF);

        auto inH =TObjectHelper::getObject<TH2F>(fH,*nH);
        RooDataHist fitDataHist((*nH+"DH").c_str(),(*nH+"DH").c_str(),varlist,&*inH);
        w.import(fitDataHist);
        w.pdf("interpPDF")->fitTo(*w.data((*nH+"DH").c_str()),RooFit::SumW2Error(kTRUE));

        TH2 * fitHist =  createTH2FromPDF(w.pdf("interpPDF"),w.var("x"),w.var("y"),
                std::string("debug_") + *nT,"",
                xAxis,yAxis);

        for(const auto& syst : systList){
            std::cout << syst <<" -> "<< w.var(syst.c_str())->getVal()<<std::endl;
        }
        plotter.add1D(fitHist);

        std::vector<double> coefList;
        for(const auto& syst : systList){ coefList.push_back(w.var(syst.c_str())->getVal());}
        std::cout <<"Doing "<< outFileName <<" Fit"<<std::endl;
        auto * hfit = makeFitHist(*nT,&*h,coefList,upHists,downHists);

        for(unsigned int iS = 0; iS < systList.size(); ++iS){
            std::cout << "Trying -> "<<systList[iS]<<std::endl;
            const bool isXSyst = (systList[iS].find("Y") == std::string::npos);
            const bool condOnX = !*xCy;
            const bool isCondOnSyst = (condOnX == isXSyst);
            if(isCondOnSyst){
                auto condOnTH1Up =  conditionalOn(*nT + "_"+systList[iS]+"Up" + "_Debug_1D",
                        hfit,&*h,&*upHists[iS],condOnX);
                auto condOnTH1Down =  conditionalOn(*nT + "_"+systList[iS]+"Down" + "_Debug_1D",
                        hfit,&*h,&*downHists[iS],condOnX);
                mergeConditionalHistos(*nT + "_"+systList[iS]+"Up",
                        hfit,condOnTH1Up,condOnX);
                mergeConditionalHistos(*nT + "_"+systList[iS]+"Down",
                        hfit,condOnTH1Down,condOnX);
            } else {
                auto condOnTH1Nom =  conditionalOn(*nT + "_"+systList[iS]+"_Nom" + "_Debug_1D",
                        hfit,&*h,&*h,condOnX);
                auto condTH2Up = conditional(*nT + "_"+systList[iS]+"Up" + "__Debug_COND2D",
                        hfit,&*h,&*upHists[iS],condOnX);
                auto condTH2Down = conditional(*nT + "_"+systList[iS]+"Down" + "__Debug_COND2D",
                        hfit,&*h,&*downHists[iS],condOnX);
                mergeConditionalHistos(*nT+"_"+systList[iS]+"Up",
                        condTH2Up,condOnTH1Nom,condOnX);
                mergeConditionalHistos(*nT+"_"+systList[iS]+"Down",
                        condTH2Down,condOnTH1Nom,condOnX);
            }
        }

        //Add in extra systematics
        std::vector<std::unique_ptr<TH2F>> downExtraHists;
        std::vector<std::unique_ptr<TH2F>> upExtraHists;
        for(const auto& syst : extraSysts){
            for(const auto& var : systVar){
                auto * histV = &(var== "Up" ? upExtraHists : downExtraHists);
                histV->emplace_back(TObjectHelper::getObject<TH2F>(fT,*nT+"_"+syst+var));
            }
        }
        for(unsigned int iS = 0; iS < extraSysts.size(); ++iS){
            const bool isXSyst = (extraSysts[iS].find("Y") == std::string::npos);
            const bool condOnX = !*xCy;
            const bool isCondOnSyst = (condOnX == isXSyst);
            if(isCondOnSyst){
                auto condOnTH1Up =  conditionalOn(*nT + "_"+extraSysts[iS]+"Up" + "_Debug_1D",
                        hfit,&*h,&*upExtraHists[iS],condOnX);
                auto condOnTH1Down =  conditionalOn(*nT + "_"+extraSysts[iS]+"Down" + "_Debug_1D",
                        hfit,&*h,&*downExtraHists[iS],condOnX);
                mergeConditionalHistos(*nT + "_"+extraSysts[iS]+"Up",
                        hfit,condOnTH1Up,condOnX);
                mergeConditionalHistos(*nT + "_"+extraSysts[iS]+"Down",
                        hfit,condOnTH1Down,condOnX);
            } else {
                auto condOnTH1Nom =  conditionalOn(*nT + "_"+extraSysts[iS]+"_Nom" + "_Debug_1D",
                        hfit,&*h,&*h,condOnX);
                auto condTH2Up = conditional(*nT + "_"+extraSysts[iS]+"Up" + "__Debug_COND2D",
                        hfit,&*h,&*upExtraHists[iS],condOnX);
                auto condTH2Down = conditional(*nT + "_"+extraSysts[iS]+"Down" + "__Debug_COND2D",
                        hfit,&*h,&*downExtraHists[iS],condOnX);
                mergeConditionalHistos(*nT + "_"+extraSysts[iS]+"Up",
                        condTH2Up,condOnTH1Nom,condOnX);
                mergeConditionalHistos(*nT + "_"+extraSysts[iS]+"Down",
                        condTH2Down,condOnTH1Nom,condOnX);
            }
        }

        plotter.write(outFileName);
        fT->Close();
        fH->Close();
    }

    std::vector<std::unique_ptr<TH2F> > pdfHistos;
    std::vector<RooRealVar> plottingVars;
    std::vector<RooRealVar> systVars;
    std::vector<RooDataHist> dataHists;
    std::vector<RooHistPdf> pdfs;


    std::vector<std::string> getList(const std::string& inList){
        std::vector<std::string> systList(
                std::sregex_token_iterator(inList.begin(), inList.end(), std::regex(","), -1),
                std::sregex_token_iterator());
        return systList;
    }

    HistGetter plotter;
    TFile * fT = 0;
    TFile * fH = 0;
    std::shared_ptr<std::vector<double>>       xb;
    std::shared_ptr<std::vector<double>>       yb;

};

#endif

void fit2DTemplate(std::string outFileName,std::string arguments){
    fit2DTemplateAnalyzer a(outFileName, arguments);
}
