
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "makeJSON.C"
#include <string.h>
#include <regex>
#include "TGraphErrors.h"





class fitBKGRatioAnalyzer {
public:

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
        bool anyNegative = false;
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
                if(newV < 0){
                         newV *= -1;
                         anyNegative = true;
                     }
                hd->SetBinContent(iX,newV);
                hd->SetBinError(iX,0);
            }
        if(anyNegative) std::cout << std::endl <<"ALERT -> Negative values!"<<std::endl;
        hd->Scale(1./hd->Integral(1,hd->GetNbinsX()) );
        plotter.add1D(hd);
        return hd;
    }


    fitBKGRatioAnalyzer(std::string outFileName,std::string arguments )
    {
        ParParser p;
        auto nF  = p.addString("nF","input file name",true);
        auto nH  = p.addString("nH","histogram base name",true);
        auto nN  = p.addString("nN" ,"Numerator prefix",true);
        auto nD  = p.addString("nD","Denominator prefix",true);
        auto b   = p.addVFloat("b","graph binning",true);
        auto lb  = p.addBool("lb","histogram list binning (explicit bin edges)");

        auto gs    = p.addString("g" ,"Function to fit like pol3 or laur3",true);
        auto fv    = p.addString("fVar","fit x var name ",false,"XX");
        auto fMin  = p.addFloat("fMin","fit minimum x",true);
        auto fMax  = p.addFloat("fMax","fit maximum x",true);


        p.parse(arguments);

        if((*lb && b-> size()  < 2) || (!*lb &&  b->size() != 3)
        )                     throw std::invalid_argument("fitBKGRatioAnalyzer::fitBKGRatioAnalyzer() -> Bad parsing");
        double* xbb = &b->at(0);
        const int nXb = b->size() -1;

        fT =  TObjectHelper::getFile(*nF);

        auto nh = TObjectHelper::getObject<TH1F>(fT,*nN + "_"+ *nH);
        auto dh = TObjectHelper::getObject<TH1F>(fT,*nD + "_"+ *nH);

        std::string xitle = std::string(";")+ dh->GetXaxis()->GetTitle();

        TH1 * numH = 0;
        TH1 * denH = 0;
        TH1 * totH = 0;
        TH1 * avgH = 0;
        if(*lb){
            numH   =  plotter.getOrMake1D("numH",xitle.c_str(),nXb,xbb);
            denH   =  plotter.getOrMake1D("denH",xitle.c_str(),nXb,xbb);
            avgH   =  plotter.getOrMake1D("avgH",xitle.c_str(),nXb,xbb);
            totH   =  plotter.getOrMake1D("totH",xitle.c_str(),nXb,xbb);
        } else {
            numH   =   plotter.getOrMake1D("numH",xitle.c_str(),xbb[0],xbb[1],xbb[2]);
            denH   =   plotter.getOrMake1D("denH",xitle.c_str(),xbb[0],xbb[1],xbb[2]);
            avgH   =   plotter.getOrMake1D("avgH",xitle.c_str(),xbb[0],xbb[1],xbb[2]);
            totH   =   plotter.getOrMake1D("totH",xitle.c_str(),xbb[0],xbb[1],xbb[2]);
        }

        auto fillHistogram =[](TH1 * inH, TH1* outH, TH1 * avgH = 0, TH1 * totH=0){
            const int nOutBins = outH->GetNbinsX();
            for(int iB = 1; iB <=inH->GetNbinsX(); ++iB ){
                const float val = inH->GetBinContent(iB);
                const float err = inH->GetBinError(iB);
                const float xval = inH->GetBinCenter(iB);
                const int oB = outH->FindFixBin(xval);
                if(oB < 1 || oB > nOutBins) continue;
                const float outV = outH->GetBinContent(oB);
                outH->SetBinContent(oB,outH->GetBinContent(oB) + val);
                (*outH->GetSumw2())[oB]+=err*err;

                if(avgH) {
                    avgH->SetBinContent(oB,avgH->GetBinContent(oB) + xval*val);
                    totH->SetBinContent(oB,totH->GetBinContent(oB) + val);
                }
            }
        };

        fillHistogram(&*nh,numH);
        fillHistogram(&*dh,denH,avgH,totH);





        std::unique_ptr<TGraphErrors> ratGraph(new TGraphErrors);
        std::unique_ptr<TGraphErrors> ratPlusOneGraph(new TGraphErrors);
        ratGraph->SetName("RATIO");
        ratPlusOneGraph->SetName("RATIOPlusOne");

        int curPt = 0;
        for(int iB = 1; iB <= denH->GetNbinsX(); ++iB){
            if(denH->GetBinContent(iB) <= 0) continue;
            if(numH->GetBinContent(iB) <= 0) continue;
            const float xV = avgH->GetBinContent(iB)/totH->GetBinContent(iB);
            const float yV = numH->GetBinContent(iB)/denH->GetBinContent(iB);
            const float yE = yV*std::sqrt(
                    numH->GetBinError(iB)*numH->GetBinError(iB)/( numH->GetBinContent(iB)* numH->GetBinContent(iB))
                    + denH->GetBinError(iB)*denH->GetBinError(iB)/( denH->GetBinContent(iB)* denH->GetBinContent(iB))
            );
            ratGraph->SetPoint(curPt,xV,yV);
            ratGraph->SetPointError(curPt,0.0,yE);

            ratPlusOneGraph->SetPoint(curPt,xV,yV+1);
            ratPlusOneGraph->SetPointError(curPt,0.0,yE);

            curPt++;
        }

        plotter.addGraph(&*ratGraph);
        plotter.addGraph(&*ratPlusOneGraph);
        plotter.write(outFileName+"_inputDebug.root");
        fT->Close();
        std::string args = "-i " + outFileName+"_inputDebug.root "
                + " -g RATIO:" + *gs +" "
                + " -var "+ *fv + " -minX "+ flt2Str(*fMin)+ " -maxX "+ flt2Str(*fMax);

        MakeJSON(outFileName,args);


    }

    HistGetter plotter;
    TFile * fT = 0;

};

#endif

void fitBKGRatio(std::string outFileName,std::string arguments){
    fitBKGRatioAnalyzer a(outFileName, arguments);
}
