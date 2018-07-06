
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "TRandom3.h"
#include <string.h>



class StatTesterAnalyzer {
public:

    StatTesterAnalyzer(const TH1D* inputModel, const TH1D* inputData, const unsigned int nToys, const bool saveToys, std::string outName = "")
{

        const auto * model = inputModel->GetArray();
        const auto * data  = inputData->GetArray();
        const int nBins    = inputModel->GetSize();

        double ts_nom_sa = calcSA(nBins,model,data);
        double ts_nom_ks = calcKS(nBins,model,data);
        std::vector<double> ts_sa; ts_sa.reserve(nToys);
        std::vector<double> ts_ks; ts_ks.reserve(nToys);

        std::unique_ptr<TRandom3> rand(new TRandom3(1234));

        for(unsigned int iT = 0; iT < nToys; ++iT){
            auto toyData = generatePseudoData(nBins,model,&*rand);
            ts_sa.push_back(calcSA(nBins,model,&toyData[0]));
            ts_ks.push_back(calcKS(nBins,model,&toyData[0]));
            if(saveToys){
                auto toyH = makeTH1(std::string("toy_")+ASTypes::int2Str(iT),nBins,&toyData[0],inputData );
                plotter.add1D(toyH);
            }
        }

        std::sort(ts_ks.begin(), ts_ks.end(), [](const double a, const double b){return a < b;});
        processTS(ts_nom_sa,ts_sa,"saturated");
        processTS(ts_nom_ks,ts_ks,"ks");
//        std::cout << outName<<std::endl;
        if(outName.size())plotter.write(outName);
}

    void processTS(const double nominal_ts,std::vector<double>& toy_ts, const std::string hName){
        std::sort(toy_ts.begin(), toy_ts.end(), [](const double a, const double b){return a < b;});

        double nToys = toy_ts.size();
        double total     = 0;
        double quant_5   = -1;
        double quant_159 = -1;
        double quant_50  = -1;
        double quant_841 = -1;
        double quant_95  = -1;

        double width = (toy_ts.back() - toy_ts.front())/nToys;
//        TH1 * h = new TH1D(hName.c_str(),";test statistic",toy_ts.size(), toy_ts[0] -width , toy_ts[toy_ts.size()-1]+width);


        for(unsigned int iT = 0; iT < toy_ts.size(); ++iT){
            const double ts = toy_ts[iT];
//            h->Fill(ts);
            double frac = double(iT)/nToys;
            total += ts/nToys;
            if(quant_5 < 0 && frac >= 0.05 ) quant_5  = ts;
            if(quant_159 < 0 && frac >= 0.159 ) quant_159  = ts;
            if(quant_50 < 0 && frac >= 0.5 ) quant_50  = ts;
            if(quant_841 < 0 && frac >= 0.841 ) quant_841  = ts;
            if(quant_95 < 0 && frac >= 0.95 ) quant_95  = ts;
        }

        std::cout <<hName <<" -> Data value: "<< nominal_ts <<" ";
        std::cout <<"Toy values mean(5%,15.9%,50%,84.1%,95%) : "<< total <<"("<<  quant_5<<","<<quant_159<<","<<quant_50<<","<<quant_841<<","<<quant_95<<")"<<std::endl;

//        plotter.add1D(h);
    }

    TH1D* makeTH1(const std::string& name, const unsigned int nBins, const double * input, const TH1D * refHist){
        TH1D* h((TH1D*)(refHist->Clone(name.c_str())));
        h->Reset("M");
        h->SetDirectory(0);
        for(unsigned int iB = 1; iB + 1  < nBins; ++iB){
            h->SetBinContent(iB,input[iB]);
        }
        return h;
    }

    std::unique_ptr<double[]> generatePseudoData(const unsigned int nBins, const double * input, TRandom* rand ){
        std::unique_ptr<double[]> pd (new double[nBins]() );
        for(unsigned int iB = 1; iB + 1  < nBins; ++iB){
            pd[iB] = rand->Poisson(input[iB]);
        }
        return pd;
    }

    double getIntegral(const unsigned int nBins, const double *input){
        double total = 0;
        for(unsigned int iB =1 ; iB + 1 < nBins; ++iB){
            total += input[iB];
        }
        return total;
    }

    double calcSA(const unsigned int nBins, const double *inputModel, const double *inputData){
        double total = 0;
        for(unsigned int iB =1 ; iB + 1  < nBins; ++iB){
            //            std::cout << total <<" :: "<< iB <<"  "<<inputModel[iB] <<"  "<< inputData[iB] <<std::endl;
            total += inputModel[iB];
            if(inputData[iB]) total += inputData[iB]*std::log(inputData[iB]/inputModel[iB]) - inputData[iB];
        }
//        std::cout << total*2<<std::endl;
        return total*2;
    }

    double calcKS(const unsigned int nBins, const double *inputModel, const double *inputData){
        auto modInt = getIntegral(nBins,inputModel);
        auto datInt = getIntegral(nBins,inputData);

        double runningModel = 0;
        double runningData = 0;
        double maxTS = 0;
        for(unsigned int iB =1 ; iB+1 < nBins; ++iB){
            runningModel += inputModel[iB]/modInt;
            runningData += inputData[iB]/datInt;
            double result = std::fabs(runningModel-runningData);
            if(result>maxTS) maxTS = result;
        }
        return maxTS;
    }


    HistGetter plotter;


};

#endif

void StatTester(const TH1D* inputModel, const TH1D* inputData, const int nToys, const bool saveToys, std::string outName = ""){
    StatTesterAnalyzer(inputModel, inputData, nToys, saveToys, outName);
}
