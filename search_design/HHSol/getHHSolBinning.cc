#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TPaletteAxis.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include <iomanip>

std::vector<int> signalMassBins =
{600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};



enum ParType {P_MEAN, P_GAUSM, P_GAUSS};

std::string flt2Str(double val, int pre = -1){
    std::stringstream stream;
    if(pre > -1) stream << std::fixed << std::setprecision(pre) << val;
    else stream << val;
    return stream.str();
}

std::string doOne(TFile *inFile, ParType type, std::string varName, std::string sumName,
        double min =-2, double max =2, double * setMean =0){
    std::string outString = sumName +"\t";
    double tot = 0;
    std::vector<TObject*> dists;
    Plotter * p = new Plotter();

    TF1 * gausFit=0;


    for(unsigned int iST = 0; iST < 2; ++iST){
        const std::string sigType = iST ? "radion":"blkgrav";
        TGraph * g = new TGraph();
        int nPts =0;
        for(unsigned int iM = 0; iM < signalMassBins.size(); ++iM ){
            const int mass = signalMassBins[iM];
            std::string pref = sigType + TString::Format("_m%i",mass).Data();
            TH1 * h = 0;
            inFile->GetObject((pref +"_"+varName).c_str(),h);
            if(h == 0) continue;
            h->SetYTitle(pref.c_str());
            dists.push_back(h);


            double sumVal = 0;
            TFitResultPtr fitres;
            switch(type){
            case P_MEAN:
                sumVal= h->GetMean();
                break;
            case P_GAUSM :
                fitres = h->Fit("gaus","S","",min,max);
                sumVal= fitres->Parameter(1);
                break;
            case P_GAUSS :
                gausFit = new TF1((sumName+"_"+pref).c_str(),"gaus",min,max);
                if(setMean)gausFit->FixParameter(1,*setMean);
                else gausFit->SetParameter(1,(min+max)/2.);
                gausFit->SetParLimits(0,10,1000);
                gausFit->SetParameter(0,100);

                gausFit->SetParLimits(2,0.1,30);
                gausFit->SetParameter(2,3);
                fitres = h->Fit(gausFit,"SB","",min,max);
                sumVal= fitres->Parameter(2);
                break;
            }

            g->SetPoint(nPts,mass,sumVal);
            nPts++;
        }
        auto sumFit = g->Fit("pol0","S","",750,3600);
        outString += flt2Str(sumFit->Parameter(0))+"\t";
        tot+= sumFit->Parameter(0);
        p->addGraph(g,sigType);
    }
    Drawing::drawAll(dists, (sumName + "_dist").c_str());

    p->draw(false,(sumName).c_str());
    new TCanvas();
    outString += flt2Str(tot/2.);
    return outString;


}

void getParams(const std::string& filename){
    TFile * f = new TFile(filename.c_str());
    std::vector<std::string> outStrings;
    double mean0 = 0;
    double meanW = 41;
    outStrings.push_back(doOne(f,P_GAUSS,"extraMetPerp","metPerp",-200,200));
    outStrings.push_back(doOne(f,P_GAUSS,"extraMetParRelhwwMag","negMETParErr",-0.5,0,&mean0));
    outStrings.push_back(doOne(f,P_GAUSS,"extraMetParRelhwwMag","posMETParErr",0,0.5,&mean0));
    outStrings.push_back(doOne(f,P_GAUSS,"wqqPTRes","jetErr",-0.5,.5));
    outStrings.push_back(doOne(f,P_GAUSM,"lt0p75_virtualSD","onWlnuMeanJet",0,120));
    outStrings.push_back(doOne(f,P_GAUSM,"lt0p75_onshellSD","offWlnuMeanJet",0,120));
    outStrings.push_back(doOne(f,P_GAUSM,"lt0p75_virtualSD_trueW","onWlnuMeanWlnu",65,95));
    outStrings.push_back(doOne(f,P_GAUSM,"lt0p75_onshellSD_trueW","offWlnuMeanWlnu",0,60));
    outStrings.push_back(doOne(f,P_GAUSS,"lt0p75_onshellSD_trueW","offWlnuPosWlnuErr",meanW,60,&meanW));
    outStrings.push_back(doOne(f,P_GAUSS,"lt0p75_onshellSD_trueW","offWnluNegWlnuErr",0,meanW,&meanW));
    outStrings.push_back(doOne(f,P_GAUSS,"lt0p75_virtualSD_trueW","onWlnuWlnuErr",65,95));
    outStrings.push_back(doOne(f,P_GAUSS,"lt0p75_onshellSD_trueH","offWnluHWWErr",100,150));
    outStrings.push_back(doOne(f,P_GAUSS,"lt0p75_virtualSD_trueH","onWlnuHWWErr",100,150));


    std::cout <<std::endl;
    for(unsigned int iP = 0; iP < outStrings.size(); ++iP ){
        std::cout << outStrings[iP]<<std::endl;
    }
    std::cout <<std::endl;



}


void getHWWParams(){
    getParams("tryNewHiggsSol_signal.root");
}







