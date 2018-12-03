#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLatex.h"
#include "../background_estimation/predTools/CutConstants.h"
#include "HistoPlotting/include/StyleInfo.h"

using namespace CutConstants;

const float minX = 700;
const float maxX = 3600;
const float minY = 0.5;
const float maxY = 10000;
const bool doLog = true;
    std::string titleX = "#it{m}_{X} [GeV]";
    std::string titleY = "#sigma#bf{#it{#Beta}}(X #rightarrow HH) [fb]";
std::vector<std::string> partLabels = {"Spin-0 X","Spin-2 X"};
std::vector<std::string> fileLabel = {"radion","blkgrav"};
std::string lumiText = "35.9 fb^{-1} (13 TeV)";
const bool prelim = false;
const float limitScale = 1000;

//https://github.com/CrossSectionsLHC/WED

std::vector<std::pair<float,float>> radionXSec = {
        {800  ,   1.36e+02},
        {850  ,   1.13e+02},
        {900  ,   9.34e+01},
        {950  ,   7.82e+01},
        {1000 ,   6.55e+01},
        {1050 ,   5.60e+01},
        {1100 ,   4.78e+01},
        {1150 ,   4.09e+01},
        {1200 ,   3.49e+01},
        {1250 ,   2.98e+01},
        {1300 ,   2.55e+01},
        {1350 ,   2.18e+01},
        {1400 ,   1.86e+01},
        {1450 ,   1.59e+01},
        {1500 ,   1.36e+01},
        {1550 ,   1.19e+01},
        {1600 ,   1.03e+01},
        {1650 ,   9.02e+00},
        {1700 ,   7.86e+00},
        {1750 ,   6.86e+00},
        {1800 ,   5.98e+00},
        {1850 ,   5.26e+00},
        {1900 ,   4.62e+00},
        {1950 ,   4.06e+00},
        {2000 ,   3.57e+00},
        {2050 ,   3.17e+00},
        {2100 ,   2.81e+00},
        {2150 ,   2.49e+00},
        {2200 ,   2.21e+00},
        {2250 ,   1.96e+00},
        {2300 ,   1.74e+00},
        {2350 ,   1.54e+00},
        {2400 ,   1.37e+00},
        {2450 ,   1.21e+00},
        {2500 ,   1.08e+00},
        {2550 ,   9.60e-01},
        {2600 ,   8.57e-01},
        {2650 ,   7.65e-01},
        {2700 ,   6.83e-01},
        {2750 ,   6.10e-01},
        {2800 ,   5.45e-01},
        {2850 ,   4.86e-01},
        {2900 ,   4.34e-01},
        {2950 ,   3.88e-01},
        {3000 ,   3.46e-01},
        {3050 ,   3.11e-01},
        {3100 ,   2.79e-01},
        {3150 ,   2.50e-01},
        {3200 ,   2.25e-01},
        {3250 ,   2.02e-01},
        {3300 ,   1.81e-01},
        {3350 ,   1.62e-01},
        {3400 ,   1.46e-01},
        {3450 ,   1.31e-01},
        {3500 ,   1.17e-01}
};

std::vector<std::pair<float,float>> bgXsec = {
{800     ,1.82e+00},
{850     ,1.34e+00},
{900     ,9.80e-01},
{950     ,7.45e-01},
{1000    ,5.67e-01},
{1050    ,4.51e-01},
{1100    ,3.58e-01},
{1150    ,2.85e-01},
{1200    ,2.27e-01},
{1250    ,1.80e-01},
{1300    ,1.43e-01},
{1350    ,1.14e-01},
{1400    ,9.06e-02},
{1450    ,7.20e-02},
{1500    ,5.73e-02},
{1550    ,4.74e-02},
{1600    ,3.92e-02},
{1650    ,3.24e-02},
{1700    ,2.68e-02},
{1750    ,2.21e-02},
{1800    ,1.83e-02},
{1850    ,1.54e-02},
{1900    ,1.29e-02},
{1950    ,1.08e-02},
{2000    ,9.06e-03},
{2050    ,7.74e-03},
{2100    ,6.61e-03},
{2150    ,5.64e-03},
{2200    ,4.82e-03},
{2250    ,4.11e-03},
{2300    ,3.51e-03},
{2350    ,3.00e-03},
{2400    ,2.56e-03},
{2450    ,2.19e-03},
{2500    ,1.87e-03},
{2550    ,1.62e-03},
{2600    ,1.40e-03},
{2650    ,1.21e-03},
{2700    ,1.05e-03},
{2750    ,9.08e-04},
{2800    ,7.86e-04},
{2850    ,6.80e-04},
{2900    ,5.89e-04},
{2950    ,5.10e-04},
{3000    ,4.41e-04},
{3050    ,3.86e-04},
{3100    ,3.37e-04},
{3150    ,2.95e-04},
{3200    ,2.58e-04},
{3250    ,2.25e-04},
{3300    ,1.97e-04},
{3350    ,1.72e-04},
{3400    ,1.51e-04},
{3450    ,1.32e-04},
{3500    ,1.15e-04}
};


TGraph * getSignalCrossSection(int sig) {
TGraph * gr = new TGraph();

const auto& sigXSec = sig == RADION ? radionXSec : bgXsec;

int pt = 0;
for(const auto & xsec : sigXSec){
    gr->SetPoint(pt,xsec.first, xsec.second);
    pt++;
}



return gr;
}




void go(const bool blind, const int sig, const std::string& inName, const std::string& outName){




    std::vector<std::pair<float,float>> masses ;
    std::vector<std::pair<float,float>> obs    ;
    std::vector<std::pair<float,float>> m2Sigma;
    std::vector<std::pair<float,float>> m1Sigma;
    std::vector<std::pair<float,float>> exp    ;
    std::vector<std::pair<float,float>> p2Sigma;
    std::vector<std::pair<float,float>> p1Sigma;

    TFile * f = new TFile(inName.c_str(),"read");
    TTree * t = 0;
    f->GetObject("limit",t);

    double limit;
    float quantileExpected;
    double mh;
    t->SetBranchAddress("limit",&limit);
    t->SetBranchAddress("quantileExpected",&quantileExpected);
    t->SetBranchAddress("mh",&mh);
    int eventNumber = 0;
    while(t->GetEntry(eventNumber)){
        if(mh < minX || mh > maxX) continue;
        if(quantileExpected<0)
            obs.emplace_back(mh,limit*limitScale);
        else if(quantileExpected>0.02 && quantileExpected<0.03)
            m2Sigma.emplace_back(mh,limit*limitScale);
        else if(quantileExpected>0.15 && quantileExpected<0.17)
            m1Sigma.emplace_back(mh,limit*limitScale);
        else if(quantileExpected>0.49 && quantileExpected<0.51)
            exp.emplace_back(mh,limit*limitScale);
        else if(quantileExpected>0.974 && quantileExpected<0.976)
            p2Sigma.emplace_back(mh,limit*limitScale);
        else if(quantileExpected>0.83  && quantileExpected<0.85)
            p1Sigma.emplace_back(mh,limit*limitScale);
        ++eventNumber;
    }
    f->Close();


    auto chk=[&](std::vector<std::pair<float,float>>& obj, std::string name){
        std::sort(obj.begin(),obj.end(),[](const std::pair<float,float> a, const std::pair<float,float>b){return a.first < b.first;} );
        if(exp.size() == obj.size()) return;
        std::cout << "bad file!!!"<<std::endl;
        throw std::range_error(name);
    };
    if(!blind) chk(obs ,"obs"   );
    chk(m2Sigma,"m2Sigma");
    chk(m1Sigma,"m1Sigma");
    chk(exp    ,"exp    ");
    chk(p2Sigma,"p2Sigma");
    chk(p1Sigma,"p1Sigma");



    auto band68= new TGraphAsymmErrors();
    auto band95= new TGraphAsymmErrors();
    auto bandObs= new TGraph();

    for(unsigned int iM = 0; iM < exp.size(); ++iM){
        float mh = exp[iM].first;
        band68->SetPoint(iM,mh,exp[iM].second);
        band95->SetPoint(iM,mh,exp[iM].second);
        if(!blind) bandObs->SetPoint(iM,mh,obs[iM].second);
        band68->SetPointError(iM,0.0,0.0,exp[iM].second-m1Sigma[iM].second,p1Sigma[iM].second-exp[iM].second);
        band95->SetPointError(iM,0.0,0.0,exp[iM].second-m2Sigma[iM].second,p2Sigma[iM].second-exp[iM].second);
        if(!blind){
            std::cout << mh <<"\t"<< obs[iM].second<<"\n";
        }

    }

    auto c = new TCanvas();
    auto frame=c->DrawFrame(minX,minY,maxX,maxY);
    frame->GetXaxis()->SetTitle(titleX.c_str());
    frame->GetXaxis()->SetTitleOffset(0.9);
    frame->GetXaxis()->SetTitleSize(0.05);

    frame->GetYaxis()->SetTitle(titleY.c_str());
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.15);

    auto sigXSec = getSignalCrossSection(sig);
    sigXSec->SetLineWidth(3);
    sigXSec->SetLineColor(kRed);
//    bandExp->SetLineStyle(7);
    sigXSec->SetMarkerStyle(0);
    sigXSec->SetMarkerColor(kRed);


    auto bandExp = (TGraphAsymmErrors*)band68->Clone();
    bandExp->SetLineWidth(3);
    bandExp->SetLineColor(kBlue);
    bandExp->SetLineStyle(7);
    bandExp->SetMarkerStyle(0);

    band68->SetFillColor(kGreen+1);
    band95->SetFillColor(kOrange);

    bandObs->SetLineWidth(3);
    bandObs->SetLineColor(kBlack);
    bandObs->SetMarkerStyle(20);

    auto leg = new TLegend(0.58,0.66,0.91,0.83,"","NDC");
    leg->SetBorderSize(0);

    leg->AddEntry(sigXSec, sig == RADION ? "Radion (#Lambda_{R}=3 TeV)" : "Bulk graviton (#tilde{k}=0.1)");
    if(!blind)leg->AddEntry(bandObs, "Observed");
    leg->AddEntry(bandExp, "Median expected","L");
    leg->AddEntry(band68, "68% expected","F");
    leg->AddEntry(band95, "95% expected","F");


    c->cd();
    frame->Draw();
    band95->Draw("3same");
    band68->Draw("3same");
    sigXSec->Draw("LX");
    bandExp->Draw("LX");
    if(!blind) bandObs->Draw("PLsame");
    leg->Draw();

    auto addText = [&](const std::string& txt, const float size, const float x, const float y) {
        auto latex = new TLatex;
        latex->SetNDC();
        latex->SetTextSize(size);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(x,y,txt.c_str());
    };

    addText(partLabels[sig],0.03,0.59,0.88);
    addText("95% CL upper limits",0.03,0.59,0.84);




    c->SetLogy(doLog);
    c->Draw();


    StyleInfo::CMS_lumi(c, 10,lumiText,prelim ? "Preliminary" : "", 1.6);
    c->Update();
    c->RedrawAxis();
    c->GetFrame()->Draw();
    c->SaveAs((outName+".pdf").c_str());
    TFile *fo= new TFile((outName+".root").c_str(),"recreate");
    c->Write();
    fo->Close();

}

#endif
void doLimits(bool blind, int sig = RADION, std::string inName = "higgsCombineTest.AsymptoticLimits.root", std::string outName = "limitPlot"){

    go(blind,sig,inName,outName +"_"+ fileLabel[sig]);
}
