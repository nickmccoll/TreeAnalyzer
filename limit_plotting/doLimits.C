#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLatex.h"
#include "../background_estimation/predTools/CutConstants.h"
#include "HistoPlotting/include/StyleInfo.h"

using namespace CutConstants;

const float minX = 700;
const float maxX = 3100;
const float minY = 0.5;
const float maxY = 10000;
const bool doLog = true;
std::string titleX = "m_{X} [GeV]";
std::string titleY = "#sigma#bf{#it{#Beta}}(X #rightarrow HH) [fb]";
std::vector<std::string> partLabels = {"Spin-0 X","Spin-2 X"};
std::vector<std::string> fileLabel = {"radion","blkgrav"};
std::string lumiText = "35.9 fb^{-1} (13 TeV)";
const bool prelim = true;
const float limitScale = 1000;




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
    }

    auto c = new TCanvas();
    auto frame=c->DrawFrame(minX,minY,maxX,maxY);
    frame->GetXaxis()->SetTitle(titleX.c_str());
    frame->GetXaxis()->SetTitleOffset(0.9);
    frame->GetXaxis()->SetTitleSize(0.05);

    frame->GetYaxis()->SetTitle(titleY.c_str());
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.15);


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

    auto leg = new TLegend(0.60,0.66,0.93,0.83,"","NDC");
    leg->SetBorderSize(0);

    if(!blind)leg->AddEntry(bandObs, "Observed");
    leg->AddEntry(bandExp, "Median expected","L");
    leg->AddEntry(band68, "68% expected","F");
    leg->AddEntry(band95, "95% expected","F");


    c->cd();
    frame->Draw();
    band95->Draw("3same");
    band68->Draw("3same");
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

    addText(partLabels[sig],0.03,0.61,0.88);
    addText("95% CL upper limits",0.03,0.61,0.84);




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
