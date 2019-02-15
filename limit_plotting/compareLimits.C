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
const float minY = 5;
const float maxY = 500;
const bool doLog = true;
    std::string titleX = "#it{m}_{X} [GeV]";
    std::string titleY = "#sigma#bf{#it{#Beta}}(X #rightarrow HH) [fb]";
std::vector<std::string> partLabels = {"Spin-0 X","Spin-2 X"};
std::vector<std::string> fileLabel = {"radion","blkgrav"};
std::string lumiText = "35.9 fb^{-1} (13 TeV)";
const bool prelim = true;
const float limitScale = 1000 * HHtobbVVBF_BUGGY/HHtobbVVBF;

//expected/observed
std::pair<TGraph*,TGraph*> getValues(const std::string& filename, bool blind){

    std::vector<std::pair<float,float>> masses ;
    std::vector<std::pair<float,float>> obs    ;
    std::vector<std::pair<float,float>> m2Sigma;
    std::vector<std::pair<float,float>> m1Sigma;
    std::vector<std::pair<float,float>> exp    ;
    std::vector<std::pair<float,float>> p2Sigma;
    std::vector<std::pair<float,float>> p1Sigma;

    TFile * f = new TFile(filename.c_str(),"read");
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
    auto bandExp= new TGraph();
    std::cout << filename<<std::endl;
    for(unsigned int iM = 0; iM < exp.size(); ++iM){

        float mh = exp[iM].first;
        std::cout << mh<<" "<< exp[iM].second<<std::endl;
        bandExp->SetPoint(iM,mh,exp[iM].second);
        band68->SetPoint(iM,mh,exp[iM].second);
        band95->SetPoint(iM,mh,exp[iM].second);
        if(!blind) bandObs->SetPoint(iM,mh,obs[iM].second);
        band68->SetPointError(iM,0.0,0.0,exp[iM].second-m1Sigma[iM].second,p1Sigma[iM].second-exp[iM].second);
        band95->SetPointError(iM,0.0,0.0,exp[iM].second-m2Sigma[iM].second,p2Sigma[iM].second-exp[iM].second);
        if(!blind){
            std::cout << mh <<"\t"<< obs[iM].second<<"\n";
        }

    }

    return std::make_pair(bandExp,bandObs);
}


void go(const bool blind, const int sig, const std::vector<std::pair<std::string,std::string>>& inFiles, const std::string& outName){




    auto c = new TCanvas();
    auto frame=c->DrawFrame(minX,minY,maxX,maxY);
    frame->GetXaxis()->SetTitle(titleX.c_str());
    frame->GetXaxis()->SetTitleOffset(0.9);
    frame->GetXaxis()->SetTitleSize(0.05);

    frame->GetYaxis()->SetTitle(titleY.c_str());
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.15);


    std::vector<EColor> colors = {kBlack, kBlue, kRed, kGreen};

    auto leg = new TLegend(0.58,0.66,0.91,0.83,"","NDC");
    leg->SetBorderSize(0);

    c->cd();
    frame->Draw();

    for(unsigned int iF = 0; iF < inFiles.size(); ++iF){
        auto limitValues = getValues(inFiles[iF].first,blind);
        auto bandExp = limitValues.first;

        bandExp->SetLineWidth(3);
        bandExp->SetLineColor(colors[iF]);
        bandExp->SetLineStyle(7);
        bandExp->SetMarkerStyle(0);
        leg->AddEntry(bandExp, (inFiles[iF].second +" expected").c_str(),"L");
        c->cd();
        if(!iF)bandExp->Draw("LX");
        else bandExp->Draw("Lsame");
        if(!blind){
            auto bandObs = limitValues.second;
            bandObs->SetLineWidth(3);
            bandObs->SetLineColor(colors[iF]);
            bandObs->SetMarkerStyle(20);
            leg->AddEntry(bandExp, (inFiles[iF].second +" observed").c_str(),"L");
            bandObs->Draw("PLsame");
        }
    }

    c->cd();
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
void compareLimits(bool blind, int sig = RADION, std::string inName = "higgsCombineTest.AsymptoticLimits.root", std::string outName = "limitPlot"){
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std::cout <<std::endl<< "WARNING! We are assuming that you made the inputs with the incorrect efficiency!"<<std::endl;
    std:: vector<std::pair<std::string,std::string>> inputs;
    inputs = {
            {"silepton_rad/higgsCombineTest.AsymptoticLimits.root","1l"},
            {"dilepton_rad/higgsCombineTest.AsymptoticLimits.root","2l"},
            {"combined_rad/higgsCombineTest.AsymptoticLimits.root","1l+2l"}
    };
//    inputs = {
//            {"higgsCombineTest.AsymptoticLimits.testBinning_std.root","Fine binning"},
            //            {"higgsCombineTest.AsymptoticLimits.testBinning_rebinX.root","rebinX"},
            //            {"higgsCombineTest.AsymptoticLimits.testBinning_rebinY.root","rebinY"},
//            {"higgsCombineTest.AsymptoticLimits.testBinning_rebinXY.root","Coarse binning"}
//    };

    go(blind,sig,inputs ,outName +"_"+ fileLabel[sig]);
}
