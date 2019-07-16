#include "../predTools/NuisPlotter.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;
std::vector<TObject*> writeables;

void doAllBKG(const std::string& outDir) {
    const auto l = lepCats[LEP_MU];
    const auto b = btagCats[BTAG_M];
    const auto p = purCats[PURE_HP];
    const auto h = hadCuts[HAD_FULL];
    const std::string cat = l +"_"+b+"_"+p +"_"+h;
    const std::string title = l.title +", "+b.title+", "+p.title;


    std::vector<TString> bkgs = {"qg","losttw","mt","mw"};

    TH2 * outH=0;

    for(unsigned int iB = 0; iB < bkgSels.size();++iB){
        TFile * fmc = new TFile("bkgInputs/HHlnujj_"+bkgs[iB]+"_distributions.root");
        TH2 * hmc = 0;
        fmc->GetObject(bkgs[iB]+"_"+cat+"_hbbMass_hhMass",hmc);
        if(hmc==0){
            std::cout << bkgs[iB]+"_"+cat+"_hbbMass_hhMass"<<std::endl;
            continue;
        }
        TFile * ftemp = new TFile("bkgInputs/HHlnujj_"+bkgs[iB]+"_2D_template_debug.root");
        TH2 * hmtemp = 0;
        ftemp->GetObject(bkgs[iB]+"_"+cat,hmtemp);
        if(hmtemp==0)continue;
        hmtemp->Scale(hmc->Integral()/hmtemp->Integral());
        if(outH == 0) outH=(TH2*)hmtemp->Clone();
        else outH->Add(hmtemp);
    }
    outH->Scale(1./outH->Integral());
    TGaxis::SetMaxDigits(3);
    TH2 * newOutH = new TH2F("newComb",";#it{m}_{b#bar{b}} [GeV];#it{m}_{HH} [TeV]; Arbitrary scale",90,30,210,132,0.7,4);
    for(unsigned int iX = 1; iX <= 90; iX++)for(unsigned int iY = 1; iY <= 132; iY++)
        newOutH->SetBinContent(iX,iY,outH->GetBinContent(iX,iY));
    newOutH->GetXaxis()->SetTitleOffset(1.05);
    newOutH->GetYaxis()->SetTitleOffset(1.29);
    newOutH->GetZaxis()->SetTitleOffset(1.29);
    TCanvas * c1 = new TCanvas();
    newOutH->Draw("COLZ");
    gPad->Update();
    c1->SetRightMargin(0.16);
       TPaletteAxis *palette = (TPaletteAxis*)newOutH->GetListOfFunctions()->FindObject("palette");
       float axSz = 0.99-0.945;
       float axSt = 0.85;
       palette->SetX1NDC(axSt);
       palette->SetX2NDC( axSt+axSz);

       gPad->Modified();
       gPad->Update();

    c1->Update();

    StyleInfo::CMS_lumi(c1, 0, "13 TeV", "Simulation Supplementary", 1.6 );
    c1->Update();
    c1->RedrawAxis();
    c1->GetFrame()->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextColor(kWhite);
    latex.DrawLatex(0.15,0.88,title.c_str());

    c1->Print( (outDir +"/allbkg_2DModel.pdf").c_str());
}

void doSig(const int mx, const std::string& mxTitle, const std::string& outDir) {
    const auto l = lepCats[LEP_MU];
    const auto b = btagCats[BTAG_M];
    const auto p = purCats[PURE_HP];
    const auto h = hadCuts[HAD_FULL];
    const std::string cat = l +"_"+b+"_"+p +"_"+h;
    const std::string title = l.title +", "+b.title+", "+p.title;
    TFile * f = new TFile(std::string("signalInputs/HHlnujj_radHH_"+cat+"_2D_fit.root").c_str());
    TH2* c=0;
    f->GetObject(("pdf_m"+int2Str(mx)+"__MJ_MR").c_str(),c);
    TH2 * newOutH = new TH2F("newComb",";#it{m}_{b#bar{b}} [GeV];#it{m}_{HH} [TeV]; Arbitrary scale",90,30,210,132,0.7,4);

    for( int iX = 1; iX <= c->GetNbinsX(); iX++)for( int iY = 1; iY <= c->GetNbinsY(); iY++){
        int oX = newOutH->GetXaxis()->FindFixBin(c->GetXaxis()->GetBinCenter(iX));
        int oY = newOutH->GetYaxis()->FindFixBin(c->GetYaxis()->GetBinCenter(iY)/1000.0);
        if(oX < 1 || oX > 90) continue;
        if(oY < 1 || oY > 132) continue;
        newOutH->SetBinContent(oX,oY,c->GetBinContent(iX,iY));
    }
    TGaxis::SetMaxDigits(3);
    newOutH->Scale(1./newOutH->Integral());
    newOutH->GetXaxis()->SetTitleOffset(1.05);
    newOutH->GetYaxis()->SetTitleOffset(1.29);
    newOutH->GetZaxis()->SetTitleOffset(1.29);
    TCanvas * c1 = new TCanvas();
    newOutH->Draw("COLZ");
    gPad->Update();
    c1->SetRightMargin(0.16);
       TPaletteAxis *palette = (TPaletteAxis*)newOutH->GetListOfFunctions()->FindObject("palette");
       float axSz = 0.99-0.945;
       float axSt = 0.85;
       palette->SetX1NDC(axSt);
       palette->SetX2NDC( axSt+axSz);

       gPad->Modified();
       gPad->Update();

    c1->Update();

    StyleInfo::CMS_lumi(c1, 0, "13 TeV", "Simulation Supplementary", 1.6 );
    c1->Update();
    c1->RedrawAxis();
    c1->GetFrame()->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextColor(kWhite);
    latex.DrawLatex(0.15,0.88,title.c_str());
    latex.DrawLatex(0.15,0.84,mxTitle.c_str());

    c1->Print( (outDir +"/radHH_m"+int2Str(mx)+ "_2DModel.pdf").c_str());
}



void doExtraPlots(){
    TFile * f = new TFile ("out_radion.root","read");
    std::vector<std::string> sigs  = {"m1000","m2000","m3000"};
    std::vector<std::string> sigTs = {"1 TeV X_{spin-0}","2 TeV X_{spin-0}","3 TeV X_{spin-0}"};
    std::vector<std::string> vars  = {"deltaRlW","deltaPhiHH"};
    std::vector<float> maxys  = {0.065,0.13};

    for(unsigned int iV = 0; iV < vars.size(); ++iV){
        Plotter * p = new Plotter();
        for(unsigned int iS = 0; iS < sigs.size(); ++iS){
            TH1 * h = 0;
            f->GetObject((sigs[iS]+"_"+vars[iV]).c_str(),h);
            if(h==0) continue;
            p->addHistLine(h,sigTs[iS].c_str());

        }
        p->normalize();
        p->setCMSLumi(10);
        p->setCMSLumiExtraText("Simulation Supplementary");
        p->setCMSLumiLumiText("13 TeV");
        p->setYTitle("Normalized events / 0.01 units");
        double xV = 0.6;
        double yV = 0.75;
        p->setLegendPos(xV,yV,xV+0.3,yV+0.15);
        p->setMinMax(0,maxys[iV]);
        if(iV==0) p->setXTitle("#DeltaR(lepton,W#rightarrowq#bar{q}')");
        auto * c = p->draw(false,vars[iV]+"_genDist.pdf");
        p->xAxis()->SetTitleOffset(1.05);
        p->yAxis()->SetTitleOffset(1.54);

        c->Print((vars[iV]+"_genDist.pdf").c_str());

    }
}




void plotSupMatPlots(int step = 0, std::string limitBaseName = ""){
//    REGION reg = REGION(inreg);
//
//    std:: string inName =  "bkgInputs" ;
//    auto srList = getSRList(reg);
//    auto srListTitles = getSRListTitles(reg);
//    std::string outName = limitBaseName +"/plots/";
//
//    if(reg == REG_TOPCR){
//        inName =  "bkgInputsTopCR";
//        hhFilename +="_TopCR";
//        limitBaseName +="_TopCR";
//        outName=limitBaseName+"/plots/TopCR_";
//    }
//    else if(reg == REG_QGCR){
//        inName =  "bkgInputsQGCR";
//        hhFilename +="_QGCR";
//        limitBaseName +="_QGCR";
//        outName=limitBaseName+"/plots/QGCR_";
//        btagCats = qgBtagCats;
//
//    }
//    std::string filename = inName +"/"+hhFilename;
//    std::string postFitFilename = limitBaseName +"/postFit.root";

    if(step==0) {
        doAllBKG("supMatPlots");

    }
    if(step==1) {
        doSig(1600,"1.6 TeV X_{spin-0}","supMatPlots");
    }
    if(step==2){
        doExtraPlots();
    }
}







