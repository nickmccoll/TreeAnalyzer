#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
#include "TRandom3.h"
#include "RooGaussian.h"
using namespace CutConstants;
using namespace ASTypes;

double minX = 0;
double maxX = 100;
double minY = 0;
double maxY = 1000;
int nBinsX = 20;
int nBinsYF = 40;

double MH = 500;

double meanX = 50;
double sigmaX = 10;

double resY = 0.25;

double nBkg = 10.*double(nBinsX*nBinsYF);
double nSig = 50.*2.*4.;

//double binX = 5;
//double binYF = 25;



//int nBinsYC = 21;
//std::vector<double> binsYC ={0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375
//        ,400,425,450,475,500,1000};

int nBinsYC = 30;
std::vector<double> binsYC ={0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000};


//int nBinsYC = 11;
//std::vector<double> binsYC ={0,50,100,150,200,250,300,350,400,450,500,1000};


void makeCard(bool fine) {
    auto card = fine  ? DataCardMaker("fine","std","13TeV",1,"std") : DataCardMaker("coarse","std","13TeV",1,"std");
    card.addVar(MOD_MJ,100,0,1000,false);
    card.addVar(MOD_MR,1000,0,10000,false);
    card.addVar(MOD_MS,2000,true);

    card.addFixedYield("sig",0,nSig);
    card.addFixedYield("bkg",1,nBkg);

    card.addSystematic("sig_norm","lnN",{{"sig",1.0391}});//lumi = 2.5 pdf= 2, PU = 0.5, btagfake=1
    card.addSystematic("bkg_norm"  ,"lnN",{{"bkg",1.25}});
    card.add2DGaussian("sig",MOD_MJ,MOD_MR,"signalJSON.json",
            {},{},{},{},MOD_MS);

    PDFAdder::InterpSysts shapeSysts;
    if(fine){
        card.addHistoShapeFromFile("bkg",{MOD_MJ,MOD_MR},"inputHist.root","fine_bkgHist",shapeSysts);
        card.importBinnedData("inputHist.root","fine_dataHist",{MOD_MJ,MOD_MR});
    } else {
        card.addHistoShapeFromFile("bkg",{MOD_MJ,MOD_MR},"inputHist.root","coarse_bkgHist",shapeSysts);
        card.importBinnedData("inputHist.root","coarse_dataHist",{MOD_MJ,MOD_MR});
    }

    card.makeCard();
}

void makeInputs() {
    TH2 * bkg = new TH2F("fine_bkgHist",";X;Y",nBinsX,minX,maxX,nBinsYF,minY,maxY);
    double bkgWeight = 1./double(nBinsYF*nBinsX);
    for(int iX = 1; iX <= nBinsX; ++iX)
        for(int iY = 1; iY <= nBinsYF; ++iY){
            bkg->SetBinContent(iX,iY,bkgWeight);
        }
    bkg->Scale(1./bkg->Integral());
    TH2 * bkgC = new TH2F("coarse_bkgHist",";X;Y",nBinsX,minX,maxX,nBinsYC,&binsYC[0]);
    for(int iX = 1; iX <= nBinsX; ++iX){
        double xV = bkg->GetXaxis()->GetBinCenter(iX);
        for(int iY = 1; iY <= nBinsYF; ++iY){
            double yV = bkg->GetYaxis()->GetBinCenter(iY);
            bkgC->Fill(xV,yV,bkg->GetBinContent(iX,iY));
        }
    }
    bkgC->Scale(1./bkgC->Integral());

    CJSON sigJSON;
    sigJSON.addEntry("meanMJ", std::string("0+(") + flt2Str(meanX)+ ")/MH^0");
    sigJSON.addEntry("sigmaMJ", std::string("0+(") + flt2Str(sigmaX)+ ")/MH^0");
    sigJSON.addEntry("meanMR", std::string("0+MH"));
    sigJSON.addEntry("sigmaMR", std::string("0+(") + flt2Str(resY)+ ")*MH");
    sigJSON.write("signalJSON.json");
    RooWorkspace * w = new RooWorkspace("w");
    w->factory((MOD_MJ + "["+ASTypes::flt2Str(100)+","+ASTypes::flt2Str(0)+","+ASTypes::flt2Str(1000)+"]").c_str());
    w->factory((MOD_MR + "["+ASTypes::flt2Str(1000)+","+ASTypes::flt2Str(0)+","+ASTypes::flt2Str(10000)+"]").c_str());
    w->factory((MOD_MS + "["+ASTypes::flt2Str(2000)+","+ASTypes::flt2Str(0)+","+ASTypes::flt2Str(10000)+"]").c_str());
    w->var(MOD_MS.c_str())->setVal(MH);
    PDFAdder::addDoubleGaus(w,"signal","PF",MOD_MJ,MOD_MR,"signalJSON.json",{},{},{},{},MOD_MS);
    auto* pdf = w->pdf("signal_PF");

    TH2 * sig = new TH2F("fine_sig",";X;Y",nBinsX,minX,maxX,nBinsYF,minY,maxY);
    for(int iX = 1; iX <= nBinsX; ++iX){
        double xV = sig->GetXaxis()->GetBinCenter(iX);
        w->var(MOD_MJ.c_str())->setVal(xV);
        for(int iY = 1; iY <= nBinsYF; ++iY){
            double yV = sig->GetYaxis()->GetBinCenter(iY);
            w->var(MOD_MR.c_str())->setVal(yV);
            sig->SetBinContent(iX,iY,pdf->getVal());
        }
    }
    sig->Scale(1./sig->Integral());

    TH2 * sigPBkg = new TH2F("fine_sigBkg",";X;Y",nBinsX,minX,maxX,nBinsYF,minY,maxY);
    sigPBkg->Add(sig,nSig);
    sigPBkg->Add(bkg,nBkg);

    TRandom3 * rand = new TRandom3(1234);

    TH2 * data = new TH2F("fine_dataHist",";X;Y",nBinsX,minX,maxX,nBinsYF,minY,maxY);
    for(int iX = 1; iX <= nBinsX; ++iX){
        double xV = data->GetXaxis()->GetBinCenter(iX);
        w->var(MOD_MJ.c_str())->setVal(xV);
        for(int iY = 1; iY <= nBinsYF; ++iY){
            double yV = data->GetYaxis()->GetBinCenter(iY);
            double avg = sigPBkg->GetBinContent(iX,iY);
            double val = rand->Poisson(avg);
            for(int iE = 0;iE < val;++iE)
                data->Fill(xV,yV);
        }
    }
    TH2 * dataC = new TH2F("coarse_dataHist",";X;Y",nBinsX,minX,maxX,nBinsYC,&binsYC[0]);
    for(int iX = 1; iX <= nBinsX; ++iX){
        double xV = bkg->GetXaxis()->GetBinCenter(iX);
        for(int iY = 1; iY <= nBinsYF; ++iY){
            double yV = bkg->GetYaxis()->GetBinCenter(iY);
            double val = data->GetBinContent(iX,iY);
            for(int iE = 0;iE < val;++iE)
                dataC->Fill(xV,yV);
        }
    }




    TFile * f = new TFile("inputHist.root","recreate");
    bkg->Write();
    bkgC->Write();
    sig->Write();
    sigPBkg->Write();
    data->Write();
    dataC->Write();
    f->Close();

}

#endif

void makeTestCard(){
    makeInputs();
    makeCard(true);
    makeCard(false);
//    std::cout <<" <<<<< "<< inreg <<" "<< condSignal <<" "<<signals[insig]<<std::endl;
//    REGION reg = REGION(inreg);
//    if(reg == REG_QGCR) btagCats = qgBtagCats;
//    std::string mainDir = "../../";
//    go(insig,hhFilename,mainDir,reg,!condSignal);
}
