#if !defined(__CINT__) || defined(__MAKECINT__)
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooProdPdf.h"
#include "TRandom3.h"
#include "TreeAnalyzer/background_estimation/predTools/DataCardMaker.h"
#include <sstream>

#include <string>


//using ASTypes::flt2Str;
//using ASTypes::flt2Str;



const std::string mR = "MR";
const std::string mJ = "MJ";
const std::string mS = "MH";
const double lowM = 700;
const double highM = 4000;
const int nBins = 132;



const double exp_a = 100;
//const double exp_b = -0.002;
const double exp_b = -0.009;
const double hh_scaleUnc = 0.02;
const double hh2_scaleUnc = 0.005*0.005;


TRandom3 * rnd;

std::string mkInput(const std::string cardname) {
    const std::string fn = std::string("cardInput_") + cardname+".root";
    TFile *f = new TFile(fn.c_str(),"recreate");

    TH1 * hD    = new TH1F("dataobs",";MR",nBins,lowM,highM);
    TH1 * hO    = new TH1F("mc_dist",";MR",nBins,lowM,highM);
//    TH1 * hO2    = new TH1F("mc_dist2",";MR",nBins,lowM,highM);
    TH1 * hM    = new TH1F("model",";MR",nBins,lowM,highM);
    TH1 * hUP   = new TH1F("model_PTUp",";MR",nBins,lowM,highM);
    TH1 * hDOWN = new TH1F("model_PTDown",";MR",nBins,lowM,highM);

    TH1 * hUP2   = new TH1F("model_PT2Up",";MR",nBins,lowM,highM);
    TH1 * hDOWN2 = new TH1F("model_PT2Down",";MR",nBins,lowM,highM);



    for(unsigned int iB = 1; iB <= nBins; ++iB){
        const double x1 = hD->GetBinLowEdge(iB);
        const double x2 = x1+ hD->GetBinWidth(iB);
        const double x = hD->GetBinCenter(iB);
        double pIntegral = exp_a*( std::exp(exp_b*x2)-std::exp(exp_b*x1))/exp_b;

//        double pIntegral = exp_a*std::exp(exp_b*x)*(x2-x1);
        double dI = pIntegral;//*3.7/(1+hh2_scaleUnc*x*x);
//        if(x>1200) dI *=  (1+1200*1200)/(1+x*x);

        double dv = rnd->Poisson(dI);

        hD->SetBinContent(iB,dv);
        hD->SetBinError(iB,std::sqrt(dv));
        hO->SetBinContent(iB,pIntegral);
        hM->SetBinContent(iB,pIntegral);


        hUP->SetBinContent(iB,pIntegral*(1. + hh_scaleUnc*x));
        hDOWN->SetBinContent(iB,pIntegral/(1. + hh_scaleUnc*x));

        hUP2->SetBinContent(iB,pIntegral*(1. + hh2_scaleUnc*x*x));
        hDOWN2->SetBinContent(iB,pIntegral/(1. + hh2_scaleUnc*x*x));
    }

    hUP->Scale(1./hUP->Integral());
    hDOWN->Scale(1./hDOWN->Integral());

    hUP2->Scale(1./hUP2->Integral());
    hDOWN2->Scale(1./hDOWN2->Integral());

    hM->Scale(1./hM->Integral());

    hD   ->Write();
    hO   ->Write();
    hM   ->Write();
    hUP  ->Write();
    hDOWN->Write();
    hUP2  ->Write();
    hDOWN2->Write();

    f->Close();
    delete f;
    return fn;
}

void mkCard (const std::string cardname) {
    const std::string iFN = mkInput(cardname);

    auto card = DataCardMaker("a","b" ,"c",1,"cardname");
    // saved card will be: datacard_CARDNAME_a_b_c.txt
    card.addVar(mR,(lowM+highM)/2,lowM,highM,false);
    card.addVar(mS,(lowM+highM)/2,true);

    card.addFixedYieldFromFile("bkg",1,iFN,"mc_dist");
    card.addSystematic("yield","lnN",{{"sig",1.1}});
    card.addSystematic("bkg_norm","lnN",{{"bkg",1.5}});
//    card.addParamSystematic("bkg_PT" ,0,1);
    card.addParamSystematic("bkg_PT2" ,0.5,1);

    PDFAdder::InterpSysts KDESysts;
//    KDESysts.addSyst("PT",{{"bkg_PT","1"  }});
    KDESysts.addSyst("PT2",{{"bkg_PT2","1"  }});
    card.addHistoShapeFromFile("bkg",{mR}, iFN,"model",KDESysts,false,0,"");

    card.addVar("sig_mean",2000,true);
    card.addVar("sig_sig",2000*.06,true);
    RooGaussian modelS("sig_cardname_a_b_c","sig_cardname_a_b_c",*card.w->var(mR.c_str()),
            *card.w->var("sig_mean"),*card.w->var("sig_sig"));
    card.w->import(modelS);
    card.addFixedYield("sig",0,10);


    card.importBinnedData(iFN,"dataobs",{mR});
    card.makeCard();

}

std::string mkInput2D(const std::string cardname) {
    const std::string fn = std::string("cardInput_") + cardname+".root";
    TFile *f = new TFile(fn.c_str(),"recreate");

    TH2 * hD    = new TH2F("dataobs",";MR"     ,nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hO    = new TH2F("mc_dist",";MR"     ,nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hM    = new TH2F("model",";MR"       ,nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hUPX   = new TH2F("model_PTXUp",";MR"  ,nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hDOWNX = new TH2F("model_PTXDown",";MR",nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hUPY   = new TH2F("model_PTYUp",";MR"  ,nBins,lowM,highM,nBins,lowM,highM);
    TH2 * hDOWNY = new TH2F("model_PTYDown",";MR",nBins,lowM,highM,nBins,lowM,highM);

    for(unsigned int iX = 1; iX <= nBins; ++iX){
        const double x1 = hD->GetXaxis()->GetBinLowEdge(iX);
        const double x2 = x1+ hD->GetXaxis()->GetBinWidth(iX);
        const double x = (x2+x1)/2.;
        const double xw = hD->GetXaxis()->GetBinWidth(iX);
        for(unsigned int iY = 1; iY <= nBins; ++iY){
            const double y1 = hD->GetYaxis()->GetBinLowEdge(iY);
            const double y2 = y1+ hD->GetYaxis()->GetBinWidth(iY);
            const double y = (y2+y1)/2.;
            const double yw = hD->GetYaxis()->GetBinWidth(iY);
            const double p = yw*xw*exp_a*exp_a*std::exp(exp_b*x)*std::exp(exp_b*y);
            const double dv = rnd->Poisson(p);

        hD->SetBinContent(iX,iY,dv);
        hD->SetBinError(iX,iY,std::sqrt(dv));
        hO->SetBinContent(iX,iY,p);
        hM->SetBinContent(iX,iY,p);
        hUPX->SetBinContent(iX,iY,p*(1. + hh_scaleUnc*x));
        hDOWNX->SetBinContent(iX,iY,p/(1. + hh_scaleUnc*x));
        hUPY->SetBinContent(iX,iY,p*(1. + hh_scaleUnc*y));
        hDOWNY->SetBinContent(iX,iY,p/(1. + hh_scaleUnc*y));
    }}

    hUPX->Scale(1./hUPX->Integral());
    hDOWNX->Scale(1./hDOWNX->Integral());
    hUPY->Scale(1./hUPY->Integral());
    hDOWNY->Scale(1./hDOWNY->Integral());
    hM->Scale(1./hM->Integral());

    hD   ->Write();
    hO   ->Write();
    hM   ->Write();
    hUPX  ->Write();
    hDOWNX->Write();
    hUPY  ->Write();
    hDOWNY->Write();

    f->Close();
    delete f;
    return fn;
}

void mkCard2D (const std::string cardname) {
    const std::string iFN = mkInput2D(cardname);

    auto card = DataCardMaker("a","b" ,"c",1,"cardname");
    // saved card will be: datacard_CARDNAME_a_b_c.txt
    card.addVar(mR,(lowM+highM)/2,lowM,highM,false);
    card.addVar(mJ,(lowM+highM)/2,lowM,highM,false);
    card.addVar(mS,(lowM+highM)/2,true);

    card.addFixedYieldFromFile("bkg",1,iFN,"mc_dist");
    card.addSystematic("yield","lnN",{{"sig",1.1}});
    card.addSystematic("bkg_norm","lnN",{{"bkg",1.5}});
    card.addParamSystematic("bkg_PTX" ,0,1);
    card.addParamSystematic("bkg_PTY" ,0,1);

    PDFAdder::InterpSysts KDESysts;
    KDESysts.addSyst("PTX",{{"bkg_PTX","1"  }});
    KDESysts.addSyst("PTY",{{"bkg_PTY","1"  }});
    card.addHistoShapeFromFile("bkg",{mJ,mR}, iFN,"model",KDESysts,false,0,"");




    card.addVar("sig_meanX",2000,true);
    card.addVar("sig_sigX",2000*.06,true);
    RooGaussian modelSX("sig_cardname_a_b_cX","sig_cardname_a_b_cX",*card.w->var(mJ.c_str()),
            *card.w->var("sig_meanX"),*card.w->var("sig_sigX"));
    card.w->import(modelSX);
    card.addVar("sig_meanY",2000,true);
    card.addVar("sig_sigY",2000*.06,true);
    RooGaussian modelSY("sig_cardname_a_b_cY","sig_cardname_a_b_cY",*card.w->var(mR.c_str()),
            *card.w->var("sig_meanY"),*card.w->var("sig_sigY"));
    card.w->import(modelSY);
    RooProdPdf prodP("sig_cardname_a_b_c", "sig_cardname_a_b_cX*sig_cardname_a_b_cY",
            RooArgSet(*card.w->pdf("sig_cardname_a_b_cX"), *card.w->pdf("sig_cardname_a_b_cY")) );
    card.w->import(prodP);
    card.addFixedYield("sig",0,10);


    card.importBinnedData(iFN,"dataobs",{mJ,mR});
    card.makeCard();

}


#endif






void testCombineGOF(){
    rnd = new TRandom3(1234);
//    mkCard2D("test");
    mkCard("test");
}
