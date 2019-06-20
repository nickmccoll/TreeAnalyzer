
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/makePlots.C"
#include "../predTools/CutConstants.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"
#include "HistoPlotting/include/PlotTools.h"

#include "TStyle.h"
#include "TGraphAsymmErrors.h"

using namespace CutConstants;

std::vector<double> bins = {700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,4000};


void makeDataDistributions(const std::string& name, const std::string inputFile,
        const std::string outDir){
    std::string cut =hhRange.cut+"&&"+hbbRange.cut ;

    std::vector<PlotVar> vars;
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,30,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,33,minHHMass,maxHHMass );
    vars.emplace_back("counts" ,std::string(";counts"),"1.0",1,0.0,2.0);

    std::vector<PlotSamp> samps;
    samps.emplace_back(name+"_"+"all","0.76");
    samps.emplace_back(name+"_"+bkgSels[BKG_QG]    ,std::string("(0.76*(")+bkgSels[BKG_QG].cut+"))");
    samps.emplace_back(name+"_"+bkgSels[BKG_LOSTTW],std::string("(0.76*(")+bkgSels[BKG_LOSTTW].cut+"))");
    samps.emplace_back(name+"_"+bkgSels[BKG_MW]    ,std::string("(0.76*(")+bkgSels[BKG_MW].cut+"))");
    samps.emplace_back(name+"_"+bkgSels[BKG_MT]    ,std::string("(0.76*(")+bkgSels[BKG_MT].cut+"))");

    std::vector<PlotSel> sels;
    sels.emplace_back("full"                ,btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);

    std::string outFileName=outDir+name+  "_systTestDistributions.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut, nomW.cut);
}
const float hh_scaleUnc = 0.13/1000;
const float hh_resUnc   = 0.28*1000;

const float hbb_scaleUnc = 0.0015;
const float hbb_resUnc   = 18;

void transform (TH1 * outH, std::function<double(double)> f) {
    for(int iB = 1; iB <= outH->GetNbinsX(); ++iB){
        outH->SetBinContent(iB,outH->GetBinContent(iB)*f(outH->GetBinCenter(iB)) );
    }
    outH->Scale(1.0/outH->Integral());
};

TH1 * getSyst(const TH1 * hN, const TH1 * htWN, bool doHbb, int type){
    TH1 * hO = (TH1*)(doHbb ? htWN : hN)->Clone();

    double scaleUnc = doHbb ? hbb_scaleUnc : hh_scaleUnc;
    double resUnc = doHbb ? hbb_resUnc : hh_resUnc;

    switch(type){
        case 1: //scale +1
            transform(hO, [&](double x){return  (1. + scaleUnc*x);})  ;
            break;
        case 2: //scale -1
            transform(hO, [&](double x){return  1./(1. + scaleUnc*x);})  ;
            break;
        case 3: //res +1
            transform(hO, [&](double x){return  (1. + resUnc/x);})  ;
            break;
        case 4: //res -1
            transform(hO, [&](double x){return  1./(1. + resUnc/x);})  ;
            break;
    }

    if(doHbb){
        hO->Scale(htWN->Integral()/hO->Integral());
        hO->Add(hN);
        hO->Add(htWN,-1);
        hO->Scale(hN->Integral()/hO->Integral());
    } else {
        hO->Scale(hN->Integral()/hO->Integral());
    }
    return hO;
}


TH1 * getNormSyst(const TH1 * hN, const TH1 * htWN,const TH1 * hMt, int type){
    TH1 * hO = (TH1*)hN->Clone();
    TH1 * hStWN = (TH1*)htWN->Clone();
    TH1 * hSMt = (TH1*)hMt->Clone();

    switch(type){
        case 1: // +1
            hStWN->Scale(1.25);
            hSMt->Scale(1.0-0.25*htWN->Integral()/hMt->Integral());
            break;
        case 2: // -1
            hStWN->Scale(0.75);
            hSMt->Scale(1.0+0.25*htWN->Integral()/hMt->Integral());
            break;
    }
    hO->Add(htWN,-1);
    hO->Add(hMt,-1);
    hO->Add(hStWN);
    hO->Add(hSMt);
    hO->Scale(hN->Integral()/hO->Integral());
    return hO;
}


TH1 * procHH( TH1* hIn){
    TH1 * hOut   = PlotTools::rebin(hIn ,bins.size()-1,&bins[0]);
    for(int iB = 1;iB <= hIn->GetNbinsX(); ++iB ){
        double width = hOut->GetBinWidth(iB);
        double count = hOut->GetBinContent(iB);
        double error = hOut->GetBinError(iB);

        hOut->SetBinContent(iB, count*100./width);
        hOut->SetBinError(iB, error*100./width);
    }
    return hOut;
}

void makePlot(TFile * f, bool doHbb,
        const std::string& nom, const std::string& nomTitle,
        const std::string& oth, const std::string& othTitle,
        const std::string& plotTitle){

    std::string var = doHbb ? "hbbMass" : "hhMass";

    Plotter * p = new Plotter();
    TH1 * hN   = 0;
    TH1 * htWN = 0;
    TH1 * hMt  = 0;
    TH1 * hO   = 0;
    f->GetObject((nom+"_all_full_"+var).c_str(),hN);
    f->GetObject((nom+"_losttw_full_"+var).c_str(),htWN);
    f->GetObject((nom+"_mt_full_"+var).c_str(),hMt);
    f->GetObject((oth+"_all_full_"+var).c_str(),hO);
    if(!hN || !htWN || !hO || !hMt) return;
    hN   = (TH1*)hN  ->Clone();
    htWN = (TH1*)htWN->Clone();
    hMt = (TH1*)hMt->Clone();
    hO   = (TH1*)hO  ->Clone();

    TH1 * hSU = getSyst(hN,htWN,doHbb,1);
    TH1 * hSD = getSyst(hN,htWN,doHbb,2);
    TH1 * hRU = getSyst(hN,htWN,doHbb,3);
    TH1 * hRD = getSyst(hN,htWN,doHbb,4);
    TH1 * hNU = getNormSyst(hN,htWN,hMt,1);
    TH1 * hND = getNormSyst(hN,htWN,hMt,2);

    if(!doHbb){
        hN   = procHH(hN );
        hO   = procHH(hO );
        hSU  = procHH(hSU);
        hSD  = procHH(hSD);
        hRU  = procHH(hRU);
        hRD  = procHH(hRD);
        hNU  = procHH(hNU);
        hND  = procHH(hND);

    }


    p->addHistLine(hN,nomTitle.c_str(),kBlack);

    if(doHbb){
        p->addHistLine(hRU,(nomTitle+", lost t/W low #it{m}_{b#bar{b}} +1#sigma").c_str(),kBlue);
        p->addHistLine(hRD,(nomTitle+", lost t/W low #it{m}_{b#bar{b}} -1#sigma").c_str(),kBlue,9);
        p->addHistLine(hNU,(nomTitle+", lost t/W fraction +1#sigma").c_str(),kRed);
        p->addHistLine(hND,(nomTitle+", lost t/W fraction -1#sigma").c_str(),kRed,9);
    } else {
        p->addHistLine(hSU,(nomTitle+", #it{m}_{HH} scale +1#sigma").c_str(),kBlue);
        p->addHistLine(hSD,(nomTitle+", #it{m}_{HH} scale -1#sigma").c_str(),kBlue,9);
        p->addHistLine(hRU,(nomTitle+", #it{m}_{HH} res. +1#sigma").c_str(),kRed);
        p->addHistLine(hRD,(nomTitle+", #it{m}_{HH} res. -1#sigma").c_str(),kRed,9);
    }


    p->addHist(hO,othTitle.c_str(),kGray);


    if(doHbb){
        p->rebin(2);
        p->setYTitle("Normalized events / 12 GeV");
        p->setMinMax(0,0.2);
    } else {
        p->setYTitle("Normalized events / 100 GeV");
//        p->rebin(bins.size()-1,&bins[0]);
        p->setMinMax(0.00001,.5);
    }
    p->setLegendPos(0.5,0.5,0.91,0.93);



    p->setBotMinMax(0.5,1.5);
    p->setYTitleBot("Alt. / " + nomTitle);
    p->normalize();
//    p->draw(true,(plotTitle+"_"+ (doHbb ? "hbb" : "hh" )+".pdf").c_str());
    auto * c = p->drawSplitRatio(0,"stack",false,true,(plotTitle+"_"+ (doHbb ? "hbb" : "hh" )+".pdf").c_str());
    if(!doHbb){
        c->GetPad(1)->SetLogy();
        c->Update();
        c->Print((plotTitle+"_"+ (doHbb ? "hbb" : "hh" )+".pdf").c_str());
    }


}

void makeComparisons(const std::string plotDir ="plots"){
    TFile * f = new TFile((plotDir +"/all.root").c_str(),"read");

    makePlot(f,true,"ttbar_1l_madgraph","LO Madgraph"
            , "ttbar_1l_madgraphNLO","NLO Madgraph","LOvsNLO");
    makePlot(f,false,"ttbar_1l_madgraph","LO Madgraph"
            , "ttbar_1l_madgraphNLO","NLO Madgraph","LOvsNLO");

    makePlot(f,true,"ttbar_madgraph","Madgraph"
            , "ttbar_powheg","Powheg","MadgraphvPowheg");
    makePlot(f,false,"ttbar_madgraph","Madgraph"
            , "ttbar_powheg","Powheg","MadgraphvPowheg");

    makePlot(f,true,"ttbar_powheg","Pythia"
            , "ttbar_powhegHerwig","Herwig","PythiavHerwig");
    makePlot(f,false,"ttbar_powheg","Pythia"
            , "ttbar_powhegHerwig","Herwig","PythiavHerwig");

}

#endif

void plotSystTests(int step=0,std::string treeDir = "trees/bkgCompLMT/",std::string out = "plots/"){

    if(step==0){
        makeDataDistributions("ttbar_1l_madgraphNLO","trees/bkgCompLMT/ttbar_1l_madgraphNLO.root",out);
        makeDataDistributions("ttbar_madgraph"      ,"trees/bkgCompLMT/ttbar_madgraph.root"      ,out);
        makeDataDistributions("ttbar_1l_madgraph"   ,"trees/bkgCompLMT/ttbar_1l_madgraph.root"   ,out);
        makeDataDistributions("ttbar_powheg"        ,"trees/bkgCompLMT/ttbar_powheg.root"        ,out);
        makeDataDistributions("ttbar_powhegHerwig"  ,"trees/bkgCompLMT/ttbar_powhegHerwig.root"  ,out);
    }
    if(step==1){
        makeComparisons();
    }

}
