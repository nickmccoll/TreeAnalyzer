
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "TStyle.h"
#include "TPRegexp.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

#include "../background_estimation/FunctionFitter.C"
#include "../background_estimation/CutConstants.h"

#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/StyleInfo.h"

using namespace TAna;

//const double pts[] = {150,200,250,300,400,500,600,700,800,950,1100,1250,1500,1750,2000,2250};
//const int nPTS  = 15;
const double pts[] = {200,250,300,400,500,600,700,800,950,1100,1250,1500,1750,2000,2250,3000};
const int nPTS  = 15;
const double etas[] = {0,0.8,1.6,2.4};
const int nETAS = 3;

const double rPTS[] = {200,400,600,800,1000,1500,2500};
const int nRPTS  = 6;
const double rEtas[] = {0,0.6,1.4,1.8,2.0, 2.4};
const int nRETAS = 5;

//const double pts[] = {200,10000};
//const int nPTS  = 1;
//const double etas[] = {0,2.4};
//const int nETAS = 1;




class Analyzer : public BaseTreeAnalyzer {
public:
    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt),
    totPT(nETAS,std::vector<float>(nPTS,0)),nEntries(nETAS,std::vector<float>(nPTS,0)),
    totPTr(nRETAS,std::vector<float>(nRPTS,0)),nEntriesr(nRETAS,std::vector<float>(nRPTS,0))
{
        mass = &jetRawSDMass;
}
    void loadVariables() override {
        setBranchAddress("","puWeight"          ,&puWeight        ,true);
        setBranchAddress("","csvCat"            ,&csvCat          ,true);
        setBranchAddress("","bosonMass"         ,&bosonMass       ,true);
        setBranchAddress("","bosonPT"           ,&bosonPT         ,true);
        setBranchAddress("","bosonETA"          ,&bosonETA        ,true);
        setBranchAddress("","drToJet"           ,&drToJet         ,true);
        setBranchAddress("","maxSJDR"           ,&maxSJDR         ,true);
        setBranchAddress("","maxSJETA"          ,&maxSJETA        ,true);
        setBranchAddress("","minSJPT"           ,&minSJPT         ,true);
        setBranchAddress("","genPT"             ,&genPT           ,true);
        setBranchAddress("","genETA"            ,&genETA          ,true);
        setBranchAddress("","jetPT"             ,&jetPT           ,true);
        setBranchAddress("","jetCorrPT"         ,&jetCorrPT       ,true);
        setBranchAddress("","jetETA"            ,&jetETA          ,true);
        setBranchAddress("","jetID"             ,&jetID           ,true);
        setBranchAddress("","jetMass"           ,&jetMass         ,true);
        setBranchAddress("","jetRawSDMass"      ,&jetRawSDMass    ,true);
        setBranchAddress("","jetL23CorrSDMass"  ,&jetL23CorrSDMass,true);
        setBranchAddress("","jetCorrSDMass"     ,&jetCorrSDMass   ,true);
        setBranchAddress("","jetWCorrSDMass"    ,&jetWCorrSDMass  ,true);

    }




    bool runEvent() override {
        if(maxSJDR > 0.8) return false;
        if(maxSJETA > 2.4) return false;
        if(minSJPT  < 20) return false;
        if(jetID < 4) return false;

        {
            static TAxis *ptAxis = new TAxis(nPTS,pts);
            static TAxis *etaAxis = new TAxis(nETAS,etas);

            int ptBin = ptAxis->FindFixBin(bosonPT);
            if(ptBin == 0) return false;
            if(ptBin > ptAxis->GetNbins()) ptBin = ptAxis->GetNbins();
            int etaBin = etaAxis->FindFixBin(std::fabs(bosonETA));
            if(etaBin > etaAxis->GetNbins()) etaBin = etaAxis->GetNbins();

            ptBin -= 1;
            etaBin -= 1;
            auto fill = [&](int ptBin, int etaBin){
                plotter.getOrMake1D(TString::Format("dist_%i_%i",etaBin,ptBin),";H(bb) mass [GeV]", 300,0,300  )->Fill(*mass,puWeight);
                nEntries[etaBin][ptBin] += puWeight;
                totPT[etaBin][ptBin] += puWeight*bosonPT;
            };

            fill(ptBin,etaBin);

            if(pts[ptBin] >= 1100 || ptBin == 0){
                if(etaBin == 1) fill(ptBin,2);
                if(etaBin == 2) fill(ptBin,1);
            }

        }

        if(bosonPT > 0)
        {
            static TAxis *rptAxis = new TAxis(nRPTS,rPTS);
            static TAxis *retaAxis = new TAxis(nRETAS,rEtas);

            int ptBin = rptAxis->FindFixBin(bosonPT);
            if(ptBin == 0) return false;
            if(ptBin > rptAxis->GetNbins()) ptBin = rptAxis->GetNbins();
            int etaBin = retaAxis->FindFixBin(std::fabs(bosonETA));
            if(etaBin > retaAxis->GetNbins()) etaBin = retaAxis->GetNbins();

            if(etaBin >= 3){
                ptBin = std::min(ptBin,4);
            }


            ptBin -= 1;
            etaBin -= 1;
            auto fill = [&](int ptBin, int etaBin){
                plotter.getOrMake1D(TString::Format("resp_%i_%i",etaBin,ptBin),";H(bb) jet response", 500,0,2  )->Fill(jetPT/bosonPT,puWeight);
                nEntriesr[etaBin][ptBin] += puWeight;
                totPTr[etaBin][ptBin] += puWeight*bosonPT;
            };

            fill(ptBin,etaBin);
            //            if(rPTS[ptBin] >= 800){
            //                if(etaBin == 3) fill(ptBin,4);
            //                if(etaBin == 4) fill(ptBin,3);
            //            }
        }



        return true;
    }

    void getResp(TFile* outFile) {
        Plotter * p = new Plotter;
        for(unsigned int iE = 0; iE <  nRETAS; ++iE){
            const float minETA = rEtas[iE];
            const float maxETA = rEtas[iE+1];
            TGraphErrors * graph = new TGraphErrors();
            int nP = 0;
            for(unsigned int iP = 0; iP <  nRPTS; ++iP){
                const float minPT = rPTS[iP];
                const float maxPT = rPTS[iP+1];

                TH1 * h = plotter.get1D(TString::Format("resp_%u_%u",iE,iP),false);
                if(h == 0) continue;
                h->SetDirectory(0);

                h = (TH1*)h->Clone();
                h->Scale(1./h->Integral(0,-1));
                h->SetYTitle(TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA));

                //
                //                new TCanvas(TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA),TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA));
                //                h->Draw();

                graph->SetPoint(nP,totPTr[iE][iP]/nEntriesr[iE][iP],h->GetMean());
                graph->SetPointError(nP,0,h->GetMeanError());
                nP++;
            }

            p->addGraph(graph, TString::Format("%.1f#geq |#eta| <%.1f",minETA,maxETA));
            //            new TCanvas();
            //            graph->Draw("APL");

        }

        p->setXTitle( "H#rightarrowbb particle #it{p}_{T} [GeV]");
        p->setYTitle("< reco #it{p}_{T} / particle #it{p}_{T} >");
        p->setMinMax(0.86,1.025);
        p->setLegendPos(.15,.65,.6,.9);
        p->setCMSLumi(33,"13 TeV","Simulation Preliminary");

        auto * respCan = p->draw(false,"HbbAverageResponse.pdf");
        p->xAxis()->SetTitleOffset(1.10);
        p->yAxis()->SetTitleOffset(1.47);

        gStyle->SetOptFit(0000);
        std::vector<float> as;
        std::vector<float> bs;
        for(Drawing::Drawable1D& h : p->getHists()){
            auto fr =  ((TGraph*) h.obj)->Fit("pol1","S EX0");
            as.push_back(fr->GetParams()[0]);
            bs.push_back(fr->GetParams()[1]);
        }
        respCan->Print("HbbAverageResponse.pdf");

        TCanvas * ptCan = new TCanvas("ptCorrection","ptCorrection");
        TLegend * legend = new TLegend(.5,.5,.9,.75);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        for(unsigned int iE = 0; iE <  nRETAS; ++iE){
            const float a = as[iE];
            const float b = bs[iE];
            const float minETA = rEtas[iE];
            const float maxETA = rEtas[iE+1];


            TF1 * f = new TF1(TString::Format("f_%u",iE),"([0] + [1]*sqrt([2] +[3]*x)/[4])/x",200, iE <= 1 ? 2000 : 1000);
            f->FixParameter(0,-a/(2.*b));
            f->FixParameter(1,1.);
            f->FixParameter(2, a*a);
            f->FixParameter(3, 4.*b);
            f->FixParameter(4, 2.*b);
            f->SetLineColor( StyleInfo::getLineColor(iE));
            f->SetLineWidth( 4);
            legend->AddEntry(f, TString::Format("%.1f#geq |#eta| <%.1f",minETA,maxETA),"L");
            f->Draw(iE ? "SAME" : "");
            f->GetYaxis()->SetRangeUser(1.0,1.2);
            f->GetYaxis()->SetTitleOffset(1.47);
            f->GetXaxis()->SetTitleOffset(1.10);
            f->GetYaxis()->SetTitle("H#rightarrowbb jet #it{p}_{T} correction");
            f->GetXaxis()->SetTitle("H#rightarrowbb reco-jet #it{p}_{T} [GeV]");

            outFile->cd();
            TString title =  TString::Format("jec_eta_%.1f_%.1f",minETA,maxETA);
            title.ReplaceAll(".","p");
            auto *f2 = (TF1*)(f->Clone());
            f2->Write(title);
        }

        legend->Draw();
        StyleInfo::CMS_lumi(ptCan,33,"13 TeV","Simulation Preliminary",1.6);
        ptCan->Print("hbbPTCorrection.pdf");



    }
    void doFits(TFile * outFile){
        auto setup1DFit = [](const TH1* hbbH, bool doExpo, std::vector<std::unique_ptr<FunctionFitter>>& fitters, double minPT, double maxPT){
            std::string pF = "SMJJ";
            auto vN=[&](std::string var)->std::string{return var+pF;};
            fitters.emplace_back(new CBFunctionFitter(hbbH,doExpo,pF,{"MJJ"}));
            auto fitter = &* fitters.back();
            fitter->setVar(vN("mean")     ,110,50,180);
            fitter->setVar(vN("sigma")       ,10,5,30);
            //             fitter->setConst(vN("sigma"),1);
            fitter->setVar(vN("alpha")     ,1 ,0.1,2);
            //            fitter->setConst(vN("alpha"),1);
            fitter->setVar(vN("alpha2")  ,2,0.1,4);
            //             fitter->setConst(vN("alpha2")  ,1);
            fitter->setVar(vN("n")   ,  5  ,1,6);
            fitter->setVar(vN("n2")  ,5,3,20);
            fitter->setConst(vN("n")  ,1);
            fitter->setConst(vN("n2")  ,1);

            if(doExpo){
                fitter->setVar(vN("slope")  ,-1,-10,0);
                fitter->setVar(vN("fE")  ,0.1,0,0.75);
            }


            fitter->w->var("MJJ")->setRange("fit",minPT,maxPT);
            fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("fit"),RooFit::Minos(0),RooFit::NumCPU(8)});
            fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("fit"),RooFit::Minos(1),RooFit::NumCPU(8)});
        };

        Plotter * p = new Plotter;
        Plotter * pm = new Plotter;
        for(unsigned int iE = 0; iE <  nETAS; ++iE){
            std::vector<std::unique_ptr<FunctionFitter>> fitters;
            FunctionParameterPlotter funcPlotter;

            const float minETA = etas[iE];
            const float maxETA = etas[iE+1];
            TGraphErrors * graph = new TGraphErrors();
            TGraphErrors * graphM = new TGraphErrors();
            int nP = 0;
            for(unsigned int iP = 0; iP <  nPTS; ++iP){
                const float minPT = pts[iP];
                const float maxPT = pts[iP+1];

                TH1 * h = plotter.get1D(TString::Format("dist_%u_%u",iE,iP),false);
                if(h == 0) continue;
                h->SetDirectory(0);

                double minFitPT = 70;
                double maxFitPT = 150;
                if(maxPT <= 250){
                    minFitPT = 30;
                }
                else if(maxPT <= 300){
                    minFitPT = 50;
                }
                else if(minPT >= 200 && iE > 0){
                    minFitPT= 50;
                    maxFitPT= 130;
                }

                setup1DFit(h,false,fitters,minFitPT,maxFitPT);
                funcPlotter.addFit(&*fitters.back(),totPT[iE][iP]/nEntries[iE][iP], ASTypes::flt2Str(minPT)+"to"+ASTypes::flt2Str(maxPT));

                h = (TH1*)h->Clone();
                h->Scale(1./h->Integral(0,-1));
                h->SetYTitle(TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA));
                new TCanvas(TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA),TString::Format("pt: %.0f-%.0f, |#eta| %.1f-%.1f",minPT,maxPT,minETA,maxETA));
                h->Draw();
                double fitMean = fitters.back()->getVal("meanSMJJ");
                double fitError = fitters.back()->getError("meanSMJJ");
                graphM->SetPoint(nP,totPT[iE][iP]/nEntries[iE][iP],fitMean);
                graphM->SetPointError(nP,0,fitError);
                graph->SetPoint(nP,totPT[iE][iP]/nEntries[iE][iP],125.0/fitMean);
                graph->SetPointError(nP,0,(125.0/(fitMean*fitMean))*fitError);
                std::cout << nP <<" "<< totPT[iE][iP]/nEntries[iE][iP] <<" "<< fitMean<<std::endl;
                nP++;
            }
            funcPlotter.write(std::string("debug_")+ASTypes::int2Str(iE)+"_MJJ_fit.root");
            p->addGraph(graph, TString::Format("%.1f#geq |#eta| <%.1f",minETA,maxETA),-1,1,4,20,1,true,true,false,"P E 0 Z L");

            outFile->cd();
            TString title =  TString::Format("sdMassCorr_eta_%.1f_%.1f",minETA,maxETA);
            title.ReplaceAll(".","p");
            auto *f2 = (TGraphErrors*)(graph->Clone());
            f2->Write(title);

            pm->addGraph(graphM, TString::Format("%.1f#geq |#eta| <%.1f",minETA,maxETA),-1,1,4,20,1,true,true,false,"P E 0 Z L");
            new TCanvas();
            graph->Draw("APL");

        }
        p->setXTitle("H#rightarrowbb particle #it{p}_{T} [GeV]");
        p->setYTitle("#it{m}_{H#rightarrowbb} correction");
        p->setMinMax(1.0,1.85);
        p->setLegendPos(.20,.6,.5,.8);
        p->setCMSLumi(33,"13 TeV","Simulation Preliminary");

        auto * can = p->draw(false,"HbbMassCorrection.pdf");
        p->xAxis()->SetTitleOffset(1.10);
        p->xAxis()->SetRangeUser(200,2500.);
        can->Print("HbbMassCorrection.pdf");


        pm->setXTitle("H#rightarrowbb particle #it{p}_{T} [GeV]");
        pm->setYTitle("#it{m}_{H#rightarrowbb} fit #mu");
        pm->setMinMax(65,125);
        pm->setLegendPos(.5,.2,.85,.5);
        pm->setCMSLumi(33,"13 TeV","Simulation Preliminary");
        auto * can2 = pm->draw(false,"HbbSoftDropMass.pdf");
        pm->xAxis()->SetTitleOffset(1.10);
        pm->yAxis()->SetTitleOffset(1.47);
        pm->xAxis()->SetRangeUser(200,2500.);
        can2->Print("HbbSoftDropMass.pdf");
    }


    void write(TString fileName){ plotter.write(fileName);}


    float puWeight         =0;
    size8 csvCat           =0;
    float bosonMass        =0;
    float bosonPT          =0;
    float bosonETA         =0;
    float drToJet          =0;
    float maxSJDR          =0;
    float maxSJETA         =0;
    float minSJPT          =0;
    float genPT            =0;
    float genETA           =0;
    float jetPT            =0;
    float jetCorrPT        =0;
    float jetETA           =0;
    size8 jetID            =0;
    float jetMass          =0;
    float jetRawSDMass     =0;
    float jetL23CorrSDMass =0;
    float jetCorrSDMass    =0;
    float jetWCorrSDMass   =0;


    float * mass;

    std::vector<std::vector<float>> totPT;
    std::vector<std::vector<float>> nEntries;

    std::vector<std::vector<float>> totPTr;
    std::vector<std::vector<float>> nEntriesr;

    HistGetter plotter;


};

#endif

void getSDMassCorr(std::string fileName, int treeInt, std::string outFileName){
    Analyzer * a = new Analyzer(fileName,"treeMaker/Events",treeInt);
    TFile * fo = new TFile(outFileName.c_str(),"recreate");
    a->analyze();
    a->doFits(fo);
    a->getResp(fo);
        fo->Close();
}
