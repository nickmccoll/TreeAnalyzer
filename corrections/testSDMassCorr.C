
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

#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/StyleInfo.h"

using namespace TAna;

double fnc_dscb(double*xx,double*pp)
{
    double x   = xx[0];
    // gaussian core
    double N   = pp[0];//norm
    double mu  = pp[1];//mean
    double sig = pp[2];//variance
    // transition parameters
    double a1  = pp[3];
    double p1  = pp[4];
    double a2  = pp[5];
    double p2  = pp[6];

    double u   = (x-mu)/sig;
    double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
    double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
    double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
    double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

    double result(N);
    if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
    else if (u<a2)  result *= TMath::Exp(-u*u/2);
    else            result *= A2*TMath::Power(B2+u,-p2);
    return result;
}




class Analyzer : public BaseTreeAnalyzer {
public:
    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt)
{
        TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1){
            signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
        }



}
    void loadVariables() override {

        setBranchAddress("","puWeight"        ,&puWeight        ,true);
        setBranchAddress("","passEvtSel"      ,&passEvtSel        ,true);
        setBranchAddress("","bosonMass"       ,&bosonMass       ,true);
        setBranchAddress("","bosonPT"         ,&bosonPT         ,true);
        setBranchAddress("","bosonETA"         ,&bosonETA         ,true);
        setBranchAddress("","drToJet"         ,&drToJet         ,true);
        setBranchAddress("","maxSJDR"         ,&maxSJDR         ,true);
        setBranchAddress("","maxSJETA"        ,&maxSJETA        ,true);
        setBranchAddress("","minSJPT"         ,&minSJPT         ,true);
        setBranchAddress("","jetPT"           ,&jetPT           ,true);
        setBranchAddress("","jetCorrPT"           ,&jetCorrPT           ,true);
        setBranchAddress("","jetETA"          ,&jetETA          ,true);
        setBranchAddress("","jetID"           ,&jetID           ,true);
        setBranchAddress("","jetMass"         ,&jetMass         ,true);
        setBranchAddress("","jetRawSDMass"    ,&jetRawSDMass    ,true);
        setBranchAddress("","jetL23CorrSDMass",&jetL23CorrSDMass,true);
        setBranchAddress("","jetCorrSDMass"   ,&jetCorrSDMass   ,true);
        setBranchAddress("","jetWCorrSDMass"   ,&jetWCorrSDMass   ,true);

    }




    bool runEvent() override {
        if(jetID < 4) return false;
        if(minSJPT  < 20) return false;
        if(maxSJETA > 2.4) return false;

        bool goodhbb = maxSJDR < 0.8;

        TString prefix = TString::Format("m%.0f",signal_mass);

        auto mkplots =[&](const TString& prefix){
            plotter.getOrMake1DPre(prefix,"rawMass",";H(bb) mass [GeV]",200,0,200 )     ->Fill(jetRawSDMass,puWeight);
            plotter.getOrMake1DPre(prefix,"L23CorrMass",";H(bb) mass [GeV]",200,0,200 ) ->Fill(jetL23CorrSDMass,puWeight);
            plotter.getOrMake1DPre(prefix,"corrMass",";H(bb) mass [GeV]",200,0,200 )    ->Fill(jetCorrSDMass,puWeight);
            plotter.getOrMake1DPre(prefix,"WCorrMass",";H(bb) mass [GeV]",200,0,200 )   ->Fill(jetWCorrSDMass,puWeight);
            plotter.getOrMake1DPre(prefix,"fatjetMass",";H(bb) mass [GeV]",200,0,200 )   ->Fill(jetMass,puWeight);
            plotter.getOrMake1DPre(prefix,"jetResp",";H(bb) jet response", 500,0,2  )    ->Fill(jetPT/bosonPT,puWeight);
            plotter.getOrMake1DPre(prefix,"jetCorrResp",";H(bb) jet response", 500,0,2  )    ->Fill(jetCorrPT/bosonPT,puWeight);

        };

        mkplots(prefix+"_incl");
        if(goodhbb) mkplots(prefix+"_goodHBB");
        if(goodhbb && passEvtSel) mkplots(prefix+"_goodHBBGoodReco");
        if(passEvtSel) mkplots(prefix+"_goodReco");
        return true;
    }




    void write(TString fileName){ plotter.write(fileName);}



    float  puWeight         = 0;
    size8  passEvtSel       = 0;
    float  bosonMass        = 0;
    float  bosonPT          = 0;
    float  bosonETA         = 0;
    float  drToJet          = 0;
    float  maxSJDR          = 0;
    float  jetPT            = 0;
    float  jetCorrPT        = 0;
    float  jetETA           = 0;
    size8  jetID            = 0;
    float  maxSJETA         = 0;
    float  minSJPT          = 0;
    float  jetMass          = 0;
    float  jetRawSDMass     = 0;
    float  jetL23CorrSDMass = 0;
    float  jetCorrSDMass    = 0;
    float  jetWCorrSDMass   = 0;





    float signal_mass;

    HistGetter plotter;


};


void doFit(TString filename){
        TF1* fdscb = new TF1("fdscb",fnc_dscb,50,150,7);
        TF1* gaus = new TF1("gaus","[0]*exp(-0.5*((x-[1])/[2])**2)",50,200);
    TFile * f = new TFile(filename,"read");
//        std::vector<TString> pres = {"incl","goodHBB"};
    std::vector<TString> pres = {"incl","goodHBB","goodHBBGoodReco"};
    std::vector<TString> vars = {"rawMass","L23CorrMass","WCorrMass","corrMass","fatjetMass"};
    std::vector<TString> varNs = {"raw mass", "L23 corr. mass","W-jet corr. mass","H(bb) corr. mass", "AK8 fat jet mass"};
    std::vector<float> masses = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};

    for(unsigned int iP = 0; iP < pres.size();++iP){
        Plotter * pM  = new Plotter();
        Plotter * pS  = new Plotter();
        Plotter * pMoS  = new Plotter();

        for(unsigned int iV = 0; iV < vars.size();++iV){
            TGraphErrors * gM = new TGraphErrors();
            TGraphErrors * gS = new TGraphErrors();
            TGraphErrors * gMoS = new TGraphErrors();
            int nP = 0;
            for(unsigned int iM = 0; iM < masses.size();++iM){
                TH1 * h = 0;
                f->GetObject(TString::Format("m%.0f_%s_%s",masses[iM],pres[iP].Data(),vars[iV].Data()),h);
                if(h == 0) continue;
                h = (TH1*)h->Clone();
                h->Scale(1./h->Integral(0,-1));
                gaus->SetParameter(0,0.03); // N
                gaus->SetParameter(1,125 ); // mean
                gaus->SetParameter(2,12);// sigma

                h->Fit(gaus,"RQB0","",50.,200.);
                const float gausM = gaus->GetParameter(1);
                TString name = TString::Format("M(%.0f), %s, %s",masses[iM],varNs[iV].Data(), pres[iP].Data());

                std::cout << name << " "<< gausM<<std::endl;
                fdscb->SetParameter(0,0.03); // N
                fdscb->SetParameter(1,gausM+5 ); // mean
                fdscb->SetParameter(2,12);// sigma
                fdscb->SetParameter(3,1.0); // a1
                fdscb->SetParameter(4,10.0); // p1
                fdscb->SetParameter(5,2.0); // a2
                fdscb->SetParameter(6,5.0); // p2

                fdscb->SetParLimits(2,0,60);
                fdscb->SetParLimits(3,0.2,3.);
                fdscb->SetParLimits(5,1,5.);

                fdscb->FixParameter(6,1 ); // p2
                fdscb->FixParameter(4,10); // p1

                if(h->Fit(fdscb,"RQB","",gausM+5 - 50,gausM+5 + 50)) continue;

//                    new TCanvas(name,name);
//                    h->Draw();

                gM->SetPoint(nP,masses[iM],fdscb->GetParameter(1));
                gM->SetPointError(nP,0,fdscb->GetParError(1));
                gS->SetPoint(nP,masses[iM],fdscb->GetParameter(2));
                gS->SetPointError(nP,0,fdscb->GetParError(2));
                gMoS->SetPoint(nP,masses[iM],fdscb->GetParameter(2)/fdscb->GetParameter(1));
                gMoS->SetPointError(nP,0,fdscb->GetParError(2)/fdscb->GetParameter(1));
                nP++;
            }

            pM->addGraph(gM, varNs[iV],-1,1,4,20,1,true,false,false,"P L");
            pS->addGraph(gS, varNs[iV],-1,1,4,20,1,true,false,false,"P L");
            pMoS->addGraph(gMoS, varNs[iV],-1,1,4,20,1,true,false,false,"P L");
        }

        pM->setXTitle("signal sample radion mass [GeV]");
        pM->setYTitle("mass from fit [GeV]");
        pM->setCMSLumi(33,"13 TeV","Simulation Preliminary");
        pM->setMinMax(100,150);
        pM->setLegendPos(.15,.65,.5,.9);
        auto * can  = pM->draw(false,TString::Format("massTrend_%s.pdf",pres[iP].Data()));
        pM->xAxis()->SetTitleOffset(1.15);
        pM->yAxis()->SetTitleOffset(1.47);
        can->Print(TString::Format("massTrend_%s.pdf",pres[iP].Data()));

        pMoS->setXTitle("signal sample radion mass [GeV]");
        pMoS->setYTitle(" #sigma(mass) / mass (from fit)");
        pMoS->setCMSLumi(33,"13 TeV","Simulation Preliminary");
        pMoS->setMinMax(0.087,0.12);
        pMoS->setLegendPos(.5,.49,.83,.75);
        TCanvas * can2 = pMoS->draw(false,TString::Format("sigmaOMassTrend_%s.pdf",pres[iP].Data()));
        pMoS->xAxis()->SetTitleOffset(1.15);
        pMoS->yAxis()->SetTitleOffset(1.60);
        can2->Print(TString::Format("sigmaOMassTrend_%s.pdf",pres[iP].Data()));

        pS->setXTitle("signal sample radion mass [GeV]");
        pS->setYTitle("#sigma(mass) from fit [GeV]");
        pS->setCMSLumi(33,"13 TeV","Simulation Preliminary");
        pS->setMinMax(9,24);
        pS->setLegendPos(.15,.65,.5,.9);
        TCanvas * can3 = pS->draw(false,TString::Format("sigmaTrend_%s.pdf",pres[iP].Data()));
        pS->xAxis()->SetTitleOffset(1.15);
        pS->yAxis()->SetTitleOffset(1.47);
        can3->Print(TString::Format("sigmaTrend_%s.pdf",pres[iP].Data()));
        new TCanvas();

    }
}

void doRespTest(TString filename){
    TFile * f = new TFile(filename,"read");
    std::vector<TString> pres = {"incl","goodHBB","goodHBBGoodReco"};
    std::vector<TString> vars = {"jetResp","jetCorrResp"};
    std::vector<TString> varNs = {"L23 corr.", "L23 corr. + H(bb) resid."  };
    std::vector<float> masses = {600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};

    for(unsigned int iP = 0; iP < pres.size();++iP){
        Plotter * pM  = new Plotter();

        for(unsigned int iV = 0; iV < vars.size();++iV){
            TGraphErrors * gM = new TGraphErrors();
            int nP = 0;
            for(unsigned int iM = 0; iM < masses.size();++iM){
                TH1 * h = 0;
                f->GetObject(TString::Format("m%.0f_%s_%s",masses[iM],pres[iP].Data(),vars[iV].Data()),h);
                if(h == 0){
                    std::cout << TString::Format("m%.0f_%s_%s",masses[iM],pres[iP].Data(),vars[iV].Data()) << std::endl;
                    continue;
                }
                h = (TH1*)h->Clone();
                h->Scale(1./h->Integral(0,-1));
                gM->SetPoint(nP,masses[iM],h->GetMean());
                gM->SetPointError(nP,0,h->GetMeanError());
                nP++;
            }

            pM->addGraph(gM, varNs[iV],-1,1,4,20,1,true,false,false,"P L");
        }

        pM->setXTitle("signal sample radion mass [GeV]");
        pM->setYTitle("<jet response>");
        pM->setCMSLumi(33,"13 TeV","Simulation Preliminary");
        pM->setMinMax(0.92,1.05);
        pM->setLegendPos(.15,.8,.6,.9);
        auto * can  = pM->draw(false,TString::Format("avgRespTrend_%s.pdf",pres[iP].Data()));
        pM->xAxis()->SetTitleOffset(1.15);
        pM->yAxis()->SetTitleOffset(1.47);
        can->Print(TString::Format("avgRespTrend_%s.pdf",pres[iP].Data()));
        new TCanvas();

    }
}

#endif

void testSDMassCorr(std::string fileName, int treeInt, std::string outFileName){
    Analyzer * a = new Analyzer(fileName,"treeMaker/Events",treeInt);
    a->analyze();
    a->write(outFileName);
}

void testSDMassCorr(std::string outFileName){
    doFit(outFileName);
    doRespTest(outFileName);
}
