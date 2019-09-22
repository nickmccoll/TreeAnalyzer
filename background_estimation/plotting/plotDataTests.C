#include "../predTools/PostFitter.h"
#include "../predTools/StatTester.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;
std::vector<TObject*> writeables;

const float signalXS = 0.2; //in pb

class Dummy {
public:
    Dummy(const std::string& outName = "") : outName(outName) {};
    ~Dummy() {
        if(outName.size()){
            TFile * f = new TFile((outName+".root").c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
                w->Print((outName +"_"+w->GetTitle() +".pdf").c_str());
            }
            f->Close();
        }
    }
    std::string outName;
};

bool preliminary = false;

enum ModelType {MOD_NONE, MOD_MC, MOD_PRE, MOD_POST};
std::string modTitles[] = {"NONE","MC","prefit","fit"};
struct DataPlotPrefs {
    bool addData = true;
    std::pair<double,double> blindRange = std::pair<double,double>(0,0);
    ModelType modelType = MOD_PRE;
    ModelType addType   = MOD_NONE;
    bool plotIntegral = false;
    std::vector<double> bins; // if bin+1 < bin, it will skip [bin,bin+1], if bin +1 == bin it will skip [bin,bin+1] and [bin+1,bin+2]
    bool binInY = true;
    std::vector<std::string> sels;
    std::vector<std::string> titles;
    std::string topTitle;
    std::vector<std::string> botTitles;

    std::string signalFilePrefix = "signalInputs/HHbb1o2l_radHH_";
    std::vector<std::string> signals;
    std::vector<std::string> signalTitles;

    int rebinFactor = 1; //-1 means go to rebins...1 means dont rebin...2+ is standard grouping
    std::vector<double> rebins;
    bool doLog = false;
    bool addRatio = false;
    bool addErrorBars = false;
    bool removeTrailingZeros = false;

    unsigned int nToys = 100;

    bool isSupp = false;

    std::vector<float> maxTops; //override maxTop
    std::vector<float> minTops; //override minTop

    double minTop = -1;
    double maxTop = -1;
    double minBot = -1;
    double maxBot = -1;
};

struct HistContainer {
    HistContainer(bool is1l) {
        toys2D.reserve(100);
        bkg2D.resize(BKG_MT+1); for(auto* b :bkg2D) b=0;
        bkg.resize(BKG_MT+1); for(auto* b :bkg) b=0;
    }
    std::vector<TH2*> bkg2D;
    TH2* tot2D=0;
    TH2* add2D=0;
    TH2* data2D=0;
    std::vector<TH2*> toys2D;
    std::vector<TH2*> sig2D;

    std::vector<TH1*> bkg;
    std::vector<TH1*> sig;
    TH1* tot=0;
    TH1* add=0;
    TH1* data=0;
    TGraphAsymmErrors* toyErr=0;
};


void doComp(REGION region, const unsigned int nToys, const std::string& postFitFilename,const std::string& postFitCompFilename, const std::vector<std::string>& signals ){
    TFile * inFile  = new TFile(postFitFilename.c_str(),"read");
    TFile * fo      = new TFile(postFitCompFilename.c_str(),"recreate");

    auto sels = getSRList(region);
    std::string compName = "emu_LMT_I_full";
    for(unsigned int iT = 1; iT <= nToys; ++iT ){
        TH2 * h = 0;
        for(const auto& s: sels){

            TH2 * th = 0;
            inFile->GetObject((std::string("toyModel_")+ASTypes::int2Str(iT)+"_"+ s+"__"+MOD_MJ+"_"+MOD_MR).c_str(),th);
            if(!th) continue;
            if(h)h->Add(th);
            else h = (TH2*)th->Clone();
        }
        fo->cd();
        h->Write((std::string("toyModel_")+ASTypes::int2Str(iT)+"_"+ compName+"__"+MOD_MJ+"_"+MOD_MR).c_str());
    }
    for(unsigned int iT = 0; iT < bkgSels.size(); ++iT){
        TH2 * h = 0;
        for(const auto& s: sels){
            TH2 * th = 0;
            inFile->GetObject((std::string("postfit_")+ bkgSels[iT]+"_"+s+"__"+MOD_MJ+"_"+MOD_MR).c_str(),th);
            if(!th) continue;
            if(h)h->Add(th);
            else h = (TH2*)th->Clone();
        }
        fo->cd();
        h->Write((std::string("postfit_")+ bkgSels[iT]+"_"+compName+"__"+MOD_MJ+"_"+MOD_MR).c_str());
    }

    fo->Close();

    TFile * fos      = new TFile(("signalInputs/HHlnujj_radHH_" + compName +"_2D_fit.root").c_str(),"recreate");
    for(unsigned int iSig = 0; iSig < signals.size(); ++iSig){
        TH2 * h = 0;
        for(const auto& s: sels){
            TFile * f = new TFile(("signalInputs/HHlnujj_radHH_" + s +"_2D_fit.root").c_str() ,"READ");
            TFile * fy = new TFile(("signalInputs/HHlnujj_radHH_" + s +"_yield.json.root").c_str() ,"READ");
            TF1 * yFunc = 0;
            TH2 * hP = 0;
            fy->GetObject("yield_func",yFunc);
            f->GetObject((std::string("pdf_m")+signals[iSig]+ "__MJ_MR").c_str(),hP);
            if(!yFunc || !hP){
                continue;
            }
            double yield = yFunc->Eval(float(std::atoi(signals[iSig].c_str())));
            //have to do this...bounds are different in the signal
            hP->Scale(yield/hP->Integral());
            //Scale so that we get 0.2pb normalization
            hP->Scale(signalXS);
            TH2 * hpNew = new TH2F(TString(hP->GetName()) + "_newC_"+ s,hP->GetTitle(),nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
            hpNew->SetDirectory(0);
            for( int iX = 1; iX <= hpNew->GetNbinsX(); ++iX ){
                int nX = hP->GetXaxis()->FindFixBin(hpNew->GetXaxis()->GetBinCenter(iX));
                for( int iY = 1; iY <= hpNew->GetNbinsY(); ++iY ){
                    int nY = hP->GetYaxis()->FindFixBin(hpNew->GetYaxis()->GetBinCenter(iY));
                    hpNew->SetBinContent(iX,iY,hP->GetBinContent(nX,nY));
                }
            }
            if(h==0) h = hpNew;
            else h->Add(hpNew);
        }
        if(h){
            fos->cd();
            h->Write((std::string("sigComp_m")+signals[iSig]+ "__MJ_MR").c_str());
        }
    }
    fos->Close();
}



class DataPlotter {
public:
    DataPlotter(const DataPlotPrefs& plotPrefs, const std::string& inputPrefix, const std::string& postFitFilename, const bool is1l ) : prefs(plotPrefs), is1lep(is1l)
{

        if(prefs.addData) dataFile  = new TFile((inputPrefix+"_data_distributions.root").c_str(),"read");
        if(prefs.modelType == MOD_POST || prefs.addType == MOD_POST || prefs.addErrorBars )
            postfitFile  = new TFile(postFitFilename.c_str(),"read");



        if(prefs.modelType == MOD_MC || prefs.addType == MOD_MC || prefs.modelType == MOD_PRE || prefs.addType == MOD_PRE ){
            for(const auto& t: bkgSels){
                TFile * fY = new TFile((inputPrefix+"_"+t+"_distributions.root").c_str(),"read");
                mcFiles.push_back(fY);
            }
        }
        if(prefs.modelType == MOD_PRE || prefs.addType == MOD_PRE ){
            for(const auto& t: bkgSels){
                TFile * fF = new TFile((inputPrefix+"_"+t+"_2D_template_debug.root").c_str(),"read");
                prefitFiles.push_back(fF);
            }
        }

        for(unsigned int iSel = 0; iSel < prefs.sels.size(); ++iSel){
            std::vector<TH2*> sHists;
            if(prefs.signals.size()){
                TFile * f = new TFile((prefs.signalFilePrefix + prefs.sels[iSel] +"_2D_fit.root").c_str() ,"READ");
                TFile * fy = new TFile((prefs.signalFilePrefix + prefs.sels[iSel] +"_yield.json.root").c_str() ,"READ");
                TF1 * yFunc = 0;
                fy->GetObject("yield_func",yFunc);
                if(yFunc==0)  for(unsigned int iSig = 0; iSig < prefs.signals.size(); ++iSig){
                    sHists.push_back(0);
                }
                for(unsigned int iSig = 0; iSig < prefs.signals.size(); ++iSig){
                    TH2 * hP = 0;
                    f->GetObject((std::string("pdf_m")+prefs.signals[iSig]+ "__MJ_MR").c_str(),hP);
                    if(hP){
                        double yield = yFunc->Eval(float(std::atoi(prefs.signals[iSig].c_str())));
                        std::cout << "WARNING!!!!! You are applying a correction to the BF!!!"<<std::endl;
                        std::cout << "WARNING!!!!! You are applying a correction to the BF!!!"<<std::endl;
                        std::cout << "WARNING!!!!! You are applying a correction to the BF!!!"<<std::endl;
                        std::cout << "WARNING!!!!! You are applying a correction to the BF!!!"<<std::endl;
                        yield *=CutConstants::HHtobbVVBF;

                        hP->Scale(yield/hP->Integral());
                        //Scale so that we get 0.2pb normalization
                        hP->Scale(signalXS);
                        TH2 * hpNew = new TH2F(TString(hP->GetName()) + "_new",hP->GetTitle(),nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
                        hpNew->SetDirectory(0);
                        for( int iX = 1; iX <= hpNew->GetNbinsX(); ++iX ){
                            int nX = hP->GetXaxis()->FindFixBin(hpNew->GetXaxis()->GetBinCenter(iX));
                            for( int iY = 1; iY <= hpNew->GetNbinsY(); ++iY ){
                                int nY = hP->GetYaxis()->FindFixBin(hpNew->GetYaxis()->GetBinCenter(iY));
                                hpNew->SetBinContent(iX,iY,hP->GetBinContent(nX,nY));
                            }
                        }
                        sHists.push_back(hpNew);


                    } else {
                        TH2 * hMC = 0;
                        f->GetObject((std::string("sigComp_m")+prefs.signals[iSig]+ "__MJ_MR").c_str(),hMC);
                        if(hMC){
                            TH2 * hpNew = (TH2*)hMC->Clone();
                            hpNew->SetDirectory(0);
                            sHists.push_back(hpNew);
                        } else {
                            sHists.push_back(0);
                            continue;
                        }
                    }
                }
            }
            signalHists.push_back(sHists);
        }
}

    void getToyModels(const std::string& sel, HistContainer& cont){
        for(unsigned int iT = 1; iT <= prefs.nToys; ++iT ){
            TH2 * th = 0;
//            postfitFile->GetObject((std::string("toyModel_")+ASTypes::int2Str(iT)+"_"+ sel+"__"+MOD_MJ+"_"+MOD_MR).c_str(),th);
            postfitFile->GetObject((std::string("toyModel_")+ASTypes::int2Str(iT)+"_"+ sel).c_str(),th);
            if(th)cont.toys2D.push_back(th);
        }
    }

    std::vector<TH2*> getComponents(const ModelType type, const std::string& sel ){
        std::vector<TH2*> bkgs;
        for(unsigned int iT = 0; iT < bkgSels.size(); ++iT){
            TH2 * h = 0;
            if(type == MOD_MC) {
                mcFiles[iT]->GetObject((bkgSels[iT]+"_"+sel+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
            }
            if(type == MOD_PRE) {
                prefitFiles[iT]->GetObject((bkgSels[iT]+"_"+sel).c_str(),h);
            }
            if(type == MOD_POST) {
//                postfitFile->GetObject((std::string("postfit_")+ bSels[iT]+"_"+sel+"_"+MOD_MJ+"_"+MOD_MR).c_str(),h);
            	postfitFile->GetObject((std::string("postfit_")+ bkgSels[iT]+"_"+sel).c_str(),h);
            }
            bkgs.push_back(h);
        }
        return bkgs;
    }
    void getBackgrounds(const std::string& sel, HistContainer& cont){
        std::vector<TH2*> postFits;
        std::vector<TH2*> preFits;
        std::vector<TH2*> mcs;
        if(prefs.modelType == MOD_POST || prefs.addType == MOD_POST ) postFits = getComponents(MOD_POST,sel);
        if(prefs.modelType == MOD_MC || prefs.addType == MOD_MC || prefs.modelType == MOD_PRE || prefs.addType == MOD_PRE )
            mcs = getComponents(MOD_MC,sel);
        if(prefs.modelType == MOD_PRE || prefs.addType == MOD_PRE ){
            preFits = getComponents(MOD_PRE,sel);
            for(unsigned int iT = 0; iT < bkgSels.size(); ++iT){
                if(preFits[iT]) preFits[iT]->Scale(mcs[iT]->Integral()/preFits[iT]->Integral());
            }
        }

        auto combine = [](std::vector<TH2*>& comps, std::vector<CutStr> bkgCats) -> TH2* {
            TH2* total = 0;
//            for(unsigned int iT = 0; iT < BkgSels.size(); ++iT){
            for(unsigned int iT = 0; iT < bkgCats.size(); ++iT){
                if(comps[iT] ==0) continue;
                if(total==0) total = (TH2*)comps[iT]->Clone();
                else total->Add(comps[iT]);
            }
            return total;
        };

        if(prefs.modelType == MOD_POST) cont.bkg2D = postFits;
        if(prefs.modelType == MOD_PRE)  cont.bkg2D = preFits ;
        if(prefs.modelType == MOD_MC)   cont.bkg2D = mcs     ;

        cont.tot2D = combine(cont.bkg2D,bkgSels);

        if(prefs.addType == MOD_POST)cont.add2D=combine(postFits,bkgSels);
        if(prefs.addType == MOD_PRE) cont.add2D=combine(preFits,bkgSels);
        if(prefs.addType == MOD_MC) cont.add2D=combine(mcs,bkgSels);
        if(prefs.addType == MOD_NONE) cont.add2D=0;
    }

    TH1* makeProjection(TH2* h, const std::string& title, const std::string s, const int iB, const int binL, const int binH){
        auto h1 = prefs.binInY ? h->ProjectionX(  (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH)
                :  h->ProjectionY( (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH);
        if(prefs.rebinFactor > 1){
            h1->Rebin(prefs.rebinFactor);
        } else if(prefs.rebinFactor < 1){
            h1 = (TH1D*)h1->Rebin(prefs.rebins.size() -1,"", &prefs.rebins[0]);
        }

        if(prefs.plotIntegral) return PlotTools::getIntegral(h1,true,false);
        return h1;
    }
    std::pair<int,int> getBins(const HistContainer& cont,const int iB){
        const TAxis * ax =0;
        for(const auto* h: cont.bkg2D){
            if(h==0) continue;
            ax = prefs.binInY ? h->GetYaxis() : h->GetXaxis();
            break;
        }
        if(!ax && cont.add2D) ax = prefs.binInY ? cont.add2D->GetYaxis() : cont.add2D->GetXaxis();
        if(!ax && cont.data2D) ax = prefs.binInY ? cont.data2D->GetYaxis() : cont.data2D->GetXaxis();
        int binL = ax->FindFixBin(prefs.bins[iB]);
        int binH = ax->FindFixBin(prefs.bins[iB+1]) -1;
        return std::make_pair(binL,binH);
    }

    void get1DHists(const std::string& s, const unsigned int iB, HistContainer& cont){
        auto bins = getBins(cont,iB);

        auto processH2 = [&](TH2* h, const std::string& title)->TH1*{
            return makeProjection(h,title,s,iB,bins.first,bins.second);
        };

        for(unsigned int iH = 0; iH < cont.bkg2D.size(); ++iH){
            if(cont.bkg2D[iH]==0) continue;
            TH1 * h = processH2(cont.bkg2D[iH],bkgSels[iH]);
            for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
            cont.bkg[iH] = h;
        }

        //add additional
        if(cont.tot2D){
            auto h = processH2(cont.tot2D,modTitles[prefs.modelType]);
            cont.tot = h;
        }
        if(cont.add2D){
            auto h = processH2(cont.add2D,modTitles[prefs.addType]);
            cont.add = h;
        }

        //add data
        if(cont.data2D){
            auto h = processH2(cont.data2D,"Data");
            //blinding
            if(prefs.blindRange.first < prefs.blindRange.second){
                for(int iB = 1; iB <= h->GetNbinsX(); ++iB){
                    if((h->GetBinLowEdge(iB) + h->GetBinWidth(iB)) > prefs.blindRange.first &&
                            h->GetBinLowEdge(iB) < prefs.blindRange.second   ){
                        h->SetBinContent(iB,-1);
                        h->SetBinError(iB,0);
                    }
                }
            }
            cont.data = h;
        }

        //do toys
        if(cont.tot&&cont.toys2D.size()){
            std::vector<TH1*> toyProj; toyProj.reserve(cont.toys2D.size());
            for(unsigned int iT = 0; iT < cont.toys2D.size(); ++iT)
                toyProj.push_back(processH2(cont.toys2D[iT],std::string("toy_")+ASTypes::int2Str(iT)) );

            const int nEntries =  toyProj.size();
            std::vector<double> toyVs(nEntries);
            auto getErrorBands = [&](const int iB,double&mean, double&  eL, double& eH, const double alpha = (1 - 0.6827)){
                double t= 0;
                for(unsigned int i = 0; i < toyProj.size(); ++i){
                    toyVs[i] = toyProj[i]->GetBinContent(iB);
                    t += toyProj[i]->GetBinContent(iB);
                }
                std::sort(toyVs.begin(), toyVs.end());
                eL = toyVs[int( double(nEntries)*alpha/2  )];
                eH = toyVs[int( double(nEntries)* (1 - alpha/2)  )];
                mean = t/double(toyProj.size());
            };

            cont.toyErr = new TGraphAsymmErrors();

            for(int iB = 1; iB <= cont.tot->GetNbinsX(); ++iB){
                double x = cont.tot->GetBinCenter(iB);
                double y = cont.tot->GetBinContent(iB);
                double m,eL,eH;
                getErrorBands(iB,m,eL,eH);
                auto setPt = [&](const int bin, float x){
                    cont.toyErr->SetPoint(bin,x,y);
                    cont.toyErr->SetPointError(bin,cont.tot->GetBinWidth(iB)/2.,cont.tot->GetBinWidth(iB)/2.,std::max(m -eL,0.),std::max(eH-m,0.));
                };
                if(iB == 1) setPt(iB-1,x-cont.tot->GetBinWidth(iB)/2.);
                setPt(iB,x);
                if(iB == cont.tot->GetNbinsX() ) setPt(iB+1,x+cont.tot->GetBinWidth(iB)/2.);
            }
            for(auto * t :toyProj) {t->SetDirectory(0); delete t;}
        }

        //do signal
        cont.sig.resize(cont.sig2D.size(),0);
        for(unsigned int iH = 0; iH < cont.sig2D.size(); ++iH){
            if(cont.sig2D[iH]==0) continue;
            TH1 * h = processH2(cont.sig2D[iH],prefs.signals[iH]);
            for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
            cont.sig[iH] = h;
        }
    }


    std::vector<TObject*> makePlots() {
        std::vector<TObject*> writeables;

        for(unsigned int iS = 0; iS < prefs.sels.size(); ++iS){
            const auto& s = prefs.sels[iS];
            const auto& st = prefs.titles[iS];
            HistContainer cont(is1lep);
            cont.sig2D = signalHists[iS];
            getBackgrounds(s,cont);
            if(prefs.modelType == MOD_POST && prefs.addErrorBars) getToyModels(s,cont);

            //add in the data to nBKS
            TH2 * h = 0;
            if(prefs.addData){
                dataFile->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
            }
            cont.data2D = h;

            //Individual bins
            for(unsigned int iB = 0; iB + 1 < prefs.bins.size(); ++iB){
                if(prefs.bins[iB+1] <= prefs.bins[iB]){
                    if(prefs.bins[iB+1] == prefs.bins[iB]) ++iB;
                    continue;
                }
                get1DHists(s,iB,cont);
                const std::string plotTitle = s +"_"+(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]);

                Plotter * p = new Plotter();
                std::vector<Drawing::TLegendEntryDef> legEntries;

                if(cont.toyErr){
                    int fillColor = kMagenta+4;
                    cont.toyErr->SetFillColor(fillColor);
                    cont.toyErr->SetFillStyle(3352);
                    gStyle->SetHatchesLineWidth(1);
                    gStyle->SetHatchesSpacing(.5);
                    //                    p->addGraph(cont.toyErr,"fit unc.",fillColor,1,1,20,1,false,true,false,"3");

                    auto * g = p->addGraph(cont.toyErr,"Fit unc.",fillColor,1,0,20,1,false,true,false,"2");
                    legEntries.push_back(std::make_tuple(2,g,"Fit unc.","F"));
                }

                if(cont.add){
                    auto * g = p->addHistLine(cont.add,modTitles[prefs.addType],kBlue);
                    legEntries.push_back(std::make_tuple(3,g,modTitles[prefs.addType],"L"));
                }

                if( int(cont.toyErr!=0) + int(cont.add!=0)  == 1  ){
                    legEntries.push_back(std::make_tuple(4,(TObject*)(0),"",""));
                }

                std::vector<int> signalColors ={kSpring+10,634};
                if(cont.sig.size() == 1 ) signalColors = {634};

                for(unsigned int iSig = 0; iSig < cont.sig.size(); ++iSig){
                    if(cont.sig[iSig]){
                        auto * g = p->addHistLine(cont.sig[iSig],prefs.signalTitles[iSig],signalColors[iSig]);
                        legEntries.push_back(std::make_tuple(100+iSig,g,prefs.signalTitles[iSig],"L"));
                    }
                }

                if(cont.data){
                    p->addHist(cont.data,"Data",kBlack,1,2,20,0.5,true,true,true);
                    TGraphAsymmErrors * g = new TGraphAsymmErrors;
                    g->SetLineColor  (kBlack);
                    g->SetLineWidth  (2);
                    g->SetLineStyle  (1);
                    g->SetMarkerStyle(20);
                    g->SetMarkerColor(kBlack);
                    g->SetMarkerSize (0.5);
                    legEntries.push_back(std::make_tuple(0,g,"Data","P E"));
                }
                legEntries.push_back(std::make_tuple(1,(TObject*)(0),"",""));

                int nB = 0;
                for(unsigned int iH = 0; iH < cont.bkg.size(); ++iH){
                    if(cont.bkg[iH]){
                        nB ++;
                        auto *g = p->addStackHist(cont.bkg[iH],bkgSels[iH].title.c_str());
                        legEntries.push_back(std::make_tuple(10+cont.bkg.size() -iH,g,bkgSels[iH].title.c_str(),"f"));
                    }
                }
                if(nB % 2)  legEntries.push_back(std::make_tuple(10+ cont.bkg.size()+1,(TObject*)(0),"",""));

                if(cont.toyErr)p->clearTotStackError();
                p->setUnderflow(false);
                p->setOverflow(false);
                double binWidth = 10;
                for(const auto* h : cont.bkg){
                    if(h==0) continue;
                    binWidth = h->GetBinWidth(1);
                    break;
                }

                double xV =0.45;
                double yV =0.63;;

                p->setXTitle( (prefs.binInY ? hbbMCS : hhMCS) .title.c_str());
                p->setYTitle((std::string("Events / ") + flt2Str(binWidth) +" GeV").c_str() );
                if(prefs.isSupp || preliminary){
                    p->setCMSLumi(0);
                    if( prefs.isSupp) p->setCMSLumiExtraText("Supplementary");
                    else  p->setCMSLumiExtraText("Preliminary");
                    p->setCMSLumiPosition(0,1.05);

                } else {
                    p->setCMSLumi(10);
                }

                float startY = (prefs.isSupp||preliminary) ? 0.83 : 0.75;
                if(prefs.topTitle.size() ){
                    p->addText(prefs.topTitle,.18,startY,0.045);
                    startY -= 0.05;
                }

                p->addText(st,.18,startY,0.045);
                startY -= 0.05;
                if(prefs.botTitles.size() && prefs.botTitles[iB].size()){
                    p->addText(prefs.botTitles[iB],.18,startY,0.045);
                }

                p->setLegendNColumns(2);

                if(prefs.signals.size()){
                    p->setLegendPos(xV,yV,xV+0.457,yV+0.245);
                    TString sigName = TString::Format("#sigma#bf{#it{#Beta}}(X #rightarrow HH) = %.1f pb",signalXS);
                    p->addText(sigName.Data(),xV+0.0075,yV-0.035,0.042);
                }
                else {
                    p->setLegendPos(xV,yV,xV+0.445,yV+0.245);
                }

                //--------------------LEGEND AND TEXT------------------------------
                p->turnOffLegend();
                TLegend * legend = prefs.signals.size() ? new TLegend(xV,yV,xV+0.457,yV+0.245) : new TLegend(xV,yV,xV+0.445,yV+0.245);
                legend->SetFillStyle(0);
                legend->SetBorderSize(0);
                legend->SetNColumns(2);
                std::sort(legEntries.begin(), legEntries.end(), [](const Drawing::TLegendEntryDef& a, const Drawing::TLegendEntryDef& b) {return  std::get<0>(a) < std::get<0>(b); }  );
                for(const auto& l : legEntries){
                    legend->AddEntry(std::get<1>(l),std::get<2>(l),std::get<3>(l));
                }

                //--------------------LEGEND AND TEXT------------------------------
                if(prefs.removeTrailingZeros == true) p->turnOffTrailingPoissonZeros();

                if(prefs.doLog) p->setMinMax(0.1,1000 );
                if(prefs.minTop != prefs.maxTop){
                    p->setMinMax(prefs.minTop,prefs.maxTop);
                } else if(prefs.maxTops.size()|| prefs.minTops.size()){
                    double tempMax = prefs.maxTops.size() ? prefs.maxTops[iS] : prefs.maxTop;
                    double tempMin = prefs.minTops.size() ? prefs.minTops[iS] : prefs.minTop;

                    p->setMinMax(tempMin,tempMax);
                }
                else if(cont.data) p->setMinMax(0,(cont.data->GetMaximum()  + std::sqrt(cont.data->GetMaximum()))*1.75);

                if(prefs.addRatio){
                    p->setBotMinMax(0.05,1.95);
                    if(prefs.minBot != prefs.maxBot){
                        p->setBotMinMax(prefs.minBot,prefs.maxBot);
                    }
                    p->setYTitleBot((std::string("Data / ") + modTitles[prefs.modelType] +"").c_str());
                    auto * c = p->drawSplitRatio(-1,"stack",false,false,plotTitle.c_str());

                    if(prefs.doLog) c->GetPad(1)->SetLogy();
                    p->botStyle.xAxis->SetTitleOffset(1.05);
                    for(unsigned int iSig = 0; iSig < cont.sig.size(); ++iSig){
                        auto prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_0" ).c_str());
                        if(prim) prim->Delete();
                        prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_1" ).c_str());
                        if(prim) prim->Delete();
                        prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_2" ).c_str());
                        if(prim) prim->Delete();
                        prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_3" ).c_str());
                        if(prim) prim->Delete();
                        prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_4" ).c_str());
                                                if(prim) prim->Delete();
                                                prim = c->GetPad(2)->GetPrimitive((s+"_"+prefs.signals[iSig]+"_5" ).c_str());
                                                                        if(prim) prim->Delete();
                    }

                    c->GetPad(1)->cd();
                    legend->Draw();
                    c->GetPad(1)->Update();

                    writeables.push_back(c);
                } else {
                    auto * c = p->draw(false,plotTitle.c_str());
                    if(prefs.doLog)c->SetLogy();
                    c->cd();
                        legend->Draw();
                    c->Update();
                    writeables.push_back(c);
                }

            }
        }

        return writeables;
    }

    std::vector<TObject*>  makeStatTest(const std::string outNamePrefix, bool useBuiltInToys = true) {
        std::vector<TH1*> dataTSs(prefs.bins.size() -1,0);
        std::vector<TGraphAsymmErrors*> toyTSs(prefs.bins.size() -1,0);
        std::vector<TObject*> writeables;

        for(unsigned int iS = 0; iS < prefs.sels.size();++iS){
            const auto& s = prefs.sels[iS];
            HistContainer cont(is1lep);
            getBackgrounds(s,cont);
            //add in the data to nBKS
            dataFile->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),cont.data2D);

            //get toys
            std::vector<TH2*> toyFits;
            std::vector<TH2*> toyData;
            if(useBuiltInToys){
                toyFits.reserve(prefs.nToys);
                toyData.reserve(prefs.nToys);
                for(unsigned int iT = 1; iT <= prefs.nToys; ++iT ){
                    TH2 *hd=0;
                    TH2* hm=0;
                    postfitFile->GetObject((std::string("toyData_")+ASTypes::int2Str(iT)+"_"+ s+"__"+MOD_MJ+"_"+MOD_MR).c_str(),hd);
                    if(prefs.modelType == MOD_POST)
                        postfitFile->GetObject((std::string("toyDataFit_")+ASTypes::int2Str(iT)+"_"+ s+"__"+MOD_MJ+"_"+MOD_MR).c_str(),hm);
                    else
                        hm = (TH2*)cont.tot2D->Clone();
                    if(hd && hm) { toyFits.push_back(hm) ,toyData.push_back(hd);}
                }
            }

            //Individual bins
            for(unsigned int iB = 0; iB + 1 < prefs.bins.size(); ++iB){
                if(prefs.bins[iB+1] <= prefs.bins[iB]){
                    if(prefs.bins[iB+1] == prefs.bins[iB]) ++iB;
                    continue;
                }

                get1DHists(s,iB,cont);

                const std::string plotTitle = s +"_"+(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]);
                StatTesterAnalyzer * a= 0;
                if(useBuiltInToys){

                    auto bins = getBins(cont,iB);
                    std::vector<StatTesterAnalyzer::ModelAndData> toys; toys.reserve(toyFits.size());
                    for(unsigned int iT = 0; iT < toyFits.size(); ++iT ){
                        toys.emplace_back((TH1D*)makeProjection(toyFits[iT],std::string("toyDataFit_")+ASTypes::int2Str(iT),s,iB,bins.first,bins.second),
                                (TH1D*)makeProjection(toyData[iT],std::string("toyData_")+ASTypes::int2Str(iT),s,iB,bins.first,bins.second));
                    }
                    a = new StatTesterAnalyzer({(TH1D*)cont.tot,(TH1D*)cont.data},toys,plotTitle,outNamePrefix+"_"+plotTitle+".root" );
                    for(unsigned int iT = 0; iT < toyFits.size(); ++iT ){
                        delete toys[iT].first;
                        delete toys[iT].second;
                    }

                } else {
                    a = new StatTesterAnalyzer((TH1D*)cont.tot,(TH1D*)cont.data,10000,false,plotTitle,outNamePrefix+"_"+plotTitle+".root" );
                }

                if(dataTSs[iB] == 0){
                    std::string sumName = (prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]);
                    dataTSs[iB] = new TH1F((sumName +"_data_TS").c_str(),";selection",prefs.sels.size(),-0.5,prefs.sels.size()-0.5);
                    toyTSs[iB] = new TGraphAsymmErrors();//(sumName +"_toy_TS",";selection",prefs.sels.size(),-0.5,prefs.sels.size()-0.5);
                }

                dataTSs[iB]->SetBinContent(iS+1,a->ts_nom_sa);
                toyTSs[iB]->SetPoint(iS,iS,a->ts_avg_sa);
                toyTSs[iB]->SetPointError(iS,0,0,a->ts_avg_sa - a->ts_down_sa,a->ts_up_sa-a->ts_avg_sa);
                delete a;
            }

            for(unsigned int iT = 0; iT < toyFits.size(); ++iT ){
                delete toyFits[iT];
                delete toyData[iT];
            }
        }

        for(unsigned int iB = 0; iB < dataTSs.size(); ++iB){
            if(dataTSs[iB] == 0) continue;
            Plotter * p = new Plotter;
            p->addGraph(toyTSs[iB],"Toy data");
            p->addHist(dataTSs[iB],"Data",-1,1,4,20,1,true,false);
            p->setXTitle(" ");
            p->setYTitle("test statistic");
            auto c = p->draw(false,(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]) +"_testStat");
            for(unsigned int iS = 0; iS < prefs.sels.size();++iS){
                p->xAxis()->SetBinLabel(iS+1,getCategoryLabel(prefs.sels[iS]).c_str());
            }
            p->xAxis()->SetTitle(" ");
            c->Update();

            writeables.push_back(c);
        }

        return writeables;
    }


    DataPlotPrefs prefs;
    std::vector<TFile*> mcFiles;
    std::vector<TFile*> prefitFiles;
    std::vector<std::vector<TH2*>> signalHists;//[selection] [signal]
    TFile * postfitFile= 0;
    TFile * dataFile = 0;
    bool is1lep = true;
};


std::vector<TObject*> doDataPlot(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename, const bool is1l){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename,is1l);
    return a.makePlots();
}

std::vector<TObject*> doStatTest(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename, const std::string& outNamePrefix, const bool is1l){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename,is1l);
    return a.makeStatTest(outNamePrefix);
}

void doGlobChi2(std::vector<TObject*>& writeables, const std::string& limitBaseName, const double mH = 2000, const std::string& toyN = "toys"){
    TFile * fd = new TFile((limitBaseName+"/higgsCombineTest.GoodnessOfFit.mH"+ASTypes::flt2Str(1000)+".root").c_str(),"read");
    if(!fd){
        std::cout <<"No data file!"<<std::endl;
    }
    TTree * dataTree = 0;
    fd->GetObject("limit",dataTree);
    if(!dataTree){
        std::cout <<"No data tree!"<<std::endl;
    }
    double dl=0;
    dataTree->SetBranchAddress("limit",&dl);
    dataTree->GetEntry(0);
    double dataLimit = dl;
    fd->Close();

    TFile * ft = new TFile((limitBaseName+"/higgsCombineTest.GoodnessOfFit.mH"+ASTypes::flt2Str(1000)+"."+toyN+".root").c_str(),"read");
    if(!ft) return ;
    TTree * toyTree = 0;
    ft->GetObject("limit",toyTree);
    if(!toyTree) return ;

    double minT = toyTree->GetMinimum("limit");
    double maxT = toyTree->GetMaximum("limit");
    int nT = toyTree->GetEntries()/25;
    double width = (maxT - minT)/float(nT);
    TH1 * toyH = new TH1F("toyGlobChi2",";test statistic; N. of toys",nT, minT -width , maxT+width);
    toyTree->Draw("limit>>+toyGlobChi2","","goff");
    toyH->SetDirectory(0);
    TH1::AddDirectory(kFALSE);

    Plotter * p = new Plotter();
    p->addStackHist(toyH,"toy data");
    auto c = p->draw(false,"globChi2");
    c->SetTitle("globChi2");
    double yV =float(toyH->GetMaximum())*.75;
    TArrow * arrow = new TArrow(dataLimit,yV,dataLimit,0);
    arrow->Draw();
    writeables.push_back(c);
    ft->Close();
}

void doBiasTest(std::vector<TObject*>& writeables, const std::string& limitBaseName ){


    auto makeBiasPlot = [&](const std::string& filename, const std::string& hName, const double rV, TFitResultPtr& fitres) ->TH1*  {
        TFile * f = new TFile(filename.c_str(),"read");
        if(!f) return 0;
        TTree * tree = 0;
        f->GetObject("tree_fit_sb",tree);
        if(!tree){return 0;}
        double fr=0;
        double frerr=0;
        tree->SetBranchAddress("r",&fr);
        tree->SetBranchAddress("rErr",&frerr);
        std::vector<double> biases;
        for(unsigned int iE = 0; tree->GetEntry(iE); ++iE ){
            if(frerr == 0) continue;
            if(std::fabs((fr - rV)/frerr ) > 5) continue;
            biases.push_back( (fr - rV)/frerr   );
        }
        f->Close();

        TH1 * h = new TH1D(hName.c_str(),";signal strength pull",10, -5,5);
        h->SetDirectory(0);
        TH1::AddDirectory(kFALSE);
        for(auto b : biases) h->Fill(b);
        auto c = new TCanvas(hName.c_str(),hName.c_str());
        h->Draw();
        fitres = h->Fit("gaus","S");;
        writeables.push_back(c);
        //        std::sort(biases.begin(),biases.end(), [](const double a, const double b){return a < b;});
        //        double nToys = biases.size();
        //        h->SetDirectory(0);
        //        TH1::AddDirectory(kFALSE);
        //        for(auto b : biases) h->Fill(b);
        //        std::cout << hName <<" "<< biases[nToys*0.5]<<std::endl;
        //        return h;

        return h;

    };

    std::vector<std::pair<int,double>> massRs = {{1000,0.0600585},{1600,0.0239257},{2500,0.0141601}};
    std::vector<std::pair<int,double>> massRx2s = {{1000,0.120117},{1600,0.0478514},{2500,0.0283202}};
    std::vector<std::pair<int,double>> massRx5s = {{1000,0.3002925},{1600,0.1196285},{2500,0.0708005}};
    Plotter * p = new Plotter();
    Plotter * pw = new Plotter();

    std::vector<std::vector<double>> means(massRs.size());
    std::vector<std::vector<double>> meanErs(massRs.size());

    auto doSet = [&] (const std::string& label,const std::string& plotlabel, const std::vector<std::pair<int,double>>& massRs ){
        TGraphErrors * gr = new TGraphErrors();
        TGraphErrors * grw = new TGraphErrors();
        int iP = 0;

        for(unsigned int iM = 0; iM < massRs.size(); ++iM){
            const auto& mr = massRs[iM];
            TFitResultPtr fitres ;
            auto dist = makeBiasPlot(limitBaseName + "/biasInput_"+label+"_"+int2Str(mr.first)+".root",
                    "biasDist_"+label+"_"+int2Str(mr.first), mr.second,fitres);
            if(dist){
                auto pars = fitres->GetParams();
                auto errs = fitres->GetErrors();
                gr->SetPoint(iP,double(mr.first),pars[1]);
                gr->SetPointError(iP,0,errs[1]);
                grw->SetPoint(iP,double(mr.first),pars[2]);
                grw->SetPointError(iP,0,errs[2]);
                means[iM].push_back(pars[1]);
                meanErs[iM].push_back(errs[1]);
            }
            ++iP;
        }
        if(iP){
            p->addGraph(gr,plotlabel);
            pw->addGraph(grw,plotlabel);
        }
    };

    //   doSet("prefit"    ,"prefit b-model: r=excluded"   ,massRs);
    doSet("postfit"   ,"postfit b-model: r=excluded"  ,massRs);
    //   doSet("prefit_t2" ,"prefit b-model: r=2*excluded" ,massRx2s);
    doSet("postfit_t2","postfit b-model: r=2*excluded",massRx2s);
    doSet("postfit_t5","postfit b-model: r=5*excluded",massRx5s);

    p->setYTitle("signal strength bias");
    p->setXTitle((sigMCS.title).c_str());
    auto c = p->draw(false,"signalInjectTest_bias");
    c->SetTitle("signalInjectTest_bias");
    writeables.push_back(c);

    pw->setYTitle("signal strength width");
    pw->setXTitle((sigMCS.title).c_str());
    auto cw = pw->draw(false,"signalInjectTest_width");
    cw->SetTitle("signalInjectTest_width");
    writeables.push_back(cw);

    std::cout << std::endl << std::endl;
    for(unsigned int iM = 0; iM < massRs.size(); ++iM){
        std::cout << int2Str(massRs[iM].first) << "&";
        for(unsigned int iR = 0; iR < means[iM].size(); ++iR ){
            std::cout <<TString::Format("$%0.2f\\pm%0.2f$",means[iM][iR],meanErs[iM][iR]);
            if(iR+1 ==  means[iM].size()) std::cout <<" \\\\ \n";
            else std::cout <<" & ";
        }

    }
    std::cout << std::endl << std::endl;



}



TCanvas * makeUncPlot(const std::vector<std::string>& uncs, const std::string& label,const std::string& title, const RooFitResult * fit_s_res,const RooFitResult * fit_b_res, const RooArgSet* prefit){
    TH1 * fit_s_uncsize = 0;
    TH1 * fit_b_uncsize = 0;
    TGraphAsymmErrors * fit_s_uncs = 0;
    TGraphAsymmErrors * fit_b_uncs = 0;
    const RooArgList * fit_s = 0;
    const RooArgList * fit_b = 0;
    int nP=0;
    TH1 * prefit_uncsize = new TH1F((label+"_prefit").c_str(),"",uncs.size(),0,uncs.size());
    if(fit_s_res){
        fit_s_uncsize = new TH1F((label+"_fit_s_uncsize").c_str(),"",uncs.size(),0,uncs.size());
        fit_s_uncs = new TGraphAsymmErrors;
        fit_s = &fit_s_res->floatParsFinal();
    }
    if(fit_b_res){
        fit_b_uncsize = new TH1F((label+"_fit_b_uncsize").c_str(),"",uncs.size(),0,uncs.size());
        fit_b_uncs = new TGraphAsymmErrors;
        fit_b = &fit_b_res->floatParsFinal();
    }

    for(const auto& unc:uncs){
        const RooRealVar * nuis_p = ((const RooRealVar*)prefit->find(unc.c_str()));
        const RooRealVar * nuis_s = fit_s ? ((const RooRealVar*)fit_s->find(unc.c_str())) : 0;
        const RooRealVar * nuis_b = fit_b ? ((const RooRealVar*)fit_b->find(unc.c_str())) : 0;
        if(nuis_p==0) continue;

        auto fillPt = [&](TH1 * h, TGraphAsymmErrors* g, const RooRealVar *  v ){
            float SF = 1./nuis_p->getError();
            float x  = SF*(v ? v->getVal() : 0 );
            float eu = SF*(v ? v->getErrorHi() : 0);
            float ed = SF*(v ? v->getErrorLo() : 0);
            float e =  SF*(v ? v->getError() : 0);
            if(ed == 0) ed = eu;
            if(g) g->SetPoint(nP,nP+0.5,x);
            if(g) g->SetPointError(nP,0,0,std::fabs(ed),eu);
            if(g) g->GetXaxis()->SetBinLabel(nP+1,unc.c_str());
            if(h) h->SetBinContent(nP+1,x);
            if(h) h->SetBinError(nP+1,e);
            if(h) h->GetXaxis()->SetBinLabel(nP+1,unc.c_str());
        };
        if(fit_b) fillPt(fit_b_uncsize,fit_b_uncs,nuis_b);
        if(fit_s) fillPt(fit_s_uncsize,fit_s_uncs,nuis_s);
        fillPt(prefit_uncsize,0,nuis_p);
        nP++;
    }

    double topScale = (1/.5);

    TCanvas * c =new TCanvas(label.c_str(),label.c_str());
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1.0);
    pad1->SetNumber(1);
    pad1->SetBottomMargin(.025); // Upper and lower plot are joined
    pad1->SetTopMargin(.1);
    pad1->SetLeftMargin(.14);
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad

    auto * leg = new TLegend(0.65,0.7,0.93,0.89);
    leg->SetHeader(title.c_str(),"C");
    leg->SetFillColor(0)                           ;
    leg->SetTextFont(42)                           ;
    leg->SetBorderSize(1);
    leg->SetNColumns(2);


    TH1 * h_fit_e_s = fit_s ? (TH1*)fit_s_uncsize->Clone() : 0;
    TH1 * h_fit_e_b = fit_b ? (TH1*)fit_b_uncsize->Clone() : 0;

    prefit_uncsize->SetLineWidth(4)                   ;
    prefit_uncsize->SetTitle("Nuisance Paramaeters")  ;
    prefit_uncsize->SetLineColor(kBlack)         ;
    prefit_uncsize->SetFillColor(kGray)          ;
    prefit_uncsize->SetMaximum(1.9)                     ;
    prefit_uncsize->SetMinimum(-1.9)                    ;
    prefit_uncsize->SetMarkerSize(0.0)           ;

    prefit_uncsize->Draw("E2")                        ;
    prefit_uncsize->Draw("HIST same")                  ;
    prefit_uncsize->GetYaxis()->SetTitle("#theta");

    prefit_uncsize->GetXaxis()->SetLabelOffset(2);
    prefit_uncsize->GetYaxis()->SetLabelSize(prefit_uncsize->GetYaxis()->GetLabelSize()*topScale );
    prefit_uncsize->GetYaxis()->SetTitleOffset(.33/.66*prefit_uncsize->GetYaxis()->GetTitleOffset());
    prefit_uncsize->GetYaxis()->SetTitleSize(prefit_uncsize->GetYaxis()->GetTitleSize()*topScale);
    leg->AddEntry(prefit_uncsize,"Prefit","FL")       ;

    if(fit_b){
        fit_b_uncs->SetLineColor(kBlue)     ;
        fit_b_uncs->SetMarkerColor(kBlue)   ;
        fit_b_uncs->SetMarkerStyle(20)           ;
        fit_b_uncs->SetMarkerSize(1.0)           ;
        fit_b_uncs->SetLineWidth(2)              ;
        fit_b_uncs->Draw("EPsame") ;
        leg->AddEntry(fit_b_uncs,"B-only fit","EPL")     ;
    }
    if(fit_s){
        fit_s_uncs->SetLineColor(kRed)      ;
        fit_s_uncs->SetMarkerColor(kRed)    ;
        fit_s_uncs->SetMarkerStyle(20)           ;
        fit_s_uncs->SetMarkerSize(1.0)           ;
        fit_s_uncs->SetLineWidth(2)              ;
        fit_s_uncs->Draw("EPsame") ;
        leg->AddEntry(fit_s_uncs,"S+B fit"   ,"EPL")     ;
    }


    pad1->SetGridx()      ;
    pad1->RedrawAxis()    ;
    pad1->RedrawAxis("g") ;
    leg->Draw()    ;

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.5);
    pad2->SetNumber(2);
    pad2->SetTopMargin(.005);
    pad2->SetLeftMargin(.14);
    pad2->SetBottomMargin(.6);
    pad2->Draw();             // Draw the upper pad: pad1
    pad2->cd();               // pad1 becomes the current pad

    for(int iB = 1; iB <= prefit_uncsize->GetNbinsX(); ++iB){
        if(fit_s){
            h_fit_e_s->SetBinContent(iB,fit_s_uncsize->GetBinError(iB)/prefit_uncsize->GetBinError(iB)) ;
            h_fit_e_s->SetBinError(iB,0)                                                      ;
        }
        if(fit_b){
            std::cout << uncs[iB-1] <<" "<< fit_b_uncsize->GetBinError(iB)/prefit_uncsize->GetBinError(iB)<<std::endl;
            h_fit_e_b->SetBinContent(iB,fit_b_uncsize->GetBinError(iB)/prefit_uncsize->GetBinError(iB)) ;
            h_fit_e_b->SetBinError(iB,0)                                                      ;
        }
    }

    if(fit_b){                                                                               ;
    h_fit_e_b->SetLineColor(kBlue)                                               ;
    h_fit_e_b->GetYaxis()->SetTitle("#sigma_{#theta}/(#sigma_{#theta} prefit)")        ;
    h_fit_e_b->SetTitle("Nuisance Parameter Uncertainty Reduction")                   ;
    h_fit_e_b->SetMaximum(1.1)                                                        ;
    h_fit_e_b->SetMinimum(0)                                                          ;
    h_fit_e_b->SetLineWidth(4)                   ;
    h_fit_e_b->Draw("hist")                                                            ;
    }
    if(fit_s){
        h_fit_e_s->SetLineWidth(4)                   ;
        h_fit_e_s->SetLineColor(kRed)   ;
        h_fit_e_s->Draw("histsame")           ;
    }
    pad2->SetGridx()      ;
    pad2->RedrawAxis();
    pad2->RedrawAxis("g");

    float botScale = (1/.5);
    TAxis * yAxis=0;
    TAxis * xAxis=0;
    if(fit_b){
        yAxis = h_fit_e_b->GetYaxis();
        xAxis = h_fit_e_b->GetXaxis();
    } else     if(fit_s){
        yAxis = h_fit_e_s->GetYaxis();
        xAxis = h_fit_e_s->GetXaxis();
    }
    if(yAxis){
        yAxis->SetTitleSize(botScale*yAxis->GetTitleSize());
        xAxis->SetTitleSize(botScale*xAxis->GetTitleSize());
        yAxis->SetLabelSize(botScale*yAxis->GetLabelSize());
        xAxis->SetLabelSize(botScale*xAxis->GetLabelSize());
        yAxis->SetTitleOffset(.33/.66*yAxis->GetTitleOffset());
        xAxis->LabelsOption("v");
    }

    return c;
}

void doUncPlots(std::vector<TObject*>& writeables, const std::string& limitBaseName, REGION reg, bool doS = false){

    TFile * fd = new TFile((limitBaseName+"/fitDiagnostics.root").c_str(),"read");
    if(!fd){
        std::cout <<"No plots file!"<<std::endl;
        return;
    }
    RooFitResult * fit_s =0;
    RooFitResult * fit_b =0;
    RooArgSet* prefit = 0;
    fd->GetObject("nuisances_prefit",prefit);
    if(doS) fd->GetObject("fit_s",fit_s);
    fd->GetObject("fit_b",fit_b);


    auto addBN = [&](std::vector<std::string>& ns, const std::string& proc,const std::string& name ){
        for(const auto& b :btagCats){
            if(b == btagCats[BTAG_LMT]) continue;
            std::string sn = proc +"_"+name+"_"+b;
            ns.push_back(sn);
        }
    };
    auto addAllN = [&](std::vector<std::string>& ns, const std::string& proc,const std::string& name ){
        for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats){
            if(l == lepCats[LEP_EMU] ) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(p == purCats[PURE_I] ) continue;
            std::string sn = proc +"_"+name+"_"+l +"_"+b+"_"+p;
            ns.push_back(sn);

        }
    };

    std::vector<std::string> sigNs = {"yield","eff_mu","eff_e","tau21_PtDependence","tau21_eff","btag_eff","unclust","jes","jer","hbb_scale","hbb_res"};
    std::vector<std::string> hbbNs = {"hbb_scale","hbb_res"};
    addBN(hbbNs,bkgSels[BKG_QG],"PTX");
    addBN(hbbNs,bkgSels[BKG_QG],"OPTX");
    addBN(hbbNs,bkgSels[BKG_LOSTTW],"PTX");
    addBN(hbbNs,bkgSels[BKG_LOSTTW],"OPTX");

    std::vector<std::string> QGhhNs;
    addAllN(QGhhNs,bkgSels[BKG_QG],"PTY");
    addAllN(QGhhNs,bkgSels[BKG_QG],"OPTY");

    std::vector<std::string> tophhNs;
    addAllN(tophhNs,"top","res");
    addAllN(tophhNs,"top","scale");

    std::vector<std::string> topNormNs;
    addAllN(topNormNs,"top","norm");
    addBN(topNormNs,"top","mt_rel_scale");
    addBN(topNormNs,"top","lostmw_rel_scale");
    if(reg == REG_NONTOPCR){
        addBN(topNormNs,"top","tFrac");
        addBN(topNormNs,"top","lostFrac");
    } else {
        addBN(topNormNs,"top","wFrac");
        addBN(topNormNs,"top","lostFrac");
    }

    std::vector<std::string> qgNormNs;
    addAllN(qgNormNs,bkgSels[BKG_QG],"norm");



    writeables.push_back(makeUncPlot(qgNormNs,"qgNormNs","Q/G bkg. norm.",fit_s,fit_b,prefit));
    writeables.push_back(makeUncPlot(topNormNs,"topNormNs","Top bkg. norm.",fit_s,fit_b,prefit));
    writeables.push_back(makeUncPlot(QGhhNs,"QGhhNs","Q/G bkg. #it{m}_{HH}",fit_s,fit_b,prefit));
    writeables.push_back(makeUncPlot(tophhNs,"tophhNs","Top bkg. #it{m}_{HH}",fit_s,fit_b,prefit));
    writeables.push_back(makeUncPlot(hbbNs,"hbbNs","Bkg. #it{m}_{b#bar{b}}",fit_s,fit_b,prefit));
    writeables.push_back(makeUncPlot(sigNs,"sigNs","Signal",fit_s,fit_b,prefit));
}




void runPostFit(const std::string& inName, const std::string& outName, double fixR=0, int channel = 0){
    PostFitter fitter(inName,fixR);
    //        fitter.addSignal("radHH");

    fitter.addBkg(bkgSels[BKG_MT]);
    fitter.addBkg(bkgSels[BKG_MW]);
    fitter.addBkg(bkgSels[BKG_LOSTTW],"_opt");
    fitter.addBkg(bkgSels[BKG_QG],"_opt");

    fitter.addVariable(MOD_MJ);
    fitter.addVariable(MOD_MR);

    if (channel == 0 || channel == 1) {
        for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats) for(const auto& h :hadCuts){
            if(l == lepCats[LEP_EMU] ) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(p == purCats[PURE_I] ) continue;
            if(h != hadCuts[HAD_FULL] ) continue;
            const std::string wsName = "std_"+ l +"_"+b+"_"+p +"_"+h+"_13TeV";
            fitter.addCategory(l +"_"+b+"_"+p +"_"+h,wsName);
        }
    }
    if (channel == 0 || channel == 2) {
        for(const auto& l :dilepCats) for(const auto& b :btagCats) for(const auto& s :selCuts){
            if(l == dilepCats[LEP_INCL] ) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(s != selCuts[SEL_FULL] ) continue;
            const std::string wsName = "std_"+ l +"_"+b+"_"+s+"_13TeV";
            fitter.addCategory(l +"_"+b+"_"+s,wsName);
        }
    }

    fitter.doDataFit();
    fitter.doToys(10);
    fitter.write(outName);

}


void plotDataTests(int step = 0, int inreg = REG_SR, bool do1lep = true, const std::string limitBaseName = ""){
    REGION reg = REGION(inreg);

    std::vector<std::string> srList, srListTitles;
    if (do1lep) {
        srList = getSRList(reg);
        srListTitles = getSRListTitles(reg);
    } else {
        srList = getDilepSRList(reg);
        srListTitles = getDilepSRListTitles(reg);
    }

    std::string inName  =  "bkgInputs" ;
    std::string outName = limitBaseName +"/plots/";

    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
        outName=limitBaseName+"/plots/TopCR_";
    }
    else if(reg == REG_NONTOPCR){
        inName = "bkgInputsNonTopCR";
        hhFilename += "_NonTopCR";
        outName=limitBaseName+"/plots/NonTopCR_";
        btagCats = qgBtagCats;
    }

    std::string filename = inName +"/"+hhFilename;
    std::string postFitFilename = limitBaseName +"/postFit.root";


    if(step==0) {//run post fit
        runPostFit(limitBaseName+"/combined.root",postFitFilename,0,0); // may want option to do one channel at a time
    }


    if(step== 1){ //prefit
        if(outName.size())         outName += "prefit_dataComp";
        DataPlotPrefs hhPlot;
        hhPlot.modelType = MOD_PRE;
        hhPlot.addRatio = true;
        hhPlot.bins = {30,210,100,150};
        hhPlot.binInY = false;
        hhPlot.sels = srList;
        hhPlot.titles = srListTitles;

        writeables = doDataPlot(hhPlot,filename,postFitFilename,do1lep);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000};
        hbbPlot.binInY = true;

        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename,do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);

    }

    if(step== 2){ //postfit for AN

        bool blind=false;
        bool doRebin = true;

        if(outName.size())         {
            outName += "postfit_dataComp";
            if(inreg == REG_SR && !blind)
                outName += "_unblind";
            if(doRebin)
                outName += "_rebinned";
        }

        DataPlotPrefs hhPlot;
        hhPlot.modelType = MOD_POST;
        hhPlot.addType = MOD_PRE;
//        if(doRebin) hhPlot.rebinFactor = 4;
        if(doRebin) hhPlot.rebinFactor = 2;
        if(inreg == REG_SR){
            if(blind) hhPlot.bins = {30,100,100,150,210};
            else {
                hhPlot.bins = {30,210};
                hhPlot.signals = {"1000","2500"};
                hhPlot.signalTitles = {"1 TeV X_{spin-0}","2.5 TeV X_{spin-0}"};
            }
            hhPlot.minTop =0.2;
            hhPlot.maxTop =3000;
            if(doRebin) hhPlot.maxTop =10000;

        } else {
            hhPlot.bins = {30,210};
            hhPlot.minTop =0.2;
            hhPlot.maxTop =10000;
        }
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.titles = srListTitles;
        hhPlot.doLog = true;

        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        //        hhPlot.addData = false;
        writeables = doDataPlot(hhPlot,filename,postFitFilename,do1lep);

        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000};
        if(doRebin) hbbPlot.rebinFactor = 3;
        hbbPlot.binInY = true;
        hbbPlot.doLog = false;
        hbbPlot.minTop =400;
        hbbPlot.maxTop =400;
        if(inreg == REG_SR && blind) hbbPlot.blindRange ={100,150};
        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename,do1lep);

        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);
    }
    if(step ==3){ //statTest
        DataPlotPrefs hhTest;
        hhTest.modelType = MOD_POST;
        if(inreg == REG_SR){
            hhTest.bins = {30,100,100,150,210};
        } else {
            //            hhTest.bins = {30,100,100,150,210};
            hhTest.bins = {30,210};
        }
        hhTest.binInY = false;
        hhTest.sels = srList;
        hhTest.titles = srListTitles;
        hhTest.addErrorBars =true;
        writeables = doStatTest(hhTest,filename,postFitFilename,outName + "statTest",do1lep);

        DataPlotPrefs hbbTest = hhTest;
        hbbTest.binInY = true;
        //        hbbTest.rebins = {30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150};
        //        hbbTest.rebinFactor = -1;
        if(inreg == REG_SR){
            //            hbbTest.bins = {30,100,100,150,210};
        } else {
            //                        hhTest.bins = {30,210,30,100,150,210};
            //            hbbTest.bins = {700,4000};
            auto writeables2 = doStatTest(hbbTest,filename,postFitFilename,outName + "statTest",do1lep);
            writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        }

        Dummy d(outName+ "statTest_summary");

    }

    if(step==4){ //summary plots
        if(outName.size())         outName += "summaryPlots";
        doGlobChi2(writeables,limitBaseName);
        //        doUncPlots(writeables,limitBaseName,reg,reg!=REG_SR);
        doUncPlots(writeables,limitBaseName,reg,true);
        Dummy d(outName);
    }
    if(step==5){ //bias test
        if(outName.size())         outName += "biasTest";
        doBiasTest(writeables,limitBaseName);
        Dummy d(outName);
    }

    if(step== 6){ //postfit CR for Paper
        //        if(reg == REG_SR) return;
        if(outName.size())         outName += "postfitPaper_dataComp";
        DataPlotPrefs hhPlot;
        //Horrible
        doComp(reg,hhPlot.nToys,  postFitFilename,limitBaseName +"/postFit_comp.root",{});
        srList ={"emu_LMT_I_full"};
        srListTitles ={ "All categories"};


        hhPlot.modelType = MOD_POST;
        hhPlot.bins = {30,210};
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.titles = srListTitles;
        hhPlot.topTitle = reg == REG_TOPCR ? "t#bar{t} CR" : "q/g CR";

        hhPlot.minTop =0.5;
        hhPlot.maxTop =100000;
        hhPlot.rebinFactor=4;
        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        hhPlot.doLog = true;
        hhPlot.rebinFactor = 4;
        hhPlot.removeTrailingZeros = true;
        writeables = doDataPlot(hhPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.binInY = true;
        hbbPlot.bins = {700,4000};
        hbbPlot.doLog = false;
        hbbPlot.minTop =0;
        hbbPlot.maxTop =reg==REG_TOPCR ?  1100 :2000;
        hbbPlot.rebinFactor = 3;
        hbbPlot.removeTrailingZeros = false;
        auto writeables2 = doDataPlot(hbbPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);
    }

    if(step== 7){ //postfit SR for Paper
        //        if(reg != REG_TOPCR) return;
        if(outName.size())         outName += "postfitPaper_dataCompSR";
        DataPlotPrefs hhPlot;
        hhPlot.modelType = MOD_POST;
        hhPlot.bins = {30,210};
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.titles = srListTitles;
        hhPlot.signals = {"1000","2500"};
        hhPlot.signalTitles = {"1 TeV X_{spin-0}","2.5 TeV X_{spin-0}"};
        hhPlot.minTop =10;
        hhPlot.maxTop =10;

        float muMax = 3600; float muMin = 0.015;
        float elMax = 3400; float elMin = 0.003;

        hhPlot.maxTops = {elMax,elMax,elMax,elMax,elMax,elMax,muMax,muMax,muMax,muMax,muMax,muMax};
        hhPlot.minTops = {elMin,elMin,elMin,elMin,elMin,elMin,muMin,muMin,muMin,muMin,muMin,muMin};


        hhPlot.rebinFactor = 4;
        hhPlot.minBot =0.05;
        hhPlot.maxBot= 3.45;
        hhPlot.removeTrailingZeros = true;
        //        hhPlot.rebinFactor = -1;
        //        hhPlot.rebins = {700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,2000,2200,2400,3000,4000};

        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        hhPlot.doLog = true;
        writeables = doDataPlot(hhPlot,filename,postFitFilename,do1lep);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000};
        hbbPlot.binInY = true;
        hbbPlot.doLog = false;
        hbbPlot.minTop =0;
        hbbPlot.maxTop =0;
        hbbPlot.minBot =-1;
        hbbPlot.maxBot= -1;
        hbbPlot.maxTops.clear();
        hbbPlot.minTops.clear();
        hbbPlot.removeTrailingZeros = false;
//        hbbPlot.maxTops = {140,90,35, 20,25,15,200,100,40,30};
        hbbPlot.rebinFactor = 3;
        //        if(inreg == REG_SR) hbbPlot.blindRange ={100,150};
        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename,do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);
    }

    if(step== 8){ //postfit SR for Supp
        //        if(reg == REG_SR) return;
        if(outName.size())         outName += "postfitSupMat_dataComp";
        DataPlotPrefs hhPlot;
        hhPlot.signals = {"1000","2500"};
        hhPlot.signalTitles = {"1 TeV X_{spin-0}","2.5 TeV X_{spin-0}"};
        //Horrible
        doComp(reg,hhPlot.nToys,  postFitFilename,limitBaseName +"/postFit_comp.root"
                ,{"800","1000","1200","2500"});
        srList ={"emu_LMT_I_full"};
        srListTitles ={ "All categories"};


        hhPlot.modelType = MOD_POST;
        hhPlot.bins = {30,210,210,100,150};
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.titles = srListTitles;
        hhPlot.isSupp = true;
        hhPlot.botTitles = {"","","","100 < #it{m}_{b#bar{b}} < 150 GeV"};

        hhPlot.minTop =0.1;
        hhPlot.maxTop =6700;

        hhPlot.minBot =0.05;
        hhPlot.maxBot= 1.95;
        hhPlot.removeTrailingZeros = true;

        hhPlot.rebinFactor=4;
        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        hhPlot.doLog = true;
        hhPlot.rebinFactor = 4;
        writeables = doDataPlot(hhPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.binInY = true;

        hbbPlot.doLog = false;
        hbbPlot.minTop =0;
        hbbPlot.maxTop =0;
        hbbPlot.minBot =-1;
        hbbPlot.maxBot= -1;
        hbbPlot.removeTrailingZeros = false;
        hbbPlot.rebinFactor = 3;
        auto writeables2 = doDataPlot(hbbPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );

        hbbPlot.bins = {700,1000};
        hbbPlot.botTitles = {"0.7 < #it{m}_{HH} < 1 TeV"};
        hbbPlot.signals = {"800"};
        hbbPlot.signalTitles = {"0.8 TeV X_{spin-0}"};
        writeables2 = doDataPlot(hbbPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );

        hbbPlot.bins = {1000,2000};
        hbbPlot.botTitles = {"1 < #it{m}_{HH} < 2 TeV"};
        hbbPlot.signals = {"1200"};
        hbbPlot.signalTitles = {"1.2 TeV X_{spin-0}"};
        writeables2 = doDataPlot(hbbPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        hbbPlot.bins = {2000,4000};
        hbbPlot.botTitles = {"2 < #it{m}_{HH} < 4 TeV"};
        hbbPlot.signals = {"2500"};
        hbbPlot.signalTitles = {"2.5 TeV X_{spin-0}"};
        hbbPlot.maxTop =20.0;
        writeables2 = doDataPlot(hbbPlot,filename,limitBaseName +"/postFit_comp.root",do1lep);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );


        Dummy d(outName);
    }




    if(step== 9){ //postfit SR for Supp 1
        if(outName.size())         outName += "postfitSupMat_dataComp";
        DataPlotPrefs hhPlot;
        hhPlot.signals = {"1000","2500"};
        hhPlot.signalTitles = {"1 TeV X_{spin-0}","2.5 TeV X_{spin-0}"};
        hhPlot.modelType = MOD_POST;
        hhPlot.bins = {100,150};
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.titles = srListTitles;
        hhPlot.isSupp = true;
        hhPlot.botTitles = {"100 < #it{m}_{b#bar{b}} < 150 GeV"};
        hhPlot.removeTrailingZeros = true;
        hhPlot.minTop =10;
        hhPlot.maxTop =10;
        float muMax = 3600; float muMin = 0.015;
        float elMax = 3400; float elMin = 0.003;

        hhPlot.maxTops = {elMax,elMax,elMax,elMax,elMax,elMax,muMax,muMax,muMax,muMax,muMax,muMax};
        hhPlot.minTops = {elMin,elMin,elMin,elMin,elMin,elMin,muMin,muMin,muMin,muMin,muMin,muMin};


        hhPlot.minBot =0.05;
        hhPlot.maxBot= 3.45;
        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        hhPlot.doLog = true;
        hhPlot.rebinFactor = 4;
        writeables = doDataPlot(hhPlot,filename,postFitFilename,do1lep);
        Dummy d(outName);


    }









}







