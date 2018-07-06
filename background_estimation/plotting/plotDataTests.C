#include "../predTools/PostFitter.h"
#include "../predTools/StatTester.h"
#include <vector>
#include "TFile.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;
std::vector<TObject*> writeables;

class Dummy {
public:
    Dummy(const std::string& outName = "") : outName(outName) {};
    ~Dummy() {
        if(outName.size()){
            TFile * f = new TFile(outName.c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
            }
            f->Close();
        }
    }
    std::string outName;
};


enum ModelType {MOD_NONE, MOD_MC, MOD_PRE, MOD_POST};
std::string modTitles[] = {"NONE","MC","prefit","postfit"};
struct DataPlotPrefs {
    bool addData = true;
    std::pair<double,double> blindRange = std::pair<double,double>(0,0);
    ModelType modelType = MOD_PRE;
    ModelType addType   = MOD_NONE;
    bool plotIntegral = false;
    std::vector<double> bins; // if bin+1 < bin, it will skip [bin,bin+1], if bin +1 == bin it will skip [bin,bin+1] and [bin+1,bin+2]
    bool binInY = true;
    std::vector<std::string> sels;


    int rebinFactor = 1; //-1 means go to rebins...1 means dont rebin...2+ is standard grouping
    std::vector<double> rebins;
    bool doLog = false;
    bool addRatio = false;
};




class DataPlotter {
public:
    DataPlotter(const DataPlotPrefs& plotPrefs, const std::string& inputPrefix, const std::string& postFitFilename ) : prefs(plotPrefs)
{
        if(prefs.addData) dataFile  = new TFile((inputPrefix+"_data_distributions.root").c_str(),"read");
        if(prefs.modelType == MOD_POST || prefs.addType == MOD_POST )
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
                postfitFile->GetObject((bkgSels[iT]+"_"+sel).c_str(),h);
            }
            bkgs.push_back(h);
        }
        return bkgs;
    }
    //0->(nBkgs-1) are the indexed backgrounds...nBKGS is the total additional...nBKGS+1 is data
    std::vector<TH2*>  getBackgrounds(const std::string& sel){
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
        std::vector<TH2*> outHists;
        if(prefs.modelType == MOD_POST) outHists = postFits;
        if(prefs.modelType == MOD_PRE) outHists = preFits;
        if(prefs.modelType == MOD_MC) outHists = mcs;

        auto combine = [](std::vector<TH2*>& comps) -> TH2* {
            TH2* total = 0;
            for(unsigned int iT = 0; iT < bkgSels.size(); ++iT){
                if(comps[iT] ==0) continue;
                if(total==0) total = (TH2*)comps[iT]->Clone();
                else total->Add(comps[iT]);
            }
            return total;
        };


        if(prefs.addType == MOD_POST) outHists.push_back(combine(postFits));
        if(prefs.addType == MOD_PRE) outHists.push_back(combine(preFits));
        if(prefs.addType == MOD_MC) outHists.push_back(combine(mcs));
        if(prefs.addType == MOD_NONE) outHists.push_back(0);

        return outHists;
    }
    std::vector<TH1*> get1DHists(std::vector<TH2*>& hists2D, const std::string& s, const unsigned int iB){
        std::vector<TH1*> outHists(hists2D.size(), 0);
        //get the binning axis
        const TAxis * ax =0;
        for(const auto* h: hists2D){
            if(h==0) continue;
            ax = prefs.binInY ? h->GetYaxis() : h->GetXaxis();
            break;
        }
        int binL = ax->FindFixBin(prefs.bins[iB]);
        int binH = ax->FindFixBin(prefs.bins[iB+1]) -1;


        auto processH2 = [&](TH2* h, const std::string& title)->TH1*{
            auto h1 = prefs.binInY ? h->ProjectionX(  (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH)
                    :  h->ProjectionY( (s + "_" + title+"_"+int2Str(iB)).c_str(),binL,binH);

            if(prefs.rebinFactor > 1){
                h1->Rebin(prefs.rebinFactor);
            } else if(prefs.rebinFactor < 1){
                h1->Rebin(prefs.rebins.size() -1,"", &prefs.rebins[0]);
            }

            if(prefs.plotIntegral) return PlotTools::getIntegral(h1,true,false);
            return h1;
        };

        for(unsigned int iH = 0; iH < hists2D.size()-2; ++iH){
            if(hists2D[iH]==0) continue;
            TH1 * h = processH2(hists2D[iH],bkgSels[iH]);
            for(int iX = 1; iX <= h->GetNbinsX(); ++iX)h->SetBinError(iX,0);
            outHists[iH] = h;
        }

        //add additional
        if(hists2D[hists2D.size()-2] ){
            auto h = processH2(hists2D[hists2D.size()-2],modTitles[prefs.addType]);
            outHists[hists2D.size()-2] = h;
        }

        //add data
        if(hists2D.back()){
            auto h = processH2(hists2D.back(),"data");
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
            outHists.back() = h;
        }
        return outHists;
    }


    std::vector<TObject*> makePlots() {
        std::vector<TObject*> writeables;
        for(const auto& s :prefs.sels){
            auto hists2D = getBackgrounds(s);
            //add in the data to nBKS
            TH2 * h = 0;
            if(prefs.addData){
                dataFile->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
            }
            hists2D.push_back(h);

            //Individual bins
            for(unsigned int iB = 0; iB + 1 < prefs.bins.size(); ++iB){
                if(prefs.bins[iB+1] <= prefs.bins[iB]){
                    if(prefs.bins[iB+1] == prefs.bins[iB]) ++iB;
                    continue;
                }
                auto hists1D = get1DHists(hists2D,s,iB);

                const std::string plotTitle = s +"_"+(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]);


                Plotter * p = new Plotter();
                if(hists1D[hists1D.size()-2])p->addHistLine(hists1D[hists1D.size()-2],modTitles[prefs.addType],kBlue);
                if(hists1D.back()) p->addHist(hists1D.back(),"data",kBlack,1,4,20,1,true,true,true);
                for(unsigned int iH = 0; iH < hists1D.size()-2; ++iH){
                    if(hists1D[iH]) p->addStackHist(hists1D[iH],bkgSels[iH].c_str());
                }
                p->setUnderflow(false);
                p->setOverflow(false);
                double binWidth = 10;
                for(const auto* h : hists1D){
                    if(h==0) continue;
                    binWidth = h->GetBinWidth(1);
                    break;
                }
                p->setXTitle( (prefs.binInY ? hbbMCS : hhMCS) .title.c_str());
                p->setYTitle((std::string("N. of events / ") + flt2Str(binWidth) +" [GeV]").c_str() );
                p->setCMSLumi();

                if(prefs.doLog) p->setMinMax(0.1,1000 );
                if(prefs.addRatio){
                    p->setBotMinMax(0,2);
                    p->setYTitleBot((std::string("N/N(") + modTitles[prefs.modelType] +")").c_str());
                    auto * c = p->drawSplitRatio(-1,"stack",false,false,plotTitle.c_str());
                    if(prefs.doLog) c->GetPad(1)->SetLogy();
                    c->GetPad(1)->Update();
                    writeables.push_back(c);
                } else {
                    auto * c = p->draw(false,plotTitle.c_str());
                    if(prefs.doLog)c->SetLogy();
                    c->Update();
                    writeables.push_back(c);
                }
            }
        }
        return writeables;
    }

    void makeStatTest(const std::string outNamePrefix, const int nToys = 100, bool saveToys = false) {
        for(const auto& s :prefs.sels){
            auto hists2D = getBackgrounds(s);
            //add in the data to nBKS
            TH2 * h = 0;
            if(prefs.addData){
                dataFile->GetObject(("data_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),h);
            }
            hists2D.push_back(h);

            //Individual bins
            for(unsigned int iB = 0; iB + 1 < prefs.bins.size(); ++iB){
                if(prefs.bins[iB+1] <= prefs.bins[iB]){
                    if(prefs.bins[iB+1] == prefs.bins[iB]) ++iB;
                    continue;
                }
                auto hists1D = get1DHists(hists2D,s,iB);
                const std::string plotTitle = s +"_"+(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]);
                if(hists1D[hists1D.size() -2]==0) continue;
                if(hists1D.back()==0) continue;
                std::cout <<plotTitle <<std::endl;
                StatTester((TH1D*)hists1D[hists1D.size() -2],(TH1D*)hists1D.back(),nToys,saveToys,outNamePrefix+"_"+plotTitle+".root" );
            }
        }
    }


    DataPlotPrefs prefs;
    std::vector<TFile*> mcFiles;
    std::vector<TFile*> prefitFiles;
    TFile * postfitFile= 0;
    TFile * dataFile = 0;
};


std::vector<TObject*> doDataPlot(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename ){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename);
    return a.makePlots();
}

void doStatTest(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename, const std::string& outNamePrefix ){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename);
    a.makeStatTest(outNamePrefix,1000,false);
}

void runPostFit(const std::string& inName, const std::string& outName, double fixR=-1){



    auto fitter = PostFitter(inName);
    fitter.fix(MOD_MS,1500);
    if(fixR >= 0) fitter.fix("r",fixR);
    fitter.doFit();
    fitter.addSignal("radHH");
    fitter.addBkg(bkgSels[BKG_MT]);
    fitter.addBkg(bkgSels[BKG_MW]);
    fitter.addBkg(bkgSels[BKG_LOSTTW],"_opt");
    fitter.addBkg(bkgSels[BKG_QG],"_opt");

    TFile * f = new TFile(outName.c_str(),"RECREATE");

    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(p == purCats[PURE_I] ) continue;
        if(h != hadCuts[HAD_FULL] ) continue;

        std::string category = "std_"+ l +"_"+b+"_"+p +"_"+h+"_13TeV";
        auto contribs = fitter.getContributions();
        for(const auto& c : contribs){
            auto* histo = fitter.get2DHistogram(MOD_MJ,MOD_MR,category,c.name);
            if(!histo) continue;
            f->cd();
            histo->Write((c.name+"_"+l +"_"+b+"_"+p +"_"+h).c_str());
        }
    }
    f->Close();
    delete f;
}


void plotDataTests(int step = 0, int inreg = REG_SR,  std::string limitBaseName = ""){
    REGION reg = REGION(inreg);

    std:: string inName =  "bkgInputs" ;
    auto srList = getSRList(reg);

    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
        limitBaseName +="_TopCR";
    }
    else if(reg == REG_QGCR){
        inName =  "bkgInputsQGCR";
        hhFilename +="_QGCR";
        limitBaseName +="_TopCR";

    }
    std::string outName = limitBaseName +"/plots/";
    std::string filename = inName +"/"+hhFilename;
        std::string postFitFilename = limitBaseName +"/postFit.root";
//    std::string postFitFilename = limitBaseName +"/postFitFake.root";

    if(step==0) {//run post fit
        runPostFit(limitBaseName+"/combined.root",postFitFilename,0);
    }


    if(step== 1){ //prefit
        if(outName.size())         outName += "prefit_dataComp.root";

        DataPlotPrefs hhPlot;
        hhPlot.modelType = MOD_PRE;
        hhPlot.bins = {30,210,100,150};
        hhPlot.binInY = false;
        hhPlot.sels = srList;
        writeables = doDataPlot(hhPlot,filename,postFitFilename);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000};
        hbbPlot.binInY = true;
        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);
    }

    if(step== 2){ //postfit
        if(outName.size())         outName += "postfit_dataComp.root";

        DataPlotPrefs hhPlot;
        hhPlot.modelType = MOD_POST;
        hhPlot.addType = MOD_PRE;
        if(inreg == REG_SR){
            hhPlot.bins = {30,100,100,150,210};
        } else {
            hhPlot.bins = {30,210,30,100,150,210};
        }
        hhPlot.binInY = false;
        hhPlot.sels = srList;
        hhPlot.sels = {"e_T_LP_full"};
        hhPlot.addRatio = true;
        //        hhPlot.doLog = true;
        //        hhPlot.addData = false;
        writeables = doDataPlot(hhPlot,filename,postFitFilename);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000,700,1000,850,1000};
        hbbPlot.binInY = true;
        if(inreg == REG_SR) hbbPlot.blindRange ={100,150};
        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename);
        writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        Dummy d(outName);
    }
    if(step ==3){ //statTest
        DataPlotPrefs hhTest;
        hhTest.addType = MOD_POST;
        hhTest.modelType = MOD_NONE;
        if(inreg == REG_SR){
            hhTest.bins = {30,100,100,150,210};
        } else {
//            hhTest.bins = {30,210,30,100,150,210};
            hhTest.bins = {30,210};
        }
        hhTest.binInY = false;
        hhTest.sels = srList;
        doStatTest(hhTest,filename,postFitFilename,outName + "statTest");
    }







}







