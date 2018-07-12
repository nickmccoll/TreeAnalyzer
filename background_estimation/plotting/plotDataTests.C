#include "../predTools/PostFitter.h"
#include "../predTools/StatTester.h"
#include <vector>
#include "TFile.h"
#include "TStyle.h"
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
    bool addErrorBars = false;
    unsigned int nToys = 100;
};

struct HistContainter {
    HistContainter() {
        bkg2D.resize(BKG_MT+1); for(auto* b :bkg2D) b=0;
        bkg.resize(BKG_MT+1); for(auto* b :bkg) b=0;
        toys2D.reserve(100);
    }
    std::vector<TH2*> bkg2D;
    TH2* tot2D=0;
    TH2* add2D=0;
    TH2* data2D=0;
    std::vector<TH2*> toys2D;

    std::vector<TH1*> bkg;
    TH1* tot=0;
    TH1* add=0;
    TH1* data=0;
    TGraphAsymmErrors* toyErr=0;
};




class DataPlotter {
public:
    DataPlotter(const DataPlotPrefs& plotPrefs, const std::string& inputPrefix, const std::string& postFitFilename ) : prefs(plotPrefs)
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
}

    void getToyModels(const std::string& sel, HistContainter& cont){
        for(unsigned int iT = 1; iT <= prefs.nToys; ++iT ){
            TH2 * th = 0;
            postfitFile->GetObject((std::string("toyModel_")+ASTypes::int2Str(iT)+"_"+ sel+"__"+MOD_MJ+"_"+MOD_MR).c_str(),th);
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
                postfitFile->GetObject((std::string("postfit_")+ bkgSels[iT]+"_"+sel+"__"+MOD_MJ+"_"+MOD_MR).c_str(),h);
            }
            bkgs.push_back(h);
        }
        return bkgs;
    }
    void getBackgrounds(const std::string& sel, HistContainter& cont){
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

        auto combine = [](std::vector<TH2*>& comps) -> TH2* {
            TH2* total = 0;
            for(unsigned int iT = 0; iT < bkgSels.size(); ++iT){
                if(comps[iT] ==0) continue;
                if(total==0) total = (TH2*)comps[iT]->Clone();
                else total->Add(comps[iT]);
            }
            return total;
        };

        if(prefs.modelType == MOD_POST) cont.bkg2D = postFits;
        if(prefs.modelType == MOD_PRE)  cont.bkg2D = preFits ;
        if(prefs.modelType == MOD_MC)   cont.bkg2D = mcs     ;

        cont.tot2D = combine(cont.bkg2D);

        if(prefs.addType == MOD_POST)cont.add2D=combine(postFits);
        if(prefs.addType == MOD_PRE) cont.add2D=combine(preFits);
        if(prefs.addType == MOD_MC) cont.add2D=combine(mcs);
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
    std::pair<int,int> getBins(const HistContainter& cont,const int iB){
        const TAxis * ax =0;
        for(const auto* h: cont.bkg2D){
            if(h==0) continue;
            ax = prefs.binInY ? h->GetYaxis() : h->GetXaxis();
            break;
        }
        if(!ax && cont.add2D) ax = prefs.binInY ? cont.add2D->GetYaxis() : cont.add2D->GetXaxis();
        if(!ax && cont.data) ax = prefs.binInY ? cont.data->GetYaxis() : cont.data->GetXaxis();
        int binL = ax->FindFixBin(prefs.bins[iB]);
        int binH = ax->FindFixBin(prefs.bins[iB+1]) -1;
        return std::make_pair(binL,binH);
    }

    void get1DHists(const std::string& s, const unsigned int iB, HistContainter& cont){
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
            auto h = processH2(cont.data2D,"data");
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
                cont.toyErr->SetPoint(iB-1,x,y);
                cont.toyErr->SetPointError(iB-1,0,0,std::max(m -eL,0.),std::max(eH-m,0.));
            }
            for(auto * t :toyProj) {t->SetDirectory(0); delete t;}
        }
    }


    std::vector<TObject*> makePlots() {
        std::vector<TObject*> writeables;
        for(const auto& s :prefs.sels){
            HistContainter cont;
            getBackgrounds(s,cont);
            if(prefs.modelType == MOD_POST && prefs.addErrorBars)
                getToyModels(s,cont);

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

                if(cont.toyErr){
                    int fillColor = 2;
                    cont.toyErr->SetFillColor(fillColor);
                    cont.toyErr->SetFillStyle(3352);
                    gStyle->SetHatchesLineWidth(1);
                    gStyle->SetHatchesSpacing(.5);
                    p->addGraph(cont.toyErr,"bkg. unc.",fillColor,1,1,20,1,false,true,false,"3");
                }
                if(cont.add)p->addHistLine(cont.add,modTitles[prefs.addType],kBlue);
                if(cont.data) p->addHist(cont.data,"data",kBlack,1,2,20,0.5,true,true,true);
                for(unsigned int iH = 0; iH < cont.bkg.size(); ++iH){
                    if(cont.bkg[iH]) p->addStackHist(cont.bkg[iH],bkgSels[iH].c_str());
                }
                p->setUnderflow(false);
                p->setOverflow(false);
                double binWidth = 10;
                for(const auto* h : cont.bkg){
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

    std::vector<TObject*>  makeStatTest(const std::string outNamePrefix, bool useBuiltInToys = true) {
        std::vector<TH1*> dataTSs(prefs.bins.size() -1,0);
        std::vector<TGraphAsymmErrors*> toyTSs(prefs.bins.size() -1,0);
        std::vector<TObject*> writeables;

        for(unsigned int iS = 0; iS < prefs.sels.size();++iS){
            const auto& s = prefs.sels[iS];
            HistContainter cont;
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
            p->addGraph(toyTSs[iB],"toy data");
            p->addHist(dataTSs[iB],"data",-1,1,4,20,1,true,false);
            p->setXTitle(" ");
            p->setYTitle("test statistic");
            auto c = p->draw(false,(prefs.binInY ? hhMCS : hbbMCS) +"_"+flt2Str(prefs.bins[iB]) +"_"+flt2Str(prefs.bins[iB+1]) +"_testStat");
            for(unsigned int iS = 0; iS < prefs.sels.size();++iS){
                p->xAxis()->SetBinLabel(iS+1,prefs.sels[iS].c_str());
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
    TFile * postfitFile= 0;
    TFile * dataFile = 0;
};


std::vector<TObject*> doDataPlot(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename ){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename);
    return a.makePlots();
}

std::vector<TObject*> doStatTest(const DataPlotPrefs& dataPlot, const std::string& inputPrefix, const std::string& postFitFilename, const std::string& outNamePrefix ){
    DataPlotter a(dataPlot,inputPrefix,postFitFilename);
    return a.makeStatTest(outNamePrefix);
}

void runPostFit(const std::string& inName, const std::string& outName, double fixR=0){
    PostFitter fitter(inName,fixR);
    //        fitter.addSignal("radHH");
    fitter.addBkg(bkgSels[BKG_MT]);
    fitter.addBkg(bkgSels[BKG_MW]);
    fitter.addBkg(bkgSels[BKG_LOSTTW],"_opt");
    fitter.addBkg(bkgSels[BKG_QG],"_opt");
    fitter.addVariable(MOD_MJ);
    fitter.addVariable(MOD_MR);
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(p == purCats[PURE_I] ) continue;
        if(h != hadCuts[HAD_FULL] ) continue;
        const std::string wsName = "std_"+ l +"_"+b+"_"+p +"_"+h+"_13TeV";
        fitter.addCategory(l +"_"+b+"_"+p +"_"+h,wsName);
    }

    fitter.doDataFit();
    fitter.doToys(100);
    fitter.write(outName);
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
        limitBaseName +="_QGCR";
        btagCats = qgBtagCats;

    }
    std::string outName = limitBaseName +"/plots/";
    std::string filename = inName +"/"+hhFilename;
    std::string postFitFilename = limitBaseName +"/postFit.root";

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
            hhPlot.bins = {30,210};
        }
        hhPlot.binInY = false;
        hhPlot.sels =  srList;
        hhPlot.addRatio = true;
        hhPlot.addErrorBars = true;
        //        hhPlot.addData = false;
        writeables = doDataPlot(hhPlot,filename,postFitFilename);
        DataPlotPrefs hbbPlot = hhPlot;
        hbbPlot.bins = {700,4000};
        hbbPlot.binInY = true;
        if(inreg == REG_SR) hbbPlot.blindRange ={100,150};
        auto writeables2 = doDataPlot(hbbPlot,filename,postFitFilename);
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
        hhTest.addErrorBars =true;
        writeables = doStatTest(hhTest,filename,postFitFilename,outName + "statTest");

        DataPlotPrefs hbbTest = hhTest;
        hbbTest.binInY = true;
//        hbbTest.rebins = {30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150};
//        hbbTest.rebinFactor = -1;
        if(inreg == REG_SR){
        //            hbbTest.bins = {30,100,100,150,210};
        } else {
            //            hhTest.bins = {30,210,30,100,150,210};
            hbbTest.bins = {700,4000};
            auto writeables2 = doStatTest(hbbTest,filename,postFitFilename,outName + "statTest");
            writeables.insert( writeables.end(), writeables2.begin(), writeables2.end() );
        }

        Dummy d(outName+ "statTest_summary.root");

    }







}







