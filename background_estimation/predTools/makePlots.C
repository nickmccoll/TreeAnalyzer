
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TTreeFormula.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooNDKeysPdf.h"
#include "RooBinning.h"
#include "TMath.h"
using namespace TAna;


struct PlotVar{

    PlotVar(std::string varName, std::string varTitle, std::string varSel,  int nBins, double min, double max) :
         varName(varName),varTitle(varTitle),varSel(varSel),nBins(nBins),min(min), max(max){}

    PlotVar(std::string varName, std::string varTitle, std::string varSel , const std::vector<double>& bins ) :
        varName(varName),varTitle(varTitle),varSel(varSel), bins(bins){}

    PlotVar(std::string varName, std::string varTitle, std::string varSel,  int nBins, double min, double max,
            std::string varYName, std::string varYTitle, std::string varYSel,  int nBinsY, double minY, double maxY
            ) :
         varName(varName),varTitle(varTitle),varSel(varSel),nBins(nBins),min(min), max(max),
         varYName(varYName),varYTitle(varYTitle),varYSel(varYSel),nBinsY(nBinsY),minY(minY), maxY(maxY){}

    std::string varName = "";
    std::string varTitle = "";
    std::string varSel  = "";
    std::vector<double> bins;
    //nBins < 1 if variable binning (use bins vector)
    int nBins =-1;
    double min = -1;
    double max = -1;


    std::string varYName = "";
    std::string varYTitle = "";
    std::string varYSel   = ""; //blank if 1D
    std::vector<double> binsY;
    //nBins < 1 if variable binning (use bins vector)
    int nBinsY =-1;
    double minY = -1;
    double maxY = -1;
};

struct PlotSamp{
    std::string name ="";
    std::string sel ="";
    bool isVec = false;
    PlotSamp(const std::string& name,const std::string& sel, bool isVec = false) :
        name(name), sel(sel), isVec(isVec){}
};

//Sample title, sample cut
typedef std::pair<std::string,std::string> PlotSel;


class MakePlots : public BaseTreeAnalyzer {
public:

    MakePlots(std::string fileName,  std::string outFileName, const std::vector<PlotSamp>& samps, const std::vector<PlotSel>& sels,
            const std::vector<PlotVar>& vars, std::string baseSel = "1.0", std::string baseWeight = "1.0", std::string prefix = "",std::string treeName = "treeMaker/Events" ) : BaseTreeAnalyzer(fileName,treeName,2),
            outFileName(outFileName), samps(samps), sels(sels), vars(vars), prefix(prefix.size() ? prefix + "_ ": prefix ){

        tree.getTree()->SetBranchStatus("*",1);
        sampSels.resize(samps.size());
        for(unsigned int iS = 0; iS < samps.size(); ++iS){
            sampSels[iS].reset(new TTreeFormula(TString::Format("sampSel_%u",iS), samps[iS].sel.c_str(),tree.getTree()));
        }
        selSels.resize(sels.size());
        for(unsigned int iS = 0; iS < sels.size(); ++iS){
            selSels[iS].reset(new TTreeFormula(TString::Format("selSel_%u",iS), sels[iS].second.c_str(),tree.getTree()));
        }
        varSels.resize(vars.size());
        varYSels.resize(vars.size());
        for(unsigned int iV = 0; iV < vars.size(); ++iV){
            varSels[iV].reset(new TTreeFormula(TString::Format("varSel_%u",iV), vars[iV].varSel.c_str(),tree.getTree()));
            if(vars[iV].varYName.size()){
                varYSels[iV].reset(new TTreeFormula(TString::Format("varYSel_%u",iV), vars[iV].varYSel.c_str(),tree.getTree()));
            }
        }

        stdSel.reset(new TTreeFormula("stdSel", TString::Format("%s*(%s)",baseWeight.c_str(),baseSel.c_str()),tree.getTree()));
        analyze();
    }
    void analyze(int reportFrequency = 10000, int numEvents = -1, int startEvent = -1) override {
        auto * evtAna = setupEventAnalyzer();
        evtAna->analyzeEvent(this,reportFrequency,numEvents,startEvent);
        delete evtAna;
        plotter.write(outFileName);
    }
    void loadVariables() override{};
    bool runEvent() override {

        double bw   = stdSel->EvalInstance();
        if(bw==0) return false;

        for(unsigned int iS = 0; iS < samps.size(); ++iS){
            if(samps[iS].isVec) sampSels[iS]->GetNdata();
            double s   = sampSels[iS]->EvalInstance();
            if(s == 0) continue;

            for(unsigned int iSe = 0; iSe < sels.size(); ++iSe){
                double sel   = selSels[iSe]->EvalInstance();
                if(sel == 0) continue;
                std::string thisPre = prefix + samps[iS].name + std::string("_") +sels[iSe].first;

                for(unsigned int iV = 0; iV < vars.size(); ++iV){
                    double v   = varSels[iV]->EvalInstance();

                    if(varYSels[iV]){
                        double vy   = varYSels[iV]->EvalInstance();
                        TString varString = vars[iV].varName + std::string("_") + vars[iV].varYName;
                        std::string varTitle = vars[iV].varTitle + vars[iV].varYTitle;
                        if(vars[iV].nBins >= 0 && vars[iV].nBinsY >= 0){
                            plotter.getOrMake2DPre(thisPre.c_str(), varString,  vars[iV].varTitle.c_str(), vars[iV].nBins, vars[iV].min, vars[iV].max, vars[iV].nBinsY, vars[iV].minY, vars[iV].maxY   )
                            ->Fill(v,vy,bw*s*sel);
                        } else if(vars[iV].nBins >= 0 && vars[iV].nBinsY < 0){
                            plotter.getOrMake2DPre(thisPre.c_str(), varString,  vars[iV].varTitle.c_str(), vars[iV].nBins, vars[iV].min, vars[iV].max, vars[iV].binsY.size(),&vars[iV].binsY[0]  )
                            ->Fill(v,vy,bw*s*sel);
                        } else if(vars[iV].nBins < 0 && vars[iV].nBinsY >= 0){
                            plotter.getOrMake2DPre(thisPre.c_str(), varString,  vars[iV].varTitle.c_str(), vars[iV].bins.size(),&vars[iV].bins[0], vars[iV].nBinsY, vars[iV].minY, vars[iV].maxY )
                            ->Fill(v,vy,bw*s*sel);
                        } else {
                            plotter.getOrMake2DPre(thisPre.c_str(), varString,  vars[iV].varTitle.c_str(), vars[iV].bins.size(),&vars[iV].bins[0], vars[iV].binsY.size(),&vars[iV].binsY[0] )
                            ->Fill(v,vy,bw*s*sel);
                        }
                    } else {
                        if(vars[iV].nBins >= 0)
                            plotter.getOrMake1DPre(thisPre.c_str(), vars[iV].varName.c_str(),  vars[iV].varTitle.c_str(), vars[iV].nBins, vars[iV].min, vars[iV].max   )
                            ->Fill(v,bw*s*sel);
                        else
                            plotter.getOrMake1DPre(thisPre.c_str(), vars[iV].varName.c_str(),  vars[iV].varTitle.c_str(), vars[iV].bins.size(),&vars[iV].bins[0]   )
                            ->Fill(v,bw*s*sel);
                    }
                }


            }
        }

        return true;
    }

    const std::vector<PlotSamp> samps;
    const std::vector<PlotSel>  sels;
    const std::vector<PlotVar>  vars;

    std::vector<std::unique_ptr<TTreeFormula> > selSels;
    std::vector<std::unique_ptr<TTreeFormula> > sampSels;
    std::vector<std::unique_ptr<TTreeFormula> > varSels;
    std::vector<std::unique_ptr<TTreeFormula> > varYSels;
    std::unique_ptr<TTreeFormula> stdSel;
    std::string prefix;
    HistGetter plotter;

    std::string outFileName;

};

#endif
