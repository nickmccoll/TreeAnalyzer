
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/makePlots.C"
#include "../predTools/CutConstants.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/Types.h"

#include "../predTools/makeJSON.C"

using namespace CutConstants;


void fitTTBarSF(std::string inputMC,std::string inputData, std::string outputName){

    double bins[]={0,500,600,700,800,900,1000,1100,1200,1400,1600,2000,2500,3000,4000,5000};
    int nBins = 15;
    HistGetter plotter;
    std::string func = "pol1";
    double minFit = 700;
    double maxFit = 2500;

    TFile * fM =  TObjectHelper::getFile(inputMC);
    TFile * fD =  TObjectHelper::getFile(inputData);
    auto hd = TObjectHelper::getObject<TH1F>(fD,"data_cr_hhMass");
    auto ho = TObjectHelper::getObject<TH1F>(fM,"other_cr_hhMass");
    auto ht = TObjectHelper::getObject<TH1F>(fM,"ttbar_cr_hhMass");

    TH1 * hd_c = plotter.getOrMake1D("data_binned",(std::string(";")+hhMCS.title).c_str(),nBins,bins  );
    TH1 * ho_c = plotter.getOrMake1D("other_binned",(std::string(";")+hhMCS.title).c_str(),nBins,bins  );
    TH1 * ht_c = plotter.getOrMake1D("ttbar_binned",(std::string(";")+hhMCS.title).c_str(),nBins,bins  );

    TH1 * h_tot = new TH1F("tot",(std::string(";")+hhMCS.title).c_str(),nBins,bins);
    TH1 * h_avg = new TH1F("avg",(std::string(";")+hhMCS.title).c_str(),nBins,bins);

    auto fillHistogram =[](TH1 * inH, TH1* outH, TH1 * avgH = 0, TH1 * totH=0){
        const int nOutBins = outH->GetNbinsX();
        for(int iB = 1; iB <=inH->GetNbinsX(); ++iB ){
            const float val = inH->GetBinContent(iB);
            const float err = inH->GetBinError(iB);
            const float xval = inH->GetBinCenter(iB);
            const int oB = outH->FindFixBin(xval);
            if(oB < 1 || oB > nOutBins) continue;
            const float outV = outH->GetBinContent(oB);
            outH->SetBinContent(oB,outH->GetBinContent(oB) + val);
            (*outH->GetSumw2())[oB]+=err*err;

            if(avgH) {
                avgH->SetBinContent(oB,avgH->GetBinContent(oB) + xval*val);
                totH->SetBinContent(oB,totH->GetBinContent(oB) + val);
            }
        }
    };

    fillHistogram(&*hd,hd_c,h_avg,h_tot);
    fillHistogram(&*ht,ht_c);
    fillHistogram(&*ho,ho_c);

    TGraphErrors * rat = new TGraphErrors(); rat->SetName("RATIO");
    TGraphErrors * ratM50 = new TGraphErrors();ratM50->SetName("RATIO_otherM50");
    TGraphErrors * ratP50 = new TGraphErrors();ratP50->SetName("RATIO_otherP50");
    int curPt = 0;
    int curPtM50 = 0;
    int curPtP50 = 0;



    for(int iB = 1; iB <= hd_c->GetNbinsX(); ++iB){
        if(ht_c->GetBinContent(iB) <= 0) continue;
        if(hd_c->GetBinContent(iB) <= 0) continue;

        double nD = hd_c->GetBinContent(iB);
        double nT = ht_c->GetBinContent(iB);
        double nO = ho_c->GetBinContent(iB);

        double nDE = hd_c->GetBinError(iB);
        double nTE = ht_c->GetBinError(iB);
        double nOE = ho_c->GetBinError(iB);

        auto addPoint =[&](TGraphErrors *g, double otherSF, int& cP ){
            double xV = h_avg->GetBinContent(iB)/h_tot->GetBinContent(iB);
            double yV = 0;
            double yE = 1;

            if(nO*otherSF<nD){
                yV = (nD-otherSF*nO)/nT;
                yE = std::sqrt(nDE*nDE+otherSF*otherSF*nOE*nOE+yV*yV*nTE*nTE)/nT;
            }

            g->SetPoint(curPt,xV,yV);
            g->SetPointError(curPt,0.0,yE);
            cP++;
        };

        addPoint(rat,1.0,curPt);
        addPoint(ratM50,0.5,curPtM50);
        addPoint(ratP50,1.5,curPtP50);
    }

    plotter.addGraph(&*rat);
    plotter.addGraph(&*ratM50);
    plotter.addGraph(&*ratP50);
    plotter.write(outputName+"_inputDebug.root");
    fD->Close();
    fM->Close();
    std::string args = "-i " + outputName+"_inputDebug.root "
            + " -g RATIO:" + func +" "
            + " -var "+ MOD_MS+ " -minX "+ ASTypes::flt2Str(minFit)+ " -maxX "+ ASTypes::flt2Str(maxFit);

    MakeJSON(outputName,args);
}

std::string getTTBarSFString(const std::string& filename){
    CJSON json(filename+"_ttbarSF.json");
    std::string qToW = json.getP(0).second;
    auto replace = [&](const std::string& vn, const std::string tf1n){
        std:size_t index = 0;
        while (true) {
            index = qToW.find(vn, index);
            if (index == std::string::npos) break;
            qToW.replace(index, vn.size(), tf1n);
            index += 1;
        }
    };
    replace(MOD_MS,hhMCS);
    return std::string("(")+processes[TTBAR].cut+"?("+qToW+"):1.0)";
}

void go(int step, std::string treeDir) {
    std::string filename = hhFilename;
    std::string mcTree = treeDir + "/betrees_mc.root";
    std::string dataTree = treeDir + "/betrees_data.root";

    std::vector<PlotVar> vars;
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);

    std::vector<PlotSamp> mcSamps =
    { {"mc","1.0"},
            {"mc_noQCD",std::string("!(")+processes[QCD].cut+")" },
            {"ttbar",processes[TTBAR].cut},
            {"other",std::string("!(")+processes[TTBAR].cut+")" },
            {"other_noQCD",std::string("!(")+processes[TTBAR].cut+"||"+processes[QCD].cut+")"}
    };
    std::vector<PlotSamp> dataSamps =
    { {"data","1.0"}
    };
    std::vector<PlotSel> sels;
    sels.emplace_back("cr",abV.cut+"&&"+hbbMCS.cut+">80"+"&&"+btagCats[BTAG_LMT].cut);
    sels.emplace_back("noHbb",abV.cut+"&&"+hbbMCS.cut+">30"+"&&"+btagCats[BTAG_LMT].cut);
    sels.emplace_back("fullHbbHH",abV.cut+"&&"+hhRange.cut+"&&"+hbbRange.cut+"&&"+btagCats[BTAG_LMT].cut);


    if(step == 0){
        MakePlots mcPlots(mcTree,filename+"_ttbarSF_mc_inputPlots.root",mcSamps,sels,vars,"1.0",nomW.cut);
        MakePlots dataPlots(dataTree,filename+"_ttbarSF_data_inputPlots.root",dataSamps,sels,vars,"1.0","1.0");

        fitTTBarSF(filename+"_ttbarSF_mc_inputPlots.root",filename+"_ttbarSF_data_inputPlots.root",filename+"_ttbarSF.json");
    }

    if(step==1){
        std::string ttbarSFStr = getTTBarSFString(filename);

        std::vector<PlotSamp> mcTestSamps =
        { {"mc",ttbarSFStr},{"mc_noSF","1.0"},
                {"ttbar",std::string("(")+processes[TTBAR].cut+")*"+ttbarSFStr},{"ttbar_noSF",processes[TTBAR].cut},{"other",std::string("!(")+processes[TTBAR].cut+")" }
        };

        MakePlots mcTestPlots(mcTree,filename+"_ttbarSF_mc_testPlots.root",mcTestSamps,sels,vars,"1.0",nomW.cut);
        MakePlots dataTestPlots(dataTree,filename+"_ttbarSF_data_testPlots.root",dataSamps,sels,vars,"1.0","1.0");
    }

}
#endif

void makeTTBarSF(int step = 0, std::string treeDir = "../trees/"){
    go(step, treeDir);
}
