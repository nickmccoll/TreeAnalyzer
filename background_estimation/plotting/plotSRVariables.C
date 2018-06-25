
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/makePlots.C"
#include "../predTools/CutConstants.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"

using namespace CutConstants;

CutStr blindCut = CutStr("blindCut",std::string("(")+hbbMCS.cut+"<100||"+hbbMCS.cut+">150)");
std::vector<PlotVar> vars;
std::vector<PlotSamp> samps;



void makeDataDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut,bool isData){
    std::vector<PlotSel> sels;
    sels.emplace_back("full"                ,btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_hhMass"        ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_hbbMass"       ,btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_nAK4Btags"     ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+nSJs.cut+"&&"+exA.cut+"&&"+wjjBC.cut);
    sels.emplace_back("loose_wjjTau2o1"     ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+nSJs.cut+"&&"+exA.cut+"&&"+bV.cut);
    sels.emplace_back("loose_hwwPT_o_hhMass",blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+nSJs.cut+"&&"+"wwDM<2"+"&&"+wjjBC.cut+"&&"+bV.cut);
    sels.emplace_back("loose_wwDM"          ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+nSJs.cut+"&&"+"(hwwPT/hhMass>0.3)"+"&&"+wjjBC.cut+"&&"+bV.cut);
    sels.emplace_back("loose_hbbCSVCat"     ,blindCut.cut +"&&"+ hadCuts[HAD_FULL].cut);



    std::string outFileName=filename+"_"+name+  "_srVarDistributions.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,isData ? "1.0" : nomW.cut);
}

void compilePlots(const std::string& prefix, const std::string& mcFile, const std::string& dataFile,  const std::vector<std::string> signalFiles,  const std::vector<std::string> signalNames){
    TFile * fd = new TFile((prefix+dataFile).c_str(),"READ");
    TFile * fm = new TFile((prefix+mcFile  ).c_str(),"READ");
    std::vector<TFile*> sFs;
    for(const auto& f :signalFiles)
        sFs.push_back(new TFile((prefix+f).c_str(),"READ"));


    for(unsigned int iV = 1; iV <= vars.size(); ++iV){
        Plotter * p = new Plotter;

        for(unsigned int iS = 0; iS < signalFiles.size(); ++iS){
            TFile * f = new TFile((prefix+signalFiles[iS]).c_str(),"READ");
            TH1 * hm = 0;
            f->GetObject((std::string("all_loose_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            p->addHistLine(hm,signalNames[iS]);
        }

        TH1 * hd = 0;
        fd->GetObject((std::string("all_loose_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hd);
        if(hd == 0) continue;
        p->addHist(hd,"data",kBlack);


//
        for(unsigned int iP = 0; iP < processes.size(); ++iP){
            TH1 * hm = 0;
            fm->GetObject((processes[iP]+"_loose_"+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            p->addStackHist(hm,processes[iP].title);
        }
        if(iV == 2) {
            p->rebin(4);
            p->setMinMax(0.1,10000);
        }
        if(iV== 7) p->setMinMax(20,1000000);

        if(iV == 4) p->setLegendPos(0.1588629,0.5895652,0.4397993,0.9043478);
        else p->setLegendPos(0.6555184,0.5530435,0.9364548,0.8678261);
        p->setCMSLumi();
        p->setYTitle("N. of events / bin width");
        auto * c = p->draw(false,prefix+vars[iV].varName+"_srvardists.pdf");
        p->yAxis()->SetTitleOffset(1.55);
        p->xAxis()->SetTitleOffset(1.0);
        if(iV == 7){
            p->xAxis()->SetBinLabel(1,"No loose");
            p->xAxis()->SetBinLabel(2,"One loose");
            p->xAxis()->SetBinLabel(3,"Two loose");
            p->xAxis()->SetBinLabel(4,"One medium, no loose");
            p->xAxis()->SetBinLabel(5,"One medium, one loose");
            p->xAxis()->SetBinLabel(6,"Two medium");
            p->xAxis()->SetTitle(" ");
            c->SetLogy();


        }
        if(iV == 2){
            c->SetLogy();
        }
        c->Update();
        c->Print((prefix+vars[iV].varName+"_srvardists.pdf").c_str());


    }
}

void categoryPlots(const std::string& prefix, const std::string& mcFile){
    TFile * fm = new TFile((prefix+mcFile  ).c_str(),"READ");

    for(unsigned int iV = 1; iV <= vars.size(); ++iV){
        Plotter * p = new Plotter;


//
        for(unsigned int iP = 0; iP < bkgSels.size(); ++iP){
            TH1 * hm = 0;
            fm->GetObject((bkgSels[iP]+"_full_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            p->addStackHist(hm,bkgSels[iP].title);
        }
        if(iV == 2) {
            p->rebin(4);
            p->setMinMax(0.1,10000);
        }
        if(iV== 7) p->setMinMax(20,1000000);

        if(iV == 4) p->setLegendPos(0.1588629,0.5895652,0.4397993,0.9043478);
        else p->setLegendPos(0.6555184,0.5530435,0.9364548,0.8678261);
        p->setCMSLumi();
        p->setYTitle("N. of events / bin width");
        auto * c = p->draw(false,prefix+vars[iV].varName+"_srvardists.pdf");
        p->yAxis()->SetTitleOffset(1.55);
        p->xAxis()->SetTitleOffset(1.0);
        if(iV == 7){
            p->xAxis()->SetBinLabel(1,"No loose");
            p->xAxis()->SetBinLabel(2,"One loose");
            p->xAxis()->SetBinLabel(3,"Two loose");
            p->xAxis()->SetBinLabel(4,"One medium, no loose");
            p->xAxis()->SetBinLabel(5,"One medium, one loose");
            p->xAxis()->SetBinLabel(6,"Two medium");
            p->xAxis()->SetTitle(" ");
            c->SetLogy();


        }
        if(iV == 2){
            c->SetLogy();
        }
        c->Update();
        c->Print((prefix+vars[iV].varName+"_catdists.pdf").c_str());


    }
}

#endif

void plotSRVariables(int step, std::string tree, std::string name){
    bool isData = ASTypes::strFind(tree,"data");
    std::string cut =hhRange.cut+"&&"+hbbRange.cut ;
    std::string out = "srVarDists/"+hhFilename;
    if(ASTypes::strFind(tree,"radion")) blindCut.cut = "(1.0)";
    if(isData) cut += "&&"+blindCut.cut;

    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,36,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back("nAK4Btags" ,"; N. of AK4 jet b tags","nAK4Btags",4,-0.5,3.5);
    vars.emplace_back("wjjTau2o1" ,"; #tau_{2}/#tau_{1}","wjjTau2o1",50,0,1);
    vars.emplace_back("hwwPT_o_hhMass" ,"; #it{p}_{T}/#it{m}","hwwPT/hhMass",50,0,1);
    vars.emplace_back("wwDM","; #it{m}_{D} [GeV]","125*wwDM/2",50,0,500);
    vars.emplace_back("hbbCSVCat","; hbbCSVCat","(hbbCSVCat-1)",6,-0.5,5.5);


    samps.emplace_back("all","1.0");
    samps.emplace_back(processes[TTBAR],processes[TTBAR].cut);
    samps.emplace_back(processes[WJETS],processes[WJETS].cut);
    samps.emplace_back(processes[QCD],processes[QCD].cut);
    samps.emplace_back(processes[OTHER],processes[OTHER].cut);
    samps.emplace_back(bkgSels[BKG_QG],bkgSels[BKG_QG].cut);
    samps.emplace_back(bkgSels[BKG_LOSTTW],bkgSels[BKG_LOSTTW].cut);
    samps.emplace_back(bkgSels[BKG_MW],bkgSels[BKG_MW].cut);
    samps.emplace_back(bkgSels[BKG_MT],bkgSels[BKG_MT].cut);

    if(step==0){
        makeDataDistributions(name,out,tree,cut ,isData);
    }
    if(step==1){
        compilePlots("srVarDists/",hhFilename+"_mc_srVarDistributions.root",hhFilename+"_data_srVarDistributions.root",
                {hhFilename+"_m1000_srVarDistributions.root",hhFilename+"_m2500_srVarDistributions.root"}, {"#it{m}_{X}=1 TeV","#it{m}_{X}=2.5 TeV"});
    }
    if(step==2) categoryPlots("srVarDists/",hhFilename+"_mc_srVarDistributions.root");
}
