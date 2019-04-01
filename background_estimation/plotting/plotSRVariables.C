
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/makePlots.C"
#include "../predTools/CutConstants.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"

#include "TStyle.h"
#include "TGraphAsymmErrors.h"

using namespace CutConstants;
//CutStr blindCut = CutStr("blindCut",std::string("(")+hbbMCS.cut+"<100||"+hbbMCS.cut+">150)");
CutStr blindCut = CutStr("blindCut",std::string("(1.0)"));
std::vector<PlotVar> vars;
std::vector<std::string> varUnits;
std::vector<PlotSamp> samps;
float SIG_CROSS = 2; //in pb
bool preliminary = false;



void makeDataDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut,bool isData){
    std::vector<PlotSel> sels;
    sels.emplace_back("full"                ,btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_hhMass"        ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_hbbMass"       ,btagCats[BTAG_LMT].cut+"&&"+hadCuts[HAD_FULL].cut);
    sels.emplace_back("loose_nAK4Btags"     ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel.cut+"&&"+exA.cut+"&&"+wjjBC.cut);
    sels.emplace_back("loose_wjjTau2o1"     ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel.cut+"&&"+exA.cut+"&&"+bV.cut);
    sels.emplace_back("loose_hwwPT_o_hhMass",blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel.cut+"&&"+"wwDM<125"+"&&"+wjjBC.cut+"&&"+bV.cut);
    sels.emplace_back("loose_wwDM"          ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel.cut+"&&"+"(hwwPT/hhMass>0.3)"+"&&"+wjjBC.cut+"&&"+bV.cut);
    sels.emplace_back("loose_hbbCSVCat"     ,blindCut.cut +"&&"+ hadCuts[HAD_FULL].cut);



    std::string outFileName=filename+"_"+name+  "_srVarDistributions.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,isData ? "1.0" : nomW.cut);
}

void compilePlots(const std::string& prefix, const std::string& mcFile, const std::string& dataFile,  const std::vector<std::string> signalFiles,  const std::vector<std::string> signalNames){
    TFile * fd = new TFile((dataFile).c_str(),"READ");
    TFile * fm = new TFile((mcFile  ).c_str(),"READ");

    for(unsigned int iV = 1; iV < vars.size(); ++iV){
        Plotter * p = new Plotter;
        std::vector<Drawing::TLegendEntryDef> legEntries;


        int rebinFactor = 0;
        switch(iV){
        case 1:
            rebinFactor = 3;
            break;
        case 2:
            rebinFactor = 4;
            break;
        }


        for(unsigned int iP = 0; iP < processes.size(); ++iP){
            TH1 * hm = 0;
            fm->GetObject((processes[iP]+"_loose_"+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            if(rebinFactor) hm->Rebin(rebinFactor);
            auto * g = p->addStackHist(hm,processes[iP].title);
            legEntries.push_back(std::make_tuple(10+processes.size() +iP,g,processes[iP].title.c_str(),"f"));

        }

        auto * tot = (TH1*)p->getTotStack()->Clone();
        auto* errBand = new TGraphAsymmErrors();

        for(int iB = 1; iB <= tot->GetNbinsX(); ++iB){
            double x = tot->GetBinCenter(iB);
            double y = tot->GetBinContent(iB);
            double err = tot->GetBinError(iB);
            if(iB == 1){
                y += tot->GetBinContent(0);
                err = std::sqrt(err*err + tot->GetBinError(0)*tot->GetBinError(0));
            }
            if(iB == tot->GetNbinsX()){
                y += tot->GetBinContent(iB+1);
                err = std::sqrt(err*err + tot->GetBinError(iB+1)*tot->GetBinError(iB+1));
            }
            auto setPt = [&](const int bin, float x){
                errBand->SetPoint(bin,x,y);
                errBand->SetPointError(bin,tot->GetBinWidth(iB)/2.,tot->GetBinWidth(iB)/2,err,err);
            };
            if(iB == 1) setPt(iB-1,x-tot->GetBinWidth(iB)/2.);
            setPt(iB,x);
            if(iV==1){
                const auto& hs =p->getStackHists();
                std::cout << iB <<" "<< x <<" "<< y << " "<< err;
                        for(const auto& h : hs){
                            std::cout << "(" <<((const TH1*)h.obj)->GetBinContent(iB)<<","<<((const TH1*)h.obj)->GetBinError(iB)<<") ";
                        }


                        std::cout <<std::endl;

            }
            if(iB == tot->GetNbinsX() ) setPt(iB+1,x+tot->GetBinWidth(iB)/2.);
        }

        p->clearTotStackError();

        int fillColor = kMagenta+4;
        errBand->SetFillColor(fillColor);
        errBand->SetFillStyle(3352);
        gStyle->SetHatchesLineWidth(1);
        gStyle->SetHatchesSpacing(.5);
        auto * g_mcunc = p->addGraph(errBand,"Sim. stat. unc.",fillColor,1,0,20,1,false,true,false,"2");
        legEntries.push_back(std::make_tuple(2,g_mcunc,"Sim. stat. unc.","F"));
        legEntries.push_back(std::make_tuple(4,(TObject*)(0),"",""));

        std::vector<int> signalColors ={kSpring+5,634};

        for(unsigned int iS = 0; iS < signalFiles.size(); ++iS){
            TFile * f = new TFile((signalFiles[iS]).c_str(),"READ");
            TH1 * hm = 0;
            f->GetObject((std::string("all_loose_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            //Scale so that we get 1pb normalization
            hm->Scale(CutConstants::HHtobbVVBF* SIG_CROSS);
            if(rebinFactor) hm->Rebin(rebinFactor);
            hm->SetName((signalNames[iS]+vars[iV].varName).c_str());
            auto * g = p->addHistLine(hm,signalNames[iS],signalColors[iS]);
            legEntries.push_back(std::make_tuple(100+iS,g,signalNames[iS],"L"));

        }

        TH1 * hd = 0;
        fd->GetObject((std::string("all_loose_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hd);
        if(hd == 0) continue;
        if(rebinFactor) hd->Rebin(rebinFactor);
         p->addHist(hd,"Data",kBlack,1,2,20,0.5,true,true,true);
        TGraphAsymmErrors * g_d = new TGraphAsymmErrors;
        g_d->SetLineColor  (kBlack);
        g_d->SetLineWidth  (2);
        g_d->SetLineStyle  (1);
        g_d->SetMarkerStyle(20);
        g_d->SetMarkerColor(kBlack);
        g_d->SetMarkerSize (0.5);
        legEntries.push_back(std::make_tuple(0,g_d,"Data","P E"));
        legEntries.push_back(std::make_tuple(1,(TObject*)(0),"",""));

        p->setCMSLumi(10);

        bool sup = false;
        switch(iV){
        case 1: //hbb
            p->setMinMax(0,575);
            break;
        case 2: //hh
            p->setMinMax(0.1,20000);
            sup = true;
            break;
        case 3://nAK4Btags
            p->setMinMax(0,10000);
            sup = true;
            break;
        case 4://wjjTau2o1
            p->setMinMax(0,1000);
            break;
        case 5://hwwPT_o_hhMass
            p->setMinMax(0,1100);
            break;
        case 6://wwDM
            p->setMinMax(0,1250);
            break;
        case 7://hbbCSVCat
            p->setMinMax(20,1000000);
            sup = true;
            break;
        }

        double xV =0.45;
        double yV =0.619;

//        p->addText("All categories",xV+0.0075,yV+0.2075,0.045);
        if(sup||preliminary){
            p->setCMSLumiPosition(0,1.05);
           if(sup) p->setCMSLumiExtraText("Supplementary");
           else p->setCMSLumiExtraText("Preliminary");
            p->addText("All categories",.18,0.83,0.045);
        } else {
            p->addText("All categories",.18,0.75,0.045);
        }

        //--------------------LEGEND AND TEXT------------------------------
        p->turnOffLegend();
        TLegend * legend = signalNames.size() ? new TLegend(xV,yV,xV+0.453,yV+0.25) : new TLegend(xV,yV,xV+0.445,yV+0.25);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetNColumns(2);
        std::sort(legEntries.begin(), legEntries.end(), [](const Drawing::TLegendEntryDef& a, const Drawing::TLegendEntryDef& b) {return  std::get<0>(a) < std::get<0>(b); }  );
        for(const auto& l : legEntries){
            legend->AddEntry(std::get<1>(l),std::get<2>(l),std::get<3>(l));
        }
        //--------------------LEGEND AND TEXT------------------------------


        if(signalNames.size()){
            p->addText(TString::Format("#sigma#bf{#it{#Beta}}(X #rightarrow HH) = %0.f pb",SIG_CROSS),xV+0.0075,yV-0.04,0.042);
        } else {
        }
        p->setLegendNColumns(2);

        std::string ytitle = "Events";
        if(varUnits[iV]!="-1"){
            ytitle += " / " + ASTypes::flt2Str(p->getTotStack()->GetBinWidth(1)) + " " + varUnits[iV];

        }
        p->setBotMinMax(0.05,1.95);
        p->setYTitle(ytitle);
        p->setYTitleBot("Data / sim.");
//        auto * c = p->draw(false,prefix+vars[iV].varName+"_srvardists.pdf");
        auto * c = p->drawSplitRatio(-1,"stack",false,false,prefix+vars[iV].varName+"_srvardists.pdf");
//        p->yAxis()->SetTitleOffset(1.55);
        p->xAxis()->SetTitleOffset(1.05);
        c->GetPad(1)->cd();
        legend->Draw();
        c->GetPad(1)->Update();

        p->botStyle.xAxis->SetTitleOffset(1.05);


        for(unsigned int iSig = 0; iSig < signalFiles.size(); ++iSig){
            auto prim = c->GetPad(2)->GetPrimitive((signalNames[iSig]+vars[iV].varName).c_str());
            if(prim) prim->Delete();
        }

        if(iV == 1){//temp hack
            p->botStyle.xAxis->SetTitle(hbbMCS.title.c_str());
        }
        if(iV == 4){//temp hack
            p->botStyle.xAxis->SetTitle("q#bar{q}' #tau_{2}/#tau_{1}");
        }
        if(iV == 7){
            p->botStyle.xAxis->SetBinLabel(1,"0 L");
            p->botStyle.xAxis->SetBinLabel(2,"1 L");
            p->botStyle.xAxis->SetBinLabel(3,"2 L");
            p->botStyle.xAxis->SetBinLabel(4,"1 M");
            p->botStyle.xAxis->SetBinLabel(5,"1 M, 1 L");
            p->botStyle.xAxis->SetBinLabel(6,"2 M");
            p->botStyle.xAxis->SetTitle("b#bar{b} jet subjet b tagging");
            c->GetPad(1)->SetLogy();


        }
        if(iV == 2){
            c->GetPad(1)->SetLogy();
        }
        c->GetPad(1)->Update();
        c->Print((prefix+vars[iV].varName+"_srvardists.pdf").c_str());


    }
}

void categoryPlots(const std::string& prefix, const std::string& mcFile){
    TFile * fm = new TFile((mcFile  ).c_str(),"READ");

    for(unsigned int iV = 0; iV < vars.size(); ++iV){
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

void plotSRVariables( int step, int reg,std::string tree, std::string name){
    std::string prefix = "srVarDists/";
    if(reg == REG_TOPCR){
        prefix+="TopCR_";
        hhFilename += "_TopCR";
        blindCut.cut = "(1.0)";
        hadCuts[HAD_NONE].cut = preSel.cut;
        hadCuts[HAD_LB].cut   = preSel.cut+"&&"+wjjBC.cut;
        hadCuts[HAD_LT].cut   = preSel.cut+ "&&"+abV.cut;
        hadCuts[HAD_LTMB].cut = preSel.cut ;
        hadCuts[HAD_FULL].cut = preSel.cut + "&&"+abV.cut+"&&"+wjjBC.cut;
    } else if(reg==REG_QGCR){
        prefix+="QGCR_";
        btagCats = qgBtagCats;
        hhFilename +="_QGCR";
        blindCut.cut = "(1.0)";
    }


    bool isData = ASTypes::strFind(tree,"data");
    std::string cut =hhRange.cut+"&&"+hbbRange.cut ;
    std::string out = "srVarDists/"+hhFilename;
    if(ASTypes::strFind(tree,"radion")) blindCut.cut = "(1.0)";
    if(isData) cut += "&&"+blindCut.cut;
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,
            minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins
            ,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    vars.emplace_back("nAK4Btags" ,"; N. of AK4 jet b tags","nAK4Btags",4,-0.5,3.5);
    vars.emplace_back("wjjTau2o1" ,"; q#bar{q}' #tau_{2}/#tau_{1}","wjjTau2o1",50,0,1);
    vars.emplace_back("hwwPT_o_hhMass" ,"; #it{p}_{T}/#it{m}","hwwPT/hhMass",50,0,1);
    vars.emplace_back("wwDM","; #it{m}_{D} [GeV]","wwDM",50,0,500);
    vars.emplace_back("hbbCSVCat","; hbbCSVCat","(hbbCSVCat-1)",6,-0.5,5.5);

    varUnits.emplace_back("GeV");
    varUnits.emplace_back("GeV");
    varUnits.emplace_back("GeV");
    varUnits.emplace_back("-1");
    varUnits.emplace_back("units");
    varUnits.emplace_back("units");
    varUnits.emplace_back("GeV");
    varUnits.emplace_back("-1");
    std::cout <<vars[1].varTitle << std::endl;

    samps.emplace_back("all","1.0");
    if(!isData){
        samps.emplace_back(processes[TTBAR],std::string("(0.72*(")+processes[TTBAR].cut+"))");
        samps.emplace_back(processes[WJETS],processes[WJETS].cut);
        samps.emplace_back(processes[QCD],processes[QCD].cut);
        samps.emplace_back(processes[OTHER],processes[OTHER].cut);
        samps.emplace_back(bkgSels[BKG_QG],bkgSels[BKG_QG].cut);
        samps.emplace_back(bkgSels[BKG_LOSTTW],std::string("(0.76*(")+bkgSels[BKG_LOSTTW].cut+"))");
        samps.emplace_back(bkgSels[BKG_MW],std::string("(0.76*(")+bkgSels[BKG_MW].cut+"))");
        samps.emplace_back(bkgSels[BKG_MT],std::string("(0.76*(")+bkgSels[BKG_MT].cut+"))");
    }


    if(step==0){
        makeDataDistributions(name,out,tree,cut ,isData);
    }
    if(step==1){

        if(reg==REG_SR){
            compilePlots(prefix,out+"_mc_srVarDistributions.root",out+"_data_srVarDistributions.root",
                    {out+"_m1000_srVarDistributions.root",out+"_m2500_srVarDistributions.root"},
                    {"1 TeV X_{spin-0}","2.5 TeV X_{spin-0}"});
        } else{
            compilePlots(prefix,out+"_mc_srVarDistributions.root",out+"_data_srVarDistributions.root"
                    ,{},{});

        }

    }
    if(step==2) categoryPlots(prefix,out+"_mc_srVarDistributions.root");
}
