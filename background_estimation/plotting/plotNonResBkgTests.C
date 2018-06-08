#include "../predTools/CutConstants.h"
#include <vector>
#include "TFile.h"
#include "HistoPlotting/include/Plotter.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;

std::vector<TObject*> writeables;

void make2DTests(std::string plotTitle, const TH2* dH, const std::vector<TH2*>& hs,const std::vector<std::string>& hNs, const std::vector<double>& bins, bool binInY, double rebin = -1) {
    bool withRatio = true;
    if(dH)for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(dH->Integral()/hs[iH]->Integral());
    const TH1 * refHist = dH ? dH : hs[0];
    const TAxis * ax =  binInY ? refHist->GetYaxis() : refHist->GetXaxis();

    for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
        int binL = ax->FindFixBin(bins[iB]);
        int binH = ax->FindFixBin(bins[iB+1]) -1;
        if(binH <= binL) continue;
        auto proj =[&](const TH2* h, const std::string& title) ->TH1*{
            return binInY ? h->ProjectionX( (title+"_"+int2Str(iB)).c_str(),binL,binH) :  h->ProjectionY( (title+"_"+int2Str(iB)).c_str(),binL,binH);
        };
        Plotter * p = new Plotter();
        if(dH){
            auto dh1 = proj(dH,"MC");
            p->addHist(dh1,"MC");
        }

        for(unsigned int iH = 0; iH < hs.size(); ++iH){
            TH1 * h = proj(hs[iH],hNs[iH]);
            if(withRatio) for(int iX = 1; iX <= h->GetNbinsX(); ++iX) h->SetBinError(iX,0);
            p->addHistLine(h,hNs[iH].c_str());
        }

        auto setupPlotter = [&](Plotter * p, std::string name){
            p->setMinMax(.0001,(rebin < 0 ? 1.0 : rebin) * refHist->Integral()/4);
            p->setUnderflow(false);
            p->setOverflow(false);
            p->setBotMinMax(0,2);
            p->setYTitle("N. of events");
            if(rebin > 0) p->rebin(rebin);
            if(withRatio){
                auto * c = p->drawSplitRatio(dH?1:0,"stack",false,false,name);
                c->GetPad(1)->SetLogy();
                c->GetPad(1)->Update();
                writeables.push_back(c);

            } else {
                auto * c = p->draw(false,name);
                c->SetLogy();
                c->Update();
                writeables.push_back(c);
            }
        };
        setupPlotter(p,plotTitle+"_"+  flt2Str(bins[iB]) +"_"+flt2Str(bins[iB+1]));
    }

}

void test2DCondTemplate(std::string name, std::string filename, std::string sel){
    TH2* dH = 0;
    std::vector<TH2*> hs;
    std::vector<std::string> hNs;
    auto addHistos =[&](std::string extraPre="", bool addSmooth=true, bool addKDE=false,  bool addData = false){
        TFile *f = new TFile((filename + "_"+name + (extraPre.size() ? std::string("_") + extraPre +"_" : std::string("_") )+sel +"_2D_cond_template.root").c_str(),"read");
//        HHlnujj_QGCR_qg_emu_AB_I_ltmb_2D_cond_template.root
        if(addData){
            f->GetObject("histo_fine_data" ,dH);
            dH->SetXTitle(hhMCS.title.c_str());
            dH->SetYTitle(hbbMCS.title.c_str());
        }
        if(addSmooth){
            TH2* h;
            f->GetObject("histo_smooth",h);
            hs.push_back(h); hNs.push_back(extraPre.size() ? extraPre  : std::string("KDE"));
        }
        if(addKDE){
            TH2* h;
            f->GetObject("histo_KDE",h);
            hs.push_back(h); hNs.push_back(extraPre.size() ? extraPre +"_NoSmear"  : std::string("KDE_NoSmear"));
        }
    };

    addHistos("",true,true,true);
    std::vector<std::string> extras = {"xs_0p75_xc_3_ys_0p5_yc_1","xs_0p75_xc_3_ys_0p5_yc_1_2"};
    //      for(const auto& s : extras ) addHistos(s,true,false);

    std::vector<double> hBBBinning = {30,40,50,60,80,100,120,140,170,210};
    std::vector<double> hhBinning  = {700,800,900,1000,1500,2000,3000,4000,5000};
    make2DTests(name + "_COND2D_HHF" ,dH,hs,hNs,hBBBinning,false,-1);
    make2DTests(name + "_COND2D_HHC" ,dH,hs,hNs,hBBBinning,false,10);
    make2DTests(name + "_COND2D_HbbF"  ,dH,hs,hNs,hhBinning,true,-1);
    make2DTests(name + "_COND2D_HbbC"  ,dH,hs,hNs,hhBinning,true,5 );
}

void test2DTemplate(const std::string& name, const std::string& filename, const std::vector<std::string>& sels){
//
    auto cutHistograms =[](const TH2* inH) -> TH2F*{ //The conditional templates have the Hbb as the y-axis..which is only done here
        TH2F * outH = new TH2F(TString(inH->GetName()) + "_cut",(std::string(";") +hbbMCS.title+";"+hhMCS.title).c_str()  ,nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
        for(int iX =1; iX <= inH->GetNbinsX(); ++iX){
            const int outIX =outH->GetXaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(iX));
            if(outIX < 1 || outIX > outH->GetNbinsX() ) continue;
            for(int iY =1; iY <= inH->GetNbinsY(); ++iY){
                const int outIY = outH->GetYaxis()->FindFixBin(inH->GetYaxis()->GetBinCenter(iY));
                if(outIY < 1 || outIY > outH->GetNbinsY() ) continue;
                outH->SetBinContent(outIX,outIY,inH->GetBinContent(iX,iY));
                outH->SetBinError(outIX,outIY,inH->GetBinError(iX,iY));
            }
        }
        return outH;
    };

    TFile *fData = new TFile((filename + "_"+name +"_distributions.root").c_str(),"read");
    for(const auto& s :sels ){
        TH2* dH = 0;
        TH2* h = 0;
        std::vector<TH2*> hs;
        std::vector<std::string> hNs;
        fData->GetObject((name+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dH); dH = cutHistograms(dH);
        TFile * fT = new TFile((filename + "_"+name+"_"+s  +"_2D_merged_template.root").c_str(),"read");
        fT->GetObject("histo",h);hs.push_back(h);hNs.push_back("Template");

//        std::vector<double> hBBBinning = {30,40,50,60,80,100,120,140,170,210};
//        std::vector<double> hhBinning  = {700,800,900,1000,1500,2000,3000,4000,5000};
        std::vector<double> hBBBinning = {30,210};
        std::vector<double> hhBinning  = {700,4000};
    //    make2DTests(name +"_"+s+ "_Temp_HbbF",dH,hs,hNs,hBBBinning,false);
        make2DTests(name +"_"+s+ "_Temp_HHC",dH,hs,hNs,hBBBinning,false,10);
    //    make2DTests(name +"_"+s+ "_Temp_HHF" ,dH,hs,hNs,hhBinning,true);
        make2DTests(name +"_"+s+ "_Temp_HbbC" ,dH,hs,hNs,hhBinning,true,5);
    }

}

void test2DFits(std::string name, std::string filename,const std::vector<std::string>& sels,const std::vector<double>& bins, bool binInY = true, double rebin = -1){


    TFile *fData = new TFile((filename + "_"+name  +"_distributions.root").c_str(),"read");
//    TFile *fTemp = new TFile((filename + "_"+name  +"_2D_template.root").c_str(),"read");
//    TFile *fTempFit = new TFile((filename + "_"+name  +"_2D_template_debug.root").c_str(),"read");


    for(const auto& s : sels){
        TH2* dH = 0;
        std::vector<TH2*> hs;
        std::vector<std::string> hNs;
        fData->GetObject((name+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dH);
        if(dH==0) continue;

        TFile *fTemp = new TFile((filename + "_"+name+"_" +s +"_2D_template.root").c_str(),"read");
        TH2 * hF = 0;
        fTemp->GetObject("histo",hF);
        if(hF == 0) continue;

        TH2 * hOT = 0;
        fTemp->GetObject("originalPDF",hOT);
        if(hOT == 0) continue;

        if(binInY){
            make2DTests(name + "_Fit_Hbb_"+s,dH,{hF,hOT},{"Search region template","Baseline template"},bins,binInY,rebin);
        } else {
            make2DTests(name + "_Fit_HH_"+s ,dH,{hF,hOT},{"Search region template","Baseline template"},bins,binInY,rebin);
        }
//        make2DTests(name + "_Fit_Hbb_"+s,dH,{hF,hOT},{"Search region template","Baseline template"},hBBBinning,false,10);
//        make2DTests(name + " Fit_HH_ "+s ,dH,{hF,hOT},{"Search region template","Baseline template"},hhBinning,true,5);

        //        make2DTests(name + " Fit HbbF "+s,dH,{hF,hOT},{"Search region template","Original template"},hBBBinning,false);
        //        make2DTests(name + " Fit HHF "+s ,dH,{hF,hOT},{"Search region template","Original template"},hhBinning,true);
        //        make2DTests(name + " Fit HHC "+s ,dH,{hF,hOT},{"Search region template","Baseline template"},hhBinning,true,5);
    }



}

void testQCDSF(std::string name, std::string filename, const std::vector<std::string>& sels){
    TFile *fData = new TFile((filename + "_"+name  +"_distributions.root").c_str(),"read");

    std::vector<double> hBBBinning = {30,210};
    std::vector<double> hhBinning  = {700,3500};
    for(const auto& s : sels){
        TH2* hMC = 0;
        fData->GetObject((name+"_wQCD_noQCDSF_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),hMC);
        TH2* hNoSF = 0;
        fData->GetObject((name+"_noQCDSF_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),hNoSF);
        TH2* hSF = 0;
        fData->GetObject((name+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),hSF);
        if(hMC == 0) continue;

            make2DTests(name + "_testQCDSF_HHC_"+s,0,{hMC,hNoSF,hSF},{"QG+QCD","QG, no SF","QG, w/ SF"},hBBBinning,false,10);
//                make2DTests(name + "testQCDSF_HbbC_"+s ,0,{hMC,hNoSF,hSF},{"QG+QCD","QG, no SF","QG, w/ SF"},hhBinning,true,5);
    }
}




std::vector<TObject*> testRatioFits(std::string name, std::string filename,  std::string fitName, const std::vector<std::string>& sels) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<TObject*> writeables;
    std::vector<TObject*> paramPads;

    for(const auto& s : sels){
        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");
        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
            TGraphErrors * can= 0;
            if(ff){
                ff->GetObject((name).c_str(),can);
            }
            if(!can) return;
            can->GetYaxis()->SetTitle(s.c_str());
            list.push_back(can);
        };

        addGraph("RATIO", paramPads);
    }



    auto * c1 =Drawing::drawAll(paramPads,"ratio_params");
    writeables.push_back(c1);
    return writeables;
}


std::vector<TObject*> testRatioUncs(std::string name, std::string filename,  std::string fitName, const std::vector<std::pair<std::string,std::string> >& targets, //target,model
        float normUnc, float scaleUnc, float resUnc) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<TObject*> writeables;
    std::vector<TObject*> paramPads;

    for(const auto& s : targets){
        TFile *ff = new TFile((filename+"_"+name+"_"+s.second+"_"+fitName+".json.root").c_str(),"read");
        TFile *ft = new TFile((filename+"_"+name+"_"+s.first+"_"+fitName+".json_inputDebug.root").c_str(),"read");

        if(ff == 0) continue;
        if(ft == 0) continue;

        TGraphErrors * graph= 0;
        ft->GetObject("RATIOPlusOne",graph);
        if(!graph) continue;
        graph->GetYaxis()->SetTitle(s.first.c_str());
        TCanvas * can = new TCanvas((std::string("syst_")+s.first).c_str());
        can->cd();
        graph->Draw();

        TF1 * defF= 0;
        ff->GetObject("RATIO_func",defF);
        if(!defF) continue;

        auto mkF = [&](std::string name,std::string factor, int color, int style){
            std::string fs = "(1+[0]+[1]*x+[2]/x)";
            if(factor.size()) fs += "*("+factor+")";
            TF1 * f = new TF1((s.first+"_"+name+"_mix").c_str(),fs.c_str(),1,13000);
            f->SetParameter(0,defF->GetParameter(0));
            f->SetParameter(1,defF->GetParameter(1));
            f->SetParameter(2,defF->GetParameter(2));
            f->SetLineColor(color);
            f->SetLineStyle(style);
            can->cd();
            f->Draw("SAME");
        };

        mkF("nominal","",StyleInfo::getLineColor(0),1 );
        mkF("PTUp","(1. + 0.0003*x)",StyleInfo::getLineColor(1),1 );
        mkF("PTDown","1/(1. + 0.0003*x)",StyleInfo::getLineColor(1),9 );
        mkF("OPTUp","(1. + 1200/x)",StyleInfo::getLineColor(2),1 );
        mkF("OPTDown","1/(1. + 1200/x)",StyleInfo::getLineColor(2),9 );
        mkF("NORMUp","1.5",StyleInfo::getLineColor(3),1 );
        mkF("NormDown","0.5",StyleInfo::getLineColor(3),9 );
        writeables.push_back(can);
    }
    //    auto * c1 =Drawing::drawAll(paramPads,"ratio_params");
    //    writeables.push_back(c1);
    return writeables;
}

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




void plotNonResBkgTests(int step = 0,bool doTW = true, int inreg = REG_SR, std::string outName = ""){
    REGION reg = REGION(inreg);

    std:: string inName =  "bkgInputs" ;
    auto srList = getSRList(reg);
    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
    }
    else if(reg == REG_QGCR){
        inName =  "bkgInputsQGCR";
        hhFilename +="_QGCR";
    }
    std::string filename = inName +"/"+hhFilename;

    CutStr mod = bkgSels [doTW ? BKG_LOSTTW : BKG_QG];
    if(outName.size()){
        outName += std::string("/") + mod;
        if(reg == REG_TOPCR) outName +=  "TopCR";
        else if(reg == REG_QGCR) outName +=  "QGCR";
    }


    std::vector<std::string> stepSels;
    switch(step){
    case 0:
        if(doTW) return;
        if(outName.size()) outName += "_QCDRatio.root";
        writeables = testRatioFits(mod,filename,"QCDSF",{"emu_I_I_ltmb","e_I_LP_ltmb","mu_I_LP_ltmb","e_I_LP_full","e_I_HP_full","mu_I_LP_full","mu_I_HP_full"});
        if(reg == REG_QGCR){
            writeables = testRatioUncs(mod,filename,"QCDSF",{ {"emu_LMT_I_ltmb","emu_I_I_ltmb"},
                    {"e_L_LP_full","e_I_LP_full"},{"e_L_HP_full","e_I_HP_full"},{"mu_L_LP_full","mu_I_LP_full"},
                    {"mu_L_HP_full","mu_I_HP_full"}},1,1,1);
        } else {
            writeables = testRatioUncs(mod,filename,"QCDSF",{ {"emu_LMT_I_ltmb","emu_I_I_ltmb"},
                    {"e_L_LP_full","e_I_LP_full"},{"e_M_LP_full","e_I_LP_full"},{"e_T_LP_full","e_I_LP_full"},
                    {"e_L_HP_full","e_I_HP_full"},{"e_M_HP_full","e_I_HP_full"},{"e_T_HP_full","e_I_HP_full"},
                    {"mu_L_LP_full","mu_I_LP_full"},{"mu_M_LP_full","mu_I_LP_full"},{"mu_T_LP_full","mu_I_LP_full"},
                    {"mu_L_HP_full","mu_I_HP_full"},{"mu_M_HP_full","mu_I_HP_full"},{"mu_T_HP_full","mu_I_HP_full"}},1,1,1);
        }
        break;
    case 1:
        if(doTW) return;
        if(outName.size()) outName += "_testQCDSF.root";
        stepSels = srList; stepSels.push_back("emu_LMT_I_ltmb");
        testQCDSF(mod,filename,stepSels);
        break;
    case 2:
        if(outName.size()) outName += "_2DCondTemp.root";
        test2DCondTemplate(mod,filename,"emu_LMT_I_ltmb");
        break;
    case 3:
        if(outName.size()) outName += "_MVVKern.root";
        if(doTW) stepSels = {"emu_LMT_I_ltmb","emu_LMT_LP_ltmb","emu_LMT_HP_ltmb"};
        else stepSels = {"emu_LMT_I_ltmb","e_LMT_I_ltmb","mu_LMT_I_ltmb"};
        writeables = test1DKern(mod,filename,"MVV",stepSels);
        break;
    case 4:
        if(outName.size()) outName += "_2DTemp.root";
        if(doTW) stepSels = {"emu_LMT_I_ltmb","emu_LMT_LP_ltmb","emu_LMT_HP_ltmb"};
        else stepSels = {"emu_LMT_I_ltmb","e_LMT_I_ltmb","mu_LMT_I_ltmb"};
        test2DTemplate(mod,filename,stepSels);
        break;
    case 5:
        if(outName.size()) outName += "_2DFits.root";
        test2DFits(mod,filename,srList,
                {30.,210, 100,150},false,10);
//        test2DFits(mod,filename,srList,
//                {800,900,1000,1500,2000,3000,4000,5000},true,5);

//        test2DFits(mod,filename,srList,
//                {30.,210},false,10);
        test2DFits(mod,filename,srList,
                {700,4000},true,5);

        break;
    case 6:
        if(outName.size()) outName += "_2DComp.root";
        writeables = test2DModel({mod},filename,
                srList,{700,4000});
        break;
    }

    Dummy d(outName);
}
