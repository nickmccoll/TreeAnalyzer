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
    for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(dH->Integral()/hs[iH]->Integral());
    const TAxis * ax =  binInY ? dH->GetYaxis() : dH->GetXaxis();

    for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
        int binL = ax->FindFixBin(bins[iB]);
        int binH = ax->FindFixBin(bins[iB+1]) -1;
        auto proj =[&](const TH2* h, const std::string& title) ->TH1*{
            return binInY ? h->ProjectionX( (title+"_"+int2Str(iB)).c_str(),binL,binH) :  h->ProjectionY( (title+"_"+int2Str(iB)).c_str(),binL,binH);
        };
        Plotter * p = new Plotter();
        auto dh1 = proj(dH,"MC");
        p->addHist(dh1,"MC");
        for(unsigned int iH = 0; iH < hs.size(); ++iH){
            TH1 * h = proj(hs[iH],hNs[iH]);
            if(withRatio) for(int iX = 1; iX <= h->GetNbinsX(); ++iX) h->SetBinError(iX,0);
            p->addHistLine(h,hNs[iH].c_str());
        }

        auto setupPlotter = [&](Plotter * p, std::string name){
            p->setMinMax(.0001,(rebin < 0 ? 1.0 : rebin) * dh1->Integral()/4);
            p->setUnderflow(false);
            p->setOverflow(false);
            p->setBotMinMax(0,2);
            p->setYTitle("N. of events");
            if(rebin > 0) p->rebin(rebin);
            if(withRatio){
                auto * c = p->drawSplitRatio(1,"stack",false,false,name);
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

void test2DCondTemplate(std::string name, std::string filename){
    TH2* dH = 0;
    std::vector<TH2*> hs;
    std::vector<std::string> hNs;
    auto addHistos =[&](std::string extraPre="", bool addSmooth=true, bool addKDE=false,  bool addData = false){
        TFile *f = new TFile((filename + "_"+name + (extraPre.size() ? std::string("_") + extraPre +"_" : std::string("_") ) +"incl_COND2D_template.root").c_str(),"read");

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
    make2DTests(name + "_COND2D_HbbF" ,dH,hs,hNs,hBBBinning,true,-1);
    make2DTests(name + "_COND2D_HbbC" ,dH,hs,hNs,hBBBinning,true,10);
    make2DTests(name + "_COND2D_HHF"  ,dH,hs,hNs,hhBinning,false,-1);
    make2DTests(name + "_COND2D_HHC"  ,dH,hs,hNs,hhBinning,false,5 );
}

void test2DTemplate(std::string name, std::string filename){
    TH2* dH = 0;
    std::vector<TH2*> hs;
    std::vector<std::string> hNs;
    auto cutAndRotateHistograms =[](const TH2* inH) -> TH2F*{ //The conditional templates have the Hbb as the y-axis..which is only done here
        TH2F * outH = new TH2F(TString(inH->GetName()) + "_cut",(std::string(";") +hbbMCS.title+";"+hhMCS.title).c_str()  ,nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
        for(int iX =1; iX <= inH->GetNbinsX(); ++iX){
            const int outIY =outH->GetYaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(iX));
            if(outIY < 1 || outIY > outH->GetNbinsY() ) continue;
            for(int iY =1; iY <= inH->GetNbinsY(); ++iY){
                const int outIX = outH->GetXaxis()->FindFixBin(inH->GetYaxis()->GetBinCenter(iY));
                if(outIX < 1 || outIX > outH->GetNbinsX() ) continue;
                outH->SetBinContent(outIX,outIY,inH->GetBinContent(iX,iY));
                outH->SetBinError(outIX,outIY,inH->GetBinError(iX,iY));
            }
        }
        return outH;
    };

    TFile *fCond = new TFile((filename + "_"+name  +"_incl_COND2D_template.root").c_str(),"read");
    fCond->GetObject("histo_fine_data",dH); dH = cutAndRotateHistograms(dH);
    TH2* h = 0;
    fCond->GetObject("histo_smooth",h);  hs.push_back(cutAndRotateHistograms(h));hNs.push_back("COND KDE");

    TFile * fT = new TFile((filename + "_"+name  +"_2D_template.root").c_str(),"read");
    fT->GetObject("histo",h);hs.push_back(h);hNs.push_back("Template");
    std::vector<double> hBBBinning = {30,40,50,60,80,100,120,140,170,210};
    std::vector<double> hhBinning  = {700,800,900,1000,1500,2000,3000,4000,5000};
    make2DTests(name + "_Temp_HbbF",dH,hs,hNs,hBBBinning,false);
    make2DTests(name + "_Temp_HbbC",dH,hs,hNs,hBBBinning,false,10);
    make2DTests(name + "_Temp_HHF" ,dH,hs,hNs,hhBinning,true);
    make2DTests(name + "_Temp_HHC" ,dH,hs,hNs,hhBinning,true,5);
}

void test2DFits(std::string name, std::string filename){


    TFile *fData = new TFile((filename + "_"+name  +"_distributions.root").c_str(),"read");
    TFile *fTemp = new TFile((filename + "_"+name  +"_2D_template.root").c_str(),"read");
    TFile *fTempFit = new TFile((filename + "_"+name  +"_2D_template_debug.root").c_str(),"read");

    TH2 * hOT = 0;
    fTemp->GetObject("histo",hOT);

    std::vector<double> hBBBinning = {30,50,100,150,210};
    std::vector<double> hhBinning  = {800,900,1000,1500,2000,3000,4000,5000};

    //    std::vector<double> hBBBinning = {30,210};
    //    std::vector<double> hhBinning  = {800,5000};

    std::vector<std::string> sels = {"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"};
    for(const auto& s : sels){
        TH2* dH = 0;
        std::vector<TH2*> hs;
        std::vector<std::string> hNs;
        fData->GetObject((name+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dH);
        if(dH==0) continue;
        TH2 * hF = 0;
        fTempFit->GetObject((name+"_"+s).c_str(),hF);
        if(hF == 0) continue;

        //        make2DTests(name + " Fit HbbF "+s,dH,{hF,hOT},{"Search region template","Original template"},hBBBinning,false);
        make2DTests(name + "_Fit_HbbC_"+s,dH,{hF,hOT},{"Search region template","Baseline template"},hBBBinning,false,10);
        //        make2DTests(name + " Fit HHF "+s ,dH,{hF,hOT},{"Search region template","Original template"},hhBinning,true);
        //        make2DTests(name + " Fit HHC "+s ,dH,{hF,hOT},{"Search region template","Baseline template"},hhBinning,true,5);
    }



}

void testMJJKern(std::string name, std::string filename) {
    bool withRatio = true;
    TFile *f = new TFile((filename + "_"+name +"_incl_template.root").c_str(),"read");
    std::vector<TH1*> hs;
    std::vector<std::string> hNs;
    TH1* h = 0;
    f->GetObject("histo_data",h);hs.push_back(h);hNs.push_back("MC");
    f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back("KDE");

    int binL = hs[0]->FindFixBin(minHbbMass);
    int binH = hs[0]->FindFixBin(maxHbbMass);
    for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral(binL,binH)/hs[iH]->Integral(binL,binH));

    Plotter * p = new Plotter();
    Plotter * pf = new Plotter();
    for(unsigned int iH = 0; iH < hs.size(); ++iH){
        hs[iH]->SetYTitle("N. of events");
        hs[iH]->SetXTitle(hbbMCS.title.c_str());
        if(iH == 0){
            p->addHist(hs[iH],hNs[iH].c_str());
            pf->addHist(hs[iH],hNs[iH].c_str());
        }
        else  {
            p->addHistLine(hs[iH],hNs[iH].c_str());
            pf->addHistLine(hs[iH],hNs[iH].c_str());
        }
    }

    auto setupPlotter = [&](Plotter * p, std::string name, double rebin =-1){
        p->setMinMax(.1,(rebin < 0 ? 1.0 : rebin) * hs[0]->Integral()/4);
        p->setUnderflow(false);
        p->setOverflow(false);
        p->setBotMinMax(0,2);
        if(rebin > 0) p->rebin(rebin);
        if(withRatio){
            auto * c = p->drawSplitRatio(1,"stack",false,false,name);
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

    setupPlotter(p,name + "_HbbKDE_C",5);
    setupPlotter(pf,name + "_HbbKDE_F");

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




void plotNonResBkgTests(int step = 0,bool doTW = true, bool doSR = true, std::string outName = ""){
    //    hhFilename +="_CR";
    std:: string inName = doSR ? "bkgInputs" : "bkgInputsCR";
    std::string filename = inName +"/"+hhFilename;

    CutStr mod = bkgSels [doTW ? BKG_LOSTTW : BKG_QG];
    if(outName.size()){
        outName += std::string("/") + mod +  (doSR ? "" : "CR");
    }


    switch(step){
    case 0:
        if(outName.size()) outName += "_2DCondTemp.root";
        test2DCondTemplate(mod,filename);
        break;
    case 1:
        if(outName.size()) outName += "_MJJKern.root";
        testMJJKern(mod,filename);
        break;
    case 2:
        if(outName.size()) outName += "_2DTemp.root";
        test2DTemplate(mod,filename);
        break;
    case 3:
        if(outName.size()) outName += "_2DFits.root";
        test2DFits(mod,filename);
        break;
    case 4:
        if(outName.size()) outName += "_2DComp.root";
        writeables = test2DModel({mod},filename,
                {"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});
        break;
    }

    Dummy d(outName);
}
