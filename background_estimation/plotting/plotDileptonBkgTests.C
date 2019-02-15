#include <vector>
#include "TFile.h"
#include "/Users/brentstone/Dropbox/Physics/GitRepos/HistoPlotting/include/Plotter.h"
#include "plotTestHelper.h"
#include "../predTools/DileptonCutConstants.h"
#include "../predTools/CutConstants.h"

using namespace DileptonCutConstants;
using namespace ASTypes;

std::vector<TObject*> writeables;
void addWr(std::vector<TObject*> nw) { writeables.insert( writeables.end(), nw.begin(), nw.end() );}


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
            hs.push_back(h); hNs.push_back(extraPre.size() ? extraPre +"_NoSmear"  : std::string("KDE w/o smoothing"));
        }
    };

    addHistos("xs_0p3_xc_3_ys_0p75_yc_2",true,true,true);
//    std::vector<std::string> extras = {"xs_0p3_xc_3_ys_0p75_yc_2","xs_0p3_xc_3_ys_0p75_yc_3","xs_0p3_xc_3_ys_0p75_yc_4"};
//    for(const auto& s : extras ) addHistos(s,false,true);

    std::vector<double> hBBBinning = {30,40,50,60,80,100,120,140,170,210};
    std::vector<double> hhBinning  = {700,4000,700,800,900,1000,1500,2000,3000,4000,5000};


    addWr(make2DTests(name + "_COND2D_HHF" ,dH,hs,hNs,hBBBinning,false,-1));
    addWr(make2DTests(name + "_COND2D_HHC" ,dH,hs,hNs,hBBBinning,false,10));
    addWr(make2DTests(name + "_COND2D_HbbF"  ,dH,hs,hNs,hhBinning,true,-1,false));
    addWr(make2DTests(name + "_COND2D_HbbC"  ,dH,hs,hNs,hhBinning,true,5 ,false));
}

void test2DTemplate(const std::string& name, const std::string& filename, const std::vector<std::string>& sels, const std::vector<std::string>& origsels){

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
    	for(unsigned int iS = 0; iS < sels.size(); ++iS){
    		const auto& s = sels[iS];
        TH2* dH = 0;
        TH2* h = 0;
        std::vector<TH2*> hs;
        std::vector<std::string> hNs;
        fData->GetObject((name+"_"+s+"_"+hbbMCS+"_"+hhMCS).c_str(),dH); dH = cutHistograms(dH);
        TFile * fT = new TFile((filename + "_"+name+"_"+s  +"_2D_merged_template.root").c_str(),"read");
        fT->GetObject("histo",h);
        hs.push_back(h);hNs.push_back("Template");

        if(origsels.size() && origsels[iS].size()){
            TH2* ho=0;
        	TFile * fTO = new TFile((filename + "_"+name+"_"+origsels[iS]  +"_2D_cond_template.root").c_str(),"read");
            fTO->GetObject("histo_smooth",ho);
        	ho = cutHistograms(ho);

            hs.push_back(ho);hNs.push_back("Initial template");
        }




//        std::vector<double> hBBBinning = {30,40,50,60,80,100,120,140,170,210};
//        std::vector<double> hhBinning  = {700,800,900,1000,1500,2000,3000,4000,5000};
        std::vector<double> hBBBinning = {30,210};
        std::vector<double> hhBinning  = {700,2250, 4000};
    //    addWr(make2DTests(name +"_"+s+ "_Temp_HbbF",dH,hs,hNs,hBBBinning,false));
        addWr(make2DTests(name +"_"+s+ "_Temp_HHC",dH,hs,hNs,hBBBinning,false,10));
    //    addWr(make2DTests(name +"_"+s+ "_Temp_HHF" ,dH,hs,hNs,hhBinning,true));
        addWr(make2DTests(name +"_"+s+ "_Temp_HbbC" ,dH,hs,hNs,hhBinning,true,5));
    }

}
//
//
//
//void testRatioFits(std::string name, std::string filename,  std::string fitName, const std::vector<std::string>& sels,std::vector<TObject*>& writeables) {
//    Plotter * p = new Plotter; //stupid CINT bugfix.....
//
//    std::vector<TObject*> paramPads;
////    gROOT->SetBatch(true);
//    for(const auto& s : sels){
//        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");
//        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
//            TGraphErrors * can= 0;
//            if(ff){
//                ff->GetObject((name).c_str(),can);
//            }
//            if(!can) return;
//            TCanvas * c = new TCanvas();
//            can->GetYaxis()->SetTitle(getCategoryLabel(s).c_str());
//            can->GetXaxis()->SetTitle(hhMCS.title.c_str());
//            can->GetYaxis()->SetRangeUser(0.05,5);
//            can->Draw();
//            c->SetLogy();
//            c->Update();
//            list.push_back(c);
//        };
//
//
//        addGraph("RATIO", paramPads);
//    }
////    gROOT->SetBatch(false);
//
//
//
//    auto * c1 =Drawing::drawAll(paramPads,"ratio_params");
//    c1->SetTitle("ratio_params");
//    writeables.push_back(c1);
//}
//
//
//void testRatioUncs(std::string name, std::string filename,  std::string fitName, const std::vector<std::pair<std::string,std::string> >& targets, //target,model
//        float normUnc, float scaleUnc, float resUnc, std::vector<TObject*>& writeable) {
//    Plotter * p = new Plotter; //stupid CINT bugfix.....
//    std::vector<TObject*> paramPads;
//
//    for(const auto& s : targets){
//        TFile *ff = new TFile((filename+"_"+name+"_"+s.second+"_"+fitName+".json.root").c_str(),"read");
//        TFile *ft = new TFile((filename+"_"+name+"_"+s.first+"_"+fitName+".json_inputDebug.root").c_str(),"read");
//
//        if(ff == 0) continue;
//        if(ft == 0) continue;
//
//        TGraphErrors * graph= 0;
//        ft->GetObject("RATIOPlusOne",graph);
//        if(!graph) continue;
//        graph->GetYaxis()->SetTitle(s.first.c_str());
//        TCanvas * can = new TCanvas((std::string("syst_")+s.first).c_str());
//        can->cd();
//        graph->Draw();
//
//        TF1 * defF= 0;
//        ff->GetObject("RATIO_func",defF);
//        if(!defF) continue;
//
//        auto mkF = [&](std::string name,std::string factor, int color, int style){
//            std::string fs = "(1+[0]+[1]*x+[2]/x)";
//            if(factor.size()) fs += "*("+factor+")";
//            TF1 * f = new TF1((s.first+"_"+name+"_mix").c_str(),fs.c_str(),1,13000);
//            f->SetParameter(0,defF->GetParameter(0));
//            f->SetParameter(1,defF->GetParameter(1));
//            f->SetParameter(2,defF->GetParameter(2));
//            f->SetLineColor(color);
//            f->SetLineStyle(style);
//            can->cd();
//            f->Draw("SAME");
//        };
//
//        mkF("nominal","",StyleInfo::getLineColor(0),1 );
//        mkF("PTUp","(1. + 0.0003*x)",StyleInfo::getLineColor(1),1 );
//        mkF("PTDown","1/(1. + 0.0003*x)",StyleInfo::getLineColor(1),9 );
//        mkF("OPTUp","(1. + 1200/x)",StyleInfo::getLineColor(2),1 );
//        mkF("OPTDown","1/(1. + 1200/x)",StyleInfo::getLineColor(2),9 );
//        mkF("NORMUp","1.5",StyleInfo::getLineColor(3),1 );
//        mkF("NormDown","0.5",StyleInfo::getLineColor(3),9 );
//        writeables.push_back(can);
//    }
//    //    auto * c1 =Drawing::drawAll(paramPads,"ratio_params");
//    //    writeables.push_back(c1);
//}

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



void plotDileptonBkgTests(int step = 0, int inreg = REG_SR, std::string outName = ""){
    REGION reg = REGION(inreg);

    std:: string inName =  "bkgInputs" ;

    auto getSRCats = [&](REGION reg) {
        std::vector<std::string> sels;
        for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& h :hadCuts){
            if(l == lepCats[LEP_INCL]) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(h != hadCuts[HAD_FULL]) continue;
            sels.emplace_back(l +"_"+b+"_"+h);
        }
        return sels;
    };

    auto srList = getSRCats(reg);
    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
    }
    else if(reg == REG_QGCR){
        inName =  "bkgInputsQGCR";
        hhFilename +="_QGCR";
    }
    std::string filename = inName +"/"+hhFilename;

    CutStr mod = bkgSels[BKG_NOM];
    if(outName.size()){
        outName += std::string("/") + mod;
        if(reg == REG_TOPCR) outName +=  "_TopCR";
        else if(reg == REG_QGCR) outName +=  "_QGCR";
    }


    std::vector<std::string> stepSels;
    std::vector<std::string> origSels;
    switch(step){

    case 2:
        if(outName.size()) outName += "_2DCondTemp";
        test2DCondTemplate(mod,filename,"inclFlav_LMT_R_phi_b");
        break;
    case 3:
        if(outName.size()) outName += "_MVVKern";
        else stepSels = {"oppFlav_LMT_R_phi_b","sameFlav_LMT_R_phi_b","inclFlav_LMT_R_phi_b"};
        writeables = test1DKern(mod,filename,"MVV",stepSels);
        break;
    case 4:
        if(outName.size()) outName += "_2DTemp";
        stepSels = {"oppFlav_LMT_R_phi_b","sameFlav_LMT_R_phi_b","inclFlav_LMT_R_phi_b"};
        origSels = {"xs_0p3_xc_3_ys_0p75_yc_2_inclFlav_LMT_R_phi_b","xs_0p3_xc_3_ys_0p75_yc_2_inclFlav_LMT_R_phi_b","xs_0p3_xc_3_ys_0p75_yc_2_inclFlav_LMT_R_phi_b"};

        test2DTemplate(mod,filename,stepSels,origSels);
        break;
    case 5:
        if(outName.size()) outName += "_2DFits";
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
    	CutConstants::CutStr mod2("nomBKG","0==0","Bkg Template");
        if(outName.size()) outName += "_2DComp";
        writeables = test2DModel({mod2},filename,
                srList,{30,210},false);
        break;
    }

    Dummy d(outName);
}
