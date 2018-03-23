
#include "plotTestHelper.h"
#include "TH1.h"
#include "TH2.h"





TCanvas* make2DTests(std::string plotTitle, int mass, const TH2* dH,TH2* pH, const std::vector<double>& bins, bool binInY, double rebin = -1) {
    const int binXmin = dH->GetXaxis()->FindFixBin(30);
    const int binXmax = dH->GetXaxis()->FindFixBin(210);
    const int binYmin = dH->GetYaxis()->FindFixBin(100);
    const int binYmax = dH->GetYaxis()->FindFixBin(6999);

    double dataINT = dH->Integral(binXmin,binXmax,binYmin,binYmax  );
    double pdfINT  = pH->Integral(binXmin,binXmax,binYmin,binYmax  );
    pH->Scale(dataINT/pdfINT);

    const TAxis * ax =  binInY ? dH->GetYaxis() : dH->GetXaxis();

    Plotter * p = new Plotter();
    int iC = 0;
    for(unsigned int iB = 0; iB + 1 < bins.size(); ++iB){
        if(bins[iB+1]<= bins[iB]) continue;
        std::string binName = std::string(binInY ? "#it{m}_{HH} " : "#it{m}_{H#rightarrowbb} ") + ASTypes::flt2Str(bins[iB]) +"-"+ASTypes::flt2Str(bins[iB+1])+" GeV";


        int binL = ax->FindFixBin(bins[iB]);
        int binH = ax->FindFixBin(bins[iB+1]) -1;
        auto proj =[&](const TH2* h, const std::string& title) ->TH1*{
            return binInY ? h->ProjectionX( (title+"_"+int2Str(iB)).c_str(),binL,binH) :  h->ProjectionY( (title+"_"+int2Str(iB)).c_str(),binL,binH);
        };

        auto dh1 = proj(dH,"MC");
        auto ph1 = proj(pH,"PDF");
        for(int iX = 0; iX <= ph1->GetNbinsX()+1; ++iX) ph1->SetBinError(iX,0);
        p->addHist(dh1,binName, StyleInfo::getLineColor(iC),1);
        p->addHistLine(ph1,binName, StyleInfo::getLineColor(iC),9);
        iC++;
    }


    auto setupPlotter = [&](Plotter * p, std::string name) ->TCanvas*{
        //        p->setMinMax(.0001,(rebin < 0 ? 1.0 : rebin) * dh1->Integral()/4);
        p->setUnderflow(false);
        p->setOverflow(false);
        p->setBotMinMax(0,2);
        p->setYTitle(std::string("m")+int2Str(mass));
        if(rebin > 0) p->rebin(rebin);
        auto * c = p->draw(false,name);
        p->xAxis()->SetRangeUser(0.6*mass,1.4*mass);
        return c;
    };
    return setupPlotter(p,plotTitle);
}



void test2DFits(std::string name, std::string filename, std::string fitName,bool binInY, const std::vector<std::string>& sels) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<double> bins = {30,210,30,115,135,210};


    for(const auto& s : sels){
        TFile * fo =0;
        fo=new TFile((filename+"_"+name+"_"+s+"_"+fitName+".root").c_str(),"read");
        if(fo == 0) continue;
        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");
        if(ff == 0) continue;
        auto addH = [&](const std::string& name,std::vector<TObject*>& list)->bool{
            TH2 * can= 0;
            fo->GetObject(name.c_str(),can);
            if(can==0) return false;
            list.push_back(can);
            return true;
        };

        std::vector<TObject*> mcPads;
        std::vector<TObject*> pdfPads;
        std::vector<TObject*> compPads;
        gROOT->SetBatch(true);
        for(const auto& sB : signalMassBins){
            if(addH(std::string("data_m") +int2Str(sB) +"__"+MOD_MJ+"_"+MOD_MR, mcPads) && addH(std::string("pdf_m") +int2Str(sB) +"__"+MOD_MJ+"_"+MOD_MR, pdfPads)){
                if(binInY){
                    std::vector<double> cbins = {700,4000};
                    cbins.push_back(700);
                    cbins.push_back(float(sB)*.95);
                    cbins.push_back(float(sB)*1.05);
                    cbins.push_back(4000);
                    compPads.push_back(make2DTests(s+"_m"+int2Str(sB),sB,(TH2*)mcPads.back(),(TH2*)pdfPads.back(),cbins,binInY,2));
                    } else
                        compPads.push_back(make2DTests(s+"_m"+int2Str(sB),sB,(TH2*)mcPads.back(),(TH2*)pdfPads.back(),bins,binInY,2));
                }
        }
        gROOT->SetBatch(false);
        if(mcPads.size()==0)continue;


        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
            TGraphErrors * can= 0;
            std::cout <<name<<" ";
            if(ff){
                ff->GetObject((name).c_str(),can);
                std::cout <<"ff: "<<can<<" ";
            }
            if(can==0) fo->GetObject(name.c_str(),can);
            std::cout <<"fo: "<<can<<" ";
            if(can == 0) return;
            std::cout <<"fl: "<<can<<" "<<std::endl;
            can->GetYaxis()->SetTitle(name.c_str());
            list.push_back(can);
        };
        std::vector<TObject*> paramPads;
        auto vX=[&](const std::string& v )->std::string{ return v +MOD_MJ;};
        auto vY=[&](const std::string& v )->std::string{ return v +MOD_MR;};
        addGraph(vX("mean"  ), paramPads);
        addGraph(vX("sigma" ), paramPads);
        addGraph(vX("alpha" ), paramPads);
        addGraph(vX("alpha2"), paramPads);
        addGraph(vX("n"     ), paramPads);
        addGraph(vX("n2"    ), paramPads);
        addGraph(vX("slope" ), paramPads);
        addGraph(vX("fE"    ), paramPads);

        addGraph(vY("mean_p1" ), paramPads);
        addGraph(vY("sigma_p1" ), paramPads);

        addGraph(vY("mean"  ), paramPads);
        addGraph(vY("sigma" ), paramPads);
        addGraph(vY("alpha" ), paramPads);
        addGraph(vY("alpha2"), paramPads);
        addGraph(vY("n"     ), paramPads);
        addGraph(vY("n2"    ), paramPads);



        Drawing::drawAll(mcPads, (s + ": MC").c_str(),"COLZ");
        Drawing::drawAll(pdfPads, (s + ": PDF").c_str(),"COLZ");
        Drawing::drawAll(compPads, (s + ": COMP").c_str());
        Drawing::drawAll(paramPads, (s + ": params").c_str());

    }

}

void plotYields(std::string name, std::string filename,std::string fitName, const std::vector<std::string>& sels) {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<TObject*> paramPads;

    for(const auto& s : sels){
        TFile * fo = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".root").c_str(),"read");
        if(fo == 0) continue;

        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");

        auto addGraph = [&](const std::string& name,std::vector<TObject*>& list){
            TGraphErrors * can= 0;
            if(ff){
                ff->GetObject((name).c_str(),can);
            }
            if(!can) fo->GetObject(name.c_str(),can);
            if(!can) return;
            can->GetYaxis()->SetTitle(s.c_str());
            list.push_back(can);
        };

        addGraph("yield", paramPads);
    }
    Drawing::drawAll(paramPads,"Yields");

}


void testSignal1DFits(std::string name, std::string filename, std::string varName, std::string fitName, const std::vector<std::string>& sels){
    std::vector<std::string> canNames;
    for(const auto& sB : signalMassBins){ canNames.push_back(std::string("can_m") +int2Str(sB));}
    test1DFits(name,filename, varName,fitName,sels,canNames);
}

void plotSignalTests(){
    std::string filename = hhFilename;
//        plotYields(radionSig,filename,"yield",{"emu_LMT_I_full","e_L_LP_full","e_M_LP_full","e_T_LP_full","e_L_HP_full","e_M_HP_full","e_T_HP_full","mu_L_LP_full","mu_M_LP_full","mu_T_LP_full","mu_L_HP_full","mu_M_HP_full","mu_T_HP_full"});

//        testSignal1DFits(radionSig,filename,MOD_MJ,"MJJ_fit1stIt",{"emu_LMT_I_ltmb","emu_L_I_ltmb","emu_M_I_ltmb","emu_T_I_ltmb"});
//    testSignal1DFits(radionSig,filename,MOD_MJ,"MJJ_fit",{"emu_LMT_I_ltmb","emu_L_I_ltmb","emu_M_I_ltmb","emu_T_I_ltmb"});
//    testSignal1DFits(radionSig,filename,MOD_MR,"MVV_fit1stIt",{"e_LMT_I_ltmb","mu_LMT_I_ltmb"});
    test2DFits(radionSig,filename,"2D_fit",false,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"});



    //OLD
    //        testSignal1DFits(radionSig,filename,"MJJ","MJJ_fit1stIt",{"emu_LMT_ltmb","emu_L_ltmb","emu_T_ltmb"});
//            testSignal1DFits(radionSig,filename,MOD_MJ,"MJJ_fit",{"emu_LMT_ltmb","emu_M_ltmb","emu_T_ltmb"});
//    testSignal1DFits(radionSig,filename,"MVV","MVV_fit",{"e_LMT_ltmb","mu_LMT_ltmb"});

//        compFitParams(radionSig,filename,"MVV","MVV_fit", "e v mu", {"e_LMT_ltmb","mu_LMT_ltmb"} );
//        compFitParams(radionSig,filename,"MVV","MVV_fit", "full v ltmb", {"mu_LMT_ltmb","mu_LMT_full","e_LMT_ltmb","e_LMT_full"} );
//        compFitParams(radionSig,filename,"MVV","MVV_fit", "L v M", {"mu_L_ltmb","mu_M_ltmb","mu_T_ltmb","e_L_ltmb","e_M_ltmb","e_T_ltmb"} );
//    test2DFits(radionSig,filename,"MVV","2D_fit1stIt");
//    test2DFits(radionSig,filename,"MVV","2D_fit");



    //Moving to MJ cond on MR
//    testSignal1DFits(radionSig,filename,MOD_MR,"MVV_fit1stIt",{"e_LMT_ltmb","mu_LMT_ltmb","emu_LMT_ltmb"});
//    testSignal1DFits(radionSig,filename,MOD_MR,"MVV_fit",{"e_LMT_ltmb","mu_LMT_ltmb","emu_LMT_ltmb"});
//    testSignal1DFits(radionSig,filename,MOD_MJ,"MJJ_fit1stIt",{"emu_LMT_ltmb","emu_L_ltmb","emu_M_ltmb","emu_T_ltmb"});
//    test2DFits(radionSig,filename,"2D_fit1stIt",true,{"emu_LMT_ltmb","e_LMT_ltmb","mu_LMT_ltmb"});


//        test2DFits(radionSig,filename,"2D_fit1stIt",true,{"e_LMT_ltmb","mu_LMT_ltmb"});
//                test2DFits(radionSig,filename,"2D_fit1stIt",false,{"e_LMT_ltmb","mu_LMT_ltmb"});

//            test2DFits(radionSig,filename,"2D_fit",true,{"e_L_full","mu_L_full","e_M_full","mu_M_full","e_T_full","mu_T_full"});
//                    test2DFits(radionSig,filename,"2D_fit",false,{"e_L_full","mu_L_full","e_M_full","mu_M_full","e_T_full","mu_T_full"});

//    test2DFits(radionSig,filename,"2D_fit",false,{"e_L_full","mu_L_full"});
//    test2DFits(radionSig,filename,"2D_fit",true,{"e_M_full","mu_M_full"});
//
//    testSignal1DFits(radionSig,filename,MOD_MJ,"MJJ_fit",{"emu_LMT_ltmb"});
//        test2DFits(radionSig,filename,"2D_fit",false,{"e_L_full","e_M_full","e_T_full"});

}
