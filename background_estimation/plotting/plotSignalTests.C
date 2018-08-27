
#include "plotTestHelper.h"
#include "TH1.h"
#include "TH2.h"

std::vector<TObject*> writeables;

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



void test2DFits(std::string name, std::string filename, const std::vector<int>& signalMassBins, std::string fitName,bool binInY, const std::vector<std::string>& sels, std::string outName = "") {
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<double> bins = {30,210,30,115,135,210};

    for(const auto& s : sels){
        TFile * fo =0;
        fo=new TFile((filename+"_"+name+"_"+s+"_"+fitName+".root").c_str(),"read");
        if(fo == 0) continue;
        TFile *ff = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".json.root").c_str(),"read");
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
            can->GetXaxis()->SetTitle(sigMCS.title.c_str());
            list.push_back(can);
        };
        std::vector<TObject*> paramPads;
        auto vX=[&](const std::string& v )->std::string{ return v +MOD_MJ;};
        auto vY=[&](const std::string& v )->std::string{ return v +MOD_MR;};
//        addGraph(vX("mean"  ), paramPads);
//        addGraph(vX("sigma" ), paramPads);
//        addGraph(vX("alpha" ), paramPads);
//        addGraph(vX("alpha2"), paramPads);
//        addGraph(vX("n"     ), paramPads);
//        addGraph(vX("n2"    ), paramPads);
//        addGraph(vX("slope" ), paramPads);
//        addGraph(vX("fE"    ), paramPads);

        addGraph(vY("mean_p1" ), paramPads);
        addGraph(vY("sigma_p1" ), paramPads);

        addGraph(vY("mean"  ), paramPads);
        addGraph(vY("sigma" ), paramPads);
        addGraph(vY("alpha" ), paramPads);
        addGraph(vY("alpha2"), paramPads);
        addGraph(vY("n"     ), paramPads);
        addGraph(vY("n2"    ), paramPads);


        TCanvas * c =Drawing::drawAll(compPads, (s + "_COMP").c_str());
        TCanvas * c1 =Drawing::drawAll(paramPads, (s + "_params").c_str());
        c->SetTitle((s + "_COMP").c_str());
        c1->SetTitle((s + "_params").c_str());
        writeables.push_back(c);
        writeables.push_back(c1);

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
            can->GetXaxis()->SetTitle(sigMCS.title.c_str());
            list.push_back(can);
        };

        addGraph("yield", paramPads);
    }
    auto * c = Drawing::drawAll(paramPads,"Yields");
    c->SetTitle((name+"_yields").c_str());
    c->SetName((name+"_yields").c_str());
    writeables.push_back(c);

}

void plotEfficiencies(std::string name, std::string filename,std::string fitName) {
//    gROOT->SetBatch(true);
    Plotter * p = new Plotter; //stupid CINT bugfix.....
    std::vector<TObject*> paramPads;
    double inclusiveN  = 1000*35.9; //1000 fb x lumi
//    inclusiveN *= 2*0.5824*(.2137+.002619); //BR to bbVV
    inclusiveN *= 2*0.5824*(.2137); //BR to bbWW
    inclusiveN *= 2*.676*(.216+.108*.3524);//BR of WW to he, hmu or htau where the tau is leptonic
    for(auto& l : lepCats){
        if(l == lepCats[LEP_EMU]) continue;
        for(auto& p : purCats){
            if(p == purCats[PURE_I]) continue;
            for(auto& h : hadCuts){
                if(h != hadCuts[HAD_FULL]) continue;
                Plotter * plot = new Plotter();
                for(auto& b : btagCats){
                    std::string s = l+"_"+b+"_"+p+"_"+h;
                    TFile * fo = new TFile((filename+"_"+name+"_"+s+"_"+fitName+".root").c_str(),"read");
                    if(fo == 0) continue;
                    TGraphErrors * gr= 0;
                    fo->GetObject("yield",gr);
                    if(gr == 0) continue;
                    gr->GetXaxis()->SetTitle(sigMCS.title.c_str());
                    gr->GetYaxis()->SetTitle("efficiency");
                    for(int iP = 0; iP < gr->GetN(); ++iP){
                        double x,y,ey;
                        gr->GetPoint(iP,x,y);
                        ey = gr->GetErrorY(iP);
                        gr->SetPoint(iP,x,y/inclusiveN);
                        gr->SetPointError(iP,0,ey/inclusiveN);
                    }
                    plot->addGraph(gr,b.title);
                }
                plot->setMinMax(0,.1);
                plot->addText(l.title+", "+p.title,0.1454849,0.866087,0.04);
                plot->setLegendPos(0.5016722,0.7286957,0.9682274,0.9165217);
                auto c = plot->draw(false,(name+"_"+l+"_"+p+"_"+h+"_sigeff") .c_str());
                writeables.push_back(c);
            }
        }

    }
//    gROOT->SetBatch(false);


}


std::vector<TObject*> testSignal1DFits(std::string name, std::string filename, const std::vector<int>& signalMassBins, std::string varName, std::string fitName, const std::vector<std::string>& sels){
    std::vector<std::string> canNames;
    for(const auto& sB : signalMassBins){ canNames.push_back(std::string("can_m") +int2Str(sB));}
    return test1DFits(name,filename, varName,fitName,sels,canNames);
}

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



void plotSignalTests(int cat = 0,int sig = RADION,  bool doCond = true, std::string outName = ""){
    std:: string inName = doCond ? "signalInputs" : "signalInputsNoCond";
    std::string filename = inName +"/"+hhFilename;
    if(outName.size()) outName += doCond ? std::string("/signal") : std::string("/signalNoCond");
    std::string name = signals[sig];


    switch(cat) {
    case 0:
        if(outName.size()) outName += "_yield";
        plotYields(name,filename,"yield",{"e_L_LP_full","e_M_LP_full","e_T_LP_full","e_L_HP_full","e_M_HP_full","e_T_HP_full","mu_L_LP_full","mu_M_LP_full","mu_T_LP_full","mu_L_HP_full","mu_M_HP_full","mu_T_HP_full"});
        plotEfficiencies(name,filename,"yield" );
        break;
    case 1:
        if(outName.size()) outName += "_MJJ_fit1stIt";
        writeables = testSignal1DFits(name,filename,signalMassBins[sig],MOD_MJ,"MJJ_fit1stIt",{"emu_LMT_I_ltmb","emu_L_I_ltmb","emu_M_I_ltmb","emu_T_I_ltmb"});
        break;
    case 2:
        if(outName.size()) outName += "_MJJ_fit";
        writeables = testSignal1DFits(name,filename,signalMassBins[sig],MOD_MJ,"MJJ_fit",{"emu_LMT_I_ltmb","emu_L_I_ltmb","emu_M_I_ltmb","emu_T_I_ltmb"});
        break;
    case 3:
        if(outName.size()) outName += "_MVV_fit";
        if(doCond) writeables =  testSignal1DFits(name,filename,signalMassBins[sig],MOD_MR,"MVV_fit1stIt",{"e_LMT_I_ltmb","mu_LMT_I_ltmb"});
        else  writeables =  testSignal1DFits(name,filename,signalMassBins[sig],MOD_MR,"MVV_fit1stIt",{"e_LMT_LP_full","e_LMT_HP_full","mu_LMT_LP_full","mu_LMT_HP_full"});
        break;
    case 4:
        if(outName.size()) outName += "_2D_fit";
        test2DFits(name,filename,signalMassBins[sig],"2D_fit",false,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},outName);
        break;
    }

    Dummy d(outName);



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
