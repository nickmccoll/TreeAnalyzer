
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



void test2DFits(std::string name, std::string filename, const std::vector<int>& signalMassBins, std::string fitName,bool binInY, const std::vector<std::string>& sels) {
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
            can->GetYaxis()->SetTitle(getCategoryLabel(s).c_str());
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

void plotNormSysts(std::string name, std::string filename,const std::vector<int>& signalMassBins,std::vector<std::string> systs) {
    const std::string inName =  filename+"_"+name + "_normSyst_distributions.root";
    TFile * f = new TFile(inName.c_str(), "read");
    std::vector<std::string> sels;
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_E]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_MU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_L]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_M]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_T]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
//    for(auto& l : lepCats)for(auto& b : btagCats)for(auto& p : purCats)for(auto& h : hadCuts){
//        if(l == lepCats[LEP_EMU]) continue;
//        if(b == btagCats[BTAG_LMT]) continue;
//        if(p == purCats[PURE_I]) continue;
//        if(h != hadCuts[HAD_FULL]) continue;
//        sels.push_back( l+"_"+b+"_"+p+"_"+h);
//    }
    for(auto& sys : systs){
        gROOT->SetBatch(true);
        std::vector<TObject*>plots;
        for(auto& sel : sels){
            TGraph * gu = new TGraph;
            TGraph * gd = new TGraph;
            TGraph * ga = new TGraph;
            int nP = 0;
            std::vector<float> vals;
            for(auto& sM : signalMassBins){
                const std::string sigName =name+"_m"+ASTypes::int2Str(sM);
                TH1 * hN = 0;
                f->GetObject((sigName+"_"+sel+"_"+hhMCS).c_str(),hN);
                TH1 * hU = 0;
                f->GetObject((sigName+"_"+sys+"Up"+"_"+sel+"_"+hhMCS).c_str(),hU);
                TH1 * hD = 0;
                f->GetObject((sigName+"_"+sys+"Down"+"_"+sel+"_"+hhMCS).c_str(),hD);
                if(hN ==0||hU==0||hD==0)continue;
                float nU = hU->Integral();
                float nN = hN->Integral();
                float nD = hD->Integral();
                gu->SetPoint(nP,float(sM),nU/nN-1.0);
                gd->SetPoint(nP,float(sM),nD/nN-1.0);
                ga->SetPoint(nP,float(sM),(nU-nD)/(2*nN));
                vals.push_back(nU/nN-1.0);
                vals.push_back(nD/nN-1.0);
                vals.push_back((nU-nD)/(2*nN));
                nP++;
            }
            Plotter * p = new Plotter;
            p->addGraph(gu, (sys+": Up").c_str());
            p->addGraph(gd, (sys+": Down").c_str());
            p->addGraph(ga, (sys+": Avg").c_str());
            p->setYTitle(sel);
            std::sort(vals.begin(),vals.end(),[](const float a, const float b){return a < b;} );
            if(vals.size()){
                float min = vals[0] < 0 ? vals[0]*1.2 : vals[0]*0.8;
                float max = vals.back() < 0 ? vals.back()*0.8 : vals.back()*1.2;
                p->setMinMax( min,max );
            }

            auto c = p->draw(false,sys+"_"+sel);
            plots.push_back(c);
        }
        gROOT->SetBatch(false);
        auto c = Drawing::drawAll(plots,sys);
        c->SetName(sys.c_str());
        c->SetTitle(sys.c_str());
        writeables.push_back(c);
    }
}

void plotPDFSysts(std::string name, std::string filename,const std::vector<int>& signalMassBins) {
    const std::string inName =  filename+"_"+name + "_pdfSyst_distributions.root";
    TFile * f = new TFile(inName.c_str(), "read");


    std::vector<std::string> sels;
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);

        gROOT->SetBatch(true);
        std::vector<TObject*>plots;
        for(auto& sel : sels){
            TGraph * gn = new TGraph;
            TGraph * gsu = new TGraph;
            TGraph * gsd = new TGraph;
            TGraph * gpu = new TGraph;
            TGraph * gpd = new TGraph;
            TGraph * gtu = new TGraph;
            TGraph * gtd = new TGraph;
            int nP = 0;
            for(auto& sM : signalMassBins){
                const std::string sigName =name+"_m"+ASTypes::int2Str(sM);


            TH1 * hIN = 0;
            f->GetObject((sigName+"_incl_yield").c_str(),hIN);
            TH1 * hN = 0;
            f->GetObject((sigName+"_"+sel+"_yield").c_str(),hN);
            if(hIN ==0) continue;
            if(hN ==0) continue;
            double nomAcc = hN->GetBinContent(1)/hIN->GetBinContent(1);


            std::vector<float> scaleAccs;
            for(unsigned int i = 0; i < 6; ++i){
                TH1 * hINS = 0;
                f->GetObject((sigName+"_s_"+ ASTypes::int2Str(i) +"_incl_yield").c_str(),hINS);
                TH1 * hNS = 0;
                f->GetObject((sigName+"_s_"+ ASTypes::int2Str(i)+"_"+sel+"_yield").c_str(),hNS);
//                std::cout <<sigName+"_s_"+sel+"_yield " << hINS <<" "<< hNS <<std::endl;
                double acc = hNS->GetBinContent(1)/hINS->GetBinContent(1);
//                double diff = std::fabs(acc-nomAcc)/nomAcc;
                scaleAccs.push_back(acc);
//                if(diff > maxDiff) maxDiff = diff;
            }
            gn->SetPoint(nP,float(sM),nomAcc);
            std::sort(scaleAccs.begin(),scaleAccs.end(),[](float a, float b){return a<b;});
            gsu->SetPoint(nP,float(sM),scaleAccs.back() /nomAcc -1.0 );
            gsd->SetPoint(nP,float(sM),scaleAccs.front()/nomAcc-1.0);

            double x = 0;
            double x2 = 0;
            for(unsigned int i = 0; i < 100; ++i){
                TH1 * hINS = 0;
                f->GetObject((sigName+"_p_"+ ASTypes::int2Str(i) +"_incl_yield").c_str(),hINS);
                TH1 * hNS = 0;
                f->GetObject((sigName+"_p_"+ ASTypes::int2Str(i)+"_"+sel+"_yield").c_str(),hNS);
                double acc = hNS->GetBinContent(1)/hINS->GetBinContent(1);
                x += acc;
                x2 += acc*acc;
            }
            double mean = x/100.;
            double sqmean = x2/100.;
            double stdDev = std::sqrt(sqmean-mean*mean);
            gpu->SetPoint(nP,float(sM),(mean+stdDev)/nomAcc-1.0);
            gpd->SetPoint(nP,float(sM),(mean-stdDev)/nomAcc-1.0);

            gtu->SetPoint(nP,float(sM),std::sqrt((mean+stdDev)*(mean+stdDev)+scaleAccs.back() *scaleAccs.back() )/nomAcc);
            gtd->SetPoint(nP,float(sM),std::sqrt((mean-stdDev)*(mean-stdDev)+scaleAccs.front()*scaleAccs.front())/nomAcc);

            nP++;
            }
            Plotter * p = new Plotter;
//            p->addGraph(gn,"nom");
            p->addGraph(gsd,"scale: Down");
            p->addGraph(gsu, "scale: Up");
            p->addGraph(gpd,"pdf: Down");
            p->addGraph(gpu, "pdf: Up");

            p->addGraph(gtd,"tot: Down");
            p->addGraph(gtu, "tot: Up");

            p->setMinMax(-0.5,0.5);
            p->setYTitle(sel);

            auto c = p->draw(false,sel);
            plots.push_back(c);
        }
        gROOT->SetBatch(false);
        auto c = Drawing::drawAll(plots,"PDFSyst");
        c->SetName("PDFSyst");
        c->SetTitle("PDFSyst");
        writeables.push_back(c);
}



void plotShapeSystNorms(std::string name, std::string filename,const std::vector<int>& signalMassBins,std::vector<std::string> systs) {
    const std::string inName =  filename+"_"+name + "_NOMSyst_distributions.root";
    TFile * fn = new TFile(inName.c_str(), "read");
    std::vector<TFile*> systUpFs;
    std::vector<TFile*> systDownFs;
    for(unsigned int iS = 0; iS < systs.size(); ++iS ){
        systUpFs.push_back(new TFile((filename+"_"+name + "_"+systs[iS]+"Up_distributions.root").c_str()));
        systDownFs.push_back(new TFile((filename+"_"+name + "_"+systs[iS]+"Down_distributions.root").c_str()));
    }
    std::vector<std::string> sels;
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_E]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_MU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_L]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_M]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_T]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_FULL]);

    TF1 *f1 = new TF1("f1","pol0",750,3600);
    for(unsigned int iS = 0; iS < systs.size(); ++iS ){
        auto sys = systs[iS];
        gROOT->SetBatch(true);
        std::vector<TObject*>plots;
        for(auto& sel : sels){
            TGraph * gu = new TGraph;
            TGraph * gd = new TGraph;
            TGraph * ga = new TGraph;
            int nP = 0;
            std::vector<float> vals;
            for(auto& sM : signalMassBins){
                const std::string sigName =name+"_m"+ASTypes::int2Str(sM);
                TH2 * hN = 0;
                fn->GetObject((sigName+"_"+sel+"_"+hbbMCS+"_"+hhMCS).c_str(),hN);
                TH2 * hU = 0;
                systUpFs[iS]->GetObject((sigName+"_"+sel+"_"+hbbMCS+"_"+hhMCS).c_str(),hU);
                TH2 * hD = 0;
                systDownFs[iS]->GetObject((sigName+"_"+sel+"_"+hbbMCS+"_"+hhMCS).c_str(),hD);
                if(hN ==0||hU==0||hD==0)continue;

                int lbx, lby, hbx, hby;
                lbx = hN->GetXaxis()->FindFixBin(minHbbMass);
                hbx = hN->GetXaxis()->FindFixBin(maxHbbMass);
                lby = hN->GetYaxis()->FindFixBin(minHHMass);
                hby = hN->GetYaxis()->FindFixBin(maxHHMass);
                std::cout << lbx <<" "<<hbx <<" "<<lby <<" "<< hby <<" "<< hN->Integral(lbx,hbx,lby,hby) <<std::endl;
                float nU = hU->Integral(lbx,hbx,lby,hby);
                float nN = hN->Integral(lbx,hbx,lby,hby);
                float nD = hD->Integral(lbx,hbx,lby,hby);
                gu->SetPoint(nP,float(sM),nU/nN-1.0);
                gd->SetPoint(nP,float(sM),nD/nN-1.0);
                ga->SetPoint(nP,float(sM),(nU-nD)/(2*nN));
                vals.push_back(nU/nN-1.0);
                vals.push_back(nD/nN-1.0);
                vals.push_back((nU-nD)/(2*nN));
                nP++;
            }
            Plotter * p = new Plotter;
            p->addGraph(gu, (sys+": Up").c_str());
            p->addGraph(gd, (sys+": Down").c_str());
            auto ng = p->addGraph(ga, (sys+": Avg").c_str());
            p->setYTitle(sel);
            std::sort(vals.begin(),vals.end(),[](const float a, const float b){return a < b;} );
            if(vals.size()){
                float min = vals[0] < 0 ? vals[0]*1.2 : vals[0]*0.8;
                float max = vals.back() < 0 ? vals.back()*0.8 : vals.back()*1.2;
                p->setMinMax( min,max );
            }

            auto c = p->draw(false,sys+"_"+sel);
            c->cd();
            ng->Fit(f1,"","",750,3600);
            f1->Draw("SAME");
            plots.push_back(c);
        }
        gROOT->SetBatch(false);
        auto c = Drawing::drawAll(plots,sys);
        c->SetName(sys.c_str());
        c->SetTitle(sys.c_str());
        writeables.push_back(c);
    }
}

void plotShapeSystParams(std::string name, std::string filename,const std::vector<int>& signalMassBins,std::vector<std::string> systs) {
    const std::string inName =  filename+"_"+name + "_NOMSyst_distributions.root";
    TFile * fn = new TFile(inName.c_str(), "read");
    std::vector<TFile*> systUpFs;
    std::vector<TFile*> systDownFs;
    for(unsigned int iS = 0; iS < systs.size(); ++iS ){
        systUpFs.push_back(new TFile((filename+"_"+name + "_"+systs[iS]+"Up_distributions.root").c_str()));
        systDownFs.push_back(new TFile((filename+"_"+name + "_"+systs[iS]+"Down_distributions.root").c_str()));
    }
    std::vector<std::string> sels;
    sels.push_back(lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB]);
//    sels.push_back(lepCats[LEP_E]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB]);/
//    sels.push_back(lepCats[LEP_MU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB]);
    std::vector<std::string> params = {"meanMR","sigmaMR"};

    for(unsigned int iS = 0; iS < systs.size(); ++iS ){
        auto sys = systs[iS];
        gROOT->SetBatch(true);
        std::vector<TObject*>plots;
        for(auto& sel : sels){

            TFile* fU = new TFile((filename+"_"+name + "_"+sel+ "_"+systs[iS]+"Up_MVV_fit1stIt.root").c_str());
            TFile* fD = new TFile((filename+"_"+name + "_"+sel+ "_"+systs[iS]+"Down_MVV_fit1stIt.root").c_str());
            TFile* fN = new TFile((filename+"_"+name + "_"+sel+ "_"+"NOMSyst_MVV_fit1stIt.root").c_str());
            if(fU ==0||fD==0||fN==0)continue;
            for(auto& par : params){
            TGraph *iGU =   0;
            fU->GetObject(par.c_str(),iGU);
            TGraph *iGD =   0;
            fD->GetObject(par.c_str(),iGD);
            TGraph *iGN =   0;
            fN->GetObject(par.c_str(),iGN);
            if(iGU ==0||iGD==0||iGN==0)continue;

            TGraph * gu = new TGraph;
            TGraph * gd = new TGraph;
            TGraph * ga = new TGraph;
            int nP = 0;
            std::vector<float> vals;
            for(auto& sM : signalMassBins){
                if(sM < 1000) continue;
                if(sM > 3500) continue;
                float nU = iGU->Eval(sM);
                float nN = iGN->Eval(sM);
                float nD = iGD->Eval(sM);
                gu->SetPoint(nP,float(sM),nU/nN-1.0);
                gd->SetPoint(nP,float(sM),nD/nN-1.0);
                ga->SetPoint(nP,float(sM),(nU-nD)/(2*nN));
                vals.push_back(nU/nN-1.0);
                vals.push_back(nD/nN-1.0);
                vals.push_back((nU-nD)/(2*nN));
                nP++;
            }
            Plotter * p = new Plotter;
            p->addGraph(gu, (sys+" "+par+": Up").c_str());
            p->addGraph(gd, (sys+" "+par+": Down").c_str());
            auto gr = p->addGraph(ga, (sys+" "+par+": Avg").c_str());
            p->setYTitle(sel);
            std::sort(vals.begin(),vals.end(),[](const float a, const float b){return a < b;} );
            if(vals.size()){
                float min = vals[0] < 0 ? vals[0]*1.2 : vals[0]*0.8;
                float max = vals.back() < 0 ? vals.back()*0.8 : vals.back()*1.2;
                p->setMinMax( min,max );
            }

            auto c = p->draw(false,sys+"_"+sel+"_"+par);
            std::cout<<std::endl <<"FIT: " <<sys+"_"+sel+"_"+par<<std::endl;
            gr->Fit("pol0");
            plots.push_back(c);
        }
        }
        gROOT->SetBatch(false);
        auto c = Drawing::drawAll(plots,sys+"Par");
        c->SetName((sys+"Par").c_str());
        c->SetTitle((sys+"Par").c_str());
        writeables.push_back(c);
    }
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
    std::string name = signals[sig];
    if(outName.size()) outName += doCond ? std::string("/") + name : std::string("/NoCond")+name;




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
        test2DFits(name,filename,signalMassBins[sig],"2D_fit",false,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"});
        //run below for an
//        test2DFits(name,filename,{1600},"2D_fit",false,{"mu_M_LP_full"}); outName +=+"_forAN";
        break;
    case 5:
        if(outName.size()) outName += "_normSyst";
        plotPDFSysts(name,filename,signalMassBins[sig]);
        plotNormSysts(name,filename,signalMassBins[sig],{"muID","muISO","elReco","elID","elISO","b_real","b_fake","pu"});
        plotShapeSystNorms(name,filename,signalMassBins[sig],{"MET","JER","JES"});
        plotShapeSystParams(name,filename,signalMassBins[sig],{"MET","JER","JES"});
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
