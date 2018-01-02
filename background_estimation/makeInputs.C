
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string>
#include <TSystem.h>
#include "make2DDetectorParam.C"
#include "make1DTemplateWithScaledKernels.C"
#include "make1DTemplateWithAdaKern.C"
#include "make2DTemplateWithAdaKern.C"
#include "fit2DTemplate.C"
#include "mergeHistosToPDF2D.C"
#include "makePlots.C"
#include "CutConstants.h"
#include <thread>
#include <chrono>
#include <future>


using namespace CutConstants;



void makeBackgroundShapesMJJScaledKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    TString resCut = TString::Format("%s&&%s",bV.c_str(),hbbBC.c_str());
    if(addCut.size()){
        resCut +="&&";
        resCut += addCut;
    }
        make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass/hbbGenMass -xb 500,0,5 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
        make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass/genhhMass -xb 500,0,5 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
        gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
        gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));

    std::string tempFile=filename+"_"+name+"_template.root";

//    TString baseSel = TString::Format("%s&&%s&&%s",bV.c_str(),wjjBC.c_str(),exA.c_str());
    TString baseSel = TString::Format("%s&&%s",bV.c_str(),wjjBC.c_str());
    if(addCut.size()){
        baseSel +="&&";
        baseSel += addCut;
    }
    for(unsigned int iP = 0; iP < purs.size() && iP < 1; ++iP){
        for(unsigned int iL = 0; iL < leps.size() && iL < 1; ++iL){
            TString sel =TString::Format("%s&&%s&&%s",baseSel.Data(),pursS[iP].c_str(),lepsS[iL].c_str());
            std::string args = TString::Format("-v -n mjj_ -x hbbMass -g hbbGenMass -xb 125,0,250 -s %s -w %s -ks 0.66 -ss 0.3 -rs 0.5 -ken -sa .1 -ra .1 -vsf %s -vsh scalexHisto -vsv hbbGenPT",sel.Data(), nomW.c_str(), resFile.c_str()).Data();

            unsigned int nP = 0;
            auto runOne = [&](std::string inputFile, std::string outputFile, std::string args){
                nP++;
                if(fork() == 0){
                    std::system(TString::Format("root -b -q 'background_estimation/make1DTemplateWithScaledKernels.C+(\"%s\",\"%s\",\"%s\")'",inputFile.c_str(),outputFile.c_str(),args.c_str()));
                    exit(0);
                }
            };

            runOne(inputFile, std::string("temp/nom_mjj_")+tempFile, args +" -t no");
            runOne(inputFile, std::string("temp/rdown_mjj_")+tempFile, args +" -t r" );
            runOne(inputFile, std::string("temp/rup_mjj_")+tempFile, args +" -t R" );
            runOne(inputFile, std::string("temp/sdown_mjj_")+tempFile, args +" -t s" );
            runOne(inputFile, std::string("temp/sup_mjj_")+tempFile, args +" -t S" );
            int status;
            for(unsigned int i = 0; i <nP; ++i) wait(&status);
            std::cout << "DONE!"<<std::endl;

            gSystem->Exec(TString::Format("hadd -f %s %s",tempFile.c_str(),(std::string("temp/*_mjj_")+tempFile).c_str()));
            gSystem->Exec(TString::Format("rm %s",(std::string("temp/*_mjj_")+tempFile).c_str()));
        }
    }

}

void makeBackgroundShapesMJJAdaKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    TString resCut = TString::Format("%s",hbbBC.c_str());
    if(addCut.size()){
        resCut +="&&";
        resCut += addCut;
    }
        make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass/hbbGenMass -xb 500,0,5 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
        make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass/genhhMass -xb 500,0,5 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
        gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
        gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));

    std::string tempFile=filename+"_"+name+"_template.root";

    TString baseSel = addCut;
    for(unsigned int iP = 0; iP < purs.size() && iP < 1; ++iP){
        for(unsigned int iL = 0; iL < leps.size() && iL < 1; ++iL){
            TString sel =TString::Format("%s&&%s&&%s",baseSel.Data(),pursS[iP].c_str(),lepsS[iL].c_str());
            std::string args = TString::Format("-v -n histo -x hbbMass -g hbbGenMass -xb 125,0,250 -s %s -w %s -khs 1.0 -kss -ks 1.5 -kr 1.5 -hs 0.00714 -hr 45 -vsf %s -vsh scalexHisto -vsv hbbGenPT ",sel.Data(), nomW.c_str(), resFile.c_str()).Data();
            make1DTemplateWithAdaKern(inputFile,tempFile, args);
        }
    }

}

void makeBackgroundShapesMVVConditional(const std::string name, const std::string filename,  const std::string inputFile, const std::string addCut="", float khxs = 1,float khxc = 5,float khys = 1,float khyc = 5) {
    TString baseSel = addCut;
    for(unsigned int iP = 0; iP < purs.size() && iP < 1; ++iP){
        for(unsigned int iL = 0; iL < leps.size() && iL < 1; ++iL){
            std::string tempFile=filename+"_"+name+"_COND2D_template.root";
            TString sel =TString::Format("%s&&%s&&%s",baseSel.Data(),pursS[iP].c_str(),lepsS[iL].c_str());
            float eMin = (name.find("T") != std::string::npos) ? 1500 : 2000;
            float eMax = (name.find("T") != std::string::npos) ? 2500 : 4500;
            std::string args = TString::Format("-v -n histo -vx hhMass -vy hbbMass -xb 200,0,5000 -yb 125,0,250 -ycb 50,60,80,100,120,140,160,180,200,220,250 -s %s -w %s -khxs %f -khxc %f -khys %f -khyc %f -kss  -hs 0.0003 -hr 1200 -emin %.0f -emax %.0f  ",
                    sel.Data(), nomW.c_str(), khxs,khxc,khys,khyc, eMin, eMax ).Data();
            make2DTemplateWithAdaKern(inputFile,tempFile, args);
        }
    }
}


void mergeBackgroundShapes(const std::string& name, const std::string& filename){
    for(unsigned int iP = 0; iP < purs.size() && iP < 1; ++iP){
        for(unsigned int iL = 0; iL < leps.size() && iL < 1; ++iL){
            std::string inFileX=filename+"_"+name+"_template.root";
            std::string inFileY=filename+"_"+name+"_COND2D_template.root";
            std::string rootFile=filename+"_"+name+"_2D.root";
            std::string args = TString::Format("-v -n histo -inX %s -inY %s -sX Scale:ScaleX,Res:ResX,PT:PTX,OPT:OPTX -sY PT:PTY,OPT:OPTY -xb 90,30,210 -yb 168,800,5000",
                    inFileX.c_str(),inFileY.c_str()).Data();
            mergeHistosToPDF2D(rootFile, args);
        }
    }

}

void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
    std::vector<PlotVar> vars;
       vars.emplace_back("hbb_mass",";Hbb soft-drop mass [GeV]","hbbMass",90,30,210,"hh_mass",";HH mass [GeV]","hhMass",168,800,5000 );
       std::vector<PlotSamp> samps = { {name,"1.0"}};
       std::string baseSel = TString::Format("hhMass>800&&hhMass<5000&&hbbMass>30&&hbbMass<210").Data();
       if(addCut.size()){
           baseSel +="&&";
           baseSel += addCut;
       }
       std::vector<PlotSel> sels;
       for(unsigned int iL = 0; iL <leps.size(); ++iL ){
           for(unsigned int iB = 0; iB <purs.size(); ++iB ){
               sels.emplace_back(TString::Format("%s_%s_lWW", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("%s && %s && %s && %s", wjjBC.c_str(), bV.c_str(),lepsS[iL].c_str(),pursS[iB].c_str()).Data());
               sels.emplace_back(TString::Format("%s_%s", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("%s && %s && %s && %s && %s", wjjBC.c_str(), bV.c_str(),lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
               sels.emplace_back(TString::Format("%s_%s_lb", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("%s && %s && %s && %s",wjjBC.c_str(),  lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
               sels.emplace_back(TString::Format("%s_%s_lt", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("wjjMass>10 && %s && %s && %s && %s", bV.c_str(),lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
               sels.emplace_back(TString::Format("%s_%s_ltb", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("wjjMass>10 && %s && %s && %s", lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
               sels.emplace_back(TString::Format("%s_%s_ltmb", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                       TString::Format("%s && %s && %s", lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
           }
       }
       std::string outFileName=filename+"_"+name+"_distributions.root";
       MakePlots a(inputFile,outFileName,samps,sels,vars,baseSel,nomW);
}

void fitBackgroundShapesMVVConditional(const std::string& name, const std::string& filename){
       std::string distFileName=filename+"_"+name+"_distributions.root";
       std::string tempFile=filename+"_"+name+"_2D.root";
       std::vector<std::string> hadSel = {"","lWW","lb","ltb","ltmb"};
       for(unsigned int iL = 0; iL <leps.size(); ++iL ){
           for(unsigned int iB = 0; iB <purs.size(); ++iB ){
               for(unsigned int iH = 0; iH <hadSel.size(); ++iH ){
                   std::string hName = TString::Format("%s_%s_%s_hbb_mass_hh_mass", name.c_str(),leps[iL].c_str(),purs[iB].c_str()).Data();
                   std::string outName = TString::Format("%s_%s_%s_explode.root", name.c_str(),leps[iL].c_str(),purs[iB].c_str()).Data();
                   if(iH){
                       hName = TString::Format("%s_%s_%s_%s_hbb_mass_hh_mass", name.c_str(),leps[iL].c_str(),purs[iB].c_str(),hadSel[iH].c_str()).Data();
                       outName = TString::Format("%s_%s_%s_%s_explode.root", name.c_str(),leps[iL].c_str(),purs[iB].c_str(),hadSel[iH].c_str()).Data();

                   }
                   std::string args = TString::Format("-v -fT %s -nT histo -s PTX,OPTX,PTY,OPTY -fH %s -nH %s",
                           tempFile.c_str(), distFileName.c_str(),hName.c_str() ).Data();
                   fit2DTemplate(outName,args);
               }
           }
       }
}


void makeBackgroundShapesMVVScaledKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
    std::string resFile=filename+"_"+name+"_template.root";

    TString baseSel = TString::Format("%s&&%s&&%s&&%s",bV.c_str(),hbbBC.c_str(),exA.c_str(),wjjBC.c_str());
    if(addCut.size()){
        baseSel +="&&";
        baseSel += addCut;
    }
    for(unsigned int iP = 0; iP < purs.size(); ++iP){
        for(unsigned int iL = 0; iL < leps.size(); ++iL){
            TString sel =TString::Format("%s&&%s&&%s",baseSel.Data(),pursS[iP].c_str(),lepsS[iL].c_str());
            make1DTemplateWithScaledKernels(inputFile,std::string("mvv_")+resFile,TString::Format("-v -n mvv_ -x hhMass -g genhhMass -xb 200,0,5000 -s %s -w %s -ken -ss 0.3 -rs 1 -sa .1 -ra .1",sel.Data(), nomW.c_str()).Data()   );
        }
    }

}

void makePlots(const std::string& name, const std::string& filename,  const std::string inputFile){
    std::vector<PlotVar> vars;
    vars.emplace_back("hbb_mass",";Hbb soft-drop mass [GeV]","hbbMass",50,0,250);
    vars.emplace_back("hh_mass",";HH mass [TeV]","hhMass/1000.0",200,0.,5.);
    std::vector<PlotSamp> samps = { {"ttbar","process==2"},{"ttbar_NR","process==2 && hbbWQuark==0"},{"ttbar_R","process==2 && hbbWQuark==1"},
            {"wjets","process==3"},{"nonRes","process!=8 && hbbWQuark==0"},{"nonResWVV","(process==3 || process == 4 || process ==6 ) && hbbWQuark==0"}
    ,{"nonResT","(process==2 || process == 5 || process ==7 ) && hbbWQuark==0"},{"res","hbbWQuark==1"}};
    std::string baseSel = TString::Format("%s && %s", bV.c_str(), wjjBC.c_str()).Data();
    std::vector<PlotSel> sels;
    for(unsigned int iL = 0; iL <leps.size(); ++iL ){
        for(unsigned int iB = 0; iB <purs.size(); ++iB ){
            sels.emplace_back(TString::Format("%s_%s_lWW", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                    TString::Format("%s && %s", lepsS[iL].c_str(),pursS[iB].c_str()).Data());
            sels.emplace_back(TString::Format("%s_%s", leps[iL].c_str(),purs[iB].c_str()).Data() ,
                    TString::Format("%s && %s && %s", lepsS[iL].c_str(),pursS[iB].c_str(),exA.c_str()).Data());
        }
    }
    std::string outFileName=filename+"_"+name+"_plots.root";

    MakePlots a(inputFile,outFileName,samps,sels,vars,baseSel,nomW);
}




void go(std::string treeDir) {
    std::string filename = "hhLnuJJ";
    std::string treeArea = treeDir + "/betrees_bkg.root";
    //do nonResW
//    {
//        std::string name = "nonResW";
//        std::string baseSel = nresS + "&&"+ nresWS + "&&" +exA;
//        std::string selMJJ = baseSel + "&&(hhMass>800&&hhMass<5000)" ;
//        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,selMJJ);
//        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,2,0.75,2);
//        mergeBackgroundShapes(name,filename);
//        makeFittingDistributions(name,filename,treeArea,nresS + "&&"+ nresWS);
//        fitBackgroundShapesMVVConditional(name,filename);

        //TestQCD
//        makeFittingDistributions(name,filename,treeArea,nresS + "&&"+ nresWQS);
//        fitBackgroundShapesMVVConditional(name,filename);

        //Test optimizations
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_2_khys_1_khyc_2"   ,filename,treeArea,baseSel,1,2,1,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_1p5_khys_1_khyc_1p5"   ,filename,treeArea,baseSel,1,1.5,1,1.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_1_khys_1_khyc_1p5"   ,filename,treeArea,baseSel,1,1,1,1.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_2_khys_1_khyc_2" ,filename,treeArea,baseSel,0.75,2,1,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_2_khys_0p75_khyc_2"   ,filename,treeArea,baseSel,0.75,2,0.75,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1_khys_0p75_khyc_1"   ,filename,treeArea,baseSel,0.75,1,0.75,1);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1_khys_1_khyc_1"   ,filename,treeArea,baseSel,0.75,1,1,1);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1_khys_0p75_khyc_2"   ,filename,treeArea,baseSel,0.75,1,0.75,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_0p5_khys_1_khyc_2" ,filename,treeArea,baseSel,1,0.5,1,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_2_khys_0p75_khyc_2",filename,treeArea,baseSel,1,2,0.75,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_2_khys_0p5_khyc_2" ,filename,treeArea,baseSel,1,2,0.5,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_2_khys_1_khyc_2",filename,treeArea,baseSel,0.75,2,1,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p5_khxc_2_khys_0p5_khyc_2" ,filename,treeArea,baseSel,0.5,2,0.5,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p25_khxc_2_khys_0p25_khyc_2" ,filename,treeArea,baseSel,0.25,2,0.25,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p5_khxc_1_khys_0p5_khyc_1" ,filename,treeArea,baseSel,0.5,1,0.5,1);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p25_khxc_2_khys_0p5_khyc_2" ,filename,treeArea,baseSel,0.25,2,0.5,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p5_khxc_2_khys_0p25_khyc_2" ,filename,treeArea,baseSel,0.5,2,0.25,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_2_khys_0p75_khyc_2" ,filename,treeArea,baseSel,0.75,2,0.75,2);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1_khys_0p75_khyc_1" ,filename,treeArea,baseSel,0.75,1,0.75,1);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_0p75_khys_0p75_khyc_0p75" ,filename,treeArea,baseSel,0.75,0.75,0.75,0.75);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_0p5_khys_0p75_khyc_0p5" ,filename,treeArea,baseSel,0.75,0.5,0.75,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p5_khxc_0p5_khys_0p5_khyc_0p5" ,filename,treeArea,baseSel,0.5,0.5,0.5,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p5_khxc_0p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.5,0.5,0.25,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p25_khxc_0p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.25,0.5,0.25,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p25_khxc_0p5_khys_0p5_khyc_0p5" ,filename,treeArea,baseSel,0.25,0.5,0.5,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_0p5_khys_0p5_khyc_0p5" ,filename,treeArea,baseSel,1,0.5,0.5,0.5);
//        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_0p5_khys_0p5_khyc_0p5" ,filename,treeArea,baseSel,0.75,0.5,0.5,0.5);
//                makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_0p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,1,0.5,0.25,0.5);
//                makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_0p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.75,0.5,0.25,0.5);
//                        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_0p75_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.75,0.75,0.25,0.5);
//                        makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.75,1,0.25,0.5);
//                                makeBackgroundShapesMVVConditional(name+"_khxs_0p75_khxc_1p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,0.75,1.5,0.25,0.5);
//                                makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_1p5_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,1,1.5,0.25,0.5);
//                                makeBackgroundShapesMVVConditional(name+"_khxs_1_khxc_2_khys_0p25_khyc_0p5" ,filename,treeArea,baseSel,1,2,0.25,0.5);
//    }
//    //do nonResT
//    {
//        std::string name = "nonResT";
//        std::string baseSel = nresS + "&&"+ nresTS +"&&" +exA;
//        std::string selMJJ = baseSel + "&&(hhMass>800&&hhMass<5000)" ;
//        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel);
//        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,2,0.75,2);
//        mergeBackgroundShapes(name,filename);
//        makeFittingDistributions(name,filename,treeArea,nresS + "&&"+ nresTS);
//        fitBackgroundShapesMVVConditional(name,filename);
//    }

    //do res
    {
        std::string name = "res";
        std::string baseSel = resS + "&&" +exA;
        std::string selMJJ = baseSel + "&&(hhMass>800&&hhMass<5000)" ;
//        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel);
//        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,2,0.75,2);
//        mergeBackgroundShapes(name,filename);
//        makeFittingDistributions(name,filename,treeArea,resS);
        makeFittingDistributions(name+"W",filename,treeArea,resS+"&&(hbbWQuark<3)");
        makeFittingDistributions(name+"T",filename,treeArea,resS+"&&(hbbWQuark>2)");
//        fitBackgroundShapesMVVConditional(name,filename);
    }
}
#endif

void makeInputs(std::string treeDir = "trees/"){
    go(treeDir);
}
