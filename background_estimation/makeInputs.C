
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string>
#include <TSystem.h>
#include "make2DDetectorParam.C"
#include "make1DTemplateWithScaledKernels.C"
#include "make1DTemplateWithAdaKern.C"
#include "make2DTemplateWithAdaKern.C"
#include "makePlots.C"
#include "CutConstants.h"
#include <thread>
#include <chrono>
#include <future>

using namespace CutConstants;

void makeBackgroundShapesMJJKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    TString resCut = TString::Format("%s&&%s",bV.c_str(),hbbBC.c_str());
    if(addCut.size()){
        resCut +="&&";
        resCut += addCut;
    }
    make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass -xb 500,0,2000 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
    make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass -xb 500,0,80000 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());

    //    make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass/hbbGenMass -xb 500,0,5 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.Data()).Data());
    //    make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass/genhhMass -xb 500,0,5 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.Data()).Data());
    gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
    gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));

    //    for(unsigned int iP = 0; iP < purs.size(); ++iP){
    //        for(unsigned int iL = 0; iL < leps.size(); ++iL){
    //
    //        }
    //    }

}

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
    TString resCut = TString::Format("%s&&%s",bV.c_str(),hbbBC.c_str());
    if(addCut.size()){
        resCut +="&&";
        resCut += addCut;
    }
//        make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass/hbbGenMass -xb 500,0,5 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
//        make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass/genhhMass -xb 500,0,5 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
//        gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
//        gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));

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
//            std::string args = TString::Format("-v -n mjj_ -x hbbMass -g hbbGenMass -xb 125,0,250 -s %s -w %s -khs 1.0 -kss -ks 1.5 -kr 1.5 -hs 1.5 -hr 1.5 -vsf %s -vsh scalexHisto -vsv hbbGenPT ",sel.Data(), nomW.c_str(), resFile.c_str()).Data();
//            make1DTemplateWithAdaKern(inputFile,tempFile, args);

            std::string args = TString::Format("-v -n mvvcond_ -vx hhMass -vy hbbMass -xb 200,0,5000 -yb 125,0,250 -ycb 50,60,80,100,120,140,160,180,200,220,250 -s %s -w %s -khxs 0.50 -khys 0.75 -kss -hs 1.5 -hr 1.5  ",sel.Data(), nomW.c_str()).Data();
            make2DTemplateWithAdaKern(inputFile,tempFile, args);
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
    vars.emplace_back("hh_mass",";HH mass [TeV]","hbbMass/1000.0",200,0.,5.);
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
//        makeBackgroundShapesMJJKernel("nonResW","hhLNuJJ",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==3||process==4||process==6)",nresS.c_str()).Data());
//        makeBackgroundShapesMJJKernel("nonResT","hhLNuJJ",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==2||process==5||process==7)",nresS.c_str()).Data());

//            makeBackgroundShapesMJJScaledKernel("nonResW","hhLNuJJ_0p66_m",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==3||process==4||process==6)",nresS.c_str()).Data());
//            makeBackgroundShapesMJJScaledKernel("nonResT","hhLNuJJ_0p66_m",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==2||process==5||process==7)",nresS.c_str()).Data());

        makeBackgroundShapesMJJAdaKernel("nonResT","hhLNuJJ_ada_0p5_0p75",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==2||process==5||process==7)",nresS.c_str()).Data());
        makeBackgroundShapesMJJAdaKernel("nonResW","hhLNuJJ_ada_0p5_0p75",treeDir + "/betrees_bkg.root",TString::Format("%s&&(process==3||process==4||process==6)",nresS.c_str()).Data());


//    makeBackgroundShapesMJJScaledKernel("nonRes","hhLNuJJ",treeDir + "/betrees_bkg.root",TString::Format("%s&&%s",aQCD.c_str(),nresS.c_str()).Data());
    //    makeBackgroundShapesMVVScaledKernel("nonRes","hhLNuJJ",treeDir + "/betrees_bkg.root",TString::Format("%s&&%s",aQCD.Data(),nresS.Data()).Data());
//        makePlots("mc","hhLNuJJ",treeDir + "/betrees_bkg.root");
}
#endif

void makeInputs(std::string treeDir = "trees/"){
    go(treeDir);
}
