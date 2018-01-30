
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <string>
#include <TSystem.h>
#include "make2DDetectorParam.C"
#include "make1DTemplateWithScaledKernels.C"
#include "make1DTemplateWithAdaKern.C"
#include "make2DTemplateWithAdaKern.C"
#include "fit2DTemplate.C"
#include "fit1DTemplate.C"
#include "makeJSON.C"
#include "mergeHistosToPDF2D.C"
#include "makePlots.C"
#include "CutConstants.h"
#include "cutHistos1D.C"
#include "FunctionFitter.C"
#include <thread>
#include <chrono>
#include <future>


using namespace CutConstants;
using namespace ASTypes;



//void makeBackgroundShapesMJJScaledKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& addCut=""){
//    std::string resFile=filename+"_"+name+"_detectorResponse.root";
//    TString resCut = TString::Format("%s&&%s",bV.c_str(),hbbBC.c_str());
//    if(addCut.size()){
//        resCut +="&&";
//        resCut += addCut;
//    }
//        make2DDetectorParam(inputFile,std::string("x_")+resFile,TString::Format("-v -n xHisto -x hbbMass/hbbGenMass -xb 500,0,5 -s %s  -w %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
//        make2DDetectorParam(inputFile,std::string("y_")+resFile,TString::Format("-v -n yHisto -x hhMass/genhhMass -xb 500,0,5 -s %s -w  %s -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000",resCut.Data(),nomW.c_str()).Data());
//        gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
//        gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));
//
//    std::string tempFile=filename+"_"+name+"_template.root";
//
////    TString baseSel = TString::Format("%s&&%s&&%s",bV.c_str(),wjjBC.c_str(),exA.c_str());
//    TString baseSel = TString::Format("%s&&%s",bV.c_str(),wjjBC.c_str());
//    if(addCut.size()){
//        baseSel +="&&";
//        baseSel += addCut;
//    }
//    for(unsigned int iP = 0; iP < purs.size() && iP < 1; ++iP){
//        for(unsigned int iL = 0; iL < leps.size() && iL < 1; ++iL){
//            TString sel =TString::Format("%s&&%s&&%s",baseSel.Data(),pursS[iP].c_str(),lepsS[iL].c_str());
//            std::string args = TString::Format("-v -n mjj_ -x hbbMass -g hbbGenMass -xb 125,0,250 -s %s -w %s -ks 0.66 -ss 0.3 -rs 0.5 -ken -sa .1 -ra .1 -vsf %s -vsh scalexHisto -vsv hbbGenPT",sel.Data(), nomW.c_str(), resFile.c_str()).Data();
//
//            unsigned int nP = 0;
//            auto runOne = [&](std::string inputFile, std::string outputFile, std::string args){
//                nP++;
//                if(fork() == 0){
//                    std::system(TString::Format("root -b -q 'background_estimation/make1DTemplateWithScaledKernels.C+(\"%s\",\"%s\",\"%s\")'",inputFile.c_str(),outputFile.c_str(),args.c_str()));
//                    exit(0);
//                }
//            };
//
//            runOne(inputFile, std::string("temp/nom_mjj_")+tempFile, args +" -t no");
//            runOne(inputFile, std::string("temp/rdown_mjj_")+tempFile, args +" -t r" );
//            runOne(inputFile, std::string("temp/rup_mjj_")+tempFile, args +" -t R" );
//            runOne(inputFile, std::string("temp/sdown_mjj_")+tempFile, args +" -t s" );
//            runOne(inputFile, std::string("temp/sup_mjj_")+tempFile, args +" -t S" );
//            int status;
//            for(unsigned int i = 0; i <nP; ++i) wait(&status);
//            std::cout << "DONE!"<<std::endl;
//
//            gSystem->Exec(TString::Format("hadd -f %s %s",tempFile.c_str(),(std::string("temp/*_mjj_")+tempFile).c_str()));
//            gSystem->Exec(TString::Format("rm %s",(std::string("temp/*_mjj_")+tempFile).c_str()));
//        }
//    }
//
//}

std::unique_ptr<TH1> proj(const TH2* inH, const std::string& newName, double min, double max, bool projX, bool noBounds){
    const TAxis * ax = projX ? inH->GetYaxis() : inH->GetXaxis();
    int binL = ax->FindFixBin(min);
    int binH = ax->FindFixBin(max) -1;
    if(noBounds){
        binL = 0;
        binH = -1;
    }
    std::cout << min <<" "<< max <<" "<<binL<<" "<<binH<<std::endl;
    TH1 * hp = projX ? inH->ProjectionX(newName.c_str(),binL,binH) : inH->ProjectionY(newName.c_str(),binL,binH);
    hp->SetDirectory(0);
    return std::unique_ptr<TH1>(hp);
}
std::unique_ptr<TH1> projX(const TH2* inH, const std::string& newName, double minY, double maxY) {return proj(inH,newName,minY,maxY,true,false);}
std::unique_ptr<TH1> projY(const TH2* inH, const std::string& newName, double minX, double maxX) {return proj(inH,newName,minX,maxX,false,false);}
std::unique_ptr<TH1> projX(const TH2* inH, const std::string& newName) {return proj(inH,newName,0,0,true,true);}
std::unique_ptr<TH1> projY(const TH2* inH, const std::string& newName) {return proj(inH,newName,0,0,false,true);}

double getMean(const TH1* h, double minX = -1, double maxX = -1){
    double mean = 0;
    double meanT = 0;
    for(int iX=0; iX <= h->GetNbinsX();++iX){
        if(h->GetBinLowEdge(iX) < minX) continue;
        if(h->GetBinLowEdge(iX) >= maxX ) continue;
        mean += h->GetBinCenter(iX) * h->GetBinContent(iX);
        meanT +=h->GetBinContent(iX);
    }
    if(meanT) mean /= meanT;
    return mean;
}

double getNT(double pt){
    double L=5;
    double R=2;
    double P = 1250;
    double D = 250;
    double S = -1;
    if(pt <= P-D) return L;
    if(pt >= P+D) return R;
    return 0.5*(L+R + (L-R)*std::sin(S*TMath::PiOver2()*(pt-P)/D ));
}


void makeDetectorParam(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0"){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    std::string xArgs = std::string("-v -n xHisto -x ")+ hbbMCS.cut +"/hbbGenMass "+ " -xb 500,0,5 "+" -s "+cut +" -w "+nomW.cut+ " -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000";
    std::string yArgs = std::string("-v -n yHisto -x ")+ hhMCS.cut +"/genhhMass "+ " -xb 500,0,5 "+" -s "+cut +" -w "+nomW.cut+ " -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000";
    make2DDetectorParam(inputFile,std::string("x_")+resFile,xArgs);
    make2DDetectorParam(inputFile,std::string("y_")+resFile,yArgs);
    gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
    gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));
}

void makeBackgroundShapesMJJAdaKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0",float khxs = 1,float khxc = 5){
    std::string tempFile=filename+"_"+name+"_incl_template.root";
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    std::string args = std::string("-v -n histo ")+" -x "+hbbMCS.cut+" -g hbbGenMass -xb "+hbbInclBinning.cut+ " -s "+cut+" -w "+nomW.cut+ " -khs "+ flt2Str(khxs) +" -khc "+ flt2Str(khxc);
    args += " -kss -ks 1.5 -kr 1.5 -hs 0.00714 -hr 45 ";
    args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT ";
    make1DTemplateWithAdaKern(inputFile,tempFile, args);
}
void makeBackgroundShapesMVVAdaKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0",float khxs = 1,float khxc = 5){
    std::string tempFile=filename+"_"+name+"_incl_template.root";
    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    float eMin = (name.find("mt") != std::string::npos) ?  2000:1500;
    float eMax = (name.find("mt") != std::string::npos) ?  4500:2500;

    std::string args = std::string("-v -n histo ")+" -x "+hhMCS.cut+" -g hbbGenMass -xb "+hhInclBinning.cut+" -s "+cut+" -w "+nomW.cut+ " -khs "+ flt2Str(khxs) +" -khc "+ flt2Str(khxc);
    args += " -kss -ks 1.5 -kr 1.5 -hs 0.0003 -hr 1200 ";
    args += std::string(" -doS -emin ") + flt2Str(eMin) + " -emax "+flt2Str(eMax) + " ";
    args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT ";
    make1DTemplateWithAdaKern(inputFile,tempFile, args);
}

void makeBackgroundShapesMVVConditional(const std::string name, const std::string filename,  const std::string inputFile, const std::string cut="1.0", float khxs = 1,float khxc = 5,float khys = 1,float khyc = 5) {
    std::string tempFile=filename+"_"+name+"_incl_COND2D_template.root";
        float eMin = (name.find("lost") != std::string::npos) ? 1500 : 2000;
        float eMax = (name.find("lost") != std::string::npos) ? 2500 : 4500;
    std::string args = std::string("-v -n histo ") + " -vx "+ hhMCS.cut+ " -vy "+hbbMCS.cut+ " -xb "+hhInclBinning.cut+" -yb "+hbbInclBinning.cut+ " -ycb 50,60,80,100,120,140,160,180,200,220,250 "+ "-s "+ cut +" -w "+nomW.cut;
    args+=std::string(" -khxs ")+ flt2Str(khxs) +" -khxc "+ flt2Str(khxc) +" -khys "+ flt2Str(khys) +" -khyc "+ flt2Str(khyc) + " -hss ";
    args+=std::string(" -hs 0.0003 -hr 1200 ") + " -emin "+ flt2Str(eMin)+ " -emax "+ flt2Str(eMax);
    make2DTemplateWithAdaKern(inputFile,tempFile, args);
}


void mergeBackgroundShapes(const std::string& name, const std::string& filename){
    std::string inFileX=filename+"_"+name+"_incl_template.root";
    std::string inFileY=filename+"_"+name+"_incl_COND2D_template.root";
    std::string rootFile=filename+"_"+name+"_2D_template.root";
    std::string args = std::string("-v -n histo ") + " -inX "+  inFileX + " -inY "+ inFileY + " -sX Scale:ScaleX,Res:ResX,PT:PTX,OPT:OPTX -sY PT:PTY,OPT:OPTY "+ " -xb "+hbbBinning.cut+" -yb "+hhBinning.cut;
    mergeHistosToPDF2D(rootFile, args);
}

void cutMVVTemplate(const std::string& name, const std::string& filename){
    std::string inFileX=filename+"_"+name+"_incl_template.root";
    std::string rootFile=filename+"_"+name+"_template.root";
    std::string args = std::string("-v -n histo ") + " -i "+  inFileX + " -s Scale:Scale,Res:Res,PT:PT,OPT:OPT "+ " -xb "+hhBinning.cut;
    cutHistos1D(rootFile, args);
}

void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    std::vector<PlotVar> vars;
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    }
    std::vector<PlotSamp> samps = { {name,"1.0"}};
    std::vector<PlotSel> sels;
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        sels.emplace_back(l +"_"+p +"_"+h,
                l.cut +"&&"+p.cut+"&&"+h.cut);
        std::cout << l +"_"+p +"_"+h <<" -->  "<<l.cut +"&&"+p.cut+"&&"+h.cut<<std::endl;
    }
    std::string outFileName=filename+"_"+name+ (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,nomW.cut);
}

void fitBackgroundShapesMVVConditional(std::string name, const std::string& filename, const std::string& fittedName =""){
    std::string tempFile=filename+"_"+name+"_2D_template.root";
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string hName = name+"_"+l+"_"+p+"_"+h+"_"+hbbMCS+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+l+"_"+p+"_"+h+"_2D_template.root";
        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PTX,OPTX,PTY,OPTY " + " -fH " + distFileName + " -nH "+ hName;
        fit2DTemplate(outName,args);
    }
}
void fitBackgroundShapesMVV(std::string name, const std::string& filename, const std::string& fittedName =""){
    std::string tempFile=filename+"_"+name+"_template.root";
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string hName = name+"_"+l+"_"+p+"_"+h+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+l+"_"+p+"_"+h+"_template.root";
        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PT,OPT " + " -fH " + distFileName + " -nH "+ hName;
        fit1DTemplate(outName,args);
    }
}


template<typename Func>
void makeResMJJShapes(const std::string& name, const std::string& filename,const std::string& jsonArgs, Func doFit) {
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        //        if(l.find("emu") == std::string::npos ) continue;
        //        if(p.find("L") == std::string::npos ) continue;
        //        if(h.find("none") == std::string::npos ) continue;
        std::string hName = name+"_"+l+"_"+p+"_"+h;
        auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
        if(hbb_hh_H==0) continue;
        auto hhH  = projY(&*hbb_hh_H,hName+"_"+hhMCS);
        FunctionParameterPlotter plotter;
        std::vector<std::unique_ptr<FunctionFitter>> fitters;
        for(unsigned int iP = 0; iP < resPTBins.size() -1; ++iP){
            std::string ptName  = flt2Str(resPTBins[iP]) +"to"+flt2Str(resPTBins[iP+1]);
            auto hbbH = projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS,resPTBins[iP],resPTBins[iP+1]);
            double mean = getMean(&*hhH,resPTBins[iP],resPTBins[iP+1] );
            doFit(&*hbbH,fitters);
            plotter.addFit(&*fitters[iP],mean,ptName);
        }
        plotter.write(filename+"_"+hName+"_fit.root");
        std::string argsP1 = std::string("-i ")+ filename+"_"+hName+"_fit.root ";
        MakeJSON(filename+"_"+hName+"_fit.json",argsP1+" "+jsonArgs);
    }

}

void makeResTopMJJShapes(const std::string& name, const std::string& filename){
    auto setupFitter =[](const TH1* hbbH, std::vector<std::unique_ptr<FunctionFitter>>& fitters ){
        fitters.emplace_back(new CBFunctionFitter(hbbH,"T"));
        auto fitter = &* fitters.back();
        fitter->setVar("meanT"     ,180,120,195);
        fitter->setVar("sigmaT"     ,15.8,14,20);
        //             fitter->setConst("sigmaT",1);
        fitter->setVar("alphaT"     ,1 ,0.1,2);
        //            fitter->setConst("alphaT",1);
        fitter->setVar("alpha2T"  ,1.5,0.5,3);
        //             fitter->setConst("alpha2T"  ,1);
        fitter->setVar("nT"   ,  5  ,1,6);
        fitter->setVar("n2T"  ,5,3,6);
        fitter->setConst("n2T"  ,1);
        fitter->setConst("nT"  ,1);
        fitter->w->var("M")->setRange("fit",70,230);
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(0)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(1)});
    };
    //        std::strin args = " -g meanT:laur3,sigmaT:laur2,alphaT:pol0,alpha2T:laur3,nT:FIXp5p2p1250p250p-1,n2T:pol0 ";
    std::string args = " -g meanT:laur3,sigmaT:laur2,alphaT:laur3,alpha2T:laur3,nT:pol0,n2T:pol0 ";
    args += " -minX 500 -maxX 3500 ";
    makeResMJJShapes(name,filename,args, setupFitter);
}

void makeResWMJJShapes(const std::string& name, const std::string& filename){
    auto setupFitter =[](const TH1* hbbH, std::vector<std::unique_ptr<FunctionFitter>>& fitters ){
        fitters.emplace_back(new CBFunctionFitter(hbbH,"W"));
        auto fitter = &* fitters.back();
        fitter->setVar("meanW"     ,90,80,100);
        fitter->setVar("sigmaW"     ,8,5,15);
        //             fitter->setConst("sigmaW",1);
        fitter->setVar("alphaW"     ,1.18,0.1,10);
        //        fitter->setConst("alphaW",1);
        fitter->setVar("alpha2W"  ,0.9,0.1,5);
        //             fitter->setConst("alpha2W"  ,1);
        fitter->setVar("nW"  ,5,1,6);
        fitter->setVar("n2W"  , 1,1,6);
        fitter->setConst("nW"  ,1);
        fitter->setConst("n2W"  ,1);
        fitter->w->var("M")->setRange("fit",30,160);
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(0)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(1)});
    };
    std::string args = " -g meanW:laur2,sigmaW:laur2,alphaW:pol0,alpha2W:laur2,nW:pol0,n2W:pol0 ";
    args += " -minX 500 -maxX 3000 ";
    makeResMJJShapes(name,filename,args,setupFitter);
}

void convertFuncFitTo2DTemplate(const std::string& name, const std::string& filename,const std::string& funcParamPostfix){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string jsonFile = filename+"_"+name+"_"+lepSels[LEP_EMU]+"_"+ (p==purSels[PUR_T] && name.find("w") != std::string::npos ? purSels[PUR_M] : p )
                                +"_"+hadSels[HAD_NONE] +"_fit.json";
        CBFunctionFitter xFit(0,funcParamPostfix);
        xFit.loadJSON(jsonFile,"MH",nHbbMassBins*10,minHbbMass,maxHbbMass);

        auto * iF =  TObjectHelper::getFile(filename + "_"+name+"_"+l+"_"+p+"_"+h+"_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH1>(iF,"histo",false,false);
        if(hh_H==0) continue;
        oF->cd();
        TH2 * h_2D = new TH2F((name+"_"+l+"_"+p+"_"+h).c_str(),(std::string(";")+hbbMCS.title+";"+hhMCS.title).c_str(),nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
        for(int iY = 1; iY <= hh_H->GetNbinsX(); ++iY){
            double yNorm = hh_H->GetBinContent(iY);
            double hhV = hh_H->GetBinCenter(iY);
            auto xHist = xFit.getHistFromJSON(name+"_"+l+"_"+p+"_"+h+ASTypes::int2Str(iY)+"_hbb",hhV);
            for(int iX = 1; iX <= h_2D->GetNbinsX(); ++iX){
                double xNorm = xHist->Integral( 10*iX -9,10*iX  );
                h_2D->SetBinContent(iX,iY,xNorm*yNorm);
            }
            delete xHist;
        }
        h_2D->Write();
        iF->Close();
    }
    oF->Close();
}

void compile2DTemplatesForDebug(const std::string& name, const std::string& filename){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_"+l+"_" +p +"_"+h+"_2D_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH2>(iF,"histo__x_y",false,false);
        if(hh_H==0) continue;
        oF->cd();
        hh_H->SetXTitle(hbbMCS.title.c_str());
        hh_H->SetYTitle(hhMCS.title.c_str());
        hh_H->Write((name+"_"+l+"_"+p+"_"+h).c_str());
        iF->Close();
    }
    oF->Close();
}

void go(std::string treeDir) {
    std::string filename = hhFilename;
    std::string treeArea = treeDir + "/betrees_bkg.root";
    //do nonResW
    {
        std::string name = bkgSels[BKG_QG];
        std::string genSel = bkgSels[BKG_QG].cut + "&&"+ aQCD.cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+ hadSels[HAD_LTMB].cut;
//                makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        //MVV

//                        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,3,0.5,1);//x = hh
        //                makeBackgroundShapesMVVConditional(name+"_xs_0p75_xc_2_ys_0p75_yc_2,filename,treeArea,baseSel,0.75,2,0.75,2);//old
//                makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hhInclRange.cut,true);

        //MJJ
//        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel+"&&"+hhRange.cut);
        mergeBackgroundShapes(name,filename);
//        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVVConditional(name,filename);

        compile2DTemplatesForDebug(name,filename);makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
    }
    //    //do lost t/w
    if(false){
        std::string name = bkgSels[BKG_LOSTTW];
        std::string genSel = bkgSels[BKG_LOSTTW].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LTMB].cut;
        makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        //MVV
        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,8,0.5,2);//x = hh
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hhInclRange.cut,true);
        //MJJ
        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel+"&&"+hhRange.cut);
        mergeBackgroundShapes(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVVConditional(name,filename);

        compile2DTemplatesForDebug(name,filename);
    }

    //    //do mt
    if(false){
        std::string name = bkgSels[BKG_MT];
        std::string genSel = bkgSels[BKG_MT].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LB].cut;
        //MVV
        makeDetectorParam(name,filename,treeArea,genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,baseSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVV(name,filename);

        //MJJ
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeResTopMJJShapes(name,filename);
        convertFuncFitTo2DTemplate(name,filename,"T");
    }
    //    //do mW
    if(false){
        std::string name = bkgSels[BKG_MW];
        std::string genSel = bkgSels[BKG_MW].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LB].cut;
        //MVV
        makeDetectorParam(name,filename,treeArea,genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,baseSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVV(name,filename);

        //MJJ
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeResWMJJShapes(name,filename);
        convertFuncFitTo2DTemplate(name,filename,"W");
    }

    //Signal
    if(false){
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);

    }

}
#endif

void makeInputs(std::string treeDir = "trees/"){
    go(treeDir);
}
