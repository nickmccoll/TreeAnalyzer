
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
TH2* cutHistogram(const TH2* inH, std::string name, double nMinX, double nMaxX, double nMinY, double nMaxY ){
    double binWX = inH->GetXaxis()->GetBinWidth(1);
    double binWY = inH->GetYaxis()->GetBinWidth(1);
    auto checkBin =[&](double val, double width){
      int binN = val/width;
      if(binN*width != val){
          std::cout <<  name <<" :: "<< nMinX <<" "<<nMaxX <<" "<<nMinY <<" "<<nMaxY <<" "<< val<<" "<<width <<" "<< binN <<" "<< binN*width<<std::endl;
          throw std::invalid_argument("cutHistogram::cutHistogram() -> Bad parsing");
      }
    };
    checkBin(nMinX,binWX);
    checkBin(nMaxX,binWX);
    checkBin(nMinY,binWY);
    checkBin(nMaxY,binWY);
  TH2 * outH = new TH2F(name.c_str(),TString(";")+inH->GetXaxis()->GetTitle()+";"+ inH->GetYaxis()->GetTitle() ,
          (nMaxX-nMinX)/binWX,nMinX,nMaxX,(nMaxY-nMinY)/binWY,nMinY,nMaxY);

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

void makeSignalFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    int nameIDX =  inputFile.find("XXX", 0);
    std::vector<PlotVar> vars;
    std::vector<PlotSel> sels;
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        sels.emplace_back(l +"_"+p +"_"+h,
                l.cut +"&&"+p.cut+"&&"+h.cut);
    }
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    }

    for(const auto& sM : signalMassBins){
        std::vector<PlotSamp> samps = { {name +"_m"+ASTypes::int2Str(sM),"1.0"}};
        std::string outFileName=filename+"_"+name+"_m"+ASTypes::int2Str(sM) + (doIncl ? "_inclM_distributions.root" : "_distributions.root");
        std::string inputName = inputFile;
        inputName.replace(nameIDX,3,ASTypes::int2Str(sM));
        MakePlots a(inputName,outFileName,samps,sels,vars,cut,nomW.cut);
    }
    std::string compiledFile =  filename+"_"+name + (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    std::string allFiles = filename+"_"+name+"_m*" + (doIncl ? "_inclM_distributions.root" : "_distributions.root");

    gSystem->Exec((std::string("hadd -f ")+ compiledFile + " " + allFiles).c_str());
    gSystem->Exec((std::string("rm ") + allFiles).c_str());
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

void makeSignalMJJShapes(const std::string& name, const std::string& filename){
    auto setup1DFit = [](const TH1* hbbH, bool doExpo, std::vector<std::unique_ptr<FunctionFitter>>& fitters ){
        std::string pF = "SMJJ";
        auto vN=[&](std::string var)->std::string{return var+pF;};
        fitters.emplace_back(new CBFunctionFitter(hbbH,doExpo,pF,{"MJJ"}));
        auto fitter = &* fitters.back();
        fitter->setVar(vN("mean")     ,125,90,180);
        fitter->setVar(vN("sigma")       ,10,5,20);
        //             fitter->setConst(vN("sigma"),1);
        fitter->setVar(vN("alpha")     ,1 ,0.1,2);
        //            fitter->setConst(vN("alpha"),1);
        fitter->setVar(vN("alpha2")  ,1.5,0.1,3);
        //             fitter->setConst(vN("alpha2")  ,1);
        fitter->setVar(vN("n")   ,  5  ,1,6);
        fitter->setVar(vN("n2")  ,5,3,6);
        fitter->setConst(vN("n")  ,1);
        fitter->setConst(vN("n2")  ,1);

        if(doExpo){
            fitter->setVar(vN("slope")  ,-1,-10,0);
            fitter->setVar(vN("fE")  ,0.1,0,0.75);
        }

        fitter->w->var("MJJ")->setRange("fit",30,210);
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("fit"),RooFit::Minos(0),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("fit"),RooFit::Minos(1),RooFit::NumCPU(8)});
    };

    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I] || p == purSels[PUR_LMT] ) continue;
        if(h != hadSels[HAD_LTMB] ) continue;
        bool doExpo = p == purSels[PUR_L];

        FunctionParameterPlotter plotter;
        std::vector<std::unique_ptr<FunctionFitter>> fitters;
        std::string catName = l+"_"+p+"_"+h;
        int iP = 0;
        for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
            std::cout<<" !!!!!!! "<< catName <<" -> "<<signalMassBins[iS]<<std::endl;
            std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
            std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
            auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
            if(hbb_hh_H ==0) continue;
            auto hbb_H = projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS);
            setup1DFit(&*hbb_H,doExpo,fitters);
            plotter.addFit(&*fitters.back(),signalMassBins[iS],ptName);
            iP++;
        }
        plotter.write(filename+"_"+name+"_"+catName+"_MJJ_fit.root");

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_MJJ_fit.root "+ " -minX 550 -maxX 4550";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0 ";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_MJJ_fit.json",argsP1+" "+  ( doExpo  ?  jsonArgsExpo : jsonArgsStd ));
        }
}


void makeSignal2DShapes(const std::string& name, const std::string& filename, std::string& catName, std::string& fitName, CJSON& prevJSON, TFile* iF, bool fitAllP, bool doExpo){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;

    auto setup2DFit = [&](const TH2* hbbH, const double HHMass){
        std::string pF = "S";
        auto vN=[&](std::string var)->std::string{return var+pF;};
        auto pn0 =[&] (std::string v) ->std::string{return v+pF+"MJJ";};
        auto pn1 =[&] (std::string v) ->std::string{return v+pF+"MVV";};

        fitters.emplace_back(new CBFunctionFitter2D(hbbH,doExpo,pF,{"MJJ","MVV"}));
        auto fitter = &* fitters.back();
        fitter->setVar(pn0("mean")  ,prevJSON.evalFunc(pn0("mean")  ,HHMass),90,180);
        fitter->setVar(pn0("sigma") ,prevJSON.evalFunc(pn0("sigma") ,HHMass),5,20);
        fitter->setVar(pn0("alpha") ,prevJSON.evalFunc(pn0("alpha") ,HHMass),0.1,3);
        fitter->setVar(pn0("alpha2"),prevJSON.evalFunc(pn0("alpha2"),HHMass),0.1,3);
        fitter->setVar(pn0("n")     ,prevJSON.evalFunc(pn0("n")     ,HHMass),1,6);
        fitter->setVar(pn0("n2")    ,prevJSON.evalFunc(pn0("n2")    ,HHMass),3,6);
        fitter->setConst(pn0("mean")  ,1);
        fitter->setConst(pn0("sigma") ,1);
        fitter->setConst(pn0("alpha") ,1);
        fitter->setConst(pn0("alpha2"),1);
        fitter->setConst(pn0("n")     ,1);
        fitter->setConst(pn0("n2")    ,1);
        if(doExpo){
            fitter->setVar(pn0("slope") ,prevJSON.evalFunc(pn0("slope")    ,HHMass),-10,0);
            fitter->setVar(pn0("fE")    ,prevJSON.evalFunc(pn0("fE")    ,HHMass),0,0.75);
            fitter->setConst(pn0("slope")     ,1);
            fitter->setConst(pn0("fE")    ,1);
        }

        fitter->setVar(pn1("mean")  ,HHMass   ,HHMass -200,HHMass+200);
        fitter->setVar(pn1("sigma") ,HHMass*0.05,HHMass*0.025,HHMass*0.1);
        fitter->setVar(pn1("alpha") ,1.5,0.75,3);
        fitter->setVar(pn1("alpha2"),1.5,0.75,3);
        fitter->setVar(pn1("n")     ,5,1,6);
        fitter->setVar(pn1("n2")    ,5,1,6);
        fitter->setVar(pn1("maxS")  ,2.5 ,0,5);
        fitter->setVar(pn1("mean_p1")  ,.044,.01,.10);
        fitter->setVar(pn1("sigma_p1") ,0,0,5);

        fitter->setConst(pn1("n")     ,1);
        fitter->setConst(pn1("n2")    ,1);
        fitter->setConst(pn1("maxS")  ,1);

        if(fitAllP){
            fitter->setVar(pn1("mean")    ,prevJSON.evalFunc(pn1("mean")    ,HHMass));
            fitter->setVar(pn1("sigma")   ,prevJSON.evalFunc(pn1("sigma")   ,HHMass));
            fitter->setVar(pn1("alpha")   ,prevJSON.evalFunc(pn1("alpha")   ,HHMass));
            fitter->setVar(pn1("alpha2")  ,prevJSON.evalFunc(pn1("alpha2")  ,HHMass));
            fitter->setVar(pn1("n")       ,prevJSON.evalFunc(pn1("n")       ,HHMass));
            fitter->setVar(pn1("n2")      ,prevJSON.evalFunc(pn1("n2")      ,HHMass));
            fitter->setVar(pn1("maxS")    ,prevJSON.evalFunc(pn1("maxS")    ,HHMass));
            fitter->setVar(pn1("mean_p1") ,prevJSON.evalFunc(pn1("mean_p1") ,HHMass));
            fitter->setVar(pn1("sigma_p1"),prevJSON.evalFunc(pn1("sigma_p1"),HHMass));
            fitter->setConst(pn1("alpha")   ,1);
            fitter->setConst(pn1("alpha2")  ,1);
            fitter->setConst(pn1("n")       ,1);
            fitter->setConst(pn1("n2")      ,1);
            fitter->setConst(pn1("maxS")    ,1);
            fitter->setConst(pn1("mean_p1") ,1);
            fitter->setConst(pn1("sigma_p1"),1);
        }

        fitter->w->var("MJJ")->setRange("fit",30,210);
        fitter->w->var("MVV")->setRange("fit",HHMass*.75,HHMass*1.25);
        fitter->w->var("MVV")->setRange("coef",minInclHHMass,maxInclHHMass);
        fitter->w->var("MJJ")->setRange("coef",30,210);
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::Minos(0),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::Minos(1),RooFit::NumCPU(8)});
    };

    int iP = 0;
    for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
        std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
        std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
        auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
        if(hbb_hh_H ==0) continue;
            setup2DFit(&*hbb_hh_H,signalMassBins[iP]);
            plotter.addFit2D(&*fitters.back(),signalMassBins[iS],ptName+"_2D");
        iP++;
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");
}


void makeSignal2DShapesFirstIteration(const std::string& name, const std::string& filename){
    std::string fitName = "fit1stIt";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I] || p == purSels[PUR_LMT] ) continue;
        if(h != hadSels[HAD_LTMB] ) continue;
        bool doExpo = p == purSels[PUR_L];

        std::string catName = l+"_"+p+"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_MJJ_fit.json");
        oldJSON.fillFunctions("MH");
        makeSignal2DShapes(name,filename,catName,fitName,oldJSON,iF,false, doExpo  );

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 550 -maxX 4550 -v ";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4";
        std::string jsonArgsVV0 = "meanSMVV:laur4,sigmaSMVV:laur4,alphaSMVV:laur4,alpha2SMVV:pol1,nSMVV:pol0,n2SMVV:pol0";
        std::string jsonArgsVV1 = "maxSSMVV:pol0,mean_p1SMVV:laur4,sigma_p1SMVV:laur4";
        std::string jsonArgs = argsP1 + (doExpo ? jsonArgsExpo : jsonArgsStd ) + ","+jsonArgsVV0+","+jsonArgsVV1;
        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        newJSON.replaceEntries(oldJSON);
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
        }
}

void makeSignal2DShapesSecondIteration(const std::string& name, const std::string& filename, bool firstIteration){
    std::string fitName = "fit";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l == lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I] || p == purSels[PUR_LMT] ) continue;
        if(h != hadSels[HAD_FULL] ) continue;
        bool doExpo = p == purSels[PUR_L];

        std::string catName = l+"_"+p+"_"+h;
        std::string oldJSONCatName = lepSels[LEP_EMU]+"_"+p+"_"+hadSels[HAD_LTMB];
        CJSON oldJSON(     filename+"_"+name+"_"+oldJSONCatName+"_fit1stIt.json");
        oldJSON.fillFunctions("MH");
        makeSignal2DShapes(name,filename,catName,fitName,oldJSON,iF,false, doExpo );

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 550 -maxX 4550";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4";
        std::string jsonArgsVV0 = "meanSMVV:laur4,sigmaSMVV:laur4,alphaSMVV:laur4,alpha2SMVV:pol1,nSMVV:pol0,n2SMVV:pol0";
        std::string jsonArgsVV1 = "maxS:pol0,mean_p1:laur4,mean_p2:laur4";
        std::string jsonArgs = argsP1 + (doExpo ? jsonArgsExpo : jsonArgsStd ) + ","+jsonArgsVV0+","+jsonArgsVV1;

        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        oldJSON.replaceEntry(std::string("meanS")+"MVV", newJSON.getP(std::string("meanS")+"MVV") );
        oldJSON.replaceEntry(std::string("sigmaS")+"MVV", newJSON.getP(std::string("sigmaS")+"MVV") );
        oldJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
        }
}

//
//void makeSignal2DShapes(const std::string& name, const std::string& filename){
//
//    auto setup1DFit = [](const TH1* hbbH, std::vector<std::unique_ptr<FunctionFitter>>& fitters ){
//        std::string pF = "SMJJ";
//        auto vN=[&](std::string var)->std::string{return var+pF;};
//        fitters.emplace_back(new CBFunctionFitter(hbbH,pF,{"MJJ"}));
//        auto fitter = &* fitters.back();
//        fitter->setVar(vN("mean")     ,125,90,180);
//        fitter->setVar(vN("sigma")       ,10,5,20);
//        //             fitter->setConst(vN("sigma"),1);
//        fitter->setVar(vN("alpha")     ,1 ,0.1,2);
//        //            fitter->setConst(vN("alpha"),1);
//        fitter->setVar(vN("alpha2")  ,1.5,0.1,3);
//        //             fitter->setConst(vN("alpha2")  ,1);
//        fitter->setVar(vN("n")   ,  5  ,1,6);
//        fitter->setVar(vN("n2")  ,5,3,6);
//        fitter->setConst(vN("n")  ,1);
//        fitter->setConst(vN("n2")  ,1);
//        fitter->w->var("MJJ")->setRange("fit",30,210);
//        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(0)});
//        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(1)});
//    };
//
//    auto setup2DFit = [](const TH2* hbbH, std::vector<std::unique_ptr<FunctionFitter>>& fitters, const double HHMass, FunctionFitter * fitter1D ){
//        std::string pF = "S";
//        auto vN=[&](std::string var)->std::string{return var+pF;};
//        auto pn0 =[&] (std::string v) ->std::string{return v+pF+"MJJ";};
//        auto pn1 =[&] (std::string v) ->std::string{return v+pF+"MVV";};
//
//        fitters.emplace_back(new CBFunctionFitter2D(hbbH,pF,{"MJJ","MVV"}));
//        auto fitter = &* fitters.back();
//        fitter->setVar(pn0("mean")  ,fitter1D->getVal(pn0("mean")  ),90,180);
//        fitter->setVar(pn0("sigma") ,fitter1D->getVal(pn0("sigma") ),5,20);
//        fitter->setVar(pn0("alpha") ,fitter1D->getVal(pn0("alpha") ),0.1,3);
//        fitter->setVar(pn0("alpha2"),fitter1D->getVal(pn0("alpha2")),0.1,3);
//        fitter->setVar(pn0("n")     ,fitter1D->getVal(pn0("n")     ),1,6);
//        fitter->setVar(pn0("n2")    ,fitter1D->getVal(pn0("n2")    ),3,6);
//        fitter->setConst(pn0("mean")  ,1);
//        fitter->setConst(pn0("sigma") ,1);
//        fitter->setConst(pn0("alpha") ,1);
//        fitter->setConst(pn0("alpha2"),1);
//        fitter->setConst(pn0("n")     ,1);
//        fitter->setConst(pn0("n2")    ,1);
//
//        fitter->setVar(pn1("p0_mean")  ,HHMass   ,HHMass -200,HHMass+200);
//        fitter->setVar(pn1("p0_sigma") ,HHMass*0.05,HHMass*0.025,HHMass*0.1);
//        fitter->setVar(pn1("p0_alpha") ,1.5,0.75,3);
//        fitter->setVar(pn1("p0_alpha2"),1.5,0.75,3);
//        fitter->setVar(pn1("p1_mean")  ,HHMass*.044,HHMass*.01,HHMass*.10);
//        fitter->setVar(pn1("p1_sigma") ,0,0,100);
//        fitter->setVar(pn1("p1_alpha") ,0,-10,10);
//        fitter->setVar(pn1("p1_alpha2"),0,-10,10);
//        fitter->setVar(pn1("maxS_mean")  ,2.5   ,0,5);
//        fitter->setVar(pn1("n")     ,5,1,6);
//        fitter->setVar(pn1("n2")    ,5,1,6);
////        fitter->setConst(pn1("p0_mean")  ,1);
////        fitter->setConst(pn1("p0_sigma") ,1);
////        fitter->setConst(pn1("p0_alpha") ,1);
////        fitter->setConst(pn1("p0_alpha2"),1);
////        fitter->setConst(pn1("p1_mean")  ,1);
////        fitter->setConst(pn1("p1_sigma") ,1);
//        fitter->setConst(pn1("p1_alpha") ,1);
//        fitter->setConst(pn1("p1_alpha2"),1);
//        fitter->setConst(pn1("n")     ,1);
//        fitter->setConst(pn1("n2")    ,1);
//        fitter->setConst(pn1("maxS_mean")    ,1);
//        fitter->w->var("MJJ")->setRange("fit",30,210);
//        fitter->w->var("MVV")->setRange("fit",HHMass*.75,HHMass*1.25);
//        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(0),RooFit::NumCPU(8)});
//        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::Minos(1),RooFit::NumCPU(8)});
////        fitter->fit({RooFit::SumW2Error(1),RooFit::Minos(0)});
////        fitter->fit({RooFit::SumW2Error(1),RooFit::Minos(1)});
//    };
//
//
//    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
//    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
//        if(l.find("emu") == std::string::npos ) continue;
//        if(p.find("L") == std::string::npos ) continue;
//        if(h.find("none") == std::string::npos ) continue;
//        FunctionParameterPlotter plotter;
//        FunctionParameterPlotter plotter2D;
//        std::vector<std::unique_ptr<FunctionFitter>> fitters;
//        std::vector<std::unique_ptr<FunctionFitter>> fitters2D;
//        std::string catName = l+"_"+p+"_"+h;
//        int iP = 0;
//        for(unsigned int iS = 0; iS < signalMassBins.size() -1; ++iS){
//            std::cout<<" !!!!!!! "<< catName <<" -> "<<signalMassBins[iS]<<std::endl;
//            std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
//            std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
//            auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
//            if(hbb_hh_H ==0) continue;
//            auto hbb_H = projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS);
////            auto small2D = cutHistogram(&*hbb_hh_H,hName+"_"+hbbMCS+"_"+hhMCS + "_cut",100,150,signalMassBins[iP]*.75,signalMassBins[iP]*1.25);
////            auto small2D = cutHistogram(&*hbb_hh_H,hName+"_"+hbbMCS+"_"+hhMCS + "_cut",100,150,0,7000);
//            setup1DFit(&*hbb_H,fitters);
//            plotter.addFit(&*fitters.back(),signalMassBins[iS],ptName);
////            setup2DFit(&*hbb_hh_H,fitters2D,signalMassBins[iP],&*fitters.back());
////            plotter2D.addFit2D(&*fitters2D.back(),signalMassBins[iS],ptName+"_2D");
//            iP++;
//        }
//        plotter.write(filename+"_"+name+"_"+catName+"_MJJ_fit.root");
//        plotter2D.write(filename+"_"+name+"_"+catName+"_fit.root");
//
////        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_MJJ_fit.root ";
////        std::string jsonArgs = " -g meanSMJJ:pol1,sigmaSMJJ:pol1,alphaSMJJ:pol1,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0 ";
////        MakeJSON(filename+"_"+name+"_"catName+"_MJJ_fit.json",argsP1+" "+jsonArgs);
//        }
//}

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
        fitters.emplace_back(new CBFunctionFitter(hbbH,false,"T"));
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
        fitters.emplace_back(new CBFunctionFitter(hbbH,false,"W"));
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
        CBFunctionFitter xFit(0,false,funcParamPostfix);
        xFit.w->var("MH")->setMin(minHbbMass);
        xFit.w->var("MH")->setMax(maxHbbMass);
        xFit.w->var("MH")->setBins(nHbbMassBins*10);
        xFit.w->var("MH")->setVal((minHbbMass+maxHbbMass)/2.);

        CJSON json(jsonFile);
        json.fillFunctions("MH");
        auto pn = [&](const std::string& v) ->std::string{return v + funcParamPostfix;};
        auto * iF =  TObjectHelper::getFile(filename + "_"+name+"_"+l+"_"+p+"_"+h+"_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH1>(iF,"histo",false,false);
        if(hh_H==0) continue;
        oF->cd();
        TH2 * h_2D = new TH2F((name+"_"+l+"_"+p+"_"+h).c_str(),(std::string(";")+hbbMCS.title+";"+hhMCS.title).c_str(),nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
        for(int iY = 1; iY <= hh_H->GetNbinsX(); ++iY){
            double yNorm = hh_H->GetBinContent(iY);
            double hhV = hh_H->GetBinCenter(iY);

            xFit.setVar(pn("mean"  ),json.evalFunc(pn("mean"  ),hhV));
            xFit.setVar(pn("sigma" ),json.evalFunc(pn("sigma" ),hhV));
            xFit.setVar(pn("alpha" ),json.evalFunc(pn("alpha" ),hhV));
            xFit.setVar(pn("alpha2"),json.evalFunc(pn("alpha2"),hhV));
            xFit.setVar(pn("n"     ),json.evalFunc(pn("n"     ),hhV));
            xFit.setVar(pn("n2"    ),json.evalFunc(pn("n2"    ),hhV));
            auto xHist = xFit.pdf1D(name+"_"+l+"_"+p+"_"+h+ASTypes::int2Str(iY)+"_hbb");
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
    std::string signalTrees = treeDir + "/out_radion_hh_bbinc_mXXX_0.root";
    //do q/g
    if(false)
    {
        std::string name = bkgSels[BKG_QG];
        std::string genSel = bkgSels[BKG_QG].cut + "&&"+ aQCD.cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+ hadSels[HAD_LTMB].cut;
        makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        //MVV
        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,3,0.5,1);//x = hh
        //                makeBackgroundShapesMVVConditional(name+"_xs_0p75_xc_2_ys_0p75_yc_2,filename,treeArea,baseSel,0.75,2,0.75,2);//old
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hhInclRange.cut,true);

        //MJJ
        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel+"&&"+hhRange.cut);
        mergeBackgroundShapes(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
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
    if(true){
        std::string name = radionSig;
//        makeSignalFittingDistributions(name,filename,signalTrees,hhRange.cut+"&&"+hbbRange.cut,false);
//        makeSignalFittingDistributions(name,filename,signalTrees,hhInclRange.cut+"&&"+hbbInclRange.cut,true);
//        makeSignal2DShapes(name,filename);
//        makeSignalMJJShapes(name,filename);
        makeSignal2DShapesFirstIteration(name,filename);

    }

}
#endif

void makeInputs(std::string treeDir = "trees/"){
    go(treeDir);
}
