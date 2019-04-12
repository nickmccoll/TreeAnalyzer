
#ifndef TREEANALYZER_BACKGROUNDESTIMATION_PREDTOOLS_INPUTSHELPER_H_
#define TREEANALYZER_BACKGROUNDESTIMATION_PREDTOOLS_INPUTSHELPER_H_
#include <string>
#include <TSystem.h>

#include "CutConstants.h"
#include <thread>
#include <chrono>
#include <future>
#include <memory>

#include "TH2.h"
#include "TMath.h"

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
void fillTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,TH2 * oH){
    RooArgSet vars(*xV) ;
    vars.add(*yV) ;
    const int nXBins = oH->GetNbinsX();
    const int nYBins = oH->GetNbinsY();

    for( int iInX = 1; iInX <= nXBins; ++iInX ){
      double centerX = oH->GetXaxis()->GetBinCenter(iInX);
      double widthX = oH->GetXaxis()->GetBinWidth(iInX);
      xV->setVal(centerX);
      for( int iInY = 1; iInY <= nYBins; ++iInY ){
        double centerY = oH->GetYaxis()->GetBinCenter(iInY);
        double widthY = oH->GetYaxis()->GetBinWidth(iInY);
        yV->setVal(centerY);
        double val = pdf->getVal(vars);
        oH->Fill(centerX,centerY,val*widthX*widthY);
      }
    }
    oH->Scale(1.0/oH->Integral());
}
TH2* createTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double xMin, const double xMax,
        const unsigned int nYBins, const double yMin, const double yMax){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xMin,xMax,nYBins,yMin,yMax);
    fillTH2FromPDF(pdf,xV,yV,oH);
    return oH;
}
TH2* createTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double* xBins,
        const unsigned int nYBins, const double yMin, const double yMax){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xBins,nYBins,yMin,yMax);
    fillTH2FromPDF(pdf,xV,yV,oH);
    return oH;
}
TH2* createTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double xMin, const double xMax,
        const unsigned int nYBins, const double* yBins){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xMin,xMax,nYBins,yBins);
    fillTH2FromPDF(pdf,xV,yV,oH);
    return oH;
}
TH2* createTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double* xBins,
        const unsigned int nYBins, const double* yBins){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xBins,nYBins,yBins);
    fillTH2FromPDF(pdf,xV,yV,oH);
    return oH;
}

TH2* createTH2FromPDF(RooAbsPdf* pdf, RooRealVar* xV,RooRealVar* yV,
        const std::string& name,const std::string& title,
        const TAxis * xAxis,const TAxis * yAxis){
    if(xAxis->GetXbins()->GetSize() == 0 && yAxis->GetXbins()->GetSize() == 0){
        return  createTH2FromPDF(pdf,xV,yV,name,title,
                xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax(),
                yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax());
    } else if(xAxis->GetXbins()->GetSize() == 0 && yAxis->GetXbins()->GetSize() != 0){
        return  createTH2FromPDF(pdf,xV,yV,name,title,
                xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax(),
                yAxis->GetNbins(),yAxis->GetXbins()->GetArray());
    } else if(xAxis->GetXbins()->GetSize() != 0 && yAxis->GetXbins()->GetSize() == 0){
        return  createTH2FromPDF(pdf,xV,yV,name,title,
                xAxis->GetNbins(),xAxis->GetXbins()->GetArray(),
                yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax());
    } else{
        return  createTH2FromPDF(pdf,xV,yV,name,title,
                xAxis->GetNbins(),xAxis->GetXbins()->GetArray(),
                yAxis->GetNbins(),yAxis->GetXbins()->GetArray());
    }
}


void fillTH1FromPDF(RooAbsPdf* pdf, RooRealVar* xV,TH1 * oH){
    RooArgSet vars(*xV) ;
    const int nXBins = oH->GetNbinsX();

    for( int iInX = 1; iInX <= nXBins; ++iInX ){
      double centerX = oH->GetXaxis()->GetBinCenter(iInX);
      double widthX = oH->GetXaxis()->GetBinWidth(iInX);
      xV->setVal(centerX);
      double val = pdf->getVal(vars);
      oH->Fill(centerX,val*widthX);
    }
    oH->Scale(1.0/oH->Integral());
}
TH1* createTH1FromPDF(RooAbsPdf* pdf, RooRealVar* xV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double xMin, const double xMax){
    TH1 * oH = new TH1F(name.c_str(),title.c_str(),nXBins,xMin,xMax);
    fillTH1FromPDF(pdf,xV,oH);
    return oH;
}
TH1* createTH1FromPDF(RooAbsPdf* pdf, RooRealVar* xV,
        const std::string& name,const std::string& title,
        const unsigned int nXBins, const double* xBins){
    TH1 * oH = new TH1F(name.c_str(),title.c_str(),nXBins,xBins);
    fillTH1FromPDF(pdf,xV,oH);
    return oH;
}

TH1* createTH1FromPDF(RooAbsPdf* pdf, RooRealVar* xV,
        const std::string& name,const std::string& title,
        const TAxis * xAxis){

    return xAxis->GetXbins()->GetSize()
            ? createTH1FromPDF(pdf,xV,name,title
                    ,xAxis->GetNbins(),xAxis->GetXbins()->GetArray())
            : createTH1FromPDF(pdf,xV,name,title
                    ,xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax());
}



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
TH2* cutHistogram(const TH2* inH, TH2* outH ){
    for(int iX =1; iX <= inH->GetNbinsX(); ++iX){
        const int outIX =outH->GetXaxis()->FindFixBin(inH->GetXaxis()->GetBinCenter(iX));
        if(outIX < 1 || outIX > outH->GetNbinsX() ) continue;
        if(inH->GetXaxis()->GetBinWidth(iX)!=outH->GetXaxis()->GetBinWidth(outIX))
            throw std::invalid_argument("cutHistogram::cutHistogram() -> Bad x-axis");
        for(int iY =1; iY <= inH->GetNbinsY(); ++iY){
            const int outIY = outH->GetYaxis()->FindFixBin(inH->GetYaxis()->GetBinCenter(iY));
            if(outIY < 1 || outIY > outH->GetNbinsY() ) continue;
            if(inH->GetYaxis()->GetBinWidth(iY)!=outH->GetYaxis()->GetBinWidth(outIY))
                throw std::invalid_argument("cutHistogram::cutHistogram() -> Bad y-axis");
            outH->SetBinContent(outIX,outIY,inH->GetBinContent(iX,iY));
            outH->SetBinError(outIX,outIY,inH->GetBinError(iX,iY));
        }
    }
    return outH;
}

TH2* cutHistogram(const TH2* inHist, const std::string& name, const std::string& title,
        const unsigned int nXBins, const double xMin, const double xMax,
        const unsigned int nYBins, const double yMin, const double yMax){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xMin,xMax,nYBins,yMin,yMax);
    cutHistogram(inHist,oH);
    return oH;
}
TH2* cutHistogram(const TH2* inHist, const std::string& name, const std::string& title,
        const unsigned int nXBins, const double* xBins,
        const unsigned int nYBins, const double yMin, const double yMax){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xBins,nYBins,yMin,yMax);
    cutHistogram(inHist,oH);
    return oH;
}
TH2* cutHistogram(const TH2* inHist, const std::string& name, const std::string& title,
        const unsigned int nXBins, const double xMin, const double xMax,
        const unsigned int nYBins, const double* yBins){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xMin,xMax,nYBins,yBins);
    cutHistogram(inHist,oH);
    return oH;
}
TH2* cutHistogram(const TH2* inHist, const std::string& name, const std::string& title,
        const unsigned int nXBins, const double* xBins,
        const unsigned int nYBins, const double* yBins){
    TH2 * oH = new TH2F(name.c_str(),title.c_str(),nXBins,xBins,nYBins,yBins);
    cutHistogram(inHist,oH);
    return oH;
}


TH2* cutHistogram(const TH2* inHist, const std::string& name, const std::string& title,
        const TAxis * xAxis,const TAxis * yAxis){
    if(xAxis->GetXbins()->GetSize() == 0 && yAxis->GetXbins()->GetSize() == 0){
        return  cutHistogram(inHist,name,title,
                xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax(),
                yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax());
    } else if(xAxis->GetXbins()->GetSize() == 0 && yAxis->GetXbins()->GetSize() != 0){
        return  cutHistogram(inHist,name,title,
                xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax(),
                yAxis->GetNbins(),yAxis->GetXbins()->GetArray());
    } else if(xAxis->GetXbins()->GetSize() != 0 && yAxis->GetXbins()->GetSize() == 0){
        return  cutHistogram(inHist,name,title,
                xAxis->GetNbins(),xAxis->GetXbins()->GetArray(),
                yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax());
    } else{
        return  cutHistogram(inHist,name,title,
                xAxis->GetNbins(),xAxis->GetXbins()->GetArray(),
                yAxis->GetNbins(),yAxis->GetXbins()->GetArray());
    }
}
#endif
