
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include <string.h>
#include <regex>





class mergeHistosToPDF2DAnalyzer {
public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    mergeHistosToPDF2DAnalyzer(std::string outFileName,std::string arguments )
{


        ParParser p;
        auto name = p.addString("n","Histogram base names",true);
        auto inX  = p.addString("inX" ,"file containing x template",true);
        auto inY  = p.addString("inY" ,"file containing 2D template");
        auto sX   = p.addString("sX"  ,"Comma seperated list of systematics-> TH1Name:SystName",true);
        auto sY   = p.addString("sY"  ,"Comma seperated list of systematics-> TH1Name:SystName",true);
        xb   = p.addVFloat("xb","x-variable binning",true);
        yb   = p.addVFloat("yb","y-variable binning",true);
        p.parse(arguments);

        SystNames xSysts;
        SystNames ySysts;
        getSystList(*sX,xSysts);
        getSystList(*sY,ySysts);

        if(xb->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");
        if(yb->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        fX =  TObjectHelper::getFile(*inX);
        fY =  TObjectHelper::getFile(*inY);

        makeHisto(*name,*name,*name);
        for(const auto& syst: xSysts){
            makeHisto(*name +"_"+syst.second+"Up"  ,*name+"_"+syst.first+"Up"  ,*name);
            makeHisto(*name +"_"+syst.second+"Down",*name+"_"+syst.first+"Down",*name);
        }
        for(const auto& syst: ySysts){
            makeHisto(*name +"_"+syst.second+"Up"  ,*name,*name+"_"+syst.first+"Up"  );
            makeHisto(*name +"_"+syst.second+"Down",*name,*name+"_"+syst.first+"Down");
        }


        plotter.write(outFileName);
        fX->Close();
        fY->Close();
}

    //the x axis of inY maps to the out y-axis
    //the y axis of inY maps to the out x-axis
    //the x axis of inX maps to the out x-axis
    void makeHisto(const std::string& outName, const std::string& xName, const std::string& yName){
//        auto inX =TObjectHelper::getObject<TH1F>(fX,xName);
//        auto inY =TObjectHelper::getObject<TH2F>(fY,yName);
        TH1 * inX = 0; TH2* inY = 0;
        fX->GetObject(xName.c_str(),inX);fY->GetObject(yName.c_str(),inY);
        if(inX == 0 || inY == 0) return;
        auto outH = new TH2F("temp","",(*xb)[0],(*xb)[1],(*xb)[2],(*yb)[0],(*yb)[1],(*yb)[2]);
        for(int iX = 1; iX <= inX->GetNbinsX(); ++iX){
            int oBinX = outH->GetXaxis()->FindFixBin(inX->GetXaxis()->GetBinCenter(iX));
            if(oBinX < 1 || oBinX > outH->GetNbinsX() ) continue;
            for(int iY = 1; iY <= inY->GetNbinsX(); ++iY){
                int oBinY = outH->GetYaxis()->FindFixBin(inY->GetXaxis()->GetBinCenter(iY));
                if(oBinY < 1 || oBinY > outH->GetNbinsY() ) continue;
                outH->SetBinContent(oBinX,oBinY,inY->GetBinContent(iY,iX));
            }
            const double SF = inX->GetBinContent(iX)/outH->Integral(oBinX,oBinX,0,-1);
            for(int oBinY = 1; oBinY <= outH->GetNbinsY(); ++oBinY){
                outH->SetBinContent(oBinX,oBinY,outH->GetBinContent(oBinX,oBinY)*SF);
            }
        }
        outH->Scale(1./outH->Integral());
        outH->SetName(outName.c_str());
        plotter.add2D(outH);
    }

    void getSystList(const std::string& inList, SystNames& outNames){
        outNames.clear();
        std::vector<std::string> systList(std::sregex_token_iterator(inList.begin(), inList.end(), std::regex(","), -1), std::sregex_token_iterator());
        for(const auto& s :systList){
            std::vector<std::string> names(std::sregex_token_iterator(s.begin(), s.end(), std::regex(":"), -1), std::sregex_token_iterator());
            if(names.size() != 2) {
                std::cout << inList<<std::endl;
                throw std::invalid_argument("mergeHistosToPDF2DAnalyzer::getSystList() -> Bad parsing");
            }
            outNames.emplace_back(names[0],names[1]);
        }
    }

    HistGetter plotter;
    TFile * fX = 0;
    TFile * fY = 0;
    std::shared_ptr<std::vector<double>>       xb;
    std::shared_ptr<std::vector<double>>       yb;

};

#endif

void mergeHistosToPDF2D(std::string outFileName,std::string arguments){
    mergeHistosToPDF2DAnalyzer a(outFileName, arguments);
}
