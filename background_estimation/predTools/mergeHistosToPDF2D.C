
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
        xIsCond   = p.addBool("xIsCond","If true, set x as the conditional variable P(X|Y)*P(Y), if false: P(Y|X)*PY(Y)");
        auto in1D      = p.addString("in1D" ,"file containing 1D template",true);
        auto in2D     = p.addString("in2D" ,"file containing 2D template");

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

        f1D =  TObjectHelper::getFile(*in1D);
        f2D =  TObjectHelper::getFile(*in2D);

        makeHisto(*name,*name,*name);
        for(const auto& syst: xSysts){
            if(*xIsCond){
                makeHisto(*name +"_"+syst.second+"Up"  ,*name,*name+"_"+syst.first+"Up"  );
                makeHisto(*name +"_"+syst.second+"Down",*name,*name+"_"+syst.first+"Down");
            } else {
                makeHisto(*name +"_"+syst.second+"Up"  ,*name+"_"+syst.first+"Up"  ,*name);
                makeHisto(*name +"_"+syst.second+"Down",*name+"_"+syst.first+"Down",*name);
            }

        }
        for(const auto& syst: ySysts){
            if(*xIsCond){
                makeHisto(*name +"_"+syst.second+"Up"  ,*name+"_"+syst.first+"Up"  ,*name);
                makeHisto(*name +"_"+syst.second+"Down",*name+"_"+syst.first+"Down",*name);
            } else {
                makeHisto(*name +"_"+syst.second+"Up"  ,*name,*name+"_"+syst.first+"Up"  );
                makeHisto(*name +"_"+syst.second+"Down",*name,*name+"_"+syst.first+"Down");
            }
        }


        plotter.write(outFileName);
        f1D->Close();
        f2D->Close();
}


    void makeHisto(const std::string& outName, const std::string& h1DName, const std::string& h2DName){
        TH1 * h1D = 0; TH2* h2D = 0;
        f1D->GetObject(h1DName.c_str(),h1D);f2D->GetObject(h2DName.c_str(),h2D);
        if(h1D == 0 || h2D == 0) return;

        auto outH = new TH2F("temp","",(*xb)[0],(*xb)[1],(*xb)[2],(*yb)[0],(*yb)[1],(*yb)[2]);

        //first cut up conditional template
        for(int iX = 1; iX <= h2D->GetNbinsX(); ++iX){
            int oBinX = outH->GetXaxis()->FindFixBin(h2D->GetXaxis()->GetBinCenter(iX));
            if(oBinX < 1 || oBinX > outH->GetNbinsX() ) continue;
            for(int iY = 1; iY <= h2D->GetNbinsY(); ++iY){
                int oBinY = outH->GetYaxis()->FindFixBin(h2D->GetYaxis()->GetBinCenter(iY));
                if(oBinY < 1 || oBinY > outH->GetNbinsY() ) continue;
                outH->SetBinContent(oBinX,oBinY,h2D->GetBinContent(iX,iY));
            }
        }

        //Now apply the 1D template
        if(*xIsCond){
            for(int iY = 1; iY <= h1D->GetNbinsX(); ++iY){
                int oBinY = outH->GetYaxis()->FindFixBin(h1D->GetXaxis()->GetBinCenter(iY));
                if(oBinY < 1 || oBinY > outH->GetNbinsY() ) continue;
                const double SF = h1D->GetBinContent(iY)/outH->Integral(0,-1,oBinY,oBinY);
                for(int oBinX = 1; oBinX <= outH->GetNbinsX(); ++oBinX){
                    outH->SetBinContent(oBinX,oBinY,outH->GetBinContent(oBinX,oBinY)*SF);
                }
            }
        } else {
            for(int iX = 1; iX <= h1D->GetNbinsX(); ++iX){
                int oBinX = outH->GetXaxis()->FindFixBin(h1D->GetXaxis()->GetBinCenter(iX));
                if(oBinX < 1 || oBinX > outH->GetNbinsX() ) continue;
                const double SF = h1D->GetBinContent(iX)/outH->Integral(oBinX,oBinX,0,-1);
                for(int oBinY = 1; oBinY <= outH->GetNbinsY(); ++oBinY){
                    outH->SetBinContent(oBinX,oBinY,outH->GetBinContent(oBinX,oBinY)*SF);
                }
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
    TFile * f1D = 0;
    TFile * f2D = 0;
    std::shared_ptr<bool>                      xIsCond;
    std::shared_ptr<std::vector<double>>       xb;
    std::shared_ptr<std::vector<double>>       yb;

};

#endif

void mergeHistosToPDF2D(std::string outFileName,std::string arguments){
    mergeHistosToPDF2DAnalyzer a(outFileName, arguments);
}
