
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include <string.h>
#include <regex>





class CutHistos1DAnalyzer {
public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    CutHistos1DAnalyzer(std::string outFileName,std::string arguments )
{


        ParParser p;
        auto name = p.addString("n","Histogram base names",true);
        auto inX  = p.addString("i" ,"file containing x template",true);
        auto sX   = p.addString("s"  ,"Comma seperated list of systematics-> TH1Name:SystName",true);
        xb   = p.addVFloat("xb","x-variable binning",true);
        p.parse(arguments);

        SystNames xSysts;
        SystNames ySysts;
        getSystList(*sX,xSysts);

        if(xb->size() != 3)                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");

        fX =  TObjectHelper::getFile(*inX);
        makeHisto(*name,*name);
        for(const auto& syst: xSysts){
            makeHisto(*name +"_"+syst.second+"Up"  ,*name+"_"+syst.first+"Up"  );
            makeHisto(*name +"_"+syst.second+"Down",*name+"_"+syst.first+"Down");
        }


        plotter.write(outFileName);
        fX->Close();
}

    void makeHisto(const std::string& outName, const std::string& xName){
        TH1 * inX = 0;
        fX->GetObject(xName.c_str(),inX);
        if(inX == 0 ) return;
        auto outH = new TH1F("temp","",(*xb)[0],(*xb)[1],(*xb)[2]);
        for(int iX = 1; iX <= inX->GetNbinsX(); ++iX){
            int oBinX = outH->GetXaxis()->FindFixBin(inX->GetXaxis()->GetBinCenter(iX));
            if(oBinX < 1 || oBinX > outH->GetNbinsX() ) continue;
            outH->SetBinContent(oBinX,inX->GetBinContent(iX));
            outH->SetBinError(oBinX,inX->GetBinError(iX));
        }
        outH->Scale(1./outH->Integral());
        outH->SetName(outName.c_str());
        plotter.add1D(outH);
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
    std::shared_ptr<std::vector<double>>       xb;

};

#endif

void cutHistos1D(std::string outFileName,std::string arguments){
    CutHistos1DAnalyzer a(outFileName, arguments);
}
