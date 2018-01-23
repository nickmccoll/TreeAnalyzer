
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include <string.h>
#include <regex>
#include <fstream>
#include <vector>
#include <utility>
#include "TGraphErrors.h"
#include "TF1.h"

class MakeJSONAnalyzer {
    typedef std::vector<std::pair<std::string,std::string> > ParamNames;
    void getParamList(const std::string& inList, ParamNames& outNames){
        outNames.clear();
        std::vector<std::string> systList(std::sregex_token_iterator(inList.begin(), inList.end(), std::regex(","), -1), std::sregex_token_iterator());
        for(const auto& s :systList){
            std::vector<std::string> names(std::sregex_token_iterator(s.begin(), s.end(), std::regex(":"), -1), std::sregex_token_iterator());
            if(names.size() != 2) {
                std::cout << inList<<std::endl;
                throw std::invalid_argument("MakeJSON::ParamNames() -> Bad parsing");
            }
            outNames.emplace_back(names[0],names[1]);
        }
    }

    std::string returnFString(const TF1* func, const std::string& var) const {
        const std::string name = func->GetName();

        auto getPower = [&](unsigned int nP) ->std::string{
          std::string out;
          for(unsigned int i = 0; i <nP; ++i) out +=std::string("*")+var;
          return out;
        };

        if(name.find("pol") != std::string::npos){
            std::string outF("0");
            for(int iP = 0; iP < func->GetNpar(); ++iP) outF = outF+"+("+ASTypes::flt2Str(func->GetParameter(iP))+")"+getPower(iP);
            return outF;
        } else if(name.find("llog") != std::string::npos){
            return ASTypes::flt2Str(func->GetParameter(0))+"+"+ASTypes::flt2Str(func->GetParameter(1))+"*log("+var+")";
        } else if(name.find("laur") != std::string::npos){
            std::string outF("0");
            for(int iP = 0; iP < func->GetNpar(); ++iP) outF = outF+"+("+ASTypes::flt2Str(func->GetParameter(iP))+")/"+var+"^"+std::to_string(iP);
            return outF;
        } else if(name.find("FIX") != std::string::npos){
            std::string fstring = var+" <= ([2]-[3]) ? [0] : ("+var+" >= ([2]+[3]) ? [1] : ( 0.5*([0]+[1] + ([0]-[1])*sin([4]*1.57*("+var+"-[2])/[3]))))";
            for(int iP = 0; iP < func->GetNpar(); ++iP){
                fstring = std::regex_replace(fstring,std::regex(std::string("\\[")+ASTypes::int2Str(iP)+"\\]"),ASTypes::flt2Str(func->GetParameter(iP)));
            }
            return fstring;
        }
        return "";
    }
    void dumpJSON(const std::vector<std::pair<std::string,std::string>>& ps,std::ofstream& f ){
        f<<"{";
        for(unsigned int iP = 0; iP < ps.size(); ++iP){
            const auto& p = ps[iP];
            std::string ostr = std::string("\"") + p.first +"\": \""+ p.second+"\"";
            if(iP + 1 < ps.size()) ostr +=", ";
            f<< ostr;
        }
        f<<"}";
    }

public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    MakeJSONAnalyzer(const std::string& outFileName,std::string arguments )
{


        ParParser p;
        auto iFn   = p.addString("i","input file name",true);
        auto gs    = p.addString("g" ,"Comma separated graphs and functions to fit  like MEAN:pol3,SIGMA:pol2",true);
        auto v     = p.addString("var","x var name ",false,"MH");
        auto minX  = p.addFloat("minX","minimum x",true);
        auto maxX  = p.addFloat("maxX","maximum x",true);
        p.parse(arguments);

        ParamNames params;
        getParamList(*gs,params);

        auto iF =  TObjectHelper::getFile(*iFn);
        TFile * oRF = new TFile((std::string("debug_") + outFileName + ".root").c_str(),"recreate" );
        oRF->cd();

        std::vector<std::pair<std::string,std::string>> fStrs;

        for(const auto& p : params){
            auto g = TObjectHelper::getObjectNoOwn<TGraphErrors>(iF,p.first);
            TF1 * func = 0;
            if(p.second.find("pol") != std::string::npos) func=new TF1(p.second.c_str(),p.second.c_str(),0,13000);
            else  if(p.second.find("llog") != std::string::npos){
                func=new  TF1("llog","[0]+[1]*log(x)",1,13000);
                func->SetParameters(1,1);
            } else if(p.second.find("laur") != std::string::npos){
                std::vector<std::string> laurPs(std::sregex_token_iterator(p.second.begin(), p.second.end(), std::regex("laur"), -1), std::sregex_token_iterator());
                int order = std::stoi(laurPs[1]);
                std::string fstr = "0";
                for(int iO = 0; iO < order; ++iO){
                    fstr=fstr+"+["+std::to_string(iO)+"]"+"/x^"+std::to_string(iO);
                }
                func=new TF1(p.second.c_str(),fstr.c_str(),1,13000);
                for(int iO = 0; iO < order; ++iO){ func->SetParameter(iO,0);}
            } else if(p.second.find("FIX")!= std::string::npos){
                std::vector<std::string> laurPs(std::sregex_token_iterator(p.second.begin(), p.second.end(), std::regex("p"), -1), std::sregex_token_iterator());
                std::string fstring = "x <= ([2]-[3]) ? [0] : (x >= ([2]+[3]) ? [1] : ( 0.5*([0]+[1] + ([0]-[1])*sin([4]*1.57*(x-[2])/[3]))))";
                func=new TF1(p.second.c_str(),fstring.c_str(),1,13000);
                for(unsigned int iP = 1; iP < laurPs.size(); ++iP){
                    func->FixParameter(iP-1,std::stof(laurPs[iP]));
                }
            }
            if(func==0) throw std::invalid_argument("MakeJSON::MakeJSON() -> Bad parsing");
            g->Fit(func,"","",*minX,*maxX);
            fStrs.emplace_back(p.first, returnFString(&*func,*v));
            g->Write(p.first.c_str());
            func->Write((p.first+"_func").c_str());
        }
        oRF->Close();

        std::ofstream outJSON(outFileName.c_str(),std::ios::out|std::ios::trunc);
        dumpJSON(fStrs,outJSON);
        outJSON.close();
}
};

#endif

void MakeJSON(std::string outFileName,std::string arguments){
    MakeJSONAnalyzer a(outFileName, arguments);
}
