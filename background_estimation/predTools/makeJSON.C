
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


class CJSON {
public:
    CJSON() {} //For filling by hand
    CJSON(const std::string& inFName){ //For reading in a file
        std::ifstream file (inFName);
        if (!file.is_open())  throw std::invalid_argument("CJSON::CJSON() -> Bad file");
        std::stringstream strStream;
        strStream << file.rdbuf();
        std::string str = strStream.str();
        str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
        str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
        auto regMtch = std::regex("\", \"");
        auto regMtch2 = std::regex("\": \"");
        std::vector<std::string> paramFits(std::sregex_token_iterator(str.begin(), str.end(),regMtch,-1), std::sregex_token_iterator());
        for(const auto& s :paramFits){
            std::vector<std::string> ps(std::sregex_token_iterator(s.begin(), s.end(), regMtch2, -1), std::sregex_token_iterator());
            if(ps.size() != 2) {
                for(auto& p :ps) std::cout << p <<" ";
                std::cout <<std::endl;
                throw std::invalid_argument("CJSON::CJSON() -> Bad parsing");
            }
            ps[0].erase(std::remove(ps[0].begin(), ps[0].end(), '"'), ps[0].end());
            ps[1].erase(std::remove(ps[1].begin(), ps[1].end(), '"'), ps[1].end());
            addEntry(ps[0],ps[1]);
        }
        file.close();
    }

    void replaceEntries(const CJSON& other) {
        for(unsigned int iP = 0; iP < other.getNP(); ++iP){
            const auto& oP = other.getP(iP);
            replaceEntry(oP.first,oP.second);
        }
    }

    void replaceEntry(const std::string& name, const std::string& value) {
        for(auto& p : parameters){
            if(p.first != name) continue;
            p.second = value;
        }
    }
    void addEntry(const std::string& name, const std::string& value) { parameters.emplace_back(name,value);}
    void write(const std::string& outFile){
        std::ofstream outJSON(outFile.c_str(),std::ios::out|std::ios::trunc);
        dumpJSON(outJSON);
        outJSON.close();
    }
    unsigned int getNP() const {return parameters.size();}

    void fillFunctions(const std::string& xVarName, const std::string& yVarName = ""){
        for(unsigned int i = 0; i< parameters.size(); ++i){
            functions.emplace_back(parameters[i].first, std::unique_ptr<TF1>(getFunction(i,xVarName,yVarName)));
        }

    }

    const std::pair<std::string,std::string>& getP(unsigned int idx) const {return parameters[idx];}

    std::string getP(const std::string& name) const {
        for(auto& p : parameters){
            if(p.first != name) continue;
            return p.second;
        }
        return "";
    }

    float evalFunc(const std::string& name, double xVal) const {
        for(auto& p : functions){
            if(p.first != name) continue;
            return p.second->Eval(xVal);
        }
        return 0;
    }

private:
    void dumpJSON(std::ofstream& f ){
        f<<"{";
        for(unsigned int iP = 0; iP < parameters.size(); ++iP){
            const auto& p = parameters[iP];
            std::string ostr = std::string("\"") + p.first +"\": \""+ p.second+"\"";
            if(iP + 1 < parameters.size()) ostr +=", ";
            f<< ostr;
        }
        f<<"}";
    }
    TF1* getFunction(unsigned int idx, const std::string& xVarName, const std::string& yVarName = ""){
        std::string pstr = parameters[idx].second;
        auto replace = [&](const std::string& vn, const std::string tf1n){
            std:size_t index = 0;
            while (true) {
                 index = pstr.find(vn, index);
                 if (index == std::string::npos) break;
                 pstr.replace(index, vn.size(), tf1n);
                 index += 1;
            }
        };
        replace(xVarName,"x");
        if(yVarName.size()) replace(yVarName,"y");
        return new TF1(parameters[idx].first.c_str(),pstr.c_str(),1,13000);
    }

    std::vector<std::pair<std::string,std::string>> parameters;
    std::vector<std::pair<std::string,std::unique_ptr<TF1>> >functions;
};

CJSON makeJSON(const std::string& outFileName,std::string& arguments){
    typedef std::vector<std::pair<std::string,std::string> > ParamNames;
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    auto getParamList =[&](const std::string& inList, ParamNames& outNames){
        outNames.clear();
        auto regMtch = std::regex(",");
        auto regMtch2 = std::regex(":");
        std::vector<std::string> systList(std::sregex_token_iterator(inList.begin(), inList.end(), regMtch, -1), std::sregex_token_iterator());
        for(const auto& s :systList){
            std::vector<std::string> names(std::sregex_token_iterator(s.begin(), s.end(), regMtch2, -1), std::sregex_token_iterator());
            if(names.size() != 2) {
                std::cout << inList<<std::endl;
                throw std::invalid_argument("MakeJSON::ParamNames() -> Bad parsing");
            }
            outNames.emplace_back(names[0],names[1]);
        }
    };

    auto returnFString =[](const TF1* func, const std::string& var) ->std::string {
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
    };



    ParParser p;
    auto iFn   = p.addString("i","input file name",true);
    auto gs    = p.addString("g" ,"Comma separated graphs and functions to fit  like MEAN:pol3,SIGMA:pol2",true);
    auto v     = p.addString("var","x var name ",false,"XX");
    auto minX  = p.addFloat("minX","minimum x",true);
    auto maxX  = p.addFloat("maxX","maximum x",true);
    p.parse(arguments);

    ParamNames params;
    getParamList(*gs,params);

    auto iF =  TObjectHelper::getFile(*iFn);
    TFile * oRF = new TFile((outFileName + ".root").c_str(),"recreate" );
    oRF->cd();
    CJSON outJSON;

    for(const auto& p : params){
        auto g = TObjectHelper::getObjectNoOwn<TGraphErrors>(iF,p.first);
        TF1 * func = 0;
        if(p.second.find("pol") != std::string::npos) func=new TF1(p.second.c_str(),p.second.c_str(),0,13000);
        else  if(p.second.find("llog") != std::string::npos){
            func=new  TF1("llog","[0]+[1]*log(x)",1,13000);
            func->SetParameters(1,1);
        } else if(p.second.find("laur") != std::string::npos){
            auto regMtch = std::regex("laur");
            std::vector<std::string> laurPs(std::sregex_token_iterator(p.second.begin(), p.second.end(), regMtch, -1), std::sregex_token_iterator());
            int order = std::stoi(laurPs[1]);
            std::string fstr = "0";
            for(int iO = 0; iO < order; ++iO){
                fstr=fstr+"+["+std::to_string(iO)+"]"+"/x^"+std::to_string(iO);
            }
            func=new TF1(p.second.c_str(),fstr.c_str(),1,13000);
            for(int iO = 0; iO < order; ++iO){ func->SetParameter(iO,0);}
        } else if(p.second.find("FIX")!= std::string::npos){
            auto regMtch = std::regex("p");
            std::vector<std::string> laurPs(std::sregex_token_iterator(p.second.begin(), p.second.end(), regMtch, -1), std::sregex_token_iterator());
            std::string fstring = "x <= ([2]-[3]) ? [0] : (x >= ([2]+[3]) ? [1] : ( 0.5*([0]+[1] + ([0]-[1])*sin([4]*1.57*(x-[2])/[3]))))";
            func=new TF1(p.second.c_str(),fstring.c_str(),1,13000);
            for(unsigned int iP = 1; iP < laurPs.size(); ++iP){
                func->FixParameter(iP-1,std::stof(laurPs[iP]));
            }
        }
        if(func==0) throw std::invalid_argument("MakeJSON::MakeJSON() -> Bad parsing");
        g->Fit(func,"","",*minX,*maxX);
        g->Write(p.first.c_str());
        func->Write((p.first+"_func").c_str());
        outJSON.addEntry(p.first,returnFString(&*func,*v));
    }
    oRF->Close();
    return outJSON;
}
#endif

void MakeJSON(std::string outFileName,std::string arguments){
    auto json = makeJSON(outFileName, arguments);
    json.write(outFileName);
}
CJSON getJSON(std::string outFileName,std::string arguments){
    return makeJSON(outFileName, arguments);
}
