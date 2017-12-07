
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TTreeFormula.h"
using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName,std::string arguments ) : BaseTreeAnalyzer(fileName,treeName,2){
        ParParser p;
        auto n   = p.addString("n","histogram name",true,"xHisto");
        auto x   = p.addString("x","x-variable",true);
        auto xb  = p.addVFloat("xb","x-variable binning",true);
        auto xlb = p.addBool("xlb","x-variable list binning (explicit bin edges)");
        auto s   = p.addString("s","selection",false,"1.0");
        auto w   = p.addString("w","weight",false,"1.0");
        auto y   = p.addString("y","y-variable",true);
        auto yb  = p.addVFloat("yb","y-variable binning",true);
        auto ylb = p.addBool("ylb","y-variable list binning (explicit bin edges)");
        p.parse(arguments);

        if((*xlb && xb-> size()  < 2) || (!*xlb &&  xb->size() != 3)
           || (*ylb && yb-> size()  < 2) || (!*ylb &&  yb->size() != 3)
        )                     throw std::invalid_argument("Analyzer::Analyzer() -> Bad parsing");


        TString title = TString::Format(";%s;%s",x->c_str(),y->c_str());
        double* ybb = &yb->at(0);
        double* xbb = &xb->at(0);
        const int nXb = xb->size() -1;
        const int nYb = yb->size()-1;

        if(*xlb && *ylb)
            histo = plotter.getOrMake2D(n->c_str(),title,nXb,xbb ,nYb, ybb  );
        else if(*xlb && !*ylb)
            histo = plotter.getOrMake2D(n->c_str(),title,nXb,xbb ,(*yb)[0],(*yb)[1],(*yb)[2] );
        else if(!*xlb && *ylb)
            histo = plotter.getOrMake2D(n->c_str(),title,(*xb)[0],(*xb)[1],(*xb)[2],nYb, ybb  );
        else
            histo = plotter.getOrMake2D(n->c_str(),title,(*xb)[0],(*xb)[1],(*xb)[2],(*yb)[0],(*yb)[1],(*yb)[2]  );

        if(*ylb){
            scaleHisto = plotter.getOrMake1D(TString::Format("scale%s",n->c_str()),y->c_str(),nYb,ybb );
            resHisto   = plotter.getOrMake1D(TString::Format("res%s",n->c_str()),y->c_str(),nYb,ybb );
        } else{
            scaleHisto = plotter.getOrMake1D(TString::Format("scale%s",n->c_str()),y->c_str(),(*yb)[0],(*yb)[1],(*yb)[2] );
            resHisto   = plotter.getOrMake1D(TString::Format("res%s",n->c_str()),y->c_str(),(*yb)[0],(*yb)[1],(*yb)[2] );
        }
        tree.getTree()->SetBranchStatus("*",1);
        tree.getTree()->Draw(TString::Format("%s:%s>>+%s",y->c_str(),x->c_str(),n->c_str()),TString::Format("%s*(%s)",w->c_str(),s->c_str()),"goff");
    }


    void loadVariables() override{};
    bool runEvent() override {
        return true;
    }

    void process(std::string outFileName) {
        const unsigned int yMax = histo->GetNbinsY()+1;
        for(unsigned int iY = 1; iY <= yMax ; ++iY){
            auto tmp = histo->ProjectionX("q",iY,iY);
            scaleHisto->SetBinContent(iY,tmp->GetMean());
            scaleHisto->SetBinError(iY,tmp->GetMeanError());
            resHisto->SetBinContent(iY,tmp->GetRMS());
            resHisto->SetBinError(iY,tmp->GetRMSError());
            delete tmp;
        }
        plotter.write(outFileName);
    }

    TString fileOutName;
    HistGetter plotter;
    TH2 * histo;
    TH1 * scaleHisto;
    TH1 * resHisto;
};

#endif

void make2DDetectorParam(std::string fileName, std::string outFileName,std::string arguments){
    Analyzer a(fileName,"treeMaker/Events",arguments);
    a.process(outFileName);
}
