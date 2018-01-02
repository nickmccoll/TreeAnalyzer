
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParParser.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "Utilities/HiggsCombineImport/interface/VerticalInterpHistPdf.h"
#include <string.h>
#include <regex>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"





class fit2DTemplateAnalyzer {
public:
    typedef std::vector<std::pair<std::string,std::string> > SystNames;

    fit2DTemplateAnalyzer(std::string outFileName,std::string arguments )
{


        ParParser p;
        auto fTN = p.addString("fT","template file name",true);
        auto nT  = p.addString("nT","template histogram base name",true);
        auto s   = p.addString("s" ,"Comma separated list of systematics",true);
        auto fHN = p.addString("fH","Fitting histogram file name",true);
        auto nH  = p.addString("nH","fitting histogram name",true);
        p.parse(arguments);

        std::vector<std::string> systList = getList(*s);

        fT =  TObjectHelper::getFile(*fTN);
        fH =  TObjectHelper::getFile(*fHN);
        auto h = TObjectHelper::getObject<TH2F>(fT,*nT);
        auto * xAxis = h->GetXaxis();
        auto * yAxis = h->GetYaxis();

        RooWorkspace w("w",false);
        RooArgSet varset;
        RooArgList varlist;
        w.factory("x[0,10000]");
        w.factory("y[0,10000]");
        w.var("x")->setBinning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax()));
        w.var("y")->setBinning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax()));
        varset.add(*w.var("x"));
        varset.add(*w.var("y"));
        varlist.add(*w.var("x"));
        varlist.add(*w.var("y"));

        RooDataHist rH("nominalHist","nominalHist",varlist,&*h);
        RooHistPdf rP("nominalPDF","nominalPDF",varlist,rH);
        w.import(rH);
        w.import(rP);
        RooArgList coeffList;
        RooArgList pdfList(*w.pdf("nominalPDF"));

        const std::vector<std::string> systVar = {"Up" , "Down"};
        for(const auto& syst : systList){
            w.factory((syst+"[-5,5]").c_str());
            coeffList.add(*w.var(syst.c_str()));
            for(const auto& var : systVar){
                auto sH = TObjectHelper::getObject<TH2F>(fT,*nT+"_"+syst+var);
                RooDataHist rSH((syst+var+"Hist").c_str(),(syst+var+"Hist").c_str(),varlist,&*sH);
                RooHistPdf Spdf((syst+var+"PDF").c_str(),(syst+var+"PDF").c_str(),varset,rSH,0);
                w.import(rSH);
                w.import(Spdf);
                pdfList.add(*w.pdf((syst+var+"PDF").c_str()));
            }

        }

        FastVerticalInterpHistPdf2D interpPDF("interpPDF","interpPDF",*w.var("x"),*w.var("y"),false,pdfList,coeffList);
        w.import(interpPDF);

        auto inH =TObjectHelper::getObject<TH2F>(fH,*nH);
        RooDataHist fitDataHist((*nH+"DH").c_str(),(*nH+"DH").c_str(),varlist,&*inH);
        w.import(fitDataHist);
        w.pdf("interpPDF")->fitTo(*w.data((*nH+"DH").c_str()),RooFit::SumW2Error(kTRUE));

                auto fitHist = w.pdf("interpPDF")->createHistogram(nT->c_str(),
                        *w.var("x"),RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())),
                        RooFit::YVar(*w.var("y"),RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax())))) ;
//                fitHist->SetName(nT->c_str());
                std::cout << yAxis->GetNbins() <<" "<< yAxis->GetXmin()<<" "<<yAxis->GetXmax() <<" "<< yAxis->GetBinLowEdge(50) <<" "<<fitHist->GetYaxis()->GetBinLowEdge(50) <<std::endl;
                for(const auto& syst : systList){
                    std::cout << syst <<" -> "<< w.var(syst.c_str())->getVal()<<std::endl;
                }

                plotter.add1D(fitHist);
//        interpretPDF.fitTo(fitDataHist);
//        auto * xAxis = pdfHistos.front()->GetXaxis();
//        auto * yAxis = pdfHistos.front()->GetYaxis();
//        plotter.add1D(interpretPDF.createHistogram(nT->c_str(),
//                plottingVars[0],RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())),
//                RooFit::YVar(plottingVars[1],RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax()))))) ;



//        plottingVars.emplace_back("x","x",1);
//        plottingVars.emplace_back("y","y",1);
////        RooArgSet varset;
//        RooArgList varlist;
//        for(auto& v : plottingVars){
////            varset.add(v);
//            varlist.add(v);
//                        varset.add(self.w.var(var))
//                        varlist.add(self.w.var(var))
//        }
//
//        pdfHistos.emplace_back(TObjectHelper::getObject<TH2F>(fT,*nT));
//        dataHists.emplace_back("nominalHist","nominalHist",varlist,&*pdfHistos.back());
//        pdfs.emplace_back("nominalPDF","nominalPDF",varlist,dataHists.back());
//
//        RooArgList coeffList;
//        RooArgList pdfList;
//        std::cout << pdfList.add(pdfs[pdfList.getSize()]) <<std::endl;
//        const std::vector<std::string> systVar = {"Up" , "Down"};
//        for(const auto& syst : systList){
//            systVars.emplace_back(syst.c_str(),syst.c_str(),1);
//            coeffList.add(systVars.back());
//            for(const auto& var : systVar){
//                pdfHistos.emplace_back(TObjectHelper::getObject<TH2F>(fT,*nT+"_"+syst+var));
//                dataHists.emplace_back((syst+var+"Hist").c_str(),(syst+var+"Hist").c_str(),varlist,&*pdfHistos.back());
//                pdfs.emplace_back((syst+var+"PDF").c_str(),(syst+var+"PDF").c_str(),varlist,dataHists.back());
////                std::cout << pdfList.add(pdfs[pdfList.getSize()]) <<std::endl;
//                for(unsigned int iP = 0; iP < pdfList.getSize(); ++iP){
//                    pdfList.at(iP)->Print();
//
//                }
//            }
//        }
//
//        for(unsigned int iP = 0; iP < pdfList.getSize(); ++iP){
//            std::cout <<"sP"<< std::endl;
//            pdfHistos.at(iP)->Print();
//            dataHists.at(iP).Print();
//            pdfs.at(iP).Print();
//            pdfList.at(iP)->Print();
//            std::cout << "eP"<<std::endl;
//        }
//        for(unsigned int iP = 0; iP < coeffList.getSize(); ++iP){
//            std::cout << "sC"<<std::endl;
//            coeffList.at(iP)->Print();
//            std::cout << "eC"<<std::endl;
//        }
//
//        FastVerticalInterpHistPdf2D interpretPDF("interpPDF","interpPDF",plottingVars[0],plottingVars[1],false,pdfList,coeffList);
//
//        auto inH =TObjectHelper::getObject<TH2F>(fH,*nH);
//        RooDataHist fitDataHist((*nH+"DH").c_str(),(*nH+"DH").c_str(),varlist,&*inH);
//        interpretPDF.fitTo(fitDataHist);
//        auto * xAxis = pdfHistos.front()->GetXaxis();
//        auto * yAxis = pdfHistos.front()->GetYaxis();
//        plotter.add1D(interpretPDF.createHistogram(nT->c_str(),
//                plottingVars[0],RooFit::Binning(RooBinning (xAxis->GetNbins(),xAxis->GetXmin(),xAxis->GetXmax())),
//                RooFit::YVar(plottingVars[1],RooFit::Binning(RooBinning (yAxis->GetNbins(),yAxis->GetXmin(),yAxis->GetXmax()))))) ;
//
//
//

        plotter.write(outFileName);
        fT->Close();
        fH->Close();
}



//        for systval in systematics:
//            splitted=systval.split(':')
//            systName=splitted[1]
//            syst=splitted[0]
//            self.w.factory(systName+"[-1,1]")
//            coeffList.add(self.w.var(systName))
//
//            for variation in ["Up","Down"]:
//                histo=FR.Get(histoname+"_"+syst+variation)
//                print 'loaded',histoname+"_"+syst+variation
//                histName="_".join([name+"_"+syst+variation+"HIST",tag])
//                roohist = ROOT.RooDataHist(histName,histName,varlist,histo)
//
//                pdfName="_".join([name+"_"+syst+variation,self.tag])
//                pdf=ROOT.RooHistPdf(pdfName,pdfName,varset,roohist,order)
//
//                getattr(self.w,'import')(roohist,ROOT.RooFit.Rename(histName))
//                getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(pdfName))
//                pdfList.add(self.w.pdf(pdfName))
//
//        pdfName="_".join([name,self.tag])
//        if len(systematics)>0:
//            if len(observables)==1:
//                total=ROOT.FastVerticalInterpHistPdf(pdfName,pdfName,self.w.var(observables[0]),pdfList, coeffList)
//            elif len(observables)==2:
//                total=ROOT.FastVerticalInterpHistPdf2D(pdfName,pdfName,self.w.var(observables[0]),self.w.var(observables[1]),conditional,pdfList, coeffList)
//            getattr(self.w,'import')(total,ROOT.RooFit.Rename(pdfName))
//
//
    std::vector<std::unique_ptr<TH2F> > pdfHistos;
    std::vector<RooRealVar> plottingVars;
    std::vector<RooRealVar> systVars;
    std::vector<RooDataHist> dataHists;
    std::vector<RooHistPdf> pdfs;


    std::vector<std::string> getList(const std::string& inList){
        std::vector<std::string> systList(std::sregex_token_iterator(inList.begin(), inList.end(), std::regex(","), -1), std::sregex_token_iterator());
        return systList;
    }

    HistGetter plotter;
    TFile * fT = 0;
    TFile * fH = 0;
    std::shared_ptr<std::vector<double>>       xb;
    std::shared_ptr<std::vector<double>>       yb;

};

#endif

void fit2DTemplate(std::string outFileName,std::string arguments){
    fit2DTemplateAnalyzer a(outFileName, arguments);
}
