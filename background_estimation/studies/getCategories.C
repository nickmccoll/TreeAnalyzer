
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../predTools/BETreeAnalyzer.h"

class Analyzer : public BETreeAnalyzer {
public:
    Analyzer(std::string fileName, std::string treeName, int treeInt) : BETreeAnalyzer(fileName,treeName,treeInt){
    }



    int getCategory(const double tauCut=0.55){
        int cat=0;

        if(hbbCSVCat==4 && wjjMass <  60 &&  wjjTau2o1 >= tauCut ) cat = 0;
        if(hbbCSVCat==4 && wjjMass >= 60 &&  wjjTau2o1 >= tauCut ) cat = 1;
        if(hbbCSVCat==4 && wjjMass <  60 &&  wjjTau2o1 <  tauCut ) cat = 2;
        if(hbbCSVCat==4 && wjjMass >= 60 &&  wjjTau2o1 <  tauCut ) cat = 3;
        if(hbbCSVCat==5 && wjjMass <  60 &&  wjjTau2o1 >= tauCut ) cat = 4;
        if(hbbCSVCat==5 && wjjMass >= 60 &&  wjjTau2o1 >= tauCut ) cat = 5;
        if(hbbCSVCat==5 && wjjMass <  60 &&  wjjTau2o1 <  tauCut ) cat = 6;
        if(hbbCSVCat==5 && wjjMass >= 60 &&  wjjTau2o1 <  tauCut ) cat = 7;
        if(hbbCSVCat==6 && wjjMass <  60 &&  wjjTau2o1 >= tauCut ) cat = 8;
        if(hbbCSVCat==6 && wjjMass >= 60 &&  wjjTau2o1 >= tauCut ) cat = 9;
        if(hbbCSVCat==6 && wjjMass <  60 &&  wjjTau2o1 <  tauCut ) cat = 10;
        if(hbbCSVCat==6 && wjjMass >= 60 &&  wjjTau2o1 <  tauCut ) cat = 11;

        return cat;
    }


    bool runEvent() override {
        if(!BETreeAnalyzer::runEvent()) return false;

        std::string prefix = isSignal() ? "signal" : "bkg";
        if(process==8) return false;
        if(hbbCSVCat<4) return false;
        if(hbbNSJs!=2||wjjNSJs!=2)return false;
        if(nAK4Btags!=0) return false;
        if(wlnuDR>=3.2) return false;
        if(wwDM>=2) return false;
        if(wjjMass<=10) return false;
        if(hbbMass<=100) return false;
        if(hbbMass>=150) return false;
        if(hhMass<=700) return false;

        if(isSignal()){
            plotter.getOrMake2DPre("signal","categories",";#it{m}_{X} [GeV];cat",40,550,4550,12,-0.5,11.5)->Fill(signal_mass,getCategory(),weight);
        } else{

            plotter.getOrMake2DPre("bkg","categories",";#it{m}_{HH} [GeV];cat",6,bins,12,-0.5,11.5)->Fill(hhMass,getCategory(),weight);
        }

        for(double tc = 0.35; tc <= 1.05; tc += 0.05  ){
            const auto cat = catNs[getCategory(tc)];
            std::string prefix = isSignal() ? "signal_" : "bkg_";
            prefix += cat;

            if(hhMass > 900 && hhMass < 1100 && (!isSignal() || signal_mass==1000))
                plotter.getOrMake1DPre(prefix.c_str(),"m1000_tauEff",";tau cut",15,0.325,1.075)->Fill(tc,weight);

            if(hhMass > 2000 && hhMass < 3000 && (!isSignal() || signal_mass==2500))
                plotter.getOrMake1DPre(prefix.c_str(),"m2500_tauEff",";tau cut",15,0.325,1.075)->Fill(tc,weight);

        }


        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    const double bins[7] = {700,900,1100,1500,2000,3000,5000};
    const std::string catNs [12] = {"bL_wML_wTH","bL_wMH_wTH","bL_wML_wTL","bL_wMH_wTL",
                              "bM_wML_wTH","bM_wMH_wTH","bM_wML_wTL","bM_wMH_wTL",
                              "bT_wML_wTH","bT_wMH_wTH","bT_wML_wTL","bT_wMH_wTL" };

    HistGetter plotter;

};

#endif
void getCategories(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);

}
void getCategories(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
