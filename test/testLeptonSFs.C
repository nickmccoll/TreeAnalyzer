
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "DataFormats/interface/Electron.h"
 #include "DataFormats/interface/Muon.h"



using namespace TAna;
using namespace std;
using namespace TAna::CorrHelp;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        turnOffCorr(CORR_XSEC);
        turnOffCorr(CORR_PU  );
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        cout <<"---------------------------------------------------"<<endl;
        if(reader_event->process == FillerConstants::SIGNAL){
            cout <<"HH decays: "<<endl;
            cout << "["<< diHiggsEvt.type; if(diHiggsEvt.w1_d1) cout << *diHiggsEvt.w1_d1; cout <<"] ";
        }
        cout <<std::endl;
        cout <<"Boson decays: ";
        for(const auto& d : smDecayEvt.bosonDecays ){
            cout << "["<< d.type; if(d.dau1) cout << *d.dau1; cout <<"] ";
        }
        cout <<std::endl;
        cout <<"Top decays: ";
        for(const auto& d : smDecayEvt.topDecays ){
            cout << "["<< d.type; if(d.W_decay.dau1) cout << *d.W_decay.dau1; cout <<"] ";
        }
        cout <<std::endl;

        const auto& pM = leptonSFProc->getPromptMuons();
        const auto& pE = leptonSFProc->getPromptElectrons();
        std::vector<bool> isPrompt(selectedLeptons.size());
        std::vector<bool> isFilledM(pM.size(),false);
        std::vector<bool> isFilledE(pE.size(),false);

        for(unsigned int iS = 0; iS < selectedLeptons.size(); ++iS){
            const auto* sl = selectedLeptons[iS];
            if(sl->isMuon()){
                for(unsigned int iP = 0; iP < pM.size(); ++iP){
                    if(sl->index() == pM[iP]->index()){
                        isFilledM[iP]= true;
                        isPrompt[iS] = true;
                        break;
                    }
                }
            } else {
                for(unsigned int iP = 0; iP < pE.size(); ++iP){
                    if(sl->index() == pE[iP]->index()){
                        isFilledE[iP]= true;
                        isPrompt[iS] = true;
                        break;
                    }
                }
            }
        }
        for(auto b : isFilledM) if(!b) std:: cout<< std::endl << "BAD M!!!!!!!"<<std::endl;
        for(auto b : isFilledE) if(!b) std:: cout<< std::endl << "BAD E!!!!!!!"<<std::endl;



        cout <<"Selected leptons: ";
        for(unsigned int iS = 0; iS < selectedLeptons.size(); ++iS){
            const auto* sl = selectedLeptons[iS];
            cout << "["<< (sl->isMuon() ? "M" : "E") << (isPrompt[iS] ? "P" : "F")  << *sl <<"] ";
        }
        cout << endl;

        cout <<"Weights: ";
        cout << "T: "<<weight <<" L: "<< leptonSFProc->getSF();
        cout << endl;
        cout <<"E ("<<leptonSFProc->getElectronSF(DOWN,NONE,NONE) <<","<< leptonSFProc->getElectronSF(NOMINAL,NONE,NONE)<<","<< leptonSFProc->getElectronSF(UP,NONE,NONE)<<") "
              <<" ("<<leptonSFProc->getElectronSF(NONE,DOWN,NONE) <<","<< leptonSFProc->getElectronSF(NONE,NOMINAL,NONE)<<","<< leptonSFProc->getElectronSF(NONE,UP,NONE)<<") "
              <<" ("<<leptonSFProc->getElectronSF(NONE,NONE,DOWN) <<","<< leptonSFProc->getElectronSF(NONE,NONE,NOMINAL)<<","<< leptonSFProc->getElectronSF(NONE,NONE,UP)<<") "<<endl;
        cout <<"M  ("<<leptonSFProc->getMuonSF(DOWN,NONE,NONE) <<","<< leptonSFProc->getMuonSF(NOMINAL,NONE,NONE)<<","<< leptonSFProc->getMuonSF(UP,NONE,NONE)<<") "
              <<" ("<<leptonSFProc->getMuonSF(NONE,DOWN,NONE) <<","<< leptonSFProc->getMuonSF(NONE,NOMINAL,NONE)<<","<< leptonSFProc->getMuonSF(NONE,UP,NONE)<<") "
              <<" ("<<leptonSFProc->getMuonSF(NONE,NONE,DOWN) <<","<< leptonSFProc->getMuonSF(NONE,NONE,NOMINAL)<<","<< leptonSFProc->getMuonSF(NONE,NONE,UP)<<") "<<endl;

        for(unsigned int iS = 0; iS < selectedLeptons.size(); ++iS){
            if(!isPrompt[iS]) ParticleInfo::printGenInfo(reader_genpart->genParticles);
        }

        return true;
    }




    void write(TString fileName){ plotter.write(fileName);}

    HistGetter plotter;

};

#endif

void testLeptonSFs(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze(1000,1000);
    a.write(outFileName);
}
void testLeptonSFs(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
