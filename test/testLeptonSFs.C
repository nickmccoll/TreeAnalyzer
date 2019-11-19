
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "Configuration/interface/FillerConstants.h"
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

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
        : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        cout <<"---------------------------------------------------"<<endl;
        if(isSignal()){
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


        if(lepChan == SINGLELEP){
            cout << "1L OFFLINE "<< endl;
        } else if(lepChan == DILEP){
            cout << "2L OFFLINE "<< endl;
        } else {
            cout << "OTHER OFFLINE "<< endl;
        }

        auto lepProc = &*(lepChan==DILEP ? dileptonSFProc : leptonSFProc);
        auto lepList = (lepChan==DILEP ? selectedDileptons : selectedLeptons);


        const auto& pM = lepProc->getPromptMuons();
        const auto& pE = lepProc->getPromptElectrons();
        std::vector<bool> isPrompt(lepList.size());
        std::vector<bool> isFilledM(pM.size(),false);
        std::vector<bool> isFilledE(pE.size(),false);

        for(unsigned int iS = 0; iS < lepList.size(); ++iS){
            const auto* sl = lepList[iS];
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
        for(unsigned int iS = 0; iS < lepList.size(); ++iS){
            const auto* sl = lepList[iS];
            cout << "["<< (sl->isMuon() ? "M" : "E") << (isPrompt[iS] ? "P" : "F")  << *sl <<
                    (sl->isElectron() ?  ASTypes::flt2Str( ((const Electron*)sl)->sc_act_o_pt()) + ","+ASTypes::flt2Str( ((const Electron*)sl)->sc_dr_act())
            : ASTypes::flt2Str(sl->lepAct_o_pt()) + ","+ASTypes::flt2Str(sl->dRnorm()) )  <<"] ";
        }
        cout << endl;

        cout <<"Weights: ";
        cout << "T: "<<weight <<" L: "<< lepProc->getSF();
        cout << endl;
        cout <<"E ("<<lepProc->getElectronSF(DOWN,NONE,NONE) <<","<< lepProc->getElectronSF(NOMINAL,NONE,NONE)<<","<< lepProc->getElectronSF(UP,NONE,NONE)<<") "
              <<" ("<<lepProc->getElectronSF(NONE,DOWN,NONE) <<","<< lepProc->getElectronSF(NONE,NOMINAL,NONE)<<","<< lepProc->getElectronSF(NONE,UP,NONE)<<") "
//              <<" ("<<lepProc->getElectronSF(NONE,NONE,DOWN) <<","<< lepProc->getElectronSF(NONE,NONE,NOMINAL)<<","<< lepProc->getElectronSF(NONE,NONE,UP)<<") "
              <<endl;
        cout <<"M"
//               <<" ("<<lepProc->getMuonSF(DOWN,NONE,NONE) <<","<< lepProc->getMuonSF(NOMINAL,NONE,NONE)<<","<< lepProc->getMuonSF(UP,NONE,NONE)<<") "
              <<" ("<<lepProc->getMuonSF(NONE,DOWN,NONE) <<","<< lepProc->getMuonSF(NONE,NOMINAL,NONE)<<","<< lepProc->getMuonSF(NONE,UP,NONE)<<") "
//              <<" ("<<lepProc->getMuonSF(NONE,NONE,DOWN) <<","<< lepProc->getMuonSF(NONE,NONE,NOMINAL)<<","<< lepProc->getMuonSF(NONE,NONE,UP)<<") "
              <<endl;

        for(unsigned int iS = 0; iS < lepList.size(); ++iS){
            if(!isPrompt[iS]) ParticleInfo::printGenInfo(reader_genpart->genParticles);
        }

        return true;
    }




    void write(TString fileName){ plotter.write(fileName);}

    HistGetter plotter;

};

#endif

void testLeptonSFs(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
