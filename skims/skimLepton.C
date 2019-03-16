
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;

class FirstAnalyzer : public BaseTreeAnalyzer {
public:

    FirstAnalyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        leptonProc .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProc);
    }

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {
        if(!isRealData())
            i_normWeight             =  outTree->add<float>  ("event","normWeight"  ,"F",0);
        i_ht                     =  outTree->add<float>  ("skim" ,"ht"          ,"F",0);
        i_selLep_pt              =  outTree->add<float>  ("skim" ,"selLep_pt"   ,"F",0);
        i_selLep_eta             =  outTree->add<float>  ("skim" ,"selLep_eta"  ,"F",0);
        i_selLep_phi             =  outTree->add<float>  ("skim" ,"selLep_phi"  ,"F",0);
        i_selLep_muon            =  outTree->add<size8>  ("skim" ,"selLep_muon" ,"b",0);
    }

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );

        if(treeType == TREE_OTHER){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }

        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData(),false);            load(reader_jetwlep     );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false);  load(reader_fatjet  );
    }

    bool runEvent() override {
        auto leptons = leptonProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
//        std::cout <<"\nE ";
//        for(const auto* l : leptons) cout << l->isMuon()<<(*l)<<" ";

        if(leptons.size() != 1.0) return false;


        float normEventWeight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        auto jets = PhysicsUtilities::selObjsMom(reader_jetwlep->jets,30);
        const float ht = JetKinematics::ht(jets);

        if(treeType != TREE_OTHER && ht < 500) return false;

        if(!isRealData())
            outTree->fill(i_normWeight ,normEventWeight);
        outTree->fill(i_ht         ,ht);
        outTree->fill(i_selLep_pt  ,float(leptons.front()->pt()));
        outTree->fill(i_selLep_eta ,float(leptons.front()->eta()));
        outTree->fill(i_selLep_phi ,float(leptons.front()->phi()));
        outTree->fill(i_selLep_muon,size8(leptons.front()->isMuon()));
        return true;

    }

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    std::shared_ptr<JetReader        > reader_jet      ;
    std::shared_ptr<JetReader        > reader_jetwlep  ;

    std::unique_ptr<LeptonProcessor> leptonProc ;

    size i_normWeight = 0;
    size i_ht         = 0;
    size i_selLep_pt  = 0;
    size i_selLep_eta = 0;
    size i_selLep_phi = 0;
    size i_selLep_muon= 0;

};

class SecondAnalyzer : public BaseTreeAnalyzer {
public:

    SecondAnalyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
    }

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {
    }

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );

        if(treeType == TREE_OTHER){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false);  load(reader_fatjet  );

       setBranchAddress("skim" ,"ht"         ,   &ht                  ,true);
       setBranchAddress("skim" ,"selLep_pt"  ,   &selLep_pt           ,true);
       setBranchAddress("skim" ,"selLep_eta" ,   &selLep_eta          ,true);
       setBranchAddress("skim" ,"selLep_phi" ,   &selLep_phi          ,true);
       setBranchAddress("skim" ,"selLep_muon",   &selLep_muon         ,true);

    }

    bool runEvent() override {
        return true;

    }

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    std::shared_ptr<JetReader        > reader_jet      ;

    float   ht         =0;
    float   selLep_pt  =0;
    float   selLep_eta =0;
    float   selLep_phi =0;
    size8   selLep_muon=0;
};

#endif

void doSecondAnalyzer(std::string fileName, int treeInt, std::string outFileName){
    SecondAnalyzer a(fileName,"treeMaker/Events",treeInt);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_LOADED);
    a.analyze();
}

void skimLepton(std::string fileName, int treeInt, std::string outFileName){
    TString tempFilename = outFileName;
    tempFilename.ReplaceAll(".root", "");
    tempFilename += "_temp.root";

    FirstAnalyzer* a = new FirstAnalyzer(fileName,"treeMaker/Events",treeInt);
    a->initializeTreeCopy(tempFilename.Data(),BaseTreeAnalyzer::COPY_LOADED);
    a->analyze();
    delete a;

    doSecondAnalyzer(tempFilename.Data(),treeInt,outFileName);
    gSystem->Exec(TString::Format("rm %s", tempFilename.Data()));
}
void skimLepton(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    TString tempFilename = outFileName;
    tempFilename.ReplaceAll(".root", "");
    tempFilename += "_temp.root";

    FirstAnalyzer* a = new FirstAnalyzer(fileName,"treeMaker/Events",treeInt);
    a->setSampleInfo(xSec,numEvent);
    a->initializeTreeCopy(tempFilename.Data(),BaseTreeAnalyzer::COPY_LOADED);
    a->analyze();
    delete a;

    doSecondAnalyzer(tempFilename.Data(),treeInt,outFileName);
    gSystem->Exec(TString::Format("rm %s", tempFilename.Data()));
}
