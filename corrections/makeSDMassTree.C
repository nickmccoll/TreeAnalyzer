
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
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
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "Processors/Variables/interface/FatJetSelection.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        resetCorr();
        turnOnCorr(CORR_PU);
        fromWCorr.reset(new SoftDropMassScaleFactorsFromW (dataDirectory));
    }

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new ManualCopierEventAnalyzer();}

    virtual void bookOutputVariables() override {
        i_puWeight        =  outTree->add<float>  ("puWeight"           ,"F",0);
        i_csvCat          =  outTree->add<size8>  ("csvCat"         ,"b",0);
        i_bosonMass       =  outTree->add<float>  ("bosonMass"          ,"F",0);
        i_bosonPT         =  outTree->add<float>  ("bosonPT"            ,"F",0);
        i_bosonETA        =  outTree->add<float>  ("bosonETA"            ,"F",0);
        i_drToJet         =  outTree->add<float>  ("drToJet"            ,"F",0);
        i_maxSJDR          =  outTree->add<float> ("maxSJDR"             ,"F",0);
        i_maxSJETA        =  outTree->add<float>  ("maxSJETA"            ,"F",0);
        i_minSJPT         =  outTree->add<float>  ("minSJPT"             ,"F",0);
        i_genPT           =  outTree->add<float>  ("genPT"              ,"F",0);
        i_genETA          =  outTree->add<float>  ("genETA"             ,"F",0);
        i_jetPT           =  outTree->add<float>  ("jetPT"              ,"F",0);
        i_jetCorrPT           =  outTree->add<float>  ("jetCorrPT"              ,"F",0);
        i_jetETA          =  outTree->add<float>  ("jetETA"             ,"F",0);
        i_jetID           =  outTree->add<size8>  ("jetID"              ,"b",0);
        i_jetMass         =  outTree->add<float>  ("jetMass"            ,"F",0);
        i_jetRawSDMass    =  outTree->add<float>  ("jetRawSDMass"       ,"F",0);
        i_jetL23CorrSDMass=  outTree->add<float>  ("jetL23CorrSDMass"   ,"F",0);
        i_jetCorrSDMass   =  outTree->add<float>  ("jetCorrSDMass"      ,"F",0);
        i_jetWCorrSDMass   =  outTree->add<float>  ("jetWCorrSDMass"      ,"F",0);
    }

    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }

        checkConfig();
    }



    void fillBoson(const FatJet* fj, const GenJet* gj) {
        double maxDR2 = -1;

        double drQ = std::sqrt( std::max(PhysicsUtilities::deltaR2(*diHiggsEvt.b1,*fj),PhysicsUtilities::deltaR2(*diHiggsEvt.b2,*fj)));
        float maxQETA = -1;
        float minQPT  = -1;
        for(const auto& sj : fj->subJets() ) {
            double dr2 = PhysicsUtilities::deltaR2(sj,*fj);
            if(maxDR2 < 0){
                maxDR2 = dr2;
                maxQETA = sj.absEta();
                minQPT  = sj.pt();
            } else {
                if(dr2 > maxDR2 ) maxDR2 = dr2;
                if(sj.absEta() > maxQETA ) maxQETA = sj.absEta();
                if(sj.pt() < minQPT ) minQPT = sj.pt();
            }
        }
        if(maxDR2 < 0) return;
        if(fj->nSubJets() != 2) return;

        resetOutData();
        outTree->fill(i_csvCat      , size8(BTagging::getCSVSJCat(fj->subJets(), fjProc->param.sj_minBTagPT, fjProc->param.sj_maxBTagETA)));
        outTree->fill(i_puWeight        , float( puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL)));
        outTree->fill(i_bosonMass       , float(diHiggsEvt.hbb->mass()));
        outTree->fill(i_bosonPT         , float(diHiggsEvt.hbb->pt()));
        outTree->fill(i_bosonETA        , float(diHiggsEvt.hbb->eta()));
        outTree->fill(i_drToJet         , float(PhysicsUtilities::deltaR(*diHiggsEvt.hbb,*fj)));
        outTree->fill(i_maxSJDR          , float(drQ));
        outTree->fill(i_jetPT           , float(fj->pt()));
        outTree->fill(i_jetCorrPT       , float( hbbFJSFProc->getJEC(fj)  * fj->pt()));
        outTree->fill(i_jetETA          , float(fj->eta()));
        outTree->fill(i_jetID           , fj->jetID() );
        outTree->fill(i_maxSJETA        , float(maxQETA));
        outTree->fill(i_minSJPT         , float(minQPT));
        outTree->fill(i_genPT           , float(gj->pt()));
        outTree->fill(i_genETA          , float(gj->eta()));
        outTree->fill(i_jetMass         , float(fj->mass()));
        outTree->fill(i_jetRawSDMass    , float(fj->rawSdMom().mass()));
        outTree->fill(i_jetL23CorrSDMass, float(fj->sdMom().mass()));
        outTree->fill(i_jetCorrSDMass   , float(hbbFJSFProc->getCorrSDMass(fj)));
        outTree->fill(i_jetWCorrSDMass  , float(fromWCorr->getCorrSDMass(fj)));
        fillOutTree();
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;

        auto genjetCands = PhysicsUtilities::selObjsMom(reader_fatjet->genJets,10,2.4);

        double dr = 999;
        int idx = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.hbb,genjetCands,dr);
        if(idx < 0) return false;

        if(PhysicsUtilities::deltaR(*diHiggsEvt.w1_d1,*genjetCands[idx]) < 1.0 ) return false;
        if(PhysicsUtilities::deltaR(*diHiggsEvt.w2_d1,*genjetCands[idx]) < 1.0 ) return false;
        if(PhysicsUtilities::deltaR(*diHiggsEvt.w2_d2,*genjetCands[idx]) < 1.0 ) return false;
        dr = 999;
        int fjIDX = PhysicsUtilities::findNearestDR(*genjetCands[idx],reader_fatjet->jets,dr,0.8);
        if(fjIDX < 0) return false;
        if(reader_fatjet->jets[fjIDX].pt() < 50 ) return false;
        if(reader_fatjet->jets[fjIDX].absEta() >= 2.4 ) return false;

        fillBoson(&reader_fatjet->jets[fjIDX],genjetCands[idx]);

        return true;

    }


    size i_puWeight        = 0;
    size i_csvCat          = 0;
    size i_bosonMass       = 0;
    size i_bosonPT         = 0;
    size i_bosonETA        = 0;
    size i_drToJet         = 0;
    size i_maxSJDR          = 0;
    size i_jetPT           = 0;
    size i_jetCorrPT           = 0;
    size i_jetETA          = 0;
    size i_jetID           = 0;
    size i_maxSJETA        = 0;
    size i_minSJPT         = 0;
    size i_genPT           = 0;
    size i_genETA          = 0;
    size i_jetMass         = 0;
    size i_jetRawSDMass    = 0;
    size i_jetL23CorrSDMass= 0;
    size i_jetCorrSDMass   = 0;
    size i_jetWCorrSDMass   = 0;


    std::unique_ptr<SoftDropMassScaleFactorsFromW>    fromWCorr ;


};

#endif

void makeSDMassTree(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}
void makeSDMassTree(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();

}
