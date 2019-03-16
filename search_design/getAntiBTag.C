
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TPRegexp.h"


using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt)
{
        fjProc->param.wjj_min_CSVSJCat    = BTagging::CSVSJ_INCL   ;
        fjProc->param.wjj_max_CSVSJCat    = BTagging::CSVSJ_INCL   ;
}

    //--------------------------------------------------------------------------------------------------
    void loadVariables()  override{
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiNoLepJet",isRealData(),false);  load(reader_jet     );
        reader_puppijetwlep =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false);  load(reader_puppijetwlep     );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }

        checkConfig();
    }


    bool runEvent() override {
        if(! DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!passWjjSel) return false;
        if(!hbbCand) return false;
        if(!wjjCand) return false;

        float wlnuDR = 999.;
        float wwDRP = 999.;


        const auto wlnu = neutrino.p4() + selectedLepton->p4();
        const auto hWW = wlnu + wjjCand->p4();
        wlnuDR = PhysicsUtilities::deltaR(neutrino,*selectedLepton);
        wwDRP = PhysicsUtilities::deltaR( wlnu,*wjjCand) * hWW.pt()/125.0;

        const bool sjL = wjjCSVCat >= BTagging::CSVSJ_LF;
        const bool sjM = wjjCSVCat >= BTagging::CSVSJ_MF;


        if(wlnuDR >= 3.2) return false;
        if(wwDRP >= 2.0) return false;
        const bool hBBTM = (hbbMass > 105 && hbbMass < 145);


        auto docat =[&](TH1* h, int& idx, int nL, int nM){
            h->Fill(idx,weight);
            if(!sjL) h->Fill(idx+1,weight);
            if(!sjM) h->Fill(idx+2,weight);
            if(!sjL && nL == 0 ) h->Fill(idx+3,weight);
            if(!sjL && nM == 0 ) h->Fill(idx+4,weight);
            if(!sjM && nL == 0 ) h->Fill(idx+5,weight);
            if(!sjM && nM == 0 ) h->Fill(idx+6,weight);
            if(nL == 0 ) h->Fill(idx+7,weight);
            if(nM == 0 ) h->Fill(idx+8,weight);
            idx +=9;
        };

        auto makeStats=[&](const TString& prefix, int nBL, int nBM ){
            auto * tots = plotter.getOrMake1DPre(prefix,"evtCounts","N. of events",100,-0.5,99.5);
            if(!passHbbSel) return;
            int idx = 0;
            docat(tots,idx,nBL,nBM);
            if(hbbCSVCat == BTagging::CSVSJ_MF) docat(tots,idx,nBL,nBM); else idx+=9;
            if(hbbCSVCat == BTagging::CSVSJ_ML) docat(tots,idx,nBL,nBM); else idx+=9;
            if(hbbCSVCat == BTagging::CSVSJ_MM) docat(tots,idx,nBL,nBM); else idx+=9;
            if(hBBTM){
            docat(tots,idx,nBL,nBM);
            if(hbbCSVCat == BTagging::CSVSJ_MF) docat(tots,idx,nBL,nBM); else idx+=9;
            if(hbbCSVCat == BTagging::CSVSJ_ML) docat(tots,idx,nBL,nBM); else idx+=9;
            if(hbbCSVCat == BTagging::CSVSJ_MM) docat(tots,idx,nBL,nBM); else idx+=9;
            }
        };


        auto pltSet = [&](const TString& prefix, int nL, int nM){
            plotter.getOrMake1DPre(prefix +"_aBI","hbb_mass",";Hbb mass [GeV]; N. events / bin width",50,0.,250)->Fill(hbbMass,weight);

            plotter.getOrMake1DPre(prefix +"_aBI","hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjL)              plotter.getOrMake1DPre(prefix+"_sjL"   ,"hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjM)              plotter.getOrMake1DPre(prefix+"_sjM"   ,"hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjL && nL == 0 )  plotter.getOrMake1DPre(prefix+"_sjL_jL","hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjL && nM == 0 )  plotter.getOrMake1DPre(prefix+"_sjL_jM","hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjM && nL == 0 )  plotter.getOrMake1DPre(prefix+"_sjL_jL","hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(!sjM && nM == 0 )  plotter.getOrMake1DPre(prefix+"_sjL_jM","hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(nL == 0 )          plotter.getOrMake1DPre(prefix+"_jL"    ,"hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
            if(nM == 0 )          plotter.getOrMake1DPre(prefix+"_jM"    ,"hh_mass",";HH mass [TeV]; N. events / bin width",100,0.,5.)->Fill(hh.mass() / 1000.,weight);
        };
        auto pltRecoSet = [&](const TString& prefix, int nBL, int nBM ){
            pltSet(prefix,nBL,nBM);
            if(passHbbSel) pltSet(prefix +"_LI",nBL,nBM);
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MF) pltSet(prefix +"_L",nBL,nBM);
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_ML) pltSet(prefix +"_M",nBL,nBM);
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MM) pltSet(prefix +"_T",nBL,nBM);

            if(hBBTM){
                pltSet(prefix+"_hbbTM",nBL,nBM);
                if(passHbbSel) pltSet(prefix +"_LI_hbbTM",nBL,nBM);
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MF) pltSet(prefix +"_L_hbbTM",nBL,nBM);
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_ML) pltSet(prefix +"_M_hbbTM",nBL,nBM);
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MM) pltSet(prefix +"_T_hbbTM",nBL,nBM);
            }

        };


        auto doTaggingSet = [&](const TString& prefix,  std::vector<const Jet*> inputJets){
            const int nBL    = inputJets.size();
            const int nBM    = PhysicsUtilities::selObjsD(inputJets,[](const Jet* j){return BTagging::isMediumCSVTagged(*j);} ).size();
            makeStats(prefix+"_emu",nBL,nBM);
            if(selectedLepton->isMuon()) makeStats(prefix+"_mu",nBL,nBM);
            else makeStats(prefix+"_e",nBL,nBM);

            pltRecoSet(prefix + "_emu",nBL,nBM);
            if(selectedLepton->isMuon()) pltRecoSet(prefix+"_mu",nBL,nBM);
            else pltRecoSet(prefix+"_e",nBL,nBM);
        };



        auto goodLooseJet_tID = [&] (const Jet* j) -> bool {
            if(!j->passTightID()) return false;
            if( PhysicsUtilities::deltaR2(*j,*hbbCand ) < 1.2*1.2) return false;
            if(!BTagging::isLooseCSVTagged(*j)) return false;
            return true;
        };

        auto goodLooseJet_lID = [&] (const Jet* j) -> bool {
            if(!j->passLooseID()) return false;
            if( PhysicsUtilities::deltaR2(*j,*hbbCand ) < 1.2*1.2) return false;
            if(!BTagging::isLooseCSVTagged(*j)) return false;
            return true;
        };

        auto goodLooseJet_nID = [&] (const Jet* j) -> bool {
            if( PhysicsUtilities::deltaR2(*j,*hbbCand ) < 1.2*1.2) return false;
            if(!BTagging::isLooseCSVTagged(*j)) return false;
            return true;
        };

        doTaggingSet(smpName+"_pnl_20_tID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_pnl_20_lID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_pnl_20_nID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4,goodLooseJet_nID));
        doTaggingSet(smpName+"_pnl_30_tID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_pnl_30_lID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_pnl_30_nID"  ,PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,goodLooseJet_nID));

        doTaggingSet(smpName+"_pwl_20_tID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,20,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_pwl_20_lID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,20,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_pwl_20_nID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,20,2.4,goodLooseJet_nID));
        doTaggingSet(smpName+"_pwl_30_tID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,30,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_pwl_30_lID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,30,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_pwl_30_nID"  ,PhysicsUtilities::selObjsMom(reader_puppijetwlep->jets,30,2.4,goodLooseJet_nID));

        doTaggingSet(smpName+"_swl_20_tID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_swl_20_lID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_swl_20_nID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20,2.4,goodLooseJet_nID));
        doTaggingSet(smpName+"_swl_30_tID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,30,2.4,goodLooseJet_tID));
        doTaggingSet(smpName+"_swl_30_lID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,30,2.4,goodLooseJet_lID));
        doTaggingSet(smpName+"_swl_30_nID"  ,PhysicsUtilities::selObjsMom(reader_jetwlep->jets,30,2.4,goodLooseJet_nID));

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::shared_ptr<JetReader        > reader_puppijetwlep  ;



};

#endif

void getAntiBTag(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getAntiBTag(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
