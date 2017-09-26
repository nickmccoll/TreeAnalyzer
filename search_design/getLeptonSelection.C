
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }


    void addPlots(const TString& prefix, const TString& varname, const TString& selName, float sel, float pt){
        plotter.getOrMake1DPre(prefix,TString::Format("%s_selection",varname.Data()),TString::Format(";%s; arbitrary units",varname.Data()),20,-0.5,19.5 )->Fill(sel,weight);
        plotter.getOrMake1DPre(prefix,TString::Format("%s_%s",selName.Data(),varname.Data()),";lepton p_{T}; arbitrary units",500,0,500 )->Fill(pt,weight);
    }
    void testMuID(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Muon* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<LepSelHelpers::muFunBool> muIDs = { &Muon::passSoftID,&Muon::passLooseID,&Muon::passMedID,&Muon::passMed16ID,&Muon::passTightID,&Muon::passHighPT};
        static const std::vector<TString> muIDNs = {"soft","loose","med","med16","tight","highPT"};

        addPlots(prefix,"muID","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","muID","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedMu","muID","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedMu_passTight","muID","incl",0,signalPT);
        }

        for(unsigned int iS = 0; iS < muIDs.size(); ++iS){
            proc.lepSelParams.mu_getID   = muIDs[iS];
            proc.lepSelParams_dataABCDEF.mu_getID   = muIDs[iS];
            auto leps = proc.getMuons(*reader_event,*reader_muon);
            if(leps.size()){
                addPlots(prefix,"muID",muIDNs[iS],iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","muID",muIDNs[iS],iS+1,isSignal ?signalPT  : leps.front()->pt());
            }

            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedMu","muID",muIDNs[iS],iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedMu_passTight","muID",muIDNs[iS],iS+1,signalPT);
                }
            }
        }

    }

    void testElID(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Electron* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<LepSelHelpers::elFunBool> elIDs = { &Electron::passVetoID_noISO,&Electron::passLooseID_noISO,&Electron::passMedID_noISO,&Electron::passTightID_noISO
                ,&Electron::passHEEPID_noISO,&Electron::passMVA80ID,&Electron::passMVA90ID};
        static const std::vector<TString> elIDNs = {"veto","loose","med","tight","heep","mva80","mva90"};

        addPlots(prefix,"elID","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","elID","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedEl","elID","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedEl_passTight","elID","incl",0,signalPT);
        }

        for(unsigned int iS = 0; iS < elIDs.size(); ++iS){
            proc.lepSelParams.el_getID   = elIDs[iS];
            auto leps = proc.getElectrons(*reader_electron);
            if(leps.size()){
                addPlots(prefix,"elID",elIDNs[iS],iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","elID",elIDNs[iS],iS+1,isSignal ?signalPT  : leps.front()->pt());
            }

            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedEl","elID",elIDNs[iS],iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedEl_passTight","elID",elIDNs[iS],iS+1,signalPT);
                }
            }
        }

    }

    void testMuISO(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Muon* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     mu_ISOVs = {0.1,0.2,0.3,0.4};
        static const std::vector<LepSelHelpers::muFunFloat> mu_ISOTs = {&Muon::miniIso,&Muon::dbRelISO};
        static const std::vector<TString>                   mu_ISOTNs = {"miniIso","relIso"};


        addPlots(prefix,"muISO","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","muISO","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedMu","muISO","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedMu_passTight","muISO","incl",0,signalPT);
        }

        for(unsigned int iT = 0; iT < mu_ISOTNs.size(); ++iT){
            proc.lepSelParams.mu_getISO   = mu_ISOTs[iT];
            proc.lepSelParams_dataABCDEF.mu_getISO   = mu_ISOTs[iT];
            for(unsigned int iS = 0; iS < mu_ISOVs.size(); ++iS){
                const size idx = iT*mu_ISOVs.size() + iS + 1;
                proc.lepSelParams.mu_maxISO   = mu_ISOVs[iS];
                proc.lepSelParams_dataABCDEF.mu_maxISO   = mu_ISOVs[iS];
                TString name = mu_ISOTNs[iT] + TString::Format("_0p%.0f",mu_ISOVs[iS]*10.0 );
                auto leps = proc.getMuons(*reader_event,*reader_muon);
                if(leps.size()){
                    addPlots(prefix,"muISO",name,idx,isSignal ?signalPT  : leps.front()->pt());
                    if(passTight)  addPlots(prefix+"_passTight","muISO",name,idx,isSignal ?signalPT  : leps.front()->pt());
                }
                if(isSignal && signalLep ){
                    bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                    if(found){
                        addPlots(prefix+"_matchedMu","muISO",name,idx,signalPT);
                        if(passTight)  addPlots(prefix+"_matchedMu_passTight","muISO",name,idx,signalPT);
                    }

                }
            }
        }


    }

    void testElISO(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Electron* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     el_ISOVs = {0.1,0.2,0.3,0.4};
        static const std::vector<LepSelHelpers::elFunFloat> el_ISOTs = {&Electron::miniIso,&Electron::eaRelISO};
        static const std::vector<TString>                   el_ISOTNs = {"miniIso","relIso"};


        addPlots(prefix,"elISO","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","elISO","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedEl","elISO","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedEl_passTight","elISO","incl",0,signalPT);
        }

        for(unsigned int iT = 0; iT < el_ISOTNs.size(); ++iT){
            proc.lepSelParams.el_getISO   = el_ISOTs[iT];
            for(unsigned int iS = 0; iS < el_ISOVs.size(); ++iS){
                const size idx = iT*el_ISOVs.size() + iS + 1;
                proc.lepSelParams.el_maxISO   = el_ISOVs[iS];
                TString name = el_ISOTNs[iT] + TString::Format("_0p%.0f",el_ISOVs[iS]*10.0 );
                auto leps = proc.getElectrons(*reader_electron);
                if(leps.size()){
                    addPlots(prefix,"elISO",name,idx,isSignal ?signalPT  : leps.front()->pt());
                    if(passTight)  addPlots(prefix+"_passTight","elISO",name,idx,isSignal ?signalPT  : leps.front()->pt());
                }
                if(isSignal && signalLep ){
                    bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                    if(found){
                        addPlots(prefix+"_matchedEl","elISO",name,idx,signalPT);
                        if(passTight)  addPlots(prefix+"_matchedEl_passTight","elISO",name,idx,signalPT);
                    }

                }
            }
        }
    }

    void testMuD0(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Muon* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     mu_D0s = {0.01,0.02,0.05,0.1,0.2};

        addPlots(prefix,"muD0","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","muD0","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedMu","muD0","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedMu_passTight","muD0","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < mu_D0s.size(); ++iS){
            proc.lepSelParams.mu_maxD0   = mu_D0s[iS];
            proc.lepSelParams_dataABCDEF.mu_maxD0   = mu_D0s[iS];
            TString name = TString::Format("%.2f",mu_D0s[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getMuons(*reader_event,*reader_muon);
            if(leps.size()){
                addPlots(prefix,"muD0",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","muD0",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
            }
            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedMu","muD0",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedMu_passTight","muD0",name,iS+1,signalPT);
                }

            }
        }
    }

    void testMuSIP3D(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Muon* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     mu_S3Ds = {2,3,4,5,6,7,8};

        addPlots(prefix,"muS3","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","muS3","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedMu","muS3","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedMu_passTight","muS3","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < mu_S3Ds.size(); ++iS){
            proc.lepSelParams.mu_maxD0   = 999.9;
            proc.lepSelParams_dataABCDEF.mu_maxD0   = 999.9;
            TString name = TString::Format("%.0f",mu_S3Ds[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getMuons(*reader_event,*reader_muon);
            float maxPT =-1;
            bool foundSig = false;
            for(const auto* l : leps ){
                if(l->sip3D() >= mu_S3Ds[iS] ) continue;
                if(l->pt() > maxPT) maxPT = l->pt();
                if(isSignal && signalLep ) if(l->index() == signalLep->index()) foundSig = true;
            }

            if(maxPT > 0){
                addPlots(prefix,"muS3",name,iS+1,isSignal ?signalPT  : maxPT);
                if(passTight)  addPlots(prefix+"_passTight","muS3",name,iS+1,isSignal ?signalPT  : maxPT);
            }
            if(isSignal && signalLep ){
                if(foundSig){
                    addPlots(prefix+"_matchedMu","muS3",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedMu_passTight","muS3",name,iS+1,signalPT);
                }
            }
        }
    }


    void testMuDZ(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Muon* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     mu_DZs = {0.1,0.2,0.3,0.4,0.5};

        addPlots(prefix,"muDZ","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","muDZ","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedMu","muDZ","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedMu_passTight","muDZ","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < mu_DZs.size(); ++iS){
            proc.lepSelParams.mu_maxDZ   = mu_DZs[iS];
            proc.lepSelParams_dataABCDEF.mu_maxDZ   = mu_DZs[iS];
            TString name = TString::Format("%.1f",mu_DZs[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getMuons(*reader_event,*reader_muon);
            if(leps.size()){
                addPlots(prefix,"muDZ",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","muDZ",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
            }
            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedMu","muDZ",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedMu_passTight","muDZ",name,iS+1,signalPT);
                }

            }
        }
    }
    void testElD0(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Electron* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     el_D0s = {0.01,0.02,0.05,0.1,0.2};

        addPlots(prefix,"elD0","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","elD0","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedEl","elD0","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedEl_passTight","elD0","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < el_D0s.size(); ++iS){
            proc.lepSelParams.el_maxD0   = el_D0s[iS];
            TString name = TString::Format("%.2f",el_D0s[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getElectrons(*reader_electron);
            if(leps.size()){
                addPlots(prefix,"elD0",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","elD0",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
            }
            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedEl","elD0",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedEl_passTight","elD0",name,iS+1,signalPT);
                }

            }
        }
    }
    void testElSIP3D(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Electron* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     el_S3Ds = {2,3,4,5,6,7,8};

        addPlots(prefix,"elS3","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","elS3","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedEl","elS3","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedEl_passTight","elS3","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < el_S3Ds.size(); ++iS){
            proc.lepSelParams.el_maxD0   = 999.9;
            TString name = TString::Format("%.0f",el_S3Ds[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getElectrons(*reader_electron);
            double maxPT =-1;
            bool foundSig = false;
            for(const auto* l : leps ){
                if(l->sip3D() >= el_S3Ds[iS] ) continue;
                if(l->pt() > maxPT) maxPT = l->pt();
                if(isSignal && signalLep ) if(l->index() == signalLep->index()) foundSig = true;
            }
            if(maxPT > 0){
                addPlots(prefix,"elS3",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","elS3",name,iS+1,isSignal ?signalPT  : maxPT);
            }
            if(isSignal && signalLep ){
                if(foundSig){
                    addPlots(prefix+"_matchedEl","elS3",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedEl_passTight","elS3",name,iS+1,signalPT);
                }

            }
        }
    }
    void testElDZ(const TString& prefix, bool passTight, bool isSignal, float signalPT, const Electron* signalLep){
        LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);

        static const std::vector<float>                     el_DZs = {0.1,0.2,0.3,0.4,0.5};

        addPlots(prefix,"elDZ","incl",0,isSignal ?signalPT :-1);
        if(passTight) addPlots(prefix+"_passTight","elDZ","incl",0,isSignal ?signalPT :-1);

        if(isSignal){
            addPlots(prefix+"_matchedEl","elDZ","incl",0,signalPT);
            if(passTight)  addPlots(prefix+"_matchedEl_passTight","elDZ","incl",0,signalPT);
        }


        for(unsigned int iS = 0; iS < el_DZs.size(); ++iS){
            proc.lepSelParams.el_maxDZ   = el_DZs[iS];
            TString name = TString::Format("%.1f",el_DZs[iS]);name.ReplaceAll(".","p");
            auto leps = proc.getElectrons(*reader_electron);
            if(leps.size()){
                addPlots(prefix,"elDZ",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
                if(passTight)  addPlots(prefix+"_passTight","elDZ",name,iS+1,isSignal ?signalPT  : leps.front()->pt());
            }
            if(isSignal && signalLep ){
                bool found = false; for(const auto* l : leps) if(l->index() == signalLep->index()) found = true;
                if(found){
                    addPlots(prefix+"_matchedEl","elDZ",name,iS+1,signalPT);
                    if(passTight)  addPlots(prefix+"_matchedEl_passTight","elDZ",name,iS+1,signalPT);
                }
            }
        }
    }

    const Lepton * getMatchedLepton(const GenParticle& genLepton,const std::vector<const Muon *> muons, const std::vector<const Electron*> electrons){
       if(genLepton.absPdgId() == ParticleInfo::p_muminus){
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,muons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return muons[idx];
       } else {
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,electrons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return electrons[idx];
       }
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_wlep < 400) return false;
        const bool passTight = ht_wlep >= 1500;
        TString prefix = reader_event->process == FillerConstants::SIGNAL ? TString::Format("m%i",signal_mass) : FillerConstants::MCProcessNames[reader_event->process];


        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type >= DiHiggsEvent::MU){
            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,20,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.4);
            const auto* recoL = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);
            if(diHiggsEvt.w1_d1->absPdgId() == ParticleInfo::p_muminus){
                testMuID (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Muon*)recoL);
                testMuISO(prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Muon*)recoL);
                testMuDZ (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Muon*)recoL);
                testMuD0 (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Muon*)recoL);
                testMuSIP3D (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Muon*)recoL);
            } else  {
                testElID (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Electron*)recoL);
                testElISO(prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Electron*)recoL);
                testElDZ (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Electron*)recoL);
                testElD0 (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Electron*)recoL);
                testElSIP3D (prefix,passTight,true,diHiggsEvt.w1_d1->pt(),(const Electron*)recoL);
            }
        }
        if(reader_event->process != FillerConstants::SIGNAL){
            testMuID (prefix,passTight,false,0,0);
            testMuISO(prefix,passTight,false,0,0);
            testMuDZ (prefix,passTight,false,0,0);
            testMuD0 (prefix,passTight,false,0,0);
            testMuSIP3D(prefix,passTight,false,0,0);
            testElID (prefix,passTight,false,0,0);
            testElISO(prefix,passTight,false,0,0);
            testElDZ (prefix,passTight,false,0,0);
            testElD0 (prefix,passTight,false,0,0);
            testElSIP3D(prefix,passTight,false,0,0);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getLeptonSelection(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getLeptonSelection(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
