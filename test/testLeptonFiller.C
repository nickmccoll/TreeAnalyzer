
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TreeReaders/interface/FillerConstants.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : BaseTreeAnalyzer(fileName,treeName,treeInt, randSeed){

    }
    void loadVariables() override {
        reader_event =  loadReader<EventReader>("event",isRealData());
        reader_electrons = loadReader<ElectronReader>("electron",true);
        reader_muons = loadReader<MuonReader>("muon");
    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& electrons = reader_electrons->electrons;
        const auto& muons = reader_muons->muons;

        auto plotEIso = [&](const Electron& l, TString pre){
            plotter.getOrMake1D(TString::Format("electron%spt",pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);

            if(l.passLooseID()         ) plotter.getOrMake1D(TString::Format("electron%spassLooseID_pt"         ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMedID  ()         ) plotter.getOrMake1D(TString::Format("electron%spassMedID_pt"           ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passTightID()         ) plotter.getOrMake1D(TString::Format("electron%spassTightID_pt"         ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passHEEPID ()         ) plotter.getOrMake1D(TString::Format("electron%spassHEEPID_pt"          ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVAHZZ()          ) plotter.getOrMake1D(TString::Format("electron%spassMVAHZZ_pt"          ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVALooseID()      ) plotter.getOrMake1D(TString::Format("electron%spassMVALooseID_pt"      ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVA80ID()         ) plotter.getOrMake1D(TString::Format("electron%spassMVA80ID_pt"         ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVA90ID()         ) plotter.getOrMake1D(TString::Format("electron%spassMVA90ID_pt"         ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passLooseID_noIso()   ) plotter.getOrMake1D(TString::Format("electron%spassLooseID_noIso_pt"   ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMedID_noIso  ()   ) plotter.getOrMake1D(TString::Format("electron%spassMedID_noIso_pt"     ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passTightID_noIso()   ) plotter.getOrMake1D(TString::Format("electron%spassTightID_noIso_pt"   ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passHEEPID_noIso ()   ) plotter.getOrMake1D(TString::Format("electron%spassHEEPID_noIso_pt"    ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVALooseID_noIso()) plotter.getOrMake1D(TString::Format("electron%spassMVALooseID_noIso_pt",pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVA80ID_noIso()   ) plotter.getOrMake1D(TString::Format("electron%spassMVA80ID_noIso_pt"   ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);
            if(l.passMVA90ID_noIso()   ) plotter.getOrMake1D(TString::Format("electron%spassMVA90ID_noIso_pt"   ,pre.Data()),",p_{T}",20,0,100)->Fill(l.pt(),weight);        };

        for(const auto& l : electrons){
            plotEIso(l,"_");

            const int idx = l.index();
            bool flag = l.passMedID() && !l.passMedID_noIso();
            if(l.passHEEPID_noIso() != l.passHEEPID()) flag = true;
            if(flag){
            std::cout << "!M! "<< l.passMedID  () <<" "<< l.passMedID_noIso() <<" "<< l.pfIso()<<" "<<reader_electrons->HoE[idx]  <<" :: ";
            for(unsigned int iC = 0; iC < 16;++iC){
                std::cout << FillerConstants::doesPass(reader_electrons->passMedCutBased[idx],iC) <<" ";
            }
            std::cout << std::endl;
            std::cout << "!T! "<< l.passTightID  () <<" "<< l.passTightID_noIso() <<" "<< l.pfIso() <<" :: ";
            for(unsigned int iC = 0; iC < 16;++iC){
                std::cout << FillerConstants::doesPass(reader_electrons->passTightCutBased[idx],iC) <<" ";
            }
            std::cout << std::endl;
            std::cout << "!H! "<< l.passHEEPID  () <<" "<< l.passHEEPID_noIso()
                    <<" "<< l.trackerIso()*reader_electrons->uncorPt[idx] <<" :: ";
             std::cout << std::endl;
            }

            plotter.getOrMake1D("electron_d0"       ,";d0"      , 500, 0,5)->Fill(l.d0(),weight);
            plotter.getOrMake1D("electron_dz"       ,";dz"      , 500, 0,5)->Fill(l.dz(),weight);
            plotter.getOrMake1D("electron_sip3D"    ,";sip3D"   , 500, 0,5)->Fill(l.sip3D(),weight);
            plotter.getOrMake1D("electron_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            plotter.getOrMake1D("electron_miniIsoFP"  ,";miniIsoFP" , 500, 0,5)->Fill(l.miniIsoFP(),weight);
            plotter.getOrMake1D("electron_eaRelISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            plotter.getOrMake1D("electron_trackerIso"  ,";miniIso" , 500, 0,5)->Fill(l.trackerIso(),weight);
            plotter.getOrMake1D("electron_mvaID"    ,";mvaID"   , 500, 0,5)->Fill(l.mvaID(),weight);
        }

        for(unsigned int iL1 = 0; iL1 < electrons.size(); ++iL1){
            const auto& l1 = electrons[iL1];
            for(unsigned int iL2 = iL1+1; iL2 < electrons.size(); ++iL2){
                const auto& l2 = electrons[iL2];
                if(l1.q() == l2.q()) continue;
                float mass = (l1.p4() + l2.p4()).mass();
                if(mass > 115 || mass < 65 ) continue;
                plotEIso(l1,"_OS_");
                plotEIso(l2,"_OS_");
        }
        }

        auto plotMIso = [&](const Muon& l, TString pre){
            plotter.getOrMake1D(TString::Format("muon%spt",pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);

            if(l.passSoftID       ()  ) plotter.getOrMake1D(TString::Format("muon%spassSoftID_pt" ,pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passLooseID      ()  ) plotter.getOrMake1D(TString::Format("muon%spassLooseID_pt",pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passMedID        ()  ) plotter.getOrMake1D(TString::Format("muon%spassMedID_pt"  ,pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passTightID      ()  ) plotter.getOrMake1D(TString::Format("muon%spassTightID_pt",pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passHighPT ()        ) plotter.getOrMake1D(TString::Format("muon%spassHighPT_pt" ,pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
        };

        for(const auto& l : muons){
            plotMIso(l,"_");

            plotter.getOrMake1D("muon_d0"       ,";d0"      , 500, 0,5)->Fill(l.d0(),weight);
            plotter.getOrMake1D("muon_dz"       ,";dz"      , 500, 0,5)->Fill(l.dz(),weight);
            plotter.getOrMake1D("muon_sip3D"    ,";sip3D"   , 500, 0,5)->Fill(l.sip3D(),weight);
            plotter.getOrMake1D("muon_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            plotter.getOrMake1D("muon_pfISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            plotter.getOrMake1D("muon_trackerISO" ,";trackerIso", 500, 0,5)->Fill(l.trackerIso(),weight);

            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_CutBasedIdMediumPrompt) ){
                plotter.getOrMake1D("muon_medPrompt_d0"       ,";d0"      , 500, 0,5)->Fill(l.d0(),weight);
                plotter.getOrMake1D("muon_medPrompt_dz"       ,";dz"      , 500, 0,5)->Fill(l.dz(),weight);
            }
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_PFIsoVeryLoose) )
                plotter.getOrMake1D("muon_pfVL_pfISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_PFIsoLoose) )
                plotter.getOrMake1D("muon_pfL_pfISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_PFIsoMedium) )
                plotter.getOrMake1D("muon_pfM_pfISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_PFIsoTight) )
                plotter.getOrMake1D("muon_pfT_pfISO" ,";pfIso", 500, 0,5)->Fill(l.pfIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_TkIsoLoose) )
                plotter.getOrMake1D("muon_tkL_trackerISO" ,";pfIso", 500, 0,5)->Fill(l.trackerIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_MiniIsoLoose) )
                plotter.getOrMake1D("muon_mL_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_MiniIsoMedium) )
                plotter.getOrMake1D("muon_mM_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            if(FillerConstants::doesPass(l.id(),FillerConstants::MUID_MiniIsoTight) )
                plotter.getOrMake1D("muon_mT_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);


        }

        for(unsigned int iL1 = 0; iL1 < muons.size(); ++iL1){
            const auto& l1 = muons[iL1];
            for(unsigned int iL2 = iL1+1; iL2 < muons.size(); ++iL2){
                const auto& l2 = muons[iL2];
                if(l1.q() == l2.q()) continue;
                float mass = (l1.p4() + l2.p4()).mass();
                if(mass > 115 || mass < 65 ) continue;
                plotMIso(l1,"_OS_");
                plotMIso(l2,"_OS_");
        }
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    bool realData = false;
    std::shared_ptr<EventReader> reader_event = 0;
    std::shared_ptr<ElectronReader> reader_electrons = 0;
    std::shared_ptr<MuonReader> reader_muons = 0;
    HistGetter plotter;

};

#endif

void testLeptonFiller(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
