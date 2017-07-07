
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, bool realData) : BaseTreeAnalyzer(fileName,treeName), realData(realData){

    }
    void loadVariables() override {
        reader_event = (EventReader*)load(new EventReader("event",realData));
        reader_electrons = (ElectronReader*)load(new ElectronReader("electron"));
        reader_muons = (MuonReader*)load(new MuonReader("muon"));
    }

    bool runEvent() override {
        const double weight = reader_event->weight;
        const auto& electrons = reader_electrons->electrons;
        const auto& muons = reader_muons->muons;

        auto plotEIso = [&](const Electron& l, TString pre){
            plotter.getOrMake1D(TString::Format("electron%spt",pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);

            if(l.passVetoID       ()  ) plotter.getOrMake1D(TString::Format("electron%spassVetoID_pt"       ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passLooseID      ()  ) plotter.getOrMake1D(TString::Format("electron%spassLooseID_pt"      ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passMedID        ()  ) plotter.getOrMake1D(TString::Format("electron%spassMedID_pt"        ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passTightID      ()  ) plotter.getOrMake1D(TString::Format("electron%spassTightID_pt"      ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passHEEPID       ()  ) plotter.getOrMake1D(TString::Format("electron%spassHEEPID_pt"       ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passVetoID_noISO ()  ) plotter.getOrMake1D(TString::Format("electron%spassVetoID_noISO_pt" ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passLooseID_noISO()  ) plotter.getOrMake1D(TString::Format("electron%spassLooseID_noISO_pt",pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passMedID_noISO  ()  ) plotter.getOrMake1D(TString::Format("electron%spassMedID_noISO_pt"  ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passTightID_noISO()  ) plotter.getOrMake1D(TString::Format("electron%spassTightID_noISO_pt",pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passHEEPID_noISO ()  ) plotter.getOrMake1D(TString::Format("electron%spassHEEPID_noISO_pt" ,pre.Data()),";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
        };

        for(const auto& l : electrons){
            plotEIso(l,"_");

            plotter.getOrMake1D("electron_d0"       ,";d0"      , 500, 0,5)->Fill(l.d0(),weight);
            plotter.getOrMake1D("electron_dz"       ,";dz"      , 500, 0,5)->Fill(l.dz(),weight);
            plotter.getOrMake1D("electron_sip3D"    ,";sip3D"   , 500, 0,5)->Fill(l.sip3D(),weight);
            plotter.getOrMake1D("electron_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            plotter.getOrMake1D("electron_eaRelISO" ,";eaRelISO", 500, 0,5)->Fill(l.eaRelISO(),weight);
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
            if(l.passMed16ID      ()  ) plotter.getOrMake1D(TString::Format("muon%spassMed16ID_pt",pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
            if(l.passHighPT ()        ) plotter.getOrMake1D(TString::Format("muon%spassHighPT_pt" ,pre.Data())      ,";p_{T}", 20, 0,100)->Fill(l.pt(),weight);
        };

        for(const auto& l : muons){
            plotMIso(l,"_");

            plotter.getOrMake1D("muon_d0"       ,";d0"      , 500, 0,5)->Fill(l.d0(),weight);
            plotter.getOrMake1D("muon_dz"       ,";dz"      , 500, 0,5)->Fill(l.dz(),weight);
            plotter.getOrMake1D("muon_sip3D"    ,";sip3D"   , 500, 0,5)->Fill(l.sip3D(),weight);
            plotter.getOrMake1D("muon_miniIso"  ,";miniIso" , 500, 0,5)->Fill(l.miniIso(),weight);
            plotter.getOrMake1D("muon_dbRelISO" ,";dbRelISO", 500, 0,5)->Fill(l.dbRelISO(),weight);
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
    EventReader * reader_event = 0;
    ElectronReader * reader_electrons = 0;
    MuonReader * reader_muons = 0;
    HistGetter plotter;

};

#endif

void testLeptonFiller(std::string fileName ="output.root",std::string outFileName = "plots.root"){
    TString path = fileName;
    Analyzer a(fileName,"treeMaker/Events",path.Contains("data",TString::kIgnoreCase));
    a.analyze();
    a.write(outFileName);
}
