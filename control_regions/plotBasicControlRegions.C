
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
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"


using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        turnOffCorr(CORR_LEP );
    }


    void plotEventVariables(const TString& prefix, const float weight){
        plotter.getOrMake1DPre(prefix,"numvtx",";# of verticies",75,0,75)->Fill(reader_event->npv,weight);
        const float puW = isRealData() ? 1.0 : puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL);
        plotter.getOrMake1DPre(prefix,"woPUW_numvtx",";# of verticies",75,0,75)->Fill(reader_event->npv,weight/puW);
        plotter.getOrMake1DPre(prefix,"ht",";#it{H}_{T} [GeV]",400,0,4000)->Fill(ht_wlep,weight);
    }

    void plotAK4Kinematics(const TString& prefix, const std::vector<Jet>& jets, const float weight) {
        std::vector<const Jet*>  filteredJets = PhysicsUtilities::selObjsMom(jets,20.0,2.4);

        plotter.getOrMake1DPre(prefix,"nJets",";# of jets",20,-0.5,19.5)->Fill(filteredJets.size(),weight);

        size nBL = 0;
        size nBM = 0;
        size nBT = 0;
        for(const auto* j : filteredJets ){
            if(BTagging::isLooseCSVTagged(*j)) nBL++;
            if(BTagging::isMediumCSVTagged(*j)) nBM++;
            if(BTagging::isTightCSVTagged(*j)) nBT++;
        }

        plotter.getOrMake1DPre(prefix,"nLooseBtags",";# of loose b-tagged jets",10,-0.5,9.5)->Fill(nBL,weight);
        plotter.getOrMake1DPre(prefix,"nMedBtags",";# of medium b-tagged jets",10,-0.5,9.5)->Fill(nBM,weight);
        plotter.getOrMake1DPre(prefix,"nTightBtags",";# of tight b-tagged jets",10,-0.5,9.5)->Fill(nBT,weight);

        if(filteredJets.size() > 0)plotter.getOrMake1DPre(prefix,"leadingJetPT",";leading jet p_{T} [GeV]",400,0,2000)->Fill(filteredJets[0]->pt(),weight);
        if(filteredJets.size() > 1)plotter.getOrMake1DPre(prefix,"secLeadingJetPT",";second jet p_{T} [GeV]",400,0,2000)->Fill(filteredJets[1]->pt(),weight);
        if(filteredJets.size() > 0)plotter.getOrMake1DPre(prefix,"leadingJetETA",";leading jet #eta",50,-2.5,2.5)->Fill(filteredJets[0]->eta(),weight);
        if(filteredJets.size() > 1)plotter.getOrMake1DPre(prefix,"secLeadingJetETA",";second jet #eta",50,-2.5,2.5)->Fill(filteredJets[1]->eta(),weight);
        if(filteredJets.size() > 0)plotter.getOrMake1DPre(prefix,"leadingJetCSV",";leading jet csv",110,-.1,1)->Fill(filteredJets[0]->csv(),weight);
        if(filteredJets.size() > 1)plotter.getOrMake1DPre(prefix,"secLeadingJetCSV",";second jet csv",110,-.1,1)->Fill(filteredJets[1]->csv(),weight);

        for(const auto* j : filteredJets ){
            plotter.getOrMake1DPre(prefix,"jetPT",";jet p_{T} [GeV]",400,0,2000)->Fill(j->pt(),weight);
            plotter.getOrMake1DPre(prefix,"jetCSV",";jet csv",110,-.1,1)->Fill(j->csv(),weight);
            plotter.getOrMake1DPre(prefix,"jetETA",";jet #eta",50,-2.5,2.5)->Fill(j->eta(),weight);
        }

    }

    void plotAK8Kinematics(const TString& prefix, const float weight) {

        plotter.getOrMake1DPre(prefix,"nJets",";# of jets",20,-0.5,19.5)->Fill(fatjetCands.size(),weight);


        auto mkPlot= [&](const FatJet* j, TString namePre, TString titlePre){

            plotter.getOrMake1DPre(prefix,TString::Format("%sPT",namePre.Data()),TString::Format(";%s p_{T} [GeV]",titlePre.Data()),400,0,2000)->Fill(j->pt(),weight);
            plotter.getOrMake1DPre(prefix,TString::Format("%sETA",namePre.Data()),TString::Format(";%s #eta",titlePre.Data()),50,-2.5,2.5)->Fill(j->eta(),weight);

            plotter.getOrMake1DPre(prefix,TString::Format("%sMass",namePre.Data()),TString::Format(";%s mass [GeV]",titlePre.Data()),500,0,500)->Fill(j->mass(),weight);
            plotter.getOrMake1DPre(prefix,TString::Format("%sSDMass",namePre.Data()),TString::Format(";%s soft-drop mass [GeV]",titlePre.Data()),500,0,500)->Fill(j->sdMom().mass(),weight);
            plotter.getOrMake1DPre(prefix,TString::Format("%sRSDMass",namePre.Data()),TString::Format(";%s raw soft-drop mass [GeV]",titlePre.Data()),500,0,500)->Fill(j->rawSdMom().mass(),weight);

            plotter.getOrMake1DPre(prefix,TString::Format("%sNSubjets",namePre.Data()),TString::Format(";%s # of subjets",titlePre.Data()),4,-0.5,3.5)->Fill(j->nSubJets(),weight);

            if(j->nSubJets() > 1){
                plotter.getOrMake1DPre(prefix,TString::Format("%sMinSJCSV",namePre.Data()),TString::Format(";%s min subjet csv",titlePre.Data()),110,-.1,1)->Fill(j->minSJCSV(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sMaxSJCSV",namePre.Data()),TString::Format(";%s max subjet csv",titlePre.Data()),110,-.1,1)->Fill(j->maxSJCSV(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sTau2oTau1",namePre.Data()),TString::Format(";%s #tau_{2}/#tau_{1}",titlePre.Data()),100,0,1)->Fill(j->tau2otau1(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sTau3oTau1",namePre.Data()),TString::Format(";%s #tau_{3}/#tau_{1}",titlePre.Data()),100,0,1)->Fill(j->tau3otau1(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sTau3oTau2",namePre.Data()),TString::Format(";%s #tau_{3}/#tau_{2}",titlePre.Data()),100,0,1)->Fill(j->tau3otau2(),weight);

                plotter.getOrMake1DPre(prefix,TString::Format("%sTau1",namePre.Data()),TString::Format(";%s #tau_{1}",titlePre.Data()),100,0,1)->Fill(j->tau1(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sTau2",namePre.Data()),TString::Format(";%s #tau_{2}",titlePre.Data()),100,0,1)->Fill(j->tau2(),weight);
                plotter.getOrMake1DPre(prefix,TString::Format("%sTau3",namePre.Data()),TString::Format(";%s #tau_{3}",titlePre.Data()),100,0,1)->Fill(j->tau3(),weight);
            }


        };

        if(fatjetCands.size() > 0)  mkPlot(fatjetCands[0],"leadingJet","leading jet");
        if(fatjetCands.size() > 1)  mkPlot(fatjetCands[1],"secLeadingJet","second jet");
        for(const auto* j : fatjetCands ){ mkPlot(j,"jet","jet");}
    }

    void plotLeptons(const TString& prefix, const float weight){
        for(const auto * l : selectedLeptons){
            plotter.getOrMake1DPre(prefix,"leptonPT",";lepton p_{T} [GeV]",200,0,1000)->Fill(l->pt(),weight);
            plotter.getOrMake1DPre(prefix,"leptonETA",";lepton #eta",50,-2.5,2.5)->Fill(l->eta(),weight);
        }

        plotter.getOrMake1DPre(prefix,"met",";#slash{E}_{T} [GeV]",200,0,1000)->Fill(reader_event->met.pt(),weight);
        auto w = reader_event->met.p4() + selectedLepton->p4();
        plotter.getOrMake1DPre(prefix,"wPT",";lepton+#slash{E}_{T}  p_{T} [GeV]",200,0,1000)->Fill(w.pt(),weight);
        plotter.getOrMake1DPre(prefix,"mt",";m_{T} [GeV]",500,0,500)->Fill(JetKinematics::transverseMass(*selectedLepton,reader_event->met),weight);
    }

    void plotHH(const TString prefix, const float weight){
        plotter.getOrMake1DPre(prefix,"hhMass",";HH mass [GeV]",400,0,2000)->Fill(hh.mass(),weight);
        plotter.getOrMake1DPre(prefix,"hbbMass",";H(bb) mass [GeV]",500,0,500)->Fill(hbbCand->mass(),weight);
        plotter.getOrMake1DPre(prefix,"hbbPT"     ,";H(bb) p_{T} [GeV]"             ,400,0,2000)->Fill(hbbCand->pt(),weight);
        plotter.getOrMake1DPre(prefix,"hbbSDMass" ,";H(bb) soft-drop mass [GeV]"    ,500,0,500)->Fill(hbbCand->sdMom().mass(),weight);
        plotter.getOrMake1DPre(prefix,"hbbRSDMass",";H(bb) raw soft-drop mass [GeV]",500,0,500)->Fill(hbbCand->rawSdMom().mass(),weight);

        if(hbbCand->sdMom().mass() >10){
            plotter.getOrMake1DPre(prefix+"_hbbM10","hhMass",";HH mass [GeV]",400,0,2000)->Fill(hh.mass(),weight);
            plotter.getOrMake1DPre(prefix+"_hbbM10","hbbMass",";H(bb) mass [GeV]",500,0,500)->Fill(hbbCand->mass(),weight);
            plotter.getOrMake1DPre(prefix+"_hbbM10","hbbPT"     ,";H(bb) p_{T} [GeV]"             ,400,0,2000)->Fill(hbbCand->pt(),weight);
            plotter.getOrMake1DPre(prefix+"_hbbM10","hbbSDMass" ,";H(bb) soft-drop mass [GeV]"    ,500,0,500)->Fill(hbbCand->sdMom().mass(),weight);
            plotter.getOrMake1DPre(prefix+"_hbbM10","hbbRSDMass",";H(bb) raw soft-drop mass [GeV]",500,0,500)->Fill(hbbCand->rawSdMom().mass(),weight);
        }
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;

        if(reader_event->process >= FillerConstants::SINGLET && reader_event->process <= FillerConstants::TTX )
            smpName = "other";


        std::vector<const Jet*>  filteredJets = PhysicsUtilities::selObjsMom(reader_jetwlep->jets,20.0,2.4);
        size nBL = 0;
        size nBM = 0;
        size nBT = 0;
        for(const auto* j : filteredJets ){
            if(BTagging::isLooseCSVTagged(*j)) nBL++;
            if(BTagging::isMediumCSVTagged(*j)) nBM++;
            if(BTagging::isTightCSVTagged(*j)) nBT++;
        }

        auto mkSTDPlots = [&](const TString& prefix, const float weight){
            plotEventVariables(prefix,weight);
            plotAK4Kinematics(prefix +"_ak4Wlep",reader_jetwlep->jets,weight);
            plotAK4Kinematics(prefix +"_ak4Nolep",reader_jet->jets,weight);
            plotAK8Kinematics(prefix +"_ak8",weight);
            plotLeptons(prefix,weight);
        };
        auto mkHHPlots = [&](const TString& prefix, const float weight){
            if(wjjCand && hbbCand){

                const bool bTaggedW = FatJetSelHelpers::passWjjSelection(wjjCand,fjProc->wjj_maxT2oT1,BTagging::CSV_INCL,fjProc->wjj_minMass,fjProc->wjj_maxMass) && wjjCand->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_T];
                const bool stdW = FatJetSelHelpers::passWjjSelection(wjjCand,fjProc->wjj_maxT2oT1,fjProc->wjj_maxCSVWP,fjProc->wjj_minMass,fjProc->wjj_maxMass);
//                const bool stdHbb =  FatJetSelHelpers::passHbbSelection(hbbCand,fjProc->hbb_maxT2oT1,fjProc->hbb_l_firMinCSVWP,fjProc->hbb_l_secMinCSVWP,fjProc->hbb_minMass,fjProc->hbb_maxMass);
                const bool antiHbb = FatJetSelHelpers::passHbbSelection(hbbCand,fjProc->hbb_maxT2oT1,BTagging::CSV_INCL,BTagging::CSV_INCL,fjProc->hbb_minMass,fjProc->hbb_maxMass)&&
                        hbbCand->maxSJCSV() < BTagging::CSVWP_VALS[BTagging::CSV_L];

                const bool tauHbb = FatJetSelHelpers::passHbbSelection(hbbCand,fjProc->hbb_maxT2oT1,BTagging::CSV_INCL,BTagging::CSV_INCL,fjProc->hbb_minMass,fjProc->hbb_maxMass);


                bool cr_btaggedWjj = bTaggedW;
                bool cr_btaggedWjj2 = bTaggedW && tauHbb;
                bool cr_antiBHbb = antiHbb && stdW;




                if(cr_btaggedWjj)
                    plotHH(prefix + "_wjjBtagCR",weight);
                if(cr_btaggedWjj2)
                    plotHH(prefix + "_wjjBtagTauCR",weight);
                if(cr_antiBHbb)
                    plotHH(prefix + "_HbbAntiBCR",weight);
            }
        };


        auto doASet = [&](TString smpName, const float weight){
            if(selectedLeptons.size() == 1){
                if(selectedLepton->isMuon()){
                    TString prefix = smpName + "_1mu";
                    mkSTDPlots(prefix,weight);
                    mkHHPlots(prefix,weight);
                    if(nBM == 0){ mkSTDPlots(prefix + "_0b",weight);}
                    else if(nBM ==1) { mkSTDPlots(prefix + "_1b",weight);}
                    else if(nBM >=2) { mkSTDPlots(prefix + "_2b",weight);}
                } else {
                    TString prefix = smpName + "_1el";
                    mkSTDPlots(prefix,weight);
                    mkHHPlots(prefix,weight);
                    if(nBM == 0){ mkSTDPlots(prefix + "_0b",weight);}
                    else if(nBM ==1) { mkSTDPlots(prefix + "_1b",weight);}
                    else if(nBM >=2) { mkSTDPlots(prefix + "_2b",weight);}
                }
            } else if(selectedLeptons.size() == 2){
                const float llMass = (selectedLeptons[0]->p4() + selectedLeptons[1]->p4()).mass();
                const bool goodMass = (llMass >= 50 && llMass < 80) || llMass >= 100;
                const bool goodCharge = (selectedLeptons[0]->q() != selectedLeptons[1]->q());
                if(goodCharge && goodMass && selectedLeptons[0]->isMuon()  == selectedLeptons[1]->isMuon()  ){
                    TString prefix = smpName + "_2lsf";
                    mkSTDPlots(prefix,weight);
                    if(nBM == 0){ mkSTDPlots(prefix + "_0b",weight);}
                    else if(nBM ==1) { mkSTDPlots(prefix + "_1b",weight);}
                    else if(nBM >=2) { mkSTDPlots(prefix + "_2b",weight);}
                } else if(goodCharge && selectedLeptons[0]->isMuon()  != selectedLeptons[1]->isMuon()  ) {
                    TString prefix = smpName + "_2lof";
                    mkSTDPlots(prefix,weight);
                    if(nBM == 0){ mkSTDPlots(prefix + "_0b",weight);}
                    else if(nBM ==1) { mkSTDPlots(prefix + "_1b",weight);}
                    else if(nBM >=2) { mkSTDPlots(prefix + "_2b",weight);}
                }


            }
        };

        doASet(smpName,weight);
        leptonSFProc->load(smDecayEvt,selectedLeptons,&jets_wlep);
        doASet(smpName+"_lepSF", isRealData() ? weight : weight*leptonSFProc->getSF());
        doASet(smpName+"_lepBTagSF", isRealData() ? weight : weight*leptonSFProc->getSF()*ak4btagSFProc->getSF(jets));
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void plotBasicControlRegions(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void plotBasicControlRegions(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
