
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
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
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"

#include "Processors/Variables/interface/HiggsSolver.h"


using namespace TAna;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){
    }



    bool runEvent() override {

        auto fill=[&](float n){
             plotter.getOrMake1DPre(smpName,"nevents",";nEvents",30,-0.5,29.5)->Fill(n,weight);
           };

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(hbbCand && wjjCand && passEventFilters && selectedLeptons.size() ==1 ){

            bool passKine = true;
            const float minHT = 400;
            const float minElePT = 30;
            const float minMuPT = 26;

            if(ht_chs < minHT) passKine = false;

            float maxElePT = 0;
            float maxMuPT = 0;
            for(const auto * l : selectedLeptons ) {
                if(l->isMuon()) maxMuPT = std::max(maxMuPT, l->pt());
                else maxElePT = std::max(maxElePT, l->pt());
            }

            if(maxElePT < minElePT && maxMuPT < minMuPT ) passKine = false;


            if(passKine) fill(12);


            if(passTriggerPreselection) fill(13);
        }
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;



        if(hbbCand && wjjCand){
            fill(9);
            if(wjjCand->tau2otau1() < 0.75)
                fill(10);
            if(wjjCand->tau2otau1() < 0.55)
                fill(11);

        }

        if(hbbCand && wjjCand && isSignal()&& diHiggsEvt.type >= DiHiggsEvent::MU){

            bool isDouble = false;
            MomentumF neutrino;
            MomentumF neutrinoAlt;
            auto getInvisible =[&](const MomentumF& met, const MomentumF& vis, const double hMass){
                const double a = hMass*hMass - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
                const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
                const double B = -4*a* vis.z();
                const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
                const double delta = B*B -4*A*C;

                double metZ = 0;
                double metZAlt = 0;
                if(delta < 0) {
                    metZ= -B/(2*A);
                } else {
                    isDouble = true;
                    const double pos = (-B + std::sqrt(delta))/(2*A);
                    const double neg = (-B - std::sqrt(delta))/(2*A);
                    if(std::fabs(pos) < std::fabs(neg)){
                        metZ = pos;
                        metZAlt = neg;
                    }
                    else{
                        metZ = neg;
                        metZAlt = pos;
                    }
                }
                ASTypes::CartLorentzVector neutrinoC(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
                ASTypes::CartLorentzVector neutrinoAltC(met.x(),met.y(),metZAlt,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZAlt*metZAlt));
                neutrino.setP4(neutrinoC);
                neutrinoAlt.setP4(neutrinoAltC);
            };
            getInvisible(reader_event->met,(selectedLepton->p4() + wjjCand->p4()),125. );
            plotter.getOrMake1DPre(smpName,"sigReco",";category",6,-0.5,5.5)->Fill(0.,weight);
            if(!isDouble)plotter.getOrMake1DPre(smpName,"sigReco",";category",6,-0.5,5.5)->Fill(1,weight);
            else {
                if(std::fabs(diHiggsEvt.w1_d2->pz() - neutrino.pz()) < std::fabs(diHiggsEvt.w1_d2->pz() - neutrinoAlt.pz())  )
                    plotter.getOrMake1DPre(smpName,"sigReco",";category",6,-0.5,5.5)->Fill(2.,weight);
                else
                    plotter.getOrMake1DPre(smpName,"sigReco",";category",6,-0.5,5.5)->Fill(3.,weight);
            }



        }


        hbbCand = 0;
        for(const auto& j : reader_fatjet->jets){
            if(j.pt() < 200) break;
            if(j.absEta() > 2) continue;
            if(!j.passTightID()) continue;
            if(PhysicsUtilities::deltaR2(j,*selectedLepton) < 1 ) continue;
            hbbCand =&j;
            break;
        }
        BTagging::CSVSJ_CAT bCat = hbbCand ? BTagging::getCSVSJCat(hbbCand->subJets(), 20, 2.4) : BTagging::CSVSJ_FF;
        hbbMass = hbbCand ? hbbFJSFProc->getCorrSDMass(hbbCand) : 0;

        if(mcProc == FillerConstants::SIGNAL && diHiggsEvt.type >= DiHiggsEvent::MU){
            double w1M = (diHiggsEvt.w1_d1->p4() +diHiggsEvt.w1_d2->p4() ).mass();
            double w2M = (diHiggsEvt.w2_d1->p4() +diHiggsEvt.w2_d2->p4() ).mass();
            double wONS, wOFFS;
            if(w1M < w2M ){
                wOFFS= w1M;wONS= w2M;
            } else {
                wOFFS= w2M;wONS= w1M;
            }

            plotter.getOrMake1DPre(smpName,"wONS",";W mass",50,0,200)->Fill(wONS,weight);
            plotter.getOrMake1DPre(smpName,"wOFFS",";W mass",50,0,200)->Fill(wOFFS,weight);
            plotter.getOrMake1DPre(smpName,"wTot",";W total mass",50,0,200)->Fill(wONS+wOFFS,weight);
            plotter.getOrMake1DPre(smpName,"wOFFSFrac",";W total mass",50,0,1)->Fill(wOFFS/(wONS+wOFFS),weight);

            double z1 = diHiggsEvt.w1_d2->pz();
            double z2 = getRndGen()->Uniform() > 0.5 ? diHiggsEvt.w2_d2->pz() : diHiggsEvt.w2_d1->pz();

            double avgZ = 2*(z1-z2)/(z1+z2);
            plotter.getOrMake1DPre(smpName,"pzDif",";pzDif",200,-500,500)->Fill(z1-z2,weight);
            plotter.getOrMake1DPre(smpName,"pzRel",";pzRel",200,-5,5)->Fill(avgZ,weight);
            plotter.getOrMake1DPre(smpName,"z1oz2",";z1/z2",200,0,5)->Fill(z1/z2,weight);


            plotter.getOrMake2DPre(smpName,"z1vz2",";z1;z2",100,-500,500,100,-500,500)->Fill(z1,z2,weight);

            plotter.getOrMake1DPre(smpName,"z1",";z1",500,0,500)->Fill(z1,weight);



            plotter.getOrMake1DPre(smpName,"qq_dr",";#DeltaR qq",120,0,3)->Fill(
                    PhysicsUtilities::deltaR(*diHiggsEvt.w2_d1,*diHiggsEvt.w2_d2),weight);
            if(bCat>=5 && hbbMass>100 && hbbMass <150 ){
                plotter.getOrMake1DPre(smpName,"goodhbb_qq_dr",";#DeltaR qq",120,0,3)->Fill(
                        PhysicsUtilities::deltaR(*diHiggsEvt.w2_d1,*diHiggsEvt.w2_d2),weight);
            }
        }

        jets.clear();
        for(const auto& j : reader_jet->jets ){
            if(j.absEta() >2.4) continue;
            if(j.pt() <20) continue;
            if(!j.passTightID()) continue;
            if(PhysicsUtilities::deltaR2(j,*selectedLepton) <.4*.4) continue;
            if(hbbCand)
                if(PhysicsUtilities::deltaR2(j,*hbbCand) < 1.4*1.4 ) continue;
            jets.push_back(&j);
        }
        int nLooseBTags = PhysicsUtilities::selObjsD(jets,[](const Jet* j){return BTagging::isLooseCSVTagged(*j);} ).size();






        weight = (EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi()))
                *(smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size() ? trigSFProc->getLeptonTriggerSF(ht_chs, (selectedLepton && selectedLepton->isMuon())) : 1.0 )
                *puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL)
                *leptonSFProc->getSF()
                *sjbtagSFProc->getSF({hbbCand})*ak4btagSFProc->getSF(jets);




        fill(0);
        if(!hbbCand) return true;
        fill(1);

        plotter.getOrMake1DPre(smpName,"beforeVeto_njets",";N. jets",10,-0.5,9.5)->Fill(jets.size(),weight);
        if(nLooseBTags == 0)
            plotter.getOrMake1DPre(smpName,"njets",";N. jets",10,-0.5,9.5)->Fill(jets.size(),weight);


        double maxCSV = -1;
        for(const auto* j : jets){if(j->csv()> maxCSV ) maxCSV = j->csv(); }
        plotter.getOrMake1DPre(smpName,"beforeVeto_maxCSV",";CSV",100,0,1)->Fill(maxCSV,weight);

        if(jets.size() < 2) return true;
        fill(2);
        if(nLooseBTags>0) return true;
        fill(3);

        MomentumF wjjCand;
        double dr2 =0;
        if(jets.size() == 2){
            wjjCand.setP4(jets[0]->p4() + jets[1]->p4());
            dr2 = PhysicsUtilities::deltaR2(*jets[0],*jets[1]);
        } else {
            double dr12 = PhysicsUtilities::deltaR2(*jets[0],*jets[1]);
            double dr13 = PhysicsUtilities::deltaR2(*jets[0],*jets[2]);
            double dr23 = PhysicsUtilities::deltaR2(*jets[1],*jets[2]);
            if(dr12 < dr13 && dr12 < dr23){
                wjjCand.setP4(jets[0]->p4() + jets[1]->p4());
                dr2=dr12;
            } else if(dr13 < dr12 && dr13 < dr23){
                wjjCand.setP4(jets[0]->p4() + jets[2]->p4());
                dr2=dr13;
            } else{
                wjjCand.setP4(jets[1]->p4() + jets[2]->p4());
                dr2=dr23;
            }
        }

        plotter.getOrMake1DPre(smpName,"dr",";#DeltaR",100,0,5)->Fill(std::sqrt(dr2),weight);
        plotter.getOrMake1DPre(smpName,"mass",";wqq mass",200,0,200)->Fill(wjjCand.mass(),weight);

        int nSJs = PhysicsUtilities::selObjsMom(hbbCand->subJets(),
                20, 2.4).size();
        if(nSJs < 2) return true;
        if(bCat < 4) return true;


        neutrino=        HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + wjjCand.p4()) );
        wlnu        =  neutrino.p4() + selectedLepton->p4();
        hWW         = wlnu.p4() + wjjCand.p4();
        hh =  (hWW.p4() + hbbCand->p4()) ;

        fill(4);
        if(hbbMass >100 && hbbMass < 150){
            fill(5);

            if(reader_event->met.pt() > 50){
                fill(6);
                plotter.getOrMake1DPre(smpName,"hhmass",";hh mass",30,0,3000)->Fill(hh.mass(),weight);
            }
        }
        if(bCat< 5) return false;
        if(hbbMass >100 && hbbMass < 150){
            fill(7);
            if(reader_event->met.pt() > 50){
                fill(8);
                plotter.getOrMake1DPre(smpName,"hhmass",";hh mass",30,0,3000)->Fill(hh.mass(),weight);
            }
        }

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void testATLASSel(std::string fileName, int treeInt,int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void testATLASSel(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
