#include <AnalysisSupport/Utilities/interface/Types.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TMath.h>
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/GenParticle.h"
#include "DataFormats/interface/Lepton.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Configuration/interface/ReaderConstants.h"
#include <TString.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include "TreeAnalyzer/framework/Processors/Variables/interface/HiggsSolver.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)


#include "HHSolTreeAnalyzer.h"



using namespace TAna;
using namespace FillerConstants;


class Analyzer : public HHSolTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed,int step) :
        HHSolTreeAnalyzer(fileName,treeName,treeInt, randSeed), step(step), HSolver(""),
        BkgHSolver(""){
        parameters = ReaderConstants::set2017Parameters();
        parameters.hww.liFileName = "hhSol_templates.root";
        parameters.hww.bkgLiFileName = "hhSol_bkgTemplates.root";
        HSolver.setParamters(parameters.hww);
        BkgHSolver.setParamters(parameters.hww);
    }


    double getCorHPt(double inPT){
        const double b = 0.9663;
        const double m = -0.00001013;
        double boundPT = inPT;
        if(inPT < 300) boundPT = 300;
        if(inPT > 2000) boundPT = 2000;
        double corPT = (-1*b + std::sqrt(b*b+4.*m*boundPT))/(2.*m);

        return inPT * corPT/boundPT; //apply as a correction due to the limted range;

    }


    void getHHPtResolution(){
        std::vector<float> ptBins {300,400,500,600,800,1000,1200,1400,1600,1800,2000};

        TString prefix = (isSignal() ? TString("signal") : smpName);
        bool found = false;
        for(unsigned int iB = 0; iB < ptBins.size()-1; ++iB){
            if(true_hwwMag >= ptBins[iB] && true_hwwMag < ptBins[iB+1]){
                prefix += TString::Format("_pt_%.0f_%.0f",ptBins[iB],ptBins[iB+1]);
                found = true;
                break;
            }
        }
        if(!found) return;

        float res = hwwMag/true_hwwMag;

        plotter.getOrMake1DPre("signal","true_Hpt",";true H pT",340,300,2000)
                ->Fill(true_hwwMag,weight);

        plotter.getOrMake1DPre(prefix,"HPt",";H pT",340,300,2000)->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(prefix,"true_Hpt",";true H pT",340,300,2000)
                ->Fill(true_hwwMag,weight);

        plotter.getOrMake1DPre(prefix,"HPtRes",";H pT: reco/true",100,0,2)
                ->Fill(res,weight);


        plotter.getOrMake1DPre(prefix,"corrHPtRes",";H pT: reco/true",100,0,2)
                ->Fill(getCorHPt(hwwMag)/true_hwwMag,weight);

        plotter.getOrMake1DPre(prefix,"mX",";mX",200,0,4000)
                ->Fill(*sampParam ,weight);

        //TString resBin = "_res_";
        //if(res < 0.6) resBin += "lt_0p6";
        //else if(res < 0.8) resBin += "0p6_to_0p8";
        //else if(res < 1.2) resBin += "0p8_to_1p2";
        //else resBin += "gt_1p2";



        //plotter.getOrMake1DPre(prefix+resBin,"HPt",";H pT",340,300,2000)->Fill(hwwMag,weight);
        //plotter.getOrMake1DPre(prefix+resBin,"true_Hpt",";true H pT",340,300,2000)
        //        ->Fill(true_hwwMag,weight);
        //plotter.getOrMake1DPre(prefix+resBin,"bbJetPT_o_HPt",";bb jet pT / true H pT",100,0,4)
        //        ->Fill(bbJet.pt()/true_hwwMag ,weight);
        //
        //plotter.getOrMake1DPre(prefix+resBin,"bbJetSDMass",";bbJetSDMass",100,0,200)
        //        ->Fill(*bbJet_SDmass ,weight);
        //
        //plotter.getOrMake1DPre(prefix+resBin,"hhMass",";HH mass",200,0,4000)
        //        ->Fill(*hh_orig ,weight);
        //
        //plotter.getOrMake1DPre(prefix+resBin,"hhMassRes",";HH mass / true HH Mass",100,0,4)
        //        ->Fill(*hh_orig/float(*sampParam) ,weight);
        //
        //
        //
        //plotter.getOrMake1DPre(prefix+resBin,"hhDR",";HH DR",100,0,4)
        //        ->Fill(PhysicsUtilities::deltaR(bbJet,true_H) ,weight);
        //
        //plotter.getOrMake1DPre(prefix+resBin,"hhDP",";HH DP",100,-4,4)
        //        ->Fill(PhysicsUtilities::deltaPhi((lepton.p4()+qqJet.p4()+met.p4())
        //,true_W) ,weight);



    }



    void getTemplatePlots(bool finalBinning) {

        double extraMetPerp = metPerp-neutPerp;
        double extraMetPar = metPar-neutPar;

        double wqqPTDif = qqJet.pt()-true_jet.pt();
        double wqqPTRes = true_jet.pt() < 10 ? 2.0 : wqqPTDif/true_jet.pt();

        double true_sf =  true_jet.pt()/qqJet.pt();

        ASTypes::CylLorentzVectorF approxJ (true_sf*qqJet.pt(),qqJet.eta(),qqJet.phi(),
                isVirtualWqq ? 31 : 80.0);


        auto approxW = lepton.p4() + true_neut.p4();

        auto approxH = approxW + approxJ;

        float corHPT = getCorHPt(hwwMag);


        TString pf = (isSignal() ? TString("signal") : smpName);


        if(approxH.mass()>135){
            std:: cout << *sampParam<<" "<<int(*dhType)<<" "<< approxH.pt() <<" "<< approxH.mass()<<" "<< approxW.pt() <<" "<< approxW.mass()
                    <<" "<< wqqPTRes <<" "<< true_jet.pt() <<" "<<qqJet.pt() <<" "<< lepton.pt()<< " "<< true_lep.pt()
                    <<" "<< PhysicsUtilities::deltaR(qqJet,true_jet)<<" "<< *chi2 <<" "<<isVirtualWqq <<std::endl;
        }
//


        auto plt =[&](const TString& prefix){

            plotter.getOrMake1DPre(prefix,"true_hwwMag",";hwwMag",300,0,3000)
                ->Fill(true_hwwMag,weight);

            plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",1000,-500,500)
                ->Fill(extraMetPerp,weight);
            plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar / hwwMag",400,-2,2)
                ->Fill(extraMetPar/hwwMag,weight);
            plotter.getOrMake1DPre(prefix,"wqqPTRes",";wqqPTRes",400,-2,2)
                ->Fill(wqqPTRes,weight);

            plotter.getOrMake1DPre(prefix,"qqSDMassCoarse",";qqSDMass",2,50,70)
                ->Fill(*qqJet_SDmass,weight);
            plotter.getOrMake1DPre(prefix
                    , "approxW",";approxW",300,0,300)->Fill(approxW.mass(),weight);
            plotter.getOrMake1DPre(prefix, "approxH",";approxH",300,0,300)
                ->Fill(approxH.mass(),weight);


            //correlations

            plotter.getOrMake2DPre(prefix,"hWW_v_Wlnu"      ,";hWW term; Wlnu term"
                    ,100,80,190,100,0,200)       ->Fill(approxH.mass(),approxW.mass(),weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_jetSF"     ,";hWW term; JetSF term"
                    ,100,80,190,50,-0.5,0.5)       ->Fill(approxH.mass(),wqqPTRes,weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_metPerp"   ,";hWW term; MET perp term"
                    ,100,80,190,100,-200,200)    ->Fill(approxH.mass(),extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_metPar"    ,";hWW term; MET par term"
                    ,100,80,190,100,-2,2)     ->Fill(approxH.mass(),extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_jetSF"    ,";Wlnu term; JetSF term"
                    ,100,0,200,100,0,2)      ->Fill(approxW.mass(),wqqPTRes,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_metPerp"  ,";Wlnu term; MET perp term"
                    ,100,0,200,100,-200,200)   ->Fill(approxW.mass(),extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_metPar"   ,";Wlnu term; MET par term"
                    ,100,0,200,100,-2,2)    ->Fill(approxW.mass(),extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"jetSF_v_metPerp" ,";Jet SF term; MET perp term"
                    ,100,-2,2,100,-200,200) ->Fill(wqqPTRes,extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"jetSF_v_metPar"  ,";Jet SF term; MET par term"
                    ,100,-2,2,100,-2,2)  ->Fill(wqqPTRes,extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"metPerp_v_metPar",";MET perp term; MET par term"
                    ,100,-200,200,100,-2,2)->Fill(extraMetPerp,extraMetPar/hwwMag,weight);

            // typedef std::pair<double,double> ROTM;
            // auto rot =[](double hM, double wM, double A) -> ROTM{
            //     double x = hM -A -wM;
            //     double y = hM -A +wM;
            //     return std::make_pair(x,y);
            // };
            // auto mkROT =[&](const TString& title, double A) {
            //     auto rotV = rot(approxH.mass(),approxW.mass(),A);
            //     auto rot1 = rot(80,0,A);
            //     auto rot2 = rot(80,200,A);
            //     auto rot3 = rot(190,0,A);
            //     auto rot4 = rot(190,200,A);
            //     double minX = std::min(std::min(std::min(rot1.first,rot2.first),rot3.first),rot4.first);
            //     double maxX = std::max(std::max(std::max(rot1.first,rot2.first),rot3.first),rot4.first);
            //     double minY = std::min(std::min(std::min(rot1.second,rot2.second),rot3.second),rot4.second);
            //     double maxY = std::max(std::max(std::max(rot1.second,rot2.second),rot3.second),rot4.second);
            //     plotter.getOrMake2DPre(prefix,title,";hWW term; Wlnu term"
            //             ,100,minX,maxX,100,minY,maxY)
            //                             ->Fill(rotV.first,rotV.second,weight);
            //
            // };
            // mkROT("hWW_v_Wlnu_rot80",80);
            // mkROT("hWW_v_Wlnu_rot31",31);
        };


        auto pltFinal =[&](const TString& prefix){



            plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",30,-150,150)
                ->Fill(extraMetPerp,weight);
            plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar/hwwMag",50,-0.6,0.4)
                ->Fill(extraMetPar/hwwMag,weight);
            plotter.getOrMake1DPre(prefix,"wqqPTRes",";wqqPTRes",50,-0.5,0.5)
                ->Fill(wqqPTRes,weight);

            plotter.getOrMake1DPre(prefix,"qqSDMassCoarse",";qqSDMass",2,50,70)
                ->Fill(*qqJet_SDmass,weight);

            if(isVirtualWqq){
                plotter.getOrMake2DPre(prefix,"hWW_v_Wlnu"      ,";hWW term; Wlnu term"
                                        ,25,90,190,15,60,90)       ->Fill(approxH.mass(),approxW.mass(),weight);
                plotter.getOrMake1DPre(prefix,"hWW"      ,";hWW term"
                        ,50,90,190)       ->Fill(approxH.mass(),weight);
                plotter.getOrMake1DPre(prefix,"Wlnu"      ,";Wlnu term"
                       ,30,60,90)       ->Fill(approxW.mass(),weight);
            }

            else {
                plotter.getOrMake2DPre(prefix,"hWW_v_Wlnu"      ,";hWW term; Wlnu term"
                        ,40,110,150,15,0,60)       ->Fill(approxH.mass(),approxW.mass(),weight);
                plotter.getOrMake1DPre(prefix,"hWW"      ,";hWW term"
                        ,80,110,150)       ->Fill(approxH.mass(),weight);
                plotter.getOrMake1DPre(prefix,"Wlnu"      ,";Wlnu term"
                       ,30,0,60)       ->Fill(approxW.mass(),weight);
            }
        };

        TString ptBin;
        if(true_hwwMag >400 && true_hwwMag< 600) ptBin = "_low";
        else if(true_hwwMag >800 && true_hwwMag< 1000) ptBin = "_test";
        else if(true_hwwMag >1400 && true_hwwMag< 2000) ptBin = "_high";


        if(!finalBinning){
            plt(pf);
            plt(pf + (isVirtualWqq ? "_vqq" : "_osqq"));
            plt(pf+ptBin);
            plt(pf +ptBin+ (isVirtualWqq ? "_vqq" : "_osqq"));
        } else {
            pltFinal(pf + (isVirtualWqq ? "_vqq" : "_osqq"));
            pltFinal(pf +ptBin+ (isVirtualWqq ? "_vqq" : "_osqq"));
        }

        plotter.getOrMake1DPre(pf +ptBin,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(true_hwwMag,weight);
        plotter.getOrMake1DPre(pf,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(true_hwwMag,weight);

    }


    void getBKGTemplatePlots(bool finalBinning) {
        if(hwwMag < 300) return;
        double extraMetPerp = metPerp-neutPerp;
        double extraMetPar = metPar-neutPar;

        ASTypes::CylLorentzVectorF approxJ (qqJet.pt(),qqJet.eta(),qqJet.phi(),
                0);
        auto approxW = lepton.p4() + true_neut.p4();
        auto approxH = approxW + approxJ;

        TString pf = (isSignal() ? TString("signal") : smpName);
        TString comb = (*process==FillerConstants::TTBAR || *process == FillerConstants::WJETS)
                            ? ("ttbarPW") : "other";
        if(*process==FillerConstants::TTBAR && hwwMag > 1500) comb = "other";


        auto plt =[&](const TString& prefix){

            plotter.getOrMake1DPre(prefix,"hwwMag",";hwwMag",300,0,3000)
                ->Fill(hwwMag,weight);

            plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",1000,-500,500)
                ->Fill(extraMetPerp,weight);
            plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar / hwwMag",400,-2,2)
                ->Fill(extraMetPar/hwwMag,weight);

            plotter.getOrMake1DPre(prefix,"qqSDMassCoarse",";qqSDMass",2,50,70)
                ->Fill(*qqJet_SDmass,weight);
            plotter.getOrMake1DPre(prefix
                    , "approxW",";approxW",300,0,300)->Fill(approxW.mass(),weight);
            plotter.getOrMake1DPre(prefix, "approxH",";approxH",600,0,600)
                ->Fill(approxH.mass(),weight);

        };


        auto pltFinal =[&](const TString& prefix){



            plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",30,-150,150)
                ->Fill(extraMetPerp,weight);
            plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar/hwwMag",50,-0.6,0.4)
                ->Fill(extraMetPar/hwwMag,weight);


            plotter.getOrMake1DPre(prefix,"hWW"      ,";hWW term"
                    ,50,100,600)       ->Fill(approxH.mass(),weight);
            plotter.getOrMake1DPre(prefix,"Wlnu"      ,";Wlnu term"
                   ,20,60,100)       ->Fill(approxW.mass(),weight);
        };

        TString ptBin;
        if(hwwMag >300 && hwwMag< 500) ptBin = "_low";
        else if(hwwMag >800 && hwwMag< 900) ptBin = "_test";
        else if(hwwMag >1000 && hwwMag < 1500) ptBin = "_high";
        else ptBin = "_other";


        if(!finalBinning){
            plt(pf);
            plt(pf+ptBin);
            plt(comb);
            plt(comb+ptBin);

        } else {
            pltFinal(pf);
            pltFinal(pf +ptBin);
            pltFinal(comb);
            pltFinal(comb +ptBin);
        }

        plotter.getOrMake1DPre(pf +ptBin,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(pf,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(comb +ptBin,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(comb,"true_hwwMag",";hwwMag",300,0,3000)
            ->Fill(hwwMag,weight);

    }



    void getQCDTemplatePlots(bool finalBinning) {
//         if(hwwMag < 300) return;

         ASTypes::CylLorentzVectorF approxJ (qqJet.pt(),qqJet.eta(),qqJet.phi(),
                 0);

         ASTypes::CylLorentzVectorF aproxN (met.pt(),lepton.eta(),met.phi(),0);

         auto approxW = metPar > 0 ? lepton.p4() + aproxN : lepton.p4();
         auto approxH = approxW + qqJet.p4();

         auto approxW2 = lepton.p4();
         auto approxH2 = approxW2 + qqJet.p4();

         auto hww2ParX      = qqJet.px() + lepton.px();
         auto hww2ParY      = qqJet.py() + lepton.py();
         auto hww2ParNormX  = hww2ParX/approxH2.pt();
         auto hww2ParNormY  = hww2ParY/approxH2.pt();
         auto getPar =[&](const MomentumF& mom)->double{
             return mom.px()*hww2ParNormX+mom.py()*hww2ParNormY;
         };
         auto metPar2       = getPar(met);
         auto approxW3 = metPar2 > 0 ? lepton.p4() + aproxN : lepton.p4();
         auto approxH3 = approxW3 + qqJet.p4();

         TString pf = smpName;


         auto mkTL = [](const ASTypes::CylLorentzVectorF& in) -> TLorentzVector{
             return TLorentzVector(in.px(),in.py(),in.pz(),in.E());
         };

         auto mkBPl = [&](const TLorentzVector& mom, const TLorentzVector&d1, const TLorentzVector& d2, const TString& name){
             TVector3 bv = -1*mom.BoostVector();
             TLorentzVector bd1 = d1;TLorentzVector bd2 = d2;
             bd1.Boost(bv);
             bd2.Boost(bv);
             plotter.getOrMake1DPre(name,"deltaPhi",";deltaPhi",100,-3.2,3.2)->Fill(PhysicsUtilities::deltaPhi(mom.Phi(),bd1.Phi()),weight);
             plotter.getOrMake1DPre(name,"deltaTheta",";deltaTheta",100,-3.2,3.2)->Fill(PhysicsUtilities::deltaPhi(mom.Theta(),bd1.Theta()),weight);
             plotter.getOrMake1DPre(name,"deltaPhi2",";deltaPhi",100,-3.2,3.2)->Fill(PhysicsUtilities::deltaPhi(mom.Phi(),bd2.Phi()),weight);
             plotter.getOrMake1DPre(name,"deltaTheta2",";deltaTheta",100,-3.2,3.2)->Fill(PhysicsUtilities::deltaPhi(mom.Theta(),bd2.Theta()),weight);

         };


         auto approxHTL = mkTL(approxH);
         auto jTL = mkTL(qqJet.p4());
         auto lTL = mkTL(lepton.p4());


         auto plt =[&](const TString& prefix){

             plotter.getOrMake1DPre(prefix,"hwwMag",";hwwMag",300,0,3000)
                 ->Fill(hwwMag,weight);

             plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",1000,-500,500)
                 ->Fill(metPerp,weight);
             plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar / hwwMag",400,-2,2)
                 ->Fill(metPar/hwwMag,weight);

             plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag2",";extraMetPar / hwwMag",400,-2,2)
                 ->Fill(metPar2/approxH2.pt(),weight);

             plotter.getOrMake1DPre(prefix,"wPTRelhwwMag",";wPTRelhwwMag / hwwMag",400,-2,2)
                 ->Fill(approxW.pt()/hwwMag,weight);

             plotter.getOrMake1DPre(prefix,"qqSDMassCoarse",";qqSDMass",2,50,70)
                 ->Fill(*qqJet_SDmass,weight);
             plotter.getOrMake1DPre(prefix
                     , "approxW",";approxW",300,0,300)->Fill(approxW.mass(),weight);
             plotter.getOrMake1DPre(prefix
                     , "approxW2",";approxW",300,0,300)->Fill(approxW2.mass(),weight);

             plotter.getOrMake1DPre(prefix
                     , "ptDR",";ptDR",300,0,300)->Fill(approxH.pt()*PhysicsUtilities::deltaR(qqJet,approxW),weight);
             plotter.getOrMake1DPre(prefix
                     , "ptDR2",";ptDR2",300,0,300)->Fill(approxH2.pt()*PhysicsUtilities::deltaR(qqJet,approxW2),weight);

             plotter.getOrMake1DPre(prefix
                     , "DR",";ptDR",100,0,5)->Fill(PhysicsUtilities::deltaR(qqJet,approxW),weight);
             plotter.getOrMake1DPre(prefix
                     , "DR2",";ptDR2",100,0,5)->Fill(PhysicsUtilities::deltaR(qqJet,approxW2),weight);

             plotter.getOrMake1DPre(prefix, "approxH",";approxH",600,0,600)
                 ->Fill(approxH.mass(),weight);
             plotter.getOrMake1DPre(prefix, "approxH2",";approxH",600,0,600)
                 ->Fill(approxH2.mass(),weight);
             plotter.getOrMake1DPre(prefix, "approxH3",";approxH",600,0,600)
                 ->Fill(approxH3.mass(),weight);

             mkBPl(approxHTL,jTL,lTL,prefix);



         };


         auto pltFinal =[&](const TString& prefix){



             plotter.getOrMake1DPre(prefix,"extraMetPerp",";extraMetPerp",30,-150,150)
                 ->Fill(extraMetPerp,weight);
             plotter.getOrMake1DPre(prefix,"extraMetParRelhwwMag",";extraMetPar/hwwMag",50,-0.6,0.4)
                 ->Fill(extraMetPar/hwwMag,weight);


             plotter.getOrMake1DPre(prefix,"hWW"      ,";hWW term"
                     ,50,100,600)       ->Fill(approxH.mass(),weight);
             plotter.getOrMake1DPre(prefix,"Wlnu"      ,";Wlnu term"
                    ,20,60,100)       ->Fill(approxW.mass(),weight);
         };

         TString ptBin;
         if(hwwMag >300 && hwwMag< 500) ptBin = "_low";
         else if(hwwMag >800 && hwwMag< 900) ptBin = "_test";
         else if(hwwMag >1000 && hwwMag < 1500) ptBin = "_high";
         else ptBin = "_other";


         if(!finalBinning){
             plt(pf);
             plt(pf+ptBin);


         } else {
             pltFinal(pf);
             pltFinal(pf +ptBin);

         }

         plotter.getOrMake1DPre(pf +ptBin,"true_hwwMag",";hwwMag",300,0,3000)
             ->Fill(hwwMag,weight);
         plotter.getOrMake1DPre(pf,"true_hwwMag",";hwwMag",300,0,3000)
             ->Fill(hwwMag,weight);

     }


    static void processTemplates(const TString& inFileName, const TString& outFileName,
            bool isSignal) {
        TFile * fIn = new TFile(inFileName,"read");
        std::vector<TH1*> outHs;
        std::vector<TString> ptCats    ;
        std::vector<TString> shells  ;
        std::vector<TString> ext1DVs ;
        std::vector<TString> std1DVs ;
        std::vector<TString> twoDVs  ;

        TString procName;

        if(isSignal){
            procName = "signal";
            ptCats    = {"low_","high_"};
            shells  = {"vqq_","osqq_"};
            ext1DVs = {"extraMetPerp","extraMetParRelhwwMag","wqqPTRes","hWW","Wlnu"};
            std1DVs = {"qqSDMassCoarse"};
            twoDVs = {"signal_vqq_hWW_v_Wlnu","signal_osqq_hWW_v_Wlnu"};
        } else {
            procName = "ttbarPW";
            ptCats   = {"low_","high_"};
            shells  = {""};
            ext1DVs = {"extraMetPerp","extraMetParRelhwwMag","hWW","Wlnu"};
        }

        // start witht the 1D histograms
        auto do1D =[&](const TString& hName, bool extrap){
            TH1 * h = 0;
            fIn->GetObject(hName,h);
            if(!h)
                throw std::invalid_argument(std::string(
                        "processTemplates() -> Couldn't find the histogram! "));
            int nBX = h->GetNbinsX();
            h = (TH1*)h->Clone();
            h->Scale(1./h->Integral(0,-1));
            if(extrap){
                h->SetBinContent(0,   h->GetBinContent(0) + 0.5*h->GetBinContent(1));
                h->SetBinContent(nBX+1,h->GetBinContent(nBX+1) +
                        0.5*h->GetBinContent(nBX));
                for(int iB = 1; iB <= nBX; ++iB){
                    h->SetBinContent( iB, h->GetBinContent(iB)/h->GetBinWidth(iB) );
                    h->SetBinError( iB, h->GetBinError(iB)/h->GetBinWidth(iB) );
                }
            } else {
                double low = h->GetBinContent(0)+ h->GetBinContent(1);
                double high = h->GetBinContent(nBX)+ h->GetBinContent(nBX+1);
                h->SetBinContent(0, low);
                h->SetBinContent(1, low);
                h->SetBinContent(nBX, high);
                h->SetBinContent(nBX+1, high);
            }
            outHs.push_back(h);
        };


        for(unsigned int iS = 0; iS < ptCats.size(); ++iS){
            for(unsigned int iSh = 0; iSh < shells.size(); ++iSh){
                for(unsigned int iV = 0; iV < ext1DVs.size(); ++iV){
                    do1D( procName + "_"+ ptCats[iS]+shells[iSh]+ext1DVs[iV],true);
                }
                for(unsigned int iV = 0; iV < std1DVs.size(); ++iV){
                    do1D(ptCats[iS]+"_"+shells[iSh]+std1DVs[iV],true);
                }
            }
        }


        //Add in the average mass plots
        TH1F * ptH = new TH1F(procName + "_avgHWWMag",";avg hWWMag",
                ptCats.size(), -0.5, float(ptCats.size())-0.5);
        for(unsigned int iS = 0; iS < ptCats.size(); ++iS){
            TH1 * h = 0;
            fIn->GetObject(procName +"_" + ptCats[iS]+"true_hwwMag",h);
            if(!h)
                throw std::invalid_argument(std::string(
                        "processTemplates() -> Couldn't find the histogram! "));
            const double mean = h->GetMean();
            ptH->SetBinContent(iS+1,mean);
        }
        outHs.push_back(ptH);

        for(unsigned int iV = 0; iV < twoDVs.size(); ++iV){
            TH2 * h = 0;
            fIn->GetObject(twoDVs[iV],h);
            if(!h)
                throw std::invalid_argument(std::string(
                        "processTemplates() -> Couldn't find the histogram! "));
            h = (TH2*)h->Clone();

            //Fill in the holes first
            for(int iY = 1; iY <= h->GetNbinsY(); ++iY){
                int fFilledBin = -1;
                for(int iX = 1; iX <= h->GetNbinsX(); ++iX){
                    double binC = h->GetBinContent(iX,iY);
                    if(fFilledBin < 0){
                        if(binC <= 0) continue;
                        fFilledBin = binC;
                        if(iX == 1) continue;
                        double val = h->GetBinContent(iX,iY);
                        h->SetBinContent(iX,iY,0.9*val);
                        double totNew = 0;
                        for(int iX2 = 1; iX2 < iX; ++iX2){
                            double newV = BASEPDF::extapDown(
                                    h->GetXaxis()->GetBinCenter(iX2),
                                    h->GetXaxis()->GetBinCenter(iX),
                                    0.9*val,0.1*5*val);
                            h->SetBinContent(iX2,iY,newV);
                            totNew += newV;

                        }

                        for(int iX2 = 1; iX2 < iX; ++iX2){
                            h->SetBinContent(iX2,iY,
                                    (0.1*val/totNew)* h->GetBinContent(iX2,iY));
                        }
                        continue;
                    }
                    if(binC <= 0){
                        double newV = 0;
                        if(int(iV) == h->GetNbinsX())
                            newV = h->GetBinContent(iX-1,iY);
                        else {
                            newV = 0.5*h->GetBinContent(iX-1,iY)
                                                    + 0.5*h->GetBinContent(iX+1,iY);
                        }
                        h->SetBinContent(iX,iY,newV);
                    }
                }
            }

            //Now do the edges
            h->Scale(1.0/h->Integral(0,-1,0,-1));
            int nBX = h->GetNbinsX();
            int nBY = h->GetNbinsY();
            for(int iX = 1; iX <= nBX; ++iX){
                h->SetBinContent(iX,0,   h->GetBinContent(iX,0) + 0.5*h->GetBinContent(iX,1));
                h->SetBinContent(iX,nBY+1,h->GetBinContent(iX,nBY+1) +
                        0.5*h->GetBinContent(iX,nBY));
            }
            for(int iY = 1; iY <= nBY; ++iY){
                h->SetBinContent(0,iY,   h->GetBinContent(0,iY) + 0.5*h->GetBinContent(1,iY));
                h->SetBinContent(nBX+1,iY,h->GetBinContent(nBX+1,iY) +
                        0.5*h->GetBinContent(nBX,iY));
            }
            for(int iX = 1; iX <= nBX; ++iX){
                for(int iY = 1; iY <= nBY; ++iY){
                    double w = h->GetXaxis()->GetBinWidth(iX)*h->GetYaxis()->GetBinWidth(iY);
                    h->SetBinContent(iX,iY, h->GetBinContent(iX,iY)/w);
                }
            }
            outHs.push_back(h);
        }


        TFile * fo = new TFile(outFileName,"recreate");
        for(auto * h : outHs) h->Write();
        fo->Close();
        fIn->Close();
    }

//    static void finalizeTemplates(const TString& inFileName, const TString& outFileName,
//            bool isSignal) {
//
//            if(isSignal){
//                HSolverLi signalSolver("");
//                signalSolver.setParamters(parameters.hww);
//
//
//            }
//
//    }


    void getInputPlots() {






        plotter.getOrMake1DPre(smpName,"hwwMag",";hwwMag",300,0,3000)->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"true_hwwMag",";hwwMag",300,0,3000)
                ->Fill(true_H.pt(),weight);

        plotter.getOrMake1DPre(smpName,"hwwMagRes",";hwwMagRes",100,-2,2)
                ->Fill((hwwMag-true_H.pt())/true_H.pt(),weight);

        plotter.getOrMake1DPre(smpName,"met",";met",500,0,1000)->Fill(met.pt(),weight);

        plotter.getOrMake1DPre(smpName,"metPerp",";metPerp",500,0,1000)
                ->Fill(std::fabs(metPerp),weight);
        plotter.getOrMake1DPre(smpName,"metPar",";metPar",500,0,1000)
                ->Fill(std::fabs(metPar),weight);

        plotter.getOrMake1DPre(smpName,"neutPerp",";neutPerp",500,0,2000)
                ->Fill(std::fabs(neutPerp),weight);
        plotter.getOrMake1DPre(smpName,"neutPar",";neutPar",500,0,2000)
                ->Fill(std::fabs(neutPar),weight);
        plotter.getOrMake1DPre(smpName,"neutZ",";neutPar",500,0,2000)
                ->Fill(std::fabs(true_neut.pz()),weight);

        plotter.getOrMake1DPre(smpName,"neutPerpRelhwwMag",";neutPerp / hwwMag",500,-2,2)
                ->Fill(neutPerp/hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"neutParRelhwwMag",";neutPar / hwwMag",500,-2,2)
                ->Fill(neutPar/hwwMag,weight);

        plotter.getOrMake1DPre(smpName,"neutZRelhwwMag",";neutZ / hwwMag",500,-10,10)
                ->Fill(true_neut.pz()/hwwMag,weight);

        double extraMetPerp = metPerp-neutPerp;
        double extraMetPar = metPar-neutPar;

        plotter.getOrMake1DPre(smpName,"extraMetPerp",";extraMetPerp",500,-500,500)
                ->Fill(extraMetPerp,weight);
        plotter.getOrMake1DPre(smpName,"extraMetPar",";extraMetPar",500,-500,500)
                ->Fill(extraMetPar,weight);

        plotter.getOrMake1DPre(smpName,"extraMetPerpRelhwwMag",";extraMetPerp / hwwMag",500,-2,2)
                ->Fill(extraMetPerp/hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"extraMetParRelhwwMag",";extraMetPar / hwwMag",500,-2,2)
                ->Fill(extraMetPar/hwwMag,weight);

        double wqqDPhi = PhysicsUtilities::deltaPhi(true_jet,qqJet);
        double wqqDR   = PhysicsUtilities::deltaR(true_jet,qqJet);
        double wqqPTDif = qqJet.pt()-true_jet.pt();


        double wqqPTRes = wqqPTDif/true_jet.pt();



        if(true_jet.pt()>10){
            plotter.getOrMake1DPre(smpName,"wqqDPhi",";wqqDPhi",500,-3,3)
                ->Fill(wqqDPhi,weight);
            plotter.getOrMake1DPre(smpName,"wqqDR",";wqqDR",500,0,3)
                    ->Fill(wqqDR,weight);
            plotter.getOrMake1DPre(smpName,"wqqPTDif",";wqqPTDif",500,-2000,2000)
                    ->Fill(wqqPTDif,weight);
            plotter.getOrMake1DPre(smpName,"wqqPTRes",";wqqPTRes",500,-2,2)
                ->Fill(wqqPTRes,weight);
            plotter.getOrMake1DPre(smpName,"wqqPTTrue",";wqqPTTrue",500,0,2000)
                ->Fill(true_jet.pt(),weight);
            plotter.getOrMake1DPre(smpName,"wqqPTReco",";wqqPTReco",500,0,2000)
                ->Fill(qqJet.pt(),weight);

        }


        plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",7,-0.5,6.5)
                ->Fill(0.0,weight);

        if(*qqJet_SDmass > 60){
            plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                ->Fill(1.0,weight);
            if(isVirtualWqq)
                plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                ->Fill(2.0,weight);
            else
                plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                ->Fill(3.0,weight);
        } else {
            plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                    ->Fill(4.0,weight);
            if(isVirtualWqq)
                plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                ->Fill(5.0,weight);
            else
                plotter.getOrMake1DPre(smpName,"shellPurity",";shellPurity",5,-0.5,6.5)
                ->Fill(6.0,weight);
        }

        auto plt =[&](const TString& prefix){
            plotter.getOrMake1DPre(prefix,"qqSDMass",";qqSDMass",300,0,300)
                ->Fill(*qqJet_SDmass,weight);
            plotter.getOrMake1DPre(prefix,"qqSDMassCoarse",";qqSDMass",2,-0.5,1.5)
                    ->Fill(*qqJet_SDmass >= 60,weight);
            auto approxW = lepton.p4() + true_neut.p4();
            plotter.getOrMake1DPre(prefix
                    , "approxW",";approxW",300,0,300)->Fill(approxW.mass(),weight);
            double sf =  true_jet.pt()/qqJet.pt();

            ASTypes::CylLorentzVectorF approxJ (sf*qqJet.pt(),qqJet.eta(),qqJet.phi(),
                    isVirtualWqq ? 31 : 80.0);
            auto approxH = approxW + approxJ;

            plotter.getOrMake1DPre(prefix, "approxH",";approxH",300,0,300)
                    ->Fill(approxH.mass(),weight);

            //correlations

            plotter.getOrMake2DPre(prefix,"hWW_v_Wlnu"      ,";hWW term; Wlnu term"
                    ,100,80,190,100,0,200)       ->Fill(approxH.mass(),approxW.mass(),weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_jetSF"     ,";hWW term; JetSF term"
                    ,100,80,190,50,-0.5,0.5)       ->Fill(approxH.mass(),wqqPTRes,weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_metPerp"   ,";hWW term; MET perp term"
                    ,100,80,190,100,-200,200)    ->Fill(approxH.mass(),extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"hWW_v_metPar"    ,";hWW term; MET par term"
                    ,100,80,190,100,-2,2)     ->Fill(approxH.mass(),extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_jetSF"    ,";Wlnu term; JetSF term"
                    ,100,0,200,100,0,2)      ->Fill(approxW.mass(),wqqPTRes,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_metPerp"  ,";Wlnu term; MET perp term"
                    ,100,0,200,100,-200,200)   ->Fill(approxW.mass(),extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"Wlnu_v_metPar"   ,";Wlnu term; MET par term"
                    ,100,0,200,100,-2,2)    ->Fill(approxW.mass(),extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"jetSF_v_metPerp" ,";Jet SF term; MET perp term"
                    ,100,-2,2,100,-200,200) ->Fill(wqqPTRes,extraMetPerp,weight);
            plotter.getOrMake2DPre(prefix,"jetSF_v_metPar"  ,";Jet SF term; MET par term"
                    ,100,-2,2,100,-2,2)  ->Fill(wqqPTRes,extraMetPar/hwwMag,weight);
            plotter.getOrMake2DPre(prefix,"metPerp_v_metPar",";MET perp term; MET par term"
                    ,100,-200,200,100,-2,2)->Fill(extraMetPerp,extraMetPar/hwwMag,weight);

        };
        plt(smpName + (isVirtualWqq ? "_virtualSD" : "_onshellSD"));
        plt(smpName);
    }


    void testSol(bool doPrintouts) {
        if(!isSignal()){
            if(*bbJet_SDmass < 30) return;
            if(*bbJet_SDmass > 210) return;
            if(*nAK4Btags > 0) return;
            if(*qqJet_t2ot1>0.75) return;
        }

        auto print = [&](const ASTypes::CylLorentzVectorF& neut,
                const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& wlnu,
                const ASTypes::CylLorentzVectorF& hWW, double mass){

            std::cout << "N: ("<<neut.pt()<<","<<neut.eta()<<","<<neut.phi()<<") J:("
                    <<jet.pt()<<","<<jet.eta()<<","<<jet.phi()<<","<<jet.mass()<<") W:("
                    <<wlnu.pt()<<","<<wlnu.eta()<<","<<wlnu.phi()<<","<<wlnu.mass()<<") H:("
                    <<hWW.pt()<<","<<hWW.eta()<<","<<hWW.phi()<<","<<hWW.mass()<<") HH:("
                    <<mass<<")"<<std::endl;
        };

        if(doPrintouts){
            std::cout <<" ---------------------------- "<<std::endl;
            std::cout <<"NEVENT : "<< *sampParam<<std::endl <<" "<< isVirtualWqq <<std::endl;
            //true
            print(true_neut.p4(),true_jet.p4(),true_W.p4(), true_H.p4(),
                    (bbJet.p4()+true_H.p4()).mass() );
            //simple
            auto simN = HSolverBasic::getInvisible(met,lepton.p4()+qqJet.p4());
            print(simN.p4(),qqJet.p4(),lepton.p4()+simN.p4(), lepton.p4()+simN.p4() + qqJet.p4(),
                    (bbJet.p4()+lepton.p4()+simN.p4() + qqJet.p4()).mass() );
            //older

            HSolverChiInfo ohwwInfo;
            double ohwwChi   = oHSolver.hSolverMinimization(lepton.p4(),qqJet.p4(),
                    met.p4(),*qqJet_SDmass <60,parameters.hww, &ohwwInfo);
            std::cout <<ohwwChi<<" -> ";
            print(ohwwInfo.neutrino,ohwwInfo.wqqjet,ohwwInfo.wlnu, ohwwInfo.hWW,
                    (bbJet.p4()+ohwwInfo.hWW).mass() );



            HSolverLiInfo hwwInfo;
            HSolverLiInfo hwwInfo_vqq;
            HSolverLiInfo hwwInfo_osqq;

            HSolver.minimize(lepton,met,qqJet,*qqJet_SDmass,hwwInfo, &hwwInfo_osqq,
                    &hwwInfo_vqq);
            print(hwwInfo.neutrino,hwwInfo.wqqjet,hwwInfo.wlnu,
                    hwwInfo.hWW,(bbJet.p4()+hwwInfo.hWW).mass() );

            std::cout << hwwInfo_osqq.likeli <<" "<< hwwInfo_osqq.noSDLikli<<" "
                    << hwwInfo_vqq.likeli<<" "<< hwwInfo_vqq.noSDLikli<<std::endl;

//            std::cout <<"("<< hwwInfo.likeli<<") ("<<hwwInfo.emetperp<<","<<
//                    hwwInfo.emetpar<<","<<hwwInfo.neutrino.pz()<<","<<hwwInfo.ptRes<<") (";

//            if(hwwInfo.min_sol.likeli == hwwInfo.vqq_sol.likeli){
//                std:: cout <<
//                        -2.*std::log(HSolver.vqq_pdfs[HiggsLi::EMET_PERP]->getProbability(hwwInfo.min_sol.emetperp))
//                <<","<<-2.*std::log(HSolver.vqq_pdfs[HiggsLi::EMET_PAR]->getProbability(hwwInfo.min_sol.emetpar))
//                <<","<<-2.*std::log(HSolver.vqq_pdfs[HiggsLi::WQQ_RES]->getProbability(hwwInfo.min_sol.ptRes))
//                <<","<<-2.*std::log(HSolver.vqq_pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hwwInfo.min_sol.hWW.mass(),hwwInfo.min_sol.wlnu.mass()))
//                <<") ";
//            } else {
//                std:: cout <<
//                        -2.*std::log(HSolver.osqq_pdfs[HiggsLi::EMET_PERP]->getProbability(hwwInfo.min_sol.emetperp))
//                <<","<<-2.*std::log(HSolver.osqq_pdfs[HiggsLi::EMET_PAR]->getProbability(hwwInfo.min_sol.emetpar))
//                <<","<<-2.*std::log(HSolver.osqq_pdfs[HiggsLi::WQQ_RES]->getProbability(hwwInfo.min_sol.ptRes))
//                <<","<<-2.*std::log(HSolver.osqq_pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hwwInfo.min_sol.hWW.mass(),hwwInfo.min_sol.wlnu.mass()))
//                <<") ";
//            }





//            print(hwwInfo.osqq_sol.neutrino,hwwInfo.osqq_sol.wqqjet,hwwInfo.osqq_sol.wlnu,
//                    hwwInfo.osqq_sol.hWW,(bbJet.p4()+hwwInfo.osqq_sol.hWW).mass() );
//
//            print(hwwInfo.vqq_sol.neutrino,hwwInfo.vqq_sol.wqqjet,hwwInfo.vqq_sol.wlnu,
//                    hwwInfo.vqq_sol.hWW,(bbJet.p4()+hwwInfo.vqq_sol.hWW).mass() );
        }
        HSolverLiInfo hwwInfo;
        HSolverLiInfo osqqHWWInfo;
        HSolverLiInfo vqqHWWInfo;
        HSolverLiInfo altHWWInfo;
        HSolver.minimize(lepton,met,qqJet,*qqJet_SDmass,hwwInfo,&osqqHWWInfo,&vqqHWWInfo,&altHWWInfo);

        HSolverLiInfo hwwBkgInfo;
        HSolverLiInfo hwwBkgAltInfo;
        BkgHSolver.minimize(lepton,met,qqJet,hwwBkgInfo,&hwwBkgAltInfo);


        auto bN = HSolverBasic::getInvisible(met, (qqJet.p4()+lepton.p4()));


        const double lHH = (hwwInfo.hWW + bbJet.p4()).mass();


        const double nBL = hwwBkgInfo.likeli/hwwBkgAltInfo.likeli;

        const double ptom = (bN.p4() + lepton.p4() + qqJet.p4()).pt() / *hh_orig;
        const double ptom2 = hwwInfo.hWW.pt() / lHH;

        const double md = PhysicsUtilities::deltaR( bN.p4() + lepton.p4() ,qqJet)
        * (bN.p4() + lepton.p4() + qqJet.p4()).pt()/2.0;

        plotter.getOrMake1DPre(smpName, "simpleHH",";HH",120,0,3000)
                ->Fill(*hh_orig,weight);
        plotter.getOrMake1DPre(smpName, "chi2HH",";HH",120,0,3000)
                ->Fill(*hh_chi2,weight);
        plotter.getOrMake1DPre(smpName, "likeliHH",";HH",120,0,3000)
                ->Fill(lHH,weight);
        plotter.getOrMake1DPre(smpName, "likeliBHH",";HH",120,0,3000)
                ->Fill((hwwBkgInfo.hWW + bbJet.p4()).mass(),weight);


        auto mkDiscPlts =[&](const TString& prefix){



            plotter.getOrMake1DPre(prefix, "likeli",";likeli",200,0,5)
                    ->Fill(hwwInfo.likeli,weight);

            plotter.getOrMake1DPre(prefix, "likeli2",";likeli",200,-50,100)
                    ->Fill(hwwInfo.rawLikeli-altHWWInfo.rawLikeli,weight);

            plotter.getOrMake1DPre(prefix, "likeli_raw",";likeli",200,0,100)
                    ->Fill(hwwInfo.rawLikeli,weight);

            plotter.getOrMake1DPre(prefix, "alt",";likeli",200,0,100)
                    ->Fill(altHWWInfo.rawLikeli,weight);

            plotter.getOrMake1DPre(prefix, "md",";MD",200,0,200)
                    ->Fill(md,weight);



            plotter.getOrMake1DPre(prefix, "Blikeli",";likeli",200,0,50)
                    ->Fill(hwwBkgInfo.likeli,weight);

            plotter.getOrMake1DPre(prefix, "Blikeli_nAlt",";likeli",200,0,5)
                    ->Fill(nBL,weight);

            plotter.getOrMake1DPre(prefix, "bAlt",";likeli",200,0,100)
                    ->Fill(hwwBkgAltInfo.likeli,weight);

            plotter.getOrMake1DPre(prefix, "SoBlikeli",";likeli",200,0,5)
                    ->Fill(nBL > 0 ? hwwInfo.likeli/nBL
                            : 50,weight);


            plotter.getOrMake1DPre(prefix, "ptom",";HH",250,0,1)
                    ->Fill( ptom,weight);
            plotter.getOrMake1DPre(prefix, "ptom2",";HH",250,0,1)
                    ->Fill( ptom2,weight);
        };




        mkDiscPlts(smpName);
        if(lHH>1000) mkDiscPlts(smpName +"_m1000");
        if(lHH>2000) mkDiscPlts(smpName +"_m2000");
        if(hwwInfo.likeli < 1.45)  mkDiscPlts(smpName +"_lllt1p45");
        if(hwwInfo.likeli < 1.45 && lHH>1000)  mkDiscPlts(smpName +"_m1000_lllt1p45");
        if(hwwInfo.likeli < 1.45 && lHH>2000)  mkDiscPlts(smpName +"_m2000_lllt1p45");

        if(md < 125 && lHH>1000)  mkDiscPlts(smpName +"_m1000_mdlt125");
        if(md < 125 && lHH>2000)  mkDiscPlts(smpName +"_m2000_mdlt125");

        plotter.getOrMake1DPre(smpName, "chi2",";chi2",100,0,20)
                ->Fill(*chi2,weight);
        if(*hh_chi2>1000)
          plotter.getOrMake1DPre(smpName, "m1000_chi2",";chi2",100,0,20)
                  ->Fill(*chi2,weight);
        if(*hh_chi2>2000)
        plotter.getOrMake1DPre(smpName, "m2000_chi2",";chi2",100,0,20)
                ->Fill(*chi2,weight);





        auto plotForDCuts = [&](const TString& selName){
            if(isSignal())
                plotter.getOrMake1DPre("signal"+selName, "sampParam",";sampParam",35,500,4000)
                    ->Fill(*sampParam,weight);
            else {
                plotter.getOrMake1DPre(smpName+selName, "chi2HH",";HH",120,0,3000)
                        ->Fill(*hh_chi2,weight);
                plotter.getOrMake1DPre(smpName+selName, "likeliHH",";HH",120,0,3000)
                        ->Fill(lHH,weight);

                plotter.getOrMake1DPre(smpName+selName, "simpleHH",";HH",120,0,3000)
                        ->Fill(*hh_orig,weight);

            }

        };
        plotForDCuts("_chiltinf");
        plotForDCuts("_llltinf");
        plotForDCuts("_mdltinf");
        if(*chi2 < 1.75) plotForDCuts("_chilt1p75");
        if(hwwInfo.likeli < 1.05) plotForDCuts("_lllt1p05");
        if(hwwInfo.likeli < 1.1) plotForDCuts("_lllt1p1");
        if(hwwInfo.likeli < 1.2) plotForDCuts("_lllt1p2");
        if(hwwInfo.likeli < 1.25) plotForDCuts("_lllt1p25");

        if(*chi2 < 11) plotForDCuts("_chilt11");
        if(hwwInfo.likeli < 1.35) plotForDCuts("_lllt1p35");
        if(*chi2 < 19) plotForDCuts("_chilt19");
        if(hwwInfo.likeli < 1.45) plotForDCuts("_lllt1p45");
        if(hwwInfo.likeli < 1.6) plotForDCuts("_lllt1p6");
        if(md < 125) plotForDCuts("_mdlt125");

    }


    bool runEvent() override {
        if(!HHSolTreeAnalyzer::runEvent()) return false;
//        if(*qqJet_t2ot1 >= 0.75) return false;
        if(isSignal() && (*bbJet_SDmass < 100 || *bbJet_SDmass>150)) return false;
        if(isSignal() && *wqqDR >= 0.7) return false;
        if(isSignal() && (PhysicsUtilities::deltaR(true_jet,qqJet) > 0.4)) return false;
//        if(!isSignal() && *process!= FillerConstants::TTBAR)return false;
        if(step==0)
            getHHPtResolution();
        if(step==1)
            getInputPlots();
        if(step==2)
            getTemplatePlots(false);
        if(step==3)
            getTemplatePlots(true);
        if(step==5)
            testSol(false);
        if(step==6)
            getBKGTemplatePlots(false);
        if(step==7)
            getBKGTemplatePlots(true);
        if(step==9)
            getQCDTemplatePlots(false);
        if(step==10)
            getQCDTemplatePlots(true);


        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    int step;

    HSolverChi oHSolver;
    HSolverLi HSolver;
    HSolverBkgLi BkgHSolver;
    ParameterSet parameters;
};

#endif

void getHHTraining(int step, std::string fileName, std::string outFileName,
        float xSec=-1, float numEvent=-1){

    if(step == 4){
        Analyzer::processTemplates("hSolTrees_temp.root","hhSol_templates.root",true);
        return;
    }

    if(step == 8){
        Analyzer::processTemplates("hSolTrees_bkgTemp_bkg.root","hhSol_bkgTemplates.root",false);
        return;
    }

    Analyzer a(fileName,"treeMaker/Events",1,1,step);
    std::string outN = "hSolTrees_";


    switch(step){
    case 0 :
        outN += "HPTRes";
        break;
    case 1:
        outN += "testVars";
        break;
    case 2:
        outN += "tempStudy";
        break;
    case 3:
        outN += "temp";
        break;
    case 5:
        outN += "test";
        break;
    case 6 :
        outN += "bkgTempStudy";
        break;
    case 7 :
        outN += "bkgTemp";
        break;

    case 9 :
        outN += "qcdTempStudy";
        break;
    case 10 :
        outN += "qcdTemp";
        break;
    }





    outN += "_"+ outFileName+".root";

    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outN);


}
