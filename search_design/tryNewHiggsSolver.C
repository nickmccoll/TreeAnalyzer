#include <AnalysisSupport/Utilities/interface/Types.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TMath.h>
#include <TreeAnalyzer/framework/DataFormats/interface/FatJet.h>
#include <TreeAnalyzer/framework/DataFormats/interface/GenParticle.h>
#include <TreeAnalyzer/framework/DataFormats/interface/Lepton.h>
#include <TreeAnalyzer/framework/DataFormats/interface/Momentum.h>
#include <TreeAnalyzer/framework/Processors/GenTools/interface/DiHiggsEvent.h>
#include <TString.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

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
#include "TFitter.h"


using namespace TAna;
using namespace FillerConstants;

class HSolverInfo {
public:
  HSolverInfo()
{};
  double chiSq = -1;
  double SF = -1;
  ASTypes::CylLorentzVector neutrino;
  ASTypes::CylLorentzVector wlnu;
  ASTypes::CylLorentzVector hWW;
  ASTypes::CylLorentzVector wqqjet;
};


class HSolver {
 public:
  static const double mW    ;
  static const double mT    ;
  static const double aW    ;
  static const double aT    ;
  static const double aCM   ;
  double wLNUMass;
  double nZP;
  TFitter *minimizer;

  HSolver();
  ~HSolver();


  // operations
  static double hSolverFunction( const double leptonX, const double leptonY, const double leptonZ, const double neutrinoX, const double neutrinoY,const double neutrinoZ,
          const double jetX,    const double jetY,    const double jetZ,    const double jetM, const double jetSF,
          const double metX,    const double metY,const double wlnuM, HSolverInfo * info = 0
                   );

  static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);
  double hSolverMinimization(const ASTypes::CylLorentzVectorF& lep, const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met, bool jetIsVirtual, HSolverInfo * info);

  double round(double num, int x){
    return ceil( ( num * pow( 10,x ) ) - 0.5 ) / pow( 10,x );
  }

}; // end of class Topness


const double HSolver::mW  = 81.;
const double HSolver::mT  = 172.;
const double HSolver::aW  = 5.;
const double HSolver::aT  = 15.;
const double HSolver::aCM = 1000.;

HSolver::HSolver() : minimizer(new TFitter(4))
  {
  // less print-outs
    cout.precision(11);
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
  }

HSolver::~HSolver(){
  delete minimizer;
}

double HSolver::hSolverFunction(
                   const double leptonX, const double leptonY, const double leptonZ, const double neutrinoX, const double neutrinoY,const double neutrinoZ,
                   const double jetX,    const double jetY,    const double jetZ,    const double jetM, const double jetSF,
                   const double metX,    const double metY,const double wlnuM, HSolverInfo * info
                   ) {

    const double hwwParX = jetX+metX+leptonX;
    const double hwwParY = jetY+metY+leptonY;
    const double hwwMag = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);

    const double hwwParNormX  = hwwParX/hwwMag;
    const double hwwParNormY  = hwwParY/hwwMag;

    const double hwwPerpNormX  = -1*hwwParNormY;
    const double hwwPerpNormY  = hwwParNormX;

    auto getPerp =[&](const double momx, const double momy)->double{ return momx*hwwPerpNormX+momy*hwwPerpNormY;};
    auto getPar =[&](const double momx, const double momy)->double{ return momx*hwwParNormX+momy*hwwParNormY;};


    const double metPerp = getPerp(metX,metY);
    const double metPar = getPar(metX,metY);
    const double neutPerp = getPerp(neutrinoX,neutrinoY);
    const double neutPar = getPar(neutrinoX,neutrinoY);


    const double leptonE = std::sqrt(leptonX*leptonX+leptonY*leptonY+leptonZ*leptonZ);
    const ASTypes::CartLorentzVector lepton(leptonX,leptonY,leptonZ,leptonE);
    const double neutrinoE = std::sqrt(neutrinoX*neutrinoX+neutrinoY*neutrinoY+neutrinoZ*neutrinoZ);
    const ASTypes::CartLorentzVector neutrino(neutrinoX,neutrinoY,neutrinoZ,neutrinoE);
    const double jetE   = std::sqrt(jetSF*jetX*jetSF*jetX+jetSF*jetY*jetSF*jetY+jetSF*jetZ*jetSF*jetZ+jetM*jetM);
    const ASTypes::CartLorentzVector jet(jetSF*jetX,jetSF*jetY,jetSF*jetZ,jetE);

    const ASTypes::CartLorentzVector wlnu = lepton + neutrino;
    const ASTypes::CartLorentzVector hww = wlnu + jet;

    const double metParOff = -0.076-0.0000165*hwwMag;
    const double extraMetPar = (metPar-neutPar)/hwwMag;// - metParOff;
    const double extraMetPerp = (metPerp-neutPerp);

    const double metParErr =   extraMetPar >= 0 ? 0.064 : 0.15;

    const double metParChiSq =  extraMetPar*extraMetPar/(metParErr*metParErr);

    const double metPerpError = 29;//26.4+0.0175*hwwMag;
    const double metPerpChiSq = extraMetPerp*extraMetPerp/(metPerpError*metPerpError);
    const double jetSFChiSq = (jetSF-1.0)*(jetSF-1.0)/(0.10*0.10);

    const double meanWlnuMass = jetM > 60 ? 43 : 80;
    const double wlnuMassError = jetM > 60 ? (wlnu.mass() > 43 ? 3 : 19) : 3;

    const double wlnuMassChiSq = std::pow( (wlnu.mass() -meanWlnuMass)/wlnuMassError,2  );
//    const double wlnuMassChiSq = std::pow( (wlnu.mass() -wlnuM)/wlnuMassError,2  );

    const double meanHWWMass = 125;
    const double hWWMassError = jetM > 60 ? 3 : 9;
    const double hWWMassChiSq = std::pow( (hww.mass() -meanHWWMass)/hWWMassError,2  );

    const double thetaErr = std::max(0.06,0.283-0.000153*hwwMag);
    const double thetaDiff =  neutrino.theta() - hww.theta();
    const double thetaChiSq = thetaDiff*thetaDiff/(thetaErr*thetaErr);


//    plotter.getOrMake1DPre(smpName,"neutThetamOT",";neutThetamOT",600,-5,5)->Fill( diHiggsEvt.w1_d2->theta() - (selectedLepton->p4() + wjjCand->p4()).theta(),weight);


//    jetZ+leptonZ

//    std::cout << neutrinoX <<","<<neutrinoY <<","<< neutrinoZ<<","<<jetSF <<","<<wlnu.mass()<<","<<hww.mass()<<": "
//            <<metPerpChiSq <<", "<<metParChiSq <<", "<<jetSFChiSq <<", "<<wlnuMassChiSq <<", "<<hWWMassChiSq << std::endl;

//    const double chisq= metPerpChiSq+metParChiSq+jetSFChiSq+wlnuMassChiSq+hWWMassChiSq;
//    const double chisq= metPerpChiSq+metParChiSq+wlnuMassChiSq+thetaChiSq   ;

    const double chisq= hWWMassChiSq +wlnuMassChiSq +metPerpChiSq +metParChiSq+jetSFChiSq ;




//

  if(info){
      info->chiSq = chisq;
      info->SF = jetSF;
      info->neutrino = neutrino;
      info->wlnu = wlnu;
      info->hWW = hww;
      info->wqqjet = jet;
  }

  return chisq;

} // ~ end of Topness function



void HSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par, int flg) {


  result = hSolverFunction(par[0],par[1],par[2],par[3],
               par[4],par[5],par[6],par[7],
               par[8],par[9],par[10],par[11],
               par[12],par[13]
               );

} // ~end of minuit function



double HSolver::hSolverMinimization(const ASTypes::CylLorentzVectorF& lep, const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met, bool jetIsVirtual, HSolverInfo * info) {
    const double jetM =jetIsVirtual ? 31 : 80;



    minimizer->SetParameter(0,"leptonX"    ,lep.x(),0,lep.x()-0.001,lep.x()+0.001);
    minimizer->SetParameter(1,"leptonY"    ,lep.y(),0,lep.y()-0.001,lep.y()+0.001);
    minimizer->SetParameter(2,"leptonZ"    ,lep.z(),0,lep.z()-0.001,lep.z()+0.001);
    minimizer->SetParameter(3,"neutrinoX"  ,met.x(),500,-3000,3000);
    minimizer->SetParameter(4,"neutrinoY"  ,met.y(),500,-3000,3000);
    minimizer->SetParameter(5,"neutrinoZ"  ,0,500,-3000,3000);
    minimizer->SetParameter(6,"jetX"       ,jet.x(),0,jet.x()-0.001,jet.x()+0.001);
    minimizer->SetParameter(7,"jetY"       ,jet.y(),0,jet.y()-0.001,jet.y()+0.001);
    minimizer->SetParameter(8,"jetZ"       ,jet.z(),0,jet.z()-0.001,jet.z()+0.001);
    minimizer->SetParameter(9,"jetM"       ,jetM,0,jetM-0.001,jetM+0.0010);
    minimizer->SetParameter(10,"jetSF"      ,1,2    ,-5,5);
    minimizer->SetParameter(11,"metX"       ,met.x(),0,met.x()-0.001,met.x()+0.001);
    minimizer->SetParameter(12,"metY"       ,met.y(),0,met.y()-0.001,met.y()+0.001);
    minimizer->SetParameter(13,"wlnuM"      ,wLNUMass,0,wLNUMass-0.001,wLNUMass+0.001);



  minimizer->FixParameter(0);
  minimizer->FixParameter(1);
  minimizer->FixParameter(2);
  minimizer->FixParameter(6);
  minimizer->FixParameter(7);
  minimizer->FixParameter(8);
  minimizer->FixParameter(9);
  minimizer->FixParameter(11);
  minimizer->FixParameter(12);
  minimizer->FixParameter(13);

//  minimizer->FixParameter(10);
//  minimizer->FixParameter(3);
//  minimizer->FixParameter(4);


//  minimizer->SetParameter(3,"neutrinoX"  ,met.x(),0,-3000,3000);
//  minimizer->SetParameter(4,"neutrinoY"  ,met.y(),0,-3000,3000);
//  minimizer->FixParameter(3);
//  minimizer->FixParameter(4);
//  minimizer->SetParameter(5,"neutrinoZ"  ,nZP,0,-3000,3000);
//    minimizer->FixParameter(5);

  // Run the simplex minimizer to get close to the minimum [no good precision]
  minimizer->ExecuteCommand("SIMPLEX",0,0);

  // Run the migrad minimizer to precisely estimate the minimum
  //    minimizer->ExecuteCommand("MIGRAD",0,0);

  return  hSolverFunction(
  minimizer->GetParameter(0 ),
  minimizer->GetParameter(1 ),
  minimizer->GetParameter(2 ),
  minimizer->GetParameter(3 ),
  minimizer->GetParameter(4 ),
  minimizer->GetParameter(5 ),
  minimizer->GetParameter(6 ),
  minimizer->GetParameter(7 ),
  minimizer->GetParameter(8 ),
  minimizer->GetParameter(9 ),
  minimizer->GetParameter(10),
  minimizer->GetParameter(11),
  minimizer->GetParameter(12),
  minimizer->GetParameter(13),
  info);


} // ~ end of Topness Minimization()






class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){
    }


    void getRes() {

        double hwwParX = wjjCand->px() + selectedLepton->px() + reader_event->met.px();
        double hwwParY = wjjCand->py() + selectedLepton->py() + reader_event->met.py();

        double hwwMag = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);

        double hwwParNormX  = hwwParX/hwwMag;
        double hwwParNormY  = hwwParY/hwwMag;

        double hwwPerpNormX  = -1*hwwParNormY;
        double hwwPerpNormY  = hwwParNormX;

        plotter.getOrMake1DPre(smpName,"hwwMag",";hwwMag",300,0,3000)->Fill(hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"met",";met",500,0,1000)->Fill(reader_event->met.pt(),weight);


        auto getPerp =[&](const MomentumF* mom)->double{
            return mom->px()*hwwPerpNormX+mom->py()*hwwPerpNormY;
        };
        auto getPar =[&](const MomentumF* mom)->double{
            return mom->px()*hwwParNormX+mom->py()*hwwParNormY;
        };

        double metPerp = getPerp(&reader_event->met);
        double metPar = getPar(&reader_event->met);


        plotter.getOrMake1DPre(smpName,"metPerp",";metPerp",500,0,1000)->Fill(std::fabs(metPerp),weight);
        plotter.getOrMake1DPre(smpName,"metPar",";metPar",500,0,1000)->Fill(std::fabs(metPar),weight);

        double neutPerp = getPerp(diHiggsEvt.w1_d2);
        double neutPar = getPar(diHiggsEvt.w1_d2);

        plotter.getOrMake1DPre(smpName,"neutPerp",";neutPerp",500,0,1000)->Fill(std::fabs(neutPerp),weight);
        plotter.getOrMake1DPre(smpName,"neutPar",";neutPar",500,0,1000)->Fill(std::fabs(neutPar),weight);

        plotter.getOrMake1DPre(smpName,"neutPerpRelhwwMag",";neutPerp / hwwMag",500,-2,2)->Fill(neutPerp/hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"neutParRelhwwMag",";neutPar / hwwMag",500,-2,2)->Fill(neutPar/hwwMag,weight);

        double extraMetPerp = metPerp-neutPerp;
        double extraMetPar = metPar-neutPar;

        plotter.getOrMake1DPre(smpName,"extraMetPerp",";extraMetPerp",500,-500,500)->Fill(extraMetPerp,weight);
        plotter.getOrMake1DPre(smpName,"extraMetPar",";extraMetPar",500,-500,500)->Fill(extraMetPar,weight);

        plotter.getOrMake1DPre(smpName,"extraMetPerpRelhwwMag",";extraMetPerp / hwwMag",500,-2,2)->Fill(extraMetPerp/hwwMag,weight);
        plotter.getOrMake1DPre(smpName,"extraMetParRelhwwMag",";extraMetPar / hwwMag",500,-2,2)->Fill(extraMetPar/hwwMag,weight);


        plotter.getOrMake1DPre(smpName,"neutZ",";neutZ",500,-1000,1000)->Fill(diHiggsEvt.w1_d2->pz(),weight);
        plotter.getOrMake1DPre(smpName,"neutZmLepZ",";neutZmLepZ",500,-1000,1000)->Fill(diHiggsEvt.w1_d2->pz() - selectedLepton->pz() ,weight);


        plotter.getOrMake1DPre(smpName,"neutZmLepZRel",";neutZmLepZRel",500,-5,5)->Fill((diHiggsEvt.w1_d2->pz() - selectedLepton->pz())/std::fabs(diHiggsEvt.w1_d2->pz() + selectedLepton->pz()) ,weight);


        plotter.getOrMake1DPre(smpName,"neutZmhWWZRel",";neutZmLepZRel",500,-5,5)->Fill(
                (diHiggsEvt.w1_d2->pz() - selectedLepton->pz() - wjjCand->pz() )/std::fabs(selectedLepton->pz()+wjjCand->pz()) ,weight);

        plotter.getOrMake1DPre(smpName,"neutZmhWWZ",";neutZmLepZ",500,-1000,1000)->Fill(
                (diHiggsEvt.w1_d2->pz() - selectedLepton->pz() - wjjCand->pz() ),weight);


        double neutCosT = diHiggsEvt.w1_d2->pz()/diHiggsEvt.w1_d2->p();
        double otherCosT = ( selectedLepton->pz() + wjjCand->pz())/( selectedLepton->p4() + wjjCand->p4()).P();

        plotter.getOrMake1DPre(smpName,"neutCosT",";neutCosT",100,-1,1)->Fill(neutCosT,weight);
        plotter.getOrMake1DPre(smpName,"neutCosTmOCT",";neutCosTmOCT",200,-2,2)->Fill(otherCosT,weight);


        plotter.getOrMake1DPre(smpName,"neutTheta",";neutTheta",100,0,3.2)->Fill( diHiggsEvt.w1_d2->theta(),weight);
        plotter.getOrMake1DPre(smpName,"neutThetamOT",";neutThetamOT",600,-5,5)->Fill( diHiggsEvt.w1_d2->theta() - (selectedLepton->p4() + wjjCand->p4()).theta(),weight);
        plotter.getOrMake1DPre(smpName,"neutThetamlT",";neutThetamlT",600,-5,5)->Fill( diHiggsEvt.w1_d2->theta() - selectedLepton->theta(),weight);
        plotter.getOrMake1DPre(smpName,"neutThetamAT",";neutThetamAT",600,-5,5)->Fill( diHiggsEvt.w1_d2->theta() -  (selectedLepton->p4() + wjjCand->p4()+diHiggsEvt.w1_d2->p4()).theta(),weight);





         auto trueWqqMom = diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4();

         double wqqDPhi = PhysicsUtilities::deltaPhi(trueWqqMom,*wjjCand);
         double wqqDR   = PhysicsUtilities::deltaR(trueWqqMom,*wjjCand);
         double wqqPTDif = wjjCand->pt()-trueWqqMom.pt();
         double wqqPDif = wjjCand->p()-trueWqqMom.P();
         double wqqPTRes = wqqPTDif/trueWqqMom.pt();
         double wqqPRes = wqqPDif/trueWqqMom.P();

         if(trueWqqMom.pt()>10){
             plotter.getOrMake1DPre(smpName,"wqqDPhi",";wqqDPhi",500,-3,3)->Fill(wqqDPhi,weight);
             plotter.getOrMake1DPre(smpName,"wqqDR",";wqqDR",500,0,3)->Fill(wqqDR,weight);
             plotter.getOrMake1DPre(smpName,"wqqPTDif",";wqqPTDif",500,-2000,2000)->Fill(wqqPTDif,weight);
             plotter.getOrMake1DPre(smpName,"wqqPDif",";wqqPDif",500,-2000,2000)->Fill(wqqPDif,weight);
             plotter.getOrMake1DPre(smpName,"wqqPTRes",";wqqPTRes",500,-2,2)->Fill(wqqPTRes,weight);
             plotter.getOrMake1DPre(smpName,"wqqPTTrue",";wqqPTTrue",500,0,2000)->Fill(trueWqqMom.pt(),weight);
             plotter.getOrMake1DPre(smpName,"wqqPTReco",";wqqPTReco",500,0,2000)->Fill(wjjCand->pt(),weight);
             plotter.getOrMake1DPre(smpName,"wqqPRes",";wqqPRes",500,-2,2)->Fill(wqqPRes,weight);
         }


         bool isVirtual = trueWqqMom.mass() < (diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).mass();
         double sdMass = wjjCand->sdMom().mass();
         auto plt =[&](const TString& prefix){
           plotter.getOrMake1DPre(prefix, isVirtual ? "virtualSD" : "onshellSD",";SD",300,0,300)->Fill(sdMass,weight);
           auto trueW = selectedLepton->p4() + diHiggsEvt.w1_d2->p4();
           plotter.getOrMake1DPre(prefix+ (isVirtual ? "_virtualSD" : "_onshellSD"), "trueW",";trueW",300,0,300)->Fill(trueW.mass(),weight);
           double sf =  (diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).P()/wjjCand->p();
//           ASTypes::CylLorentzVectorF approxJ (sf*wjjCand->pt(),wjjCand->eta(),wjjCand->phi(),(diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).mass());

           ASTypes::CylLorentzVectorF approxJ (sf*wjjCand->pt(),wjjCand->eta(),wjjCand->phi(),isVirtual ? 31 : 80.0);

           auto trueH = trueW + approxJ;
           plotter.getOrMake1DPre(prefix+ (isVirtual ? "_virtualSD" : "_onshellSD"), "trueH",";trueH",300,0,300)->Fill(trueH.mass(),weight);


           plotter.getOrMake1DPre(prefix+ (isVirtual ? "_virtualSD" : "_onshellSD"), "totalM","totalM",300,0,300)->Fill((diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).mass() +(diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).mass() ,weight);




         };
         plt(smpName);
         if(wjjCand->tau2otau1() < 0.55) plt(smpName+"_lt0p55");
         if(wjjCand->tau2otau1() < 0.75) plt(smpName+"_lt0p75");









    }

    struct HiggSol {
        void setNeutrino(const MomentumF& lepton, const double nET, const float nPhi, const float nZ){
            ASTypes::CartLorentzVector neutrinoVec(nET*std::cos(nPhi),nET*std::sin(nPhi),nZ, std::sqrt(nZ*nZ+nET*nET));
            neutrino.setP4<float>(nET,neutrinoVec.eta(),nPhi,0.0);
//            neutrino.setP4(nET,PhysicsUtilities::thetaToEta(nZ/std::sqrt(nET*nET+nZ*nZ)),nPhi,0.0);
            wlnu.setP4( neutrino.p4()+lepton.p4() );
        }
        void setJet(const MomentumF& jet, const double sf,const double jetMass){
            wqqjet.setP4(jet.pt()*float(sf),jet.eta(),jet.phi(),float(jetMass));
        }
        //solution
        double chiSq = -1;
        MomentumF neutrino;
        MomentumF wlnu;
        MomentumF wqqjet;
    };


    int getNeutrinoZ(const double wMass, const MomentumF& vis, const double invx, const double invy,  double& sol1, double& sol2 ){
        const double a = wMass*wMass - vis.mass()*vis.mass() +2*vis.x()*invx +2*vis.y()*invy;
        const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
        const double B = -4*a* vis.z();
        const double C = 4*vis.E()*vis.E()*(invx*invx + invy*invy) - a*a;
        const double delta = B*B -4*A*C;
        if(delta < 0) {
            return 0;
        } else {
           sol1 = (-B + std::sqrt(delta))/(2*A);
           sol2 = (-B - std::sqrt(delta))/(2*A);
           return 2;
        }
    }


    int getMomentumJ(const double jMass, const MomentumF& wlnu, const MomentumF& wqq,  double& sol1, double& sol2 ){
        const double alpha = 4162.5;// (125*125-30*30-80*80)/2
        const double beta = (wqq.px()*wlnu.px()+wqq.py()*wlnu.py()+wqq.pz()*wlnu.pz())/wqq.p();
        const double A = beta*beta-wlnu.E()*wlnu.E();
        const double B = 2*alpha*beta;
        const double C = alpha*alpha-wlnu.E()*wlnu.E()*jMass*jMass;
        const double delta = B*B -4*A*C;
        if(delta < 0) {
            return 0;
        } else {
           sol1 = (-B + std::sqrt(delta))/(2*A);
           sol2 = (-B - std::sqrt(delta))/(2*A);
           return 2;
        }
    }





    HiggSol solve(const MomentumF& met, const MomentumF& lepton, const MomentumF& jet, const double sdMass){
//        const double lepW = sdMass < 60 ? 80 : 30;
//        const double hadW = sdMass < 60 ? 30 : 80;

        const double lepW = (diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).mass() ;
        const double hadW = (diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()).mass();

        const double hwwParX = jet.px() + lepton.px() +met.px();
        const double hwwParY = jet.py() + lepton.py() +met.py();

        const double hwwMag = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);



        const double hwwParNormX  = hwwParX/hwwMag;
        const double hwwParNormY  = hwwParY/hwwMag;

        const double hwwPerpNormX  = -1*hwwParNormY;
        const double hwwPerpNormY  = hwwParNormX;


        auto getPerp =[&](const MomentumF& mom)->double{
            return mom.px()*hwwPerpNormX+mom.py()*hwwPerpNormY;
        };
        auto getPar =[&](const MomentumF& mom)->double{
            return mom.px()*hwwParNormX+mom.py()*hwwParNormY;
        };

        const double metPerp = getPerp(met);
        const double metPar = getPar(met);


        const double minNeutrinoET = 5;
        const double maxNeutrinoET = std::max(10., hwwMag);
        const double neutrinoStepSize =  (maxNeutrinoET-minNeutrinoET)/100;



        auto getChiSq = [&](const MomentumF& neutrino, const MomentumF& wqqjet){
            const double neutPerp = getPerp(neutrino);
            const double neutPar = getPar(neutrino);
            const double extraMetPar = (metPar-neutPar)/hwwMag + 0.089;
            const double extraMetPerp = (metPerp-neutPerp);
            const double wqqPRes = (wqqjet.p()/jet.p() - 1.0);
            return extraMetPerp*extraMetPerp/(30*30) + extraMetPar*extraMetPar/(0.16*0.16)+ wqqPRes*wqqPRes/(0.14*0.14);

        };

        HiggSol bestSol;
        HiggSol curSol;

        auto testJet = [&](const double jetsol){
            curSol.setJet(jet, jetsol/jet.p(), hadW);
            curSol.chiSq = getChiSq(curSol.neutrino,curSol.wqqjet);
            if(bestSol.chiSq < 0 || curSol.chiSq  < bestSol.chiSq){
                bestSol = curSol;
            }
//            std::cout <<"("<<jetsol<<","<<curSol.wqqjet.pt()<<","<<curSol.chiSq<<","<< (curSol.wqqjet.p4()+curSol.wlnu.p4()).mass()  <<")";
        };
        auto testNeutrino = [&](const double neutrinoET, const double neutrinoPhi, const double neutrinoZ){
            curSol.setNeutrino(lepton,neutrinoET,neutrinoPhi,neutrinoZ);
            double jetsol1,jetsol2;
//            std::cout <<" "<<neutrinoZ <<"("<<curSol.wlnu.mass()<<")"<<": ";
            if(getMomentumJ(hadW,curSol.wlnu,jet,jetsol1,jetsol2)){
                if(jetsol1>0) testJet(jetsol1);
                if(jetsol2>0) testJet(jetsol2);
            }
        };


//        std::cout << "Startig >> MET: "<< met <<  "lepton: "<< lepton <<" jet: "<< jet <<" sdmass:  "<< sdMass<<std::endl;
//        std::cout << "hwwMag: "<< hwwMag << " hwwParX: "<< hwwParX <<" hwwParY: "<< hwwParY <<" hwwParNormX:  "<< hwwParNormX<<" hwwParNormY:  "<< hwwParNormY<<" hwwPerpNormX:  "<< hwwPerpNormX<<" hwwPerpNormY:  "<< hwwPerpNormY<<std::endl;
//        std::cout << "minNeutrinoET: "<< minNeutrinoET << " maxNeutrinoET: "<< maxNeutrinoET <<" neutrinoStepSize: "<< neutrinoStepSize <<std::endl;


        for(double neutrinoET = minNeutrinoET; neutrinoET <= maxNeutrinoET ; neutrinoET+=neutrinoStepSize ){
            const double minCosPhi = 1 - lepW*lepW/(2*lepton.pt()*neutrinoET );
            const double maxDeltaPhi = minCosPhi <= -1 ? TMath::Pi() : std::acos(minCosPhi);
            const double phiStepSize = maxDeltaPhi/50;

//            std::cout <<"nET: " << neutrinoET <<" minCosPhi: "<<minCosPhi<<" maxDeltaPhi: "<< maxDeltaPhi<<" phiStepSize: "<<phiStepSize<<std::endl;
            for(double deltaPhi = -maxDeltaPhi; deltaPhi < maxDeltaPhi; deltaPhi+=phiStepSize){
                const double neutrinoPhi =  TVector2::Phi_mpi_pi(lepton.phi() + deltaPhi);
//                std::cout <<"\tneutrinoPhi->" << neutrinoPhi <<":: ";
                double neutrinoZ1, neutrinoZ2;
                if(!getNeutrinoZ(lepW,lepton,neutrinoET*std::cos(neutrinoPhi),neutrinoET*std::sin(neutrinoPhi),neutrinoZ1,neutrinoZ2)){
//                    std::cout <<std::endl;
                    continue;
                }


                testNeutrino(neutrinoET,neutrinoPhi,neutrinoZ1);
                testNeutrino(neutrinoET,neutrinoPhi,neutrinoZ2);
//                std::cout << std::endl;
          }
        }

        return bestSol;
    }

    void oldGo(){


        auto wqqjet = diHiggsEvt.w2_d1->p4() + diHiggsEvt.w2_d2->p4();
        auto wlnu = diHiggsEvt.w1_d1->p4() + diHiggsEvt.w1_d2->p4();

            HiggsSolverInfo info;
            double chisq = solver.hSolverMinimization(selectedLepton->p4(),wjjCand->p4(), reader_event->met.p4(),wjjCand->sdMom().mass() < 60,&info);
            double masSol = (info.hWW + hbbCand->p4()).mass();

            auto stupidNeutrino = neutrino;
            stupidNeutrino.setP4(reader_event->met.pt(),float(0.0),reader_event->met.mass(),float(0.0));



        plotter.getOrMake1DPre(smpName,"orig_mass",";mass",500,0,8000)->Fill(hh.mass(),weight);
        plotter.getOrMake1DPre(smpName,"new_mass",";mass",500,0,8000)->Fill( ( info.wlnu + wjjCand->p4()+  hbbCand->p4()).mass(),weight);
        plotter.getOrMake1DPre(smpName,"perfect_mass",";mass",500,0,8000)->Fill( ( diHiggsEvt.w1_d2->p4()+  selectedLepton->p4() + wjjCand->p4()+  hbbCand->p4()).mass(),weight);
        plotter.getOrMake1DPre(smpName,"stupid_mass",";mass",500,0,8000)->Fill( ( stupidNeutrino.p4()+  selectedLepton->p4() + wjjCand->p4()+  hbbCand->p4()).mass(),weight);

        auto testR = [&](TString prefix, const ASTypes::CylLorentzVectorF val){
            plotter.getOrMake1DPre(prefix,"nzR",";nzR",200,-5,5)->Fill(val.pz()/diHiggsEvt.w1_d2->pz() -1,weight);
            plotter.getOrMake1DPre(prefix,"nPTR",";nPTR",200,-5,5)->Fill(val.pt()/diHiggsEvt.w1_d2->pt() -1,weight);
            plotter.getOrMake1DPre(prefix,"nDPhi",";nDPhi",200,-5,5)->Fill(PhysicsUtilities::deltaPhi(val,*diHiggsEvt.w1_d2),weight);

        };
        testR(smpName+"_oldSol",neutrino.p4());
        testR(smpName+"_newSol",ASTypes::CylLorentzVectorF(info.neutrino));


//
//        plotter.getOrMake1DPre(smpName,"massRbb",";massRbb",500,0,2000)->Fill((hWW.p4() + diHiggsEvt.hbb->p4()).mass(),weight);
//        plotter.getOrMake1DPre(smpName,"massRN",";massRN",500,0,2000)->Fill((diHiggsEvt.w1_d2->p4()+ wjjCand->p4() + selectedLepton->p4() + hbbCand->p4()).mass(),weight);
//        plotter.getOrMake1DPre(smpName,"massRJ",";massRJ",500,0,2000)->Fill((neutrino.p4()+ wqqjet + selectedLepton->p4() + hbbCand->p4()).mass(),weight);
//        plotter.getOrMake1DPre(smpName,"massRJN",";massRJN",500,0,2000)->Fill((diHiggsEvt.w1_d2->p4()+ wqqjet + selectedLepton->p4() + hbbCand->p4()).mass(),weight);
    }

    void goTest() {
        HiggsSolverInfo info;
        double chisq = solver.hSolverMinimization(selectedLepton->p4(),wjjCand->p4(), reader_event->met.p4(),wjjCand->sdMom().mass() <60 ,&info);

        auto plotgrp = [&](TString prefix, const ASTypes::CylLorentzVectorF& hh,const ASTypes::CylLorentzVectorF& hWW,const ASTypes::CylLorentzVectorF& wlnu, const ASTypes::CylLorentzVectorF& neutrino){
            plotter.getOrMake1DPre(prefix,"hh_mass" ,";hh_mass",500,0,4000)->Fill(hh.mass(),weight);
            plotter.getOrMake1DPre(prefix,"hWW_mass",";hWW_mass",500,0,2000)->Fill(hWW.mass(),weight);
            plotter.getOrMake1DPre(prefix,"hWW_pt"  ,";hWW_pt",500,0,2000)->Fill(hWW.pt(),weight);
            plotter.getOrMake1DPre(prefix,"hWW_eta" ,";hWW_eta",100,-5,5)  ->Fill(hWW.eta(),weight);
            plotter.getOrMake1DPre(prefix,"wlnu_mass",";wlnu_mass",500,0,2000)->Fill(wlnu.mass(),weight);
            plotter.getOrMake1DPre(prefix,"wlnu_pt"  ,";wlnu_pt",500,0,2000)->Fill(wlnu.pt()  ,weight);
            plotter.getOrMake1DPre(prefix,"wlnu_eta" ,";wlnu_eta",100,-5,5  )->Fill(wlnu.eta() ,weight);
            plotter.getOrMake1DPre(prefix,"neutrino_pt"  ,";neutrino_pt",500,0,2000)->Fill(neutrino.pt()  ,weight);
            plotter.getOrMake1DPre(prefix,"neutrino_eta" ,";neutrino_eta",100,-5,5  )->Fill(neutrino.eta() ,weight);


            plotter.getOrMake1DPre(prefix,"hom" ,";hom",100,0,1  )->Fill(hWW.pt()/hh.mass() ,weight);
            if(hWW.pt()/hh.mass()>0.3){
                plotter.getOrMake1DPre(prefix,"hom_hh_mass" ,";hh_mass",500,0,4000)->Fill(hh.mass(),weight);
            }

        };
        plotgrp(smpName +"_orig",hh.p4(),hWW.p4(),wlnu.p4(), neutrino.p4());
        plotgrp(smpName +"_pN",(hbbCand->p4() + selectedLepton->p4() +diHiggsEvt.w1_d2->p4() + wjjCand->p4() ),( selectedLepton->p4() +diHiggsEvt.w1_d2->p4() + wjjCand->p4()),(selectedLepton->p4() +diHiggsEvt.w1_d2->p4()), diHiggsEvt.w1_d2->p4());
        ASTypes::CylLorentzVectorF perfectZN(reader_event->met.pt(),diHiggsEvt.w1_d2->eta(),  reader_event->met.phi(),  0);
        plotgrp(smpName +"_pNZ",(hbbCand->p4() + selectedLepton->p4() +perfectZN + wjjCand->p4() ),( selectedLepton->p4() +perfectZN + wjjCand->p4()),(selectedLepton->p4() +perfectZN), perfectZN);

        auto newhh = (info.hWW + hbbCand->p4());

        plotgrp(smpName +"_new" ,ASTypes::CylLorentzVectorF(newhh),ASTypes::CylLorentzVectorF(info.hWW),ASTypes::CylLorentzVectorF(info.wlnu), ASTypes::CylLorentzVectorF(info.neutrino));

        if(chisq<10)
            plotgrp(smpName +"_newLowChisq" ,ASTypes::CylLorentzVectorF(newhh),ASTypes::CylLorentzVectorF(info.hWW),ASTypes::CylLorentzVectorF(info.wlnu), ASTypes::CylLorentzVectorF(info.neutrino));
        if(wwDM <125)
            plotgrp(smpName +"_origLowMD",hh.p4(),hWW.p4(),wlnu.p4(), neutrino.p4());

        plotter.getOrMake1DPre(smpName,"chisq",";chisq",500,0,100)->Fill(chisq,weight);
        if(wwDM <125)
            plotter.getOrMake1DPre(smpName,"lowmd_chisq",";chisq",500,0,100)->Fill(chisq,weight);
        if(chisq <10)
            plotter.getOrMake1DPre(smpName,"lowchisq_md",";md",500,0,500)->Fill(wwDM,weight);
        plotter.getOrMake1DPre(smpName,"md",";md",500,0,500)->Fill(wwDM,weight);

        plotter.getOrMake1DPre(smpName,"sf",";sf",100,-2,2)->Fill(info.SF,weight);

        auto plotJetGroup = [&](TString prefix){

            plotter.getOrMake1DPre(prefix,"jetMass",";jetMass",500,0,500)->Fill(wjjCand->mass(),weight);
            plotter.getOrMake1DPre(prefix,"jetSDMass",";jetSDMass",500,0,500)->Fill(wjjCand->sdMom().mass(),weight);

            double dr = PhysicsUtilities::deltaR(wjjCand->subJet(0),wjjCand->subJet(1));
            double pt  = wjjCand->sdMom().pt();


            plotter.getOrMake1DPre(prefix,"jetDM",";jetDM",500,0,500)->Fill(dr*wjjCand->pt()/2,weight);
            plotter.getOrMake1DPre(prefix,"jetSJDM",";jetSJDM",500,0,500)->Fill(dr*pt/2,weight);
        };

        plotJetGroup(smpName);
        if(hh.mass()>2000)plotJetGroup(smpName+"_hh2TeV");

        auto mkTL = [](const ASTypes::CylLorentzVectorF& in) -> TLorentzVector{
            return TLorentzVector(in.px(),in.py(),in.pz(),in.E());
        };
        auto hhV = mkTL(info.hWW + hbbCand->p4());
        auto hWWV = mkTL(info.hWW);
        auto hbbV = mkTL(hbbCand->p4());
        auto wlnuVect = mkTL(info.wlnu);
        auto wjjVect = mkTL(info.wqqjet);
        auto neutVect = mkTL(info.neutrino);
        auto lepVect = mkTL(selectedLepton->p4());

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
        mkBPl(hhV,hWWV,hbbV,smpName+"_hhB");
        mkBPl(hWWV,wlnuVect,wjjVect,smpName+"_hwwB");
        mkBPl(wlnuVect,neutVect,lepVect,smpName+"_wlnuB");

        if(newhh.mass()>2000){
            mkBPl(hhV,hWWV,hbbV,smpName+"_hh2TeV_hhB");
            mkBPl(hWWV,wlnuVect,wjjVect,smpName+"_hh2TeV_hwwB");
            mkBPl(wlnuVect,neutVect,lepVect,smpName+"_hh2TeV_wlnuB");
        }



    }




    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!hbbCand || !wjjCand) return false;
        if(wjjCand->tau2otau1() >0.75) return false;
        if(hbbCSVCat <4) return false;
        if(isSignal() && diHiggsEvt.type < DiHiggsEvent::MU) return false;
//        if(isSignal())getRes();
//        if(isSignal()) oldGo();

        goTest();


        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    HiggsSolver solver;

};

#endif

void tryNewHiggsSolver(std::string fileName, int treeInt,int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void tryNewHiggsSolver(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
