#include "Processors/Variables/interface/HiggsSolver.h"

namespace TAna {


HiggsSolver::HiggsSolver() : minimizer(new TFitter(4))
  {
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
  }

HiggsSolver::~HiggsSolver(){
  delete minimizer;
}

double HiggsSolver::hSolverFunction(
                   const double leptonX, const double leptonY, const double leptonZ, const double neutrinoX, const double neutrinoY,const double neutrinoZ,
                   const double jetX,    const double jetY,    const double jetZ,    const double jetM, const double jetSF,
                   const double metX,    const double metY, HiggsSolverInfo * info
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

    const double extraMetPar = (metPar-neutPar)/hwwMag;
    const double metParErr =   extraMetPar >= 0 ? 0.064 : 0.15;
    const double metParChiSq =  extraMetPar*extraMetPar/(metParErr*metParErr);


    const double extraMetPerp = (metPerp-neutPerp);
    const double metPerpError = 29;
    const double metPerpChiSq = extraMetPerp*extraMetPerp/(metPerpError*metPerpError);

    const double jetSFChiSq = (jetSF-1.0)*(jetSF-1.0)/(0.10*0.10);

    const double meanWlnuMass = jetM > 60 ? 43 : 80;
    const double wlnuMassError = jetM > 60 ? (wlnu.mass() > 43 ? 3 : 19) : 3;
    const double wlnuMassChiSq = std::pow( (wlnu.mass() -meanWlnuMass)/wlnuMassError,2  );

    const double hWWMassError = jetM > 60 ? 3 : 9;
    const double hWWMassChiSq = std::pow( (hww.mass() -125)/hWWMassError,2  );


    const double chisq= hWWMassChiSq +wlnuMassChiSq +metPerpChiSq +metParChiSq+jetSFChiSq ;


  if(info){
      info->chiSq = chisq;
      info->SF = jetSF;
      info->neutrino = neutrino;
      info->wlnu = wlnu;
      info->hWW = hww;
      info->wqqjet = jet;
  }

  return chisq;

}

void HiggsSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par, int flg) {
    if(nDim){}
        if(gout){}
            if(flg){}



  result = hSolverFunction(par[0],par[1],par[2],par[3],
               par[4],par[5],par[6],par[7],
               par[8],par[9],par[10],par[11],
               par[12]
               );

} // ~end of minuit function

double HiggsSolver::hSolverMinimization(const ASTypes::CylLorentzVectorF& lep, const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met, bool jetIsVirtual, HiggsSolverInfo * info) {
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



  minimizer->FixParameter(0);
  minimizer->FixParameter(1);
  minimizer->FixParameter(2);
  minimizer->FixParameter(6);
  minimizer->FixParameter(7);
  minimizer->FixParameter(8);
  minimizer->FixParameter(9);
  minimizer->FixParameter(11);
  minimizer->FixParameter(12);

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
  info);


}



MomentumF HiggsSolver::getInvisible(const MomentumF& met, const MomentumF& vis, const double hMass){
    const double a = hMass*hMass - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
    const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
    const double B = -4*a* vis.z();
    const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
    const double delta = B*B -4*A*C;

    double metZ = 0;
    if(delta < 0) {
        metZ= -B/(2*A);
    } else {
        const double pos = (-B + std::sqrt(delta))/(2*A);
        const double neg = (-B - std::sqrt(delta))/(2*A);
        if(std::fabs(pos) < std::fabs(neg)) metZ = pos;
        else metZ = neg;
    }
    ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
    return MomentumF(neutrino);
}

}
