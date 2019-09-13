#include "Processors/Variables/interface/HiggsSolver.h"
#include "Configuration/interface/ReaderConstants.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
namespace TAna {




//--------------------------------------------------------------------------------------------------
// HSolverBasic
//--------------------------------------------------------------------------------------------------
MomentumF HSolverBasic::getInvisible(const MomentumF& met, const MomentumF& vis, const double hMass){
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
    ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,
            std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
    return MomentumF(neutrino);
}



//--------------------------------------------------------------------------------------------------
// HSolverChi
//--------------------------------------------------------------------------------------------------
HSolverChi::HSolverChi() : minimizer(new TFitter(4))
{
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
}
//--------------------------------------------------------------------------------------------------
HSolverChi::~HSolverChi(){
    delete minimizer;
}
//--------------------------------------------------------------------------------------------------
double HSolverChi::hSolverFunction(
        const double leptonX, const double leptonY, const double leptonZ,
        const double neutrinoX, const double neutrinoY,const double neutrinoZ,
        const double jetX, const double jetY, const double jetZ, const double jetM,
        const double jetSF, const double metX, const double metY, HSolverChiInfo * info
) {

    const double hwwParX = jetX+metX+leptonX;
    const double hwwParY = jetY+metY+leptonY;
    const double hwwMag = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);

    const double hwwParNormX  = hwwParX/hwwMag;
    const double hwwParNormY  = hwwParY/hwwMag;

    const double hwwPerpNormX  = -1*hwwParNormY;
    const double hwwPerpNormY  = hwwParNormX;

    auto getPerp =[&](const double momx, const double momy)->double{
        return momx*hwwPerpNormX+momy*hwwPerpNormY;};
    auto getPar =[&](const double momx, const double momy)->double{
        return momx*hwwParNormX+momy*hwwParNormY;};

    const double metPerp = getPerp(metX,metY);
    const double metPar = getPar(metX,metY);
    const double neutPerp = getPerp(neutrinoX,neutrinoY);
    const double neutPar = getPar(neutrinoX,neutrinoY);


    const double leptonE = std::sqrt(leptonX*leptonX+leptonY*leptonY+leptonZ*leptonZ);
    const ASTypes::CartLorentzVector lepton(leptonX,leptonY,leptonZ,leptonE);
    const double neutrinoE = std::sqrt(neutrinoX*neutrinoX+neutrinoY*neutrinoY+neutrinoZ*neutrinoZ);
    const ASTypes::CartLorentzVector neutrino(neutrinoX,neutrinoY,neutrinoZ,neutrinoE);
    const double jetE   = std::sqrt(jetSF*jetSF*(jetX*jetX+jetY*jetY+jetZ*jetZ)+jetM*jetM);
    const ASTypes::CartLorentzVector jet(jetSF*jetX,jetSF*jetY,jetSF*jetZ,jetE);

    const ASTypes::CartLorentzVector wlnu = lepton + neutrino;
    const ASTypes::CartLorentzVector hww = wlnu + jet;

    const double extraMetPar = (metPar-neutPar)/hwwMag;
    const double metParErr =   extraMetPar >= 0 ? posMETParErr : negMETParErr;
    const double metParChiSq =  extraMetPar*extraMetPar/(metParErr*metParErr);


    const double extraMetPerp = (metPerp-neutPerp);
    const double metPerpError = metPerpErr;
    const double metPerpChiSq = extraMetPerp*extraMetPerp/(metPerpError*metPerpError);

    const double jetSFChiSq = (jetSF-1.0)*(jetSF-1.0)/(jetErr*jetErr);

    const double meanWlnuMass = jetM > 60 ? offWlnuMeanWlnu : onWlnuMeanWlnu;
    const double wlnuMassError = jetM > 60
            ? (wlnu.mass() > offWlnuMeanWlnu ? offWlnuPosWlnuErr : offWnluNegWlnuErr)
                    : onWlnuWlnuErr;
    const double wlnuMassChiSq = std::pow( (wlnu.mass() -meanWlnuMass)/wlnuMassError,2  );

    const double hWWMassError = jetM > 60 ? offWlnuHWWErr : onWlnuHWWErr;
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

//    std::cout <<" A "<<chisq <<" "<< neutrinoX<<" "<<neutrinoY<<" "<< neutrinoZ <<" "<< jetSF
//    <<" "<< extraMetPar <<" "<< extraMetPerp<<" "<<" "<< jet.pt() <<" "<< jetM<<" "
//    <<wlnu.mass() <<" "<< hww.mass()  <<std::endl;
    return chisq;

}
//--------------------------------------------------------------------------------------------------
void HSolverChi::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par,
        int flg) {
    if(nDim){}
    if(gout){}
    if(flg){}

    result = hSolverFunction(par[0],par[1],par[2],par[3],
            par[4],par[5],par[6],par[7],
            par[8],par[9],par[10],par[11],
            par[12]
    );

} // ~end of minuit function
//--------------------------------------------------------------------------------------------------
double HSolverChi::hSolverMinimization(const ASTypes::CylLorentzVectorF& lep,
        const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met,
        bool jetIsVirtual, const HWWParameters& pars, HSolverChiInfo * info) {

    posMETParErr        = pars.posMETParErr    ;
    negMETParErr        = pars.negMETParErr    ;
    metPerpErr          = pars.metPerpErr      ;
    jetErr              = pars.jetErr          ;
    onWlnuMeanJet       = pars.onWlnuMeanJet   ;
    offWlnuMeanJet      = pars.offWlnuMeanJet  ;
    onWlnuMeanWlnu      = pars.onWlnuMeanWlnu  ;
    offWlnuMeanWlnu     = pars.offWlnuMeanWlnu ;
    offWlnuPosWlnuErr   = pars.offWlnuPosWlnuErr;
    offWnluNegWlnuErr   = pars.offWnluNegWlnuErr;
    onWlnuWlnuErr       = pars.onWlnuWlnuErr   ;
    onWlnuHWWErr        = pars.onWlnuHWWErr    ;
    offWlnuHWWErr       = pars.offWlnuHWWErr   ;

    const double jetM =jetIsVirtual ? onWlnuMeanJet : offWlnuMeanJet;

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
//--------------------------------------------------------------------------------------------------
double HSolverChi::posMETParErr     = 0.067;
double HSolverChi::negMETParErr     = 0.16;
double HSolverChi::metPerpErr       = 31;
double HSolverChi::jetErr           = 0.11;
double HSolverChi::onWlnuMeanJet    = 30;
double HSolverChi::offWlnuMeanJet   = 80;
double HSolverChi::onWlnuMeanWlnu   = 80;
double HSolverChi::offWlnuMeanWlnu  = 41;
double HSolverChi::offWlnuPosWlnuErr= 5;
double HSolverChi::offWnluNegWlnuErr= 16;
double HSolverChi::onWlnuWlnuErr    = 2.3;
double HSolverChi::onWlnuHWWErr     = 8.3;
double HSolverChi::offWlnuHWWErr    = 3;

//--------------------------------------------------------------------------------------------------
// BASEPDF
//--------------------------------------------------------------------------------------------------
double BASEPDF::extapDown(const double x, const double bound, const double probAtBound,
        const double integBound)  {
    return probAtBound*std::exp((x-bound)*probAtBound/integBound);
}
//--------------------------------------------------------------------------------------------------
double BASEPDF::extapUp(const double x, const double bound, const double probAtBound,
        const double integBound)  {
    return probAtBound*std::exp((bound-x)*probAtBound/integBound);
}
//--------------------------------------------------------------------------------------------------
double BASEPDF::interpolate(const double x, const double x1, const double x2,
        const double y1, const double y2) {
    return y1 + (x-x1)*((y2-y1)/(x2-x1));
}

//--------------------------------------------------------------------------------------------------
// OneDimPDFWInterp
//--------------------------------------------------------------------------------------------------
void OneDimPDFWInterp::setup(TFile * inFile, const std::string& hLowName,
        const std::string& hHighName, const double lowValue, const double highValue, bool verbose ){
    hLow = TObjectHelper::getObject<TH1>(inFile,hLowName,verbose);
    hHigh = TObjectHelper::getObject<TH1>(inFile,hHighName,verbose);
    hCur.reset((TH1*)hLow->Clone());
    hM.reset((TH1*)hLow->Clone());
    hB.reset((TH1*)hLow->Clone());

    hCur->SetDirectory(0);
    hM->SetDirectory(0);
    hB->SetDirectory(0);

    nB = hCur->GetNbinsX();
    iLow= lowValue;
    iHigh= highValue;
    bW = hCur->GetBinWidth(1);

    if(iHigh <= iLow)
        throw std::invalid_argument(std::string(
                "OneDimPDFWInterp::setup() -> high value must be higher than the low value "));

    const double vD = iHigh - iLow;
    for(int iB = 0; iB <= nB+1; ++iB){
        const double yL = hLow->GetBinContent(iB);
        const double yH = hHigh->GetBinContent(iB);
        const double m = (yH-yL)/vD;
        const double b = (yL*iHigh-yH*iLow)/vD;
        hB->SetBinContent(iB, b);
        hM->SetBinContent(iB, m);
    }

}
//--------------------------------------------------------------------------------------------------
void OneDimPDFWInterp::setInterp(double i){
    if(i < iLow) i = iLow;
    else if(i > iHigh) i = iHigh;
    for(int iB = 0; iB <= nB+1; ++iB){
        hCur->SetBinContent(iB,
                hM->GetBinContent(iB)*i
                +hB->GetBinContent(iB));
    }
}
//--------------------------------------------------------------------------------------------------
double OneDimPDFWInterp::getProbability(const double x) const {
    return hCur->GetBinContent(hCur->FindFixBin(x));
}

//--------------------------------------------------------------------------------------------------
// OneDimPDFWExtrap
//--------------------------------------------------------------------------------------------------
void OneDimPDFWExtrap::setup(TFile * inFile, const std::string& hName, bool verbose ){
    hCur = TObjectHelper::getObject<TH1>(inFile,hName,verbose);
    hCur->SetDirectory(0);

    nB = hCur->GetNbinsX();
    bW = hCur->GetBinWidth(1);
    minX = hCur->GetBinCenter(1);
    maxX = hCur->GetBinCenter(nB);
}
//--------------------------------------------------------------------------------------------------
double OneDimPDFWExtrap::getProbability(const double x) const{
    if(x <= minX){
        return bW*extapDown(x,minX,
                hCur->GetBinContent(1),hCur->GetBinContent(0));
    }
    if( x >= maxX){
        return bW*extapUp(x,maxX,
                hCur->GetBinContent(nB),hCur->GetBinContent(nB+1));
    }
    return bW*hCur->Interpolate(x);
}
//--------------------------------------------------------------------------------------------------
// OneDimPDFWInterpAndExtrap
//--------------------------------------------------------------------------------------------------
void OneDimPDFWInterpAndExtrap::setup(TFile * inFile, const std::string& hLowName,
        const std::string& hHighName, const double lowValue, const double highValue, bool verbose ){
    OneDimPDFWInterp::setup(inFile,hLowName,hHighName,lowValue,highValue,verbose);
    minX = hCur->GetBinCenter(1);
    maxX = hCur->GetBinCenter(nB);
}
//--------------------------------------------------------------------------------------------------
double OneDimPDFWInterpAndExtrap::getProbability(const double x) const{
    if(x <= minX){
        return bW*extapDown(x,minX,
                hCur->GetBinContent(1),hCur->GetBinContent(0));
    }
    if( x >= maxX){
        return bW*extapUp(x,maxX,
                hCur->GetBinContent(nB),hCur->GetBinContent(nB+1));
    }
    return bW*hCur->Interpolate(x);
}
//--------------------------------------------------------------------------------------------------
// TwoDimPDF
//--------------------------------------------------------------------------------------------------
void TwoDimPDF::setup(TFile * inFile, const std::string& hName, bool verbose ){
    h = TObjectHelper::getObject<TH1>(inFile,hName,verbose);

    nBX = h->GetNbinsX();
    minX = h->GetXaxis()->GetBinCenter(1);
    maxX = h->GetXaxis()->GetBinCenter(nBX);

    nBY = h->GetNbinsY();
    minY = h->GetYaxis()->GetBinCenter(1);
    maxY = h->GetYaxis()->GetBinCenter(nBY);

    bW = h->GetXaxis()->GetBinWidth(1)*h->GetYaxis()->GetBinWidth(1);
}
//--------------------------------------------------------------------------------------------------
double TwoDimPDF::getProbability(const double x, const double y) const {
    if(std::isinf(x)||std::isinf(y)) return 0;
    if(std::isnan(x)||std::isnan(y)) return 0;

    //Do the corners first
    if(x <= minX){
        if(y<= minY){
            return bW*extapDown(x,minX,
                    h->GetBinContent(1,1),h->GetBinContent(0,1))
            *  extapDown(y,minY,
                    h->GetBinContent(1,1),h->GetBinContent(1,0));

        }
        if(y>= maxY){
            return bW*extapDown(x,minX,
                    h->GetBinContent(1,nBY),h->GetBinContent(0,nBY))
            *  extapUp(y,maxY,
                    h->GetBinContent(1,nBY),h->GetBinContent(1,nBY+1));


        }
    }
    if(x >= maxX){
        if(y<= minY){
                            return bW*extapUp(x,maxX,
                    h->GetBinContent(nBX,1),h->GetBinContent(nBX+1,1))
            *  extapDown(y,minY,
                    h->GetBinContent(nBX,1),h->GetBinContent(nBX,0));
        }
        if(y>= maxY){
            return bW*extapUp(x,maxX,
                    h->GetBinContent(nBX,nBY),h->GetBinContent(nBX+1,nBY))
            *  extapUp(y,maxY,
                    h->GetBinContent(nBX,nBY),h->GetBinContent(nBX,nBY+1));
        }
    }

    //Now the partial overflows
    if(x <= minX || x >= maxX){
        int xBP, xBI;
        if(x <= minX){
            xBP = 1;
            xBI = 0;
        } else {
            xBP = nBX;
            xBI = nBX+1;
        }
        const int    yB = h->GetYaxis()->FindFixBin(y);
        const double yBC = h->GetYaxis()->GetBinCenter(yB);
        const double yP = h->GetBinContent(xBP,yB);
        const double yI = h->GetBinContent(xBI,yB);
        double yEP, yEI;
        if(y < yBC){
            const double yM1BC = h->GetYaxis()->GetBinCenter(yB-1);
            const double yM1P = h->GetBinContent(xBP,yB-1);
            const double yM1I = h->GetBinContent(xBI,yB-1);
            yEP = interpolate(y,yM1BC,yBC,yM1P,yP);
            yEI = interpolate(y,yM1BC,yBC,yM1I,yI);
        } else {
            const double yP1BC = h->GetYaxis()->GetBinCenter(yB+1);
            const double yP1P = h->GetBinContent(xBP,yB+1);
            const double yP1I = h->GetBinContent(xBI,yB+1);
            yEP = interpolate(y,yBC,yP1BC,yP,yP1P);
            yEI = interpolate(y,yBC,yP1BC,yI,yP1I);
        }
        if(x <= minX){
            return bW*extapDown(x,minX,yEP,yEI);
        }

        return bW*extapUp(x,maxX,yEP,yEI);
    }
    if(y <= minY || y >= maxY){
        int yBP, yBI;
        if(y <= minY){
            yBP = 1;
            yBI = 0;
        } else {
            yBP = nBY;
            yBI = nBY+1;
        }
        const int    xB  = h->GetXaxis()->FindFixBin(x);
        const double xBC = h->GetXaxis()->GetBinCenter(xB);
        const double xP  = h->GetBinContent(xB,yBP);
        const double xI  = h->GetBinContent(xB,yBI);
        double xEP, xEI;
        if(x < xBC){
            const double xM1BC = h->GetXaxis()->GetBinCenter(xB-1);
            const double xM1P = h->GetBinContent(xB-1,yBP);
            const double xM1I = h->GetBinContent(xB-1,yBI);
            xEP = interpolate(x,xM1BC,xBC,xM1P,xP);
            xEI = interpolate(x,xM1BC,xBC,xM1I,xI);
        } else {
            const double xP1BC = h->GetXaxis()->GetBinCenter(xB+1);
            const double xP1P = h->GetBinContent(xB+1,yBP);
            const double xP1I = h->GetBinContent(xB+1,yBI);
            xEP = interpolate(x,xBC,xP1BC,xP,xP1P);
            xEI = interpolate(x,xBC,xP1BC,xI,xP1I);
        }
        if(y <= minY){
            return bW*extapDown(y,minY,xEP,xEP);
        }
        return bW*extapUp(y,maxY,xEP,xEI);
    }

    //And everything else
    return bW*h->Interpolate(x,y);
}
//--------------------------------------------------------------------------------------------------
// HSolverLiFunction
//--------------------------------------------------------------------------------------------------
void HSolverLiFunction::setObservables(const MomentumF& inL,
        const MomentumF& inM, const MomentumF& inJ){


    lepton = inL.p4();met = inM.p4();qqJet = inJ.p4();
    origQQPT = inJ.pt();

    hwwParX = qqJet.px()+ lepton.px() + met.px();
    hwwParY = qqJet.py() + lepton.py() + met.py();
    hwwMag  = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);

    hwwParNormX = hwwParX/hwwMag;
    hwwParNormY = hwwParY/hwwMag;
    hwwPerpNormX = -1*hwwParNormY;
    hwwPerpNormY = hwwParNormX;

    metPerp = met.px()*hwwPerpNormX+met.py()*hwwPerpNormY;
    metPar  = met.px()*hwwParNormX+met.py()*hwwParNormY;
}
//--------------------------------------------------------------------------------------------------
void HSolverLiFunction::setIterationStorage(const double * p){
    wqqSF = 1.0/(p[WQQ_RES]+1.0);
    jetE   = std::sqrt(
            wqqSF*wqqSF*(qqJet.px()*qqJet.px()+qqJet.py()*qqJet.py()+qqJet.pz()*qqJet.pz())
            +jetM*jetM);

    scaledQQJet.SetPxPyPzE(wqqSF*qqJet.px(),wqqSF*qqJet.py(),wqqSF*qqJet.pz(),jetE);

    neutPerp = metPerp- p[EMET_PERP];
    neutPar  = metPar- p[EMET_PAR]*hwwMag;
    neutE = std::sqrt(neutPerp*neutPerp+neutPar*neutPar+p[NEUT_Z]*p[NEUT_Z]);
    neutX = neutPerp*hwwPerpNormX + neutPar*hwwParNormX;
    neutY = neutPerp*hwwPerpNormY + neutPar*hwwParNormY;
    neutE = std::sqrt(neutPerp*neutPerp+neutPar*neutPar+p[NEUT_Z]*p[NEUT_Z]);
    neutrino.SetPxPyPzE(neutX,neutY,p[NEUT_Z],neutE);


    wlnu = lepton + neutrino;
    hww = wlnu + scaledQQJet;


}
//--------------------------------------------------------------------------------------------------
double HSolverLiFunction::operator()(const double* p){

    setIterationStorage(p);

    LL = 0;
    LL += std::log(pdfs[HSolverLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[HSolverLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[HSolverLi::WQQ_RES]->getProbability(p[WQQ_RES]));
//    LL += std::log(pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hww.mass(),wlnu.mass()));
    LL += std::log(pdfs[HSolverLi::WLNU_MASS]->getProbability(wlnu.mass()));
    LL += std::log(pdfs[HSolverLi::HWW_MASS]->getProbability(hww.mass()));
    return -2.0*LL;
}
//--------------------------------------------------------------------------------------------------
// HSolverLiAltFunction
//--------------------------------------------------------------------------------------------------
double HSolverLiAltFunction::operator()(const double* p){
    LL = 0;
    LL += std::log(pdfs[HSolverLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[HSolverLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[HSolverLi::WQQ_RES]->getProbability(p[WQQ_RES]));
    LL += std::log(pdfs[HSolverLi::WLNU_MASS]->getProbability(p[WLNU_MASS]));
    LL += std::log(pdfs[HSolverLi::HWW_MASS]->getProbability(p[HWW_MASS]));
    return -2.0*LL;
}
//--------------------------------------------------------------------------------------------------
// HSolverLi
//--------------------------------------------------------------------------------------------------
HSolverLi::HSolverLi(const std::string& dataDir) : dataDir(dataDir),
        osqq_functor(&osqq_function,&HSolverLiFunction::operator () ,HSolverLiFunction::NUMPARAMS),
        vqq_functor(&vqq_function,&HSolverLiFunction::operator (),HSolverLiFunction::NUMPARAMS),
        osqq_alt_functor(&osqq_alt_function,&HSolverLiAltFunction::operator (),HSolverLiFunction::NUMPARAMS),
        vqq_alt_functor(&vqq_alt_function,&HSolverLiAltFunction::operator (),HSolverLiFunction::NUMPARAMS),
        minimizer( ROOT::Minuit::kSimplex,HSolverLiFunction::NUMPARAMS)
{

    if(int(HSolverLiFunction::NUMPARAMS) != int(HSolverLiAltFunction::NUMPARAMS) ||
            int(HSolverLiFunction::NUMPARAMS) != int(HSolverBkgLiFunction::NUMPARAMS))
        throw std::invalid_argument(
                "HSolverLi::HSolverLi -> The number of all function parameters must match!");

    auto setupMinimizer=[&](Minimizer& min, ROOT::Math::Functor& f) {
        min.SetFunction(f);
        resetParameters(min);

    };
    setupMinimizer(minimizer,vqq_functor);


    auto initializePDFs = [&] (std::vector<std::shared_ptr<BASEPDF>>* pdfs){
        pdfs->resize(NPDFS);
        (*pdfs)[EMET_PERP]    .reset(new OneDimPDFWInterpAndExtrap());
        (*pdfs)[EMET_PAR]     .reset(new OneDimPDFWInterpAndExtrap());
        (*pdfs)[WQQ_RES]      .reset(new OneDimPDFWInterpAndExtrap());
        (*pdfs)[WQQ_SDMASS]   .reset(new OneDimPDFWInterp());
//        (*pdfs)[HWW_WLNU_MASS].reset(new TwoDimPDF());
        (*pdfs)[WLNU_MASS]    .reset(new OneDimPDFWInterpAndExtrap());
        (*pdfs)[HWW_MASS]     .reset(new OneDimPDFWInterpAndExtrap());
    };

    auto setupPDFS = [&](HSolverFunction * func, std::vector<std::shared_ptr<BASEPDF>>* pdfs) {
        func->pdfs.resize(NPDFS);
        (*func).pdfs[EMET_PERP]     = (*pdfs)[EMET_PERP]    ;
        (*func).pdfs[EMET_PAR]      = (*pdfs)[EMET_PAR]     ;
        (*func).pdfs[WQQ_RES]       = (*pdfs)[WQQ_RES]      ;
        (*func).pdfs[WQQ_SDMASS]    = (*pdfs)[WQQ_SDMASS]   ;
//        (*func).pdfs[HWW_WLNU_MASS] = (*pdfs)[HWW_WLNU_MASS];
        (*func).pdfs[WLNU_MASS]     = (*pdfs)[WLNU_MASS];
        (*func).pdfs[HWW_MASS]      = (*pdfs)[HWW_MASS];

    };


    initializePDFs(&osqq_pdfs);
    initializePDFs(&vqq_pdfs);

    setupPDFS(&osqq_function,&osqq_pdfs);
    setupPDFS(&vqq_function,&vqq_pdfs);
    osqq_function.jetM = 80;
    vqq_function.jetM = 31;
    setupPDFS(&osqq_alt_function,&osqq_pdfs);
    setupPDFS(&vqq_alt_function,&vqq_pdfs);
}
//--------------------------------------------------------------------------------------------------
void HSolverLi::resetParameters(Minimizer& min,const double neutZ) {

    min.SetVariable(HSolverLiFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(HSolverLiFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetVariable(HSolverLiFunction::NEUT_Z  ,"NEUT_Z",neutZ,500);
    min.SetLimitedVariable(HSolverLiFunction::WQQ_RES ,"QQJET_RES",0,0.5,-0.9,5);
    min.SetFixedVariable(HSolverLiFunction::DUMMY,"DUMMY",1);
}
//--------------------------------------------------------------------------------------------------
void HSolverLi::resetAltParameters(Minimizer& min) {

    min.SetVariable(HSolverLiAltFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(HSolverLiAltFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetLimitedVariable(HSolverLiAltFunction::WQQ_RES  ,"WQQ_RES",0,0.5,-0.9,5);
    min.SetVariable(HSolverLiAltFunction::WLNU_MASS  ,"WLNU_MASS",55,500);
    min.SetVariable(HSolverLiAltFunction::HWW_MASS  ,"HWW_MASS",125,500);
}
//--------------------------------------------------------------------------------------------------
void HSolverLi::setParamters(const HWWParameters& hwwParam, bool verbose) {
    TFile * inFile = TObjectHelper::getFile(dataDir+hwwParam.liFileName,"read",verbose);
    ptCorB = hwwParam.ptCorB;
    ptCorM = hwwParam.ptCorM;

    hwwPT = TObjectHelper::getObject<TH1>(inFile,"signal_avgHWWMag",verbose);
    const double lowM = hwwPT->GetBinContent(1);
    const double highM = hwwPT->GetBinContent(2);

    auto setup1D = [&] (int entry, const std::string& vName){
        ((OneDimPDFWInterp*)(vqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_low_vqq_")+vName,
                std::string("signal_high_vqq_")+vName,lowM,highM,verbose);
        ((OneDimPDFWInterp*)(osqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_low_osqq_")+vName,
                std::string("signal_high_osqq_")+vName,lowM,highM,verbose);
    };
//    auto setup2D = [&] (int entry, const std::string& vName){
//        ((TwoDimPDF*)(vqq_pdfs[entry].get()))-> setup(inFile,
//                std::string("signal_vqq_")+vName,verbose);
//        ((TwoDimPDF*)(osqq_pdfs[entry].get()))-> setup(inFile,
//                std::string("signal_osqq_")+vName,verbose);
//    };

    setup1D(EMET_PERP,"extraMetPerp");
    setup1D(EMET_PAR,"extraMetParRelhwwMag");
    setup1D(WQQ_RES,"wqqPTRes");
    setup1D(WQQ_SDMASS,"qqSDMassCoarse");
    setup1D(WLNU_MASS,"Wlnu");
    setup1D(HWW_MASS,"hWW");
//    setup2D(HWW_WLNU_MASS,"hWW_v_Wlnu");

    delete inFile;
};
//--------------------------------------------------------------------------------------------------
double HSolverLi::getCorrHWWPT(const double recoPT) const{
    double boundPT = recoPT;
    if(recoPT < 300) boundPT = 300;
    if(recoPT > 2000) boundPT = 2000;
    double corPT = (-1*ptCorB + std::sqrt(ptCorB*ptCorB+4.*ptCorM*boundPT))/(2.*ptCorM);
    return recoPT * corPT/boundPT; //apply as a correction due to the limted range;
}
//--------------------------------------------------------------------------------------------------
double HSolverLi::fitNormalization(double qqSDMass,HSolverLiInfo * out){


    double minLikeli = -1;
    auto doNomFit = [&](Minimizer& min, HSolverLiAltFunction& f) {
        resetAltParameters(min);

        const int minOut = min.Minimize();
        const double noSDLikli = min.MinValue();
        const double likeli =
                noSDLikli - 2.0*std::log(f.pdfs[WQQ_SDMASS]->getProbability(qqSDMass));
        const double *xs = min.X();

        auto fill = [&] (HSolverLiInfo& o){
            o.minOut=minOut;
            o.noSDLikli = noSDLikli;
            o.rawLikeli = likeli;
            for(unsigned int iP = 0; iP <HSolverLiAltFunction::NUMPARAMS; ++iP )
                o.params[iP] = xs[iP];
        };
        if(minLikeli < 0 || likeli < minLikeli){
            minLikeli = likeli;
            if(out) fill(*out);
        }
    };

    minimizer.SetFunction(osqq_alt_functor);
    doNomFit(minimizer,osqq_alt_function);
    minimizer.SetFunction(vqq_alt_functor);
    doNomFit(minimizer,vqq_alt_function);
    return minLikeli;
}
//--------------------------------------------------------------------------------------------------
double HSolverLi::fit(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
        double qqSDMass, HSolverLiInfo& out, HSolverLiInfo * osqq_sol, HSolverLiInfo * vqq_sol){

    double minLikeli = -1;
    auto doFit = [&](
            Minimizer& min, HSolverLiFunction& f, HSolverLiInfo* extraInfo){
        f.setObservables(lepton,met,qqJet);
        resetParameters(min,0);
        const int minOut = min.Minimize();
        const double noSDLikli = min.MinValue();
        const double likeli =
                noSDLikli - 2.0*std::log(f.pdfs[WQQ_SDMASS]->getProbability(qqSDMass));
        const double *xs = min.X();
        f.setIterationStorage(xs);

        auto fill = [&] (HSolverLiInfo& o){
            o.minOut=minOut;
            o.rawLikeli = likeli;
            o.noSDLikli = noSDLikli;
            for(unsigned int iP = 0; iP <HSolverLiAltFunction::NUMPARAMS; ++iP )
                o.params[iP] = xs[iP];
            o.neutrino = f.neutrino;
            o.wqqjet = f.scaledQQJet;
            o.wlnu = f.wlnu;
            o.hWW = f.hww;
        };

        if(minLikeli < 0 || likeli < minLikeli){
            minLikeli = likeli;
            fill(out);
        }
        if(extraInfo)fill(*extraInfo);
    };


    minimizer.SetFunction(osqq_functor);
    doFit(minimizer,osqq_function,osqq_sol);
    minimizer.SetFunction(vqq_functor);
    doFit(minimizer,vqq_function,vqq_sol);

    return minLikeli;
}

//--------------------------------------------------------------------------------------------------
void HSolverLi::interpolatePDFs(const MomentumF& lepton, const MomentumF& met,
        const MomentumF& qqJet){

    const double oHWWMag = (lepton.p4()+met.p4()+qqJet.p4()).pt();
    const double cHWWMag = getCorrHWWPT(oHWWMag);

    auto interp = [&](PDFList p) {
      ((OneDimPDFWInterp*)(osqq_pdfs[p].get()))->setInterp(cHWWMag);
      ((OneDimPDFWInterp*)(vqq_pdfs[p].get()))->setInterp(cHWWMag);
    };
    interp(EMET_PERP);
    interp(EMET_PAR);
    interp(WQQ_RES);
    interp(WQQ_SDMASS);
    interp(WLNU_MASS);
    interp(HWW_MASS);

}

//--------------------------------------------------------------------------------------------------
double HSolverLi::minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
        double qqSDMass, HSolverLiInfo& out, HSolverLiInfo * osqq_sol, HSolverLiInfo * vqq_sol,
        HSolverLiInfo * nom_sol){

    interpolatePDFs(lepton,met,qqJet);
    const double rawLikeli = fit(lepton,met,qqJet,qqSDMass,out, osqq_sol,vqq_sol);
    const double normLikeli = fitNormalization(qqSDMass,nom_sol);

    out.likeli = rawLikeli / normLikeli;
    return out.likeli;
}
//--------------------------------------------------------------------------------------------------
// HSolverBkgLiFunction
//--------------------------------------------------------------------------------------------------
void HSolverBkgLiFunction::setObservables(const MomentumF& inL,
        const MomentumF& inM, const MomentumF& inJ){

    lepton = inL.p4();met = inM.p4();qqJet = inJ.p4();


    hwwParX = qqJet.px()+ lepton.px() + met.px();
    hwwParY = qqJet.py() + lepton.py() + met.py();
    hwwMag  = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);

    hwwParNormX = hwwParX/hwwMag;
    hwwParNormY = hwwParY/hwwMag;
    hwwPerpNormX = -1*hwwParNormY;
    hwwPerpNormY = hwwParNormX;

    metPerp = met.px()*hwwPerpNormX+met.py()*hwwPerpNormY;
    metPar  = met.px()*hwwParNormX+met.py()*hwwParNormY;

    double jetE   = std::sqrt(
            (qqJet.px()*qqJet.px()+qqJet.py()*qqJet.py()+qqJet.pz()*qqJet.pz()));

    qqJet.SetPxPyPzE(qqJet.px(),qqJet.py(),qqJet.pz(),jetE);
}
//--------------------------------------------------------------------------------------------------
void HSolverBkgLiFunction::setIterationStorage(const double * p){

    neutPerp = metPerp- p[EMET_PERP];
    neutPar  = metPar- p[EMET_PAR]*hwwMag;
    neutE = std::sqrt(neutPerp*neutPerp+neutPar*neutPar+p[NEUT_Z]*p[NEUT_Z]);
    neutX = neutPerp*hwwPerpNormX + neutPar*hwwParNormX;
    neutY = neutPerp*hwwPerpNormY + neutPar*hwwParNormY;
    neutE = std::sqrt(neutPerp*neutPerp+neutPar*neutPar+p[NEUT_Z]*p[NEUT_Z]);
    neutrino.SetPxPyPzE(neutX,neutY,p[NEUT_Z],neutE);


    wlnu = lepton + neutrino;
    hww = wlnu + qqJet;
}
//--------------------------------------------------------------------------------------------------
double HSolverBkgLiFunction::operator()(const double* p){
    setIterationStorage(p);
    LL = 0;
    LL += std::log(pdfs[HSolverBkgLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[HSolverBkgLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[HSolverBkgLi::WLNU_MASS]->getProbability(wlnu.mass()));
    LL += std::log(pdfs[HSolverBkgLi::HWW_MASS]->getProbability(hww.mass()));
    return -2.0*LL;
}
//--------------------------------------------------------------------------------------------------
// HSolverBkgLiAltFunction
//--------------------------------------------------------------------------------------------------
double HSolverBkgLiAltFunction::operator()(const double* p){
    LL = 0;
    LL += std::log(pdfs[HSolverBkgLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[HSolverBkgLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[HSolverBkgLi::WLNU_MASS]->getProbability(p[WLNU_MASS]));
    LL += std::log(pdfs[HSolverBkgLi::HWW_MASS]->getProbability(p[HWW_MASS]));
    return -2.0*LL;
}

//--------------------------------------------------------------------------------------------------
// HSolverBkgLi
//--------------------------------------------------------------------------------------------------
HSolverBkgLi::HSolverBkgLi(const std::string& dataDir) : dataDir(dataDir),
        functor(&function,&HSolverBkgLiFunction::operator (),HSolverBkgLiFunction::NUMPARAMS),
        altFunctor(&altFunction,&HSolverBkgLiAltFunction::operator (),
                HSolverBkgLiAltFunction::NUMPARAMS),
        minimizer( ROOT::Minuit::kSimplex,HSolverBkgLiFunction::NUMPARAMS )
{
    auto setupMinimizer=[&](Minimizer& min, ROOT::Math::Functor& f) {
        min.SetFunction(f);
        resetParameters(min);
    };
    setupMinimizer(minimizer,functor);

    pdfs.resize(NPDFS);
    pdfs[EMET_PERP]    .reset(new OneDimPDFWInterpAndExtrap());
    pdfs[EMET_PAR]     .reset(new OneDimPDFWInterpAndExtrap());
    pdfs[WLNU_MASS].reset(new OneDimPDFWInterpAndExtrap());
    pdfs[HWW_MASS].reset(new OneDimPDFWInterpAndExtrap());

    auto setPDFS = [&](HSolverFunction* f){
        f->pdfs.resize(NPDFS);
        f->pdfs[EMET_PERP]     = pdfs[EMET_PERP]    ;
        f->pdfs[EMET_PAR]      = pdfs[EMET_PAR]     ;
        f->pdfs[WLNU_MASS]     = pdfs[WLNU_MASS]      ;
        f->pdfs[HWW_MASS]      = pdfs[HWW_MASS]   ;
    };
    setPDFS(&function);
    setPDFS(&altFunction);


}
//--------------------------------------------------------------------------------------------------
void HSolverBkgLi::setParamters(const HWWParameters& hwwParam, bool verbose){
    TFile * inFile = TObjectHelper::getFile(dataDir+hwwParam.bkgLiFileName,"read",verbose);
    hwwPT = TObjectHelper::getObject<TH1>(inFile,"ttbarPW_avgHWWMag",verbose);
    const double lowM = hwwPT->GetBinContent(1);
    const double highM = hwwPT->GetBinContent(2);
    auto setup1D = [&] (int entry, const std::string& vName){
        ((OneDimPDFWInterp*)(pdfs[entry].get()))-> setup(inFile,
                std::string("ttbarPW_low_")+vName,
                std::string("ttbarPW_high_")+vName,lowM,highM,verbose);
    };
    setup1D(EMET_PERP,"extraMetPerp");
    setup1D(EMET_PAR,"extraMetParRelhwwMag");
    setup1D(WLNU_MASS,"Wlnu");
    setup1D(HWW_MASS,"hWW");

    delete inFile;
}
//--------------------------------------------------------------------------------------------------
void HSolverBkgLi::resetParameters(Minimizer& min) {
    min.SetVariable(HSolverBkgLiFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(HSolverBkgLiFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetVariable(HSolverBkgLiFunction::NEUT_Z  ,"NEUT_Z",0,500);
    min.SetFixedVariable(HSolverBkgLiFunction::DUMMY ,"DUMMY",1);
    min.SetFixedVariable(HSolverBkgLiFunction::DUMMY2 ,"DUMMY2",1);
}
//--------------------------------------------------------------------------------------------------
void HSolverBkgLi::resetAltParameters(Minimizer& min) {
    min.SetVariable(HSolverBkgLiAltFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(HSolverBkgLiAltFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetVariable(HSolverBkgLiAltFunction::WLNU_MASS  ,"WLNU_MASS",80,500);
    min.SetVariable(HSolverBkgLiAltFunction::HWW_MASS  ,"HWW_MASS",125,500);
    min.SetFixedVariable(HSolverBkgLiAltFunction::DUMMY ,"DUMMY",1);

}
//--------------------------------------------------------------------------------------------------
double HSolverBkgLi::minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
        HSolverLiInfo& out, HSolverLiInfo* altOut){

    const double oHWWMag = (lepton.p4()+met.p4()+qqJet.p4()).pt();
    auto interp = [&](PDFList p) {
      ((OneDimPDFWInterp*)(pdfs[p].get()))->setInterp(oHWWMag);

    };
    interp(EMET_PERP);
    interp(EMET_PAR);
    interp(WLNU_MASS);
    interp(HWW_MASS);

    auto doFit = [&](
            Minimizer& min, HSolverBkgLiFunction& f,HSolverLiInfo& out  ){

        f.setObservables(lepton,met,qqJet);
        resetParameters(min);

        out.minOut=min.Minimize();

        out.likeli = min.MinValue();
        const double *xs = min.X();
        f.setIterationStorage(xs);
        for(unsigned int iP = 0; iP <HSolverLiAltFunction::NUMPARAMS; ++iP )
            out.params[iP] = xs[iP];
        out.neutrino = f.neutrino;
        out.wqqjet = f.qqJet;
        out.wlnu = f.wlnu;
        out.hWW = f.hww;

    };

    auto doNomFit = [&](
            Minimizer& min,HSolverLiInfo& out  ){
        resetAltParameters(min);
        out.minOut=min.Minimize();
        out.likeli = min.MinValue();

        const double *xs = min.X();
        for(unsigned int iP = 0; iP <HSolverLiAltFunction::NUMPARAMS; ++iP )
            out.params[iP] = xs[iP];
    };

    minimizer.SetFunction(functor);
    doFit(minimizer,function,out);


    if(altOut){
        minimizer.SetFunction(altFunctor);
        doNomFit(minimizer,*altOut);
    }
    return out.likeli;

}
}




