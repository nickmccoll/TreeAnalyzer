#include "Processors/Variables/interface/HiggsSolver.h"
#include "Configuration/interface/ReaderConstants.h"
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
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
        const double leptonX, const double leptonY, const double leptonZ,
        const double neutrinoX, const double neutrinoY,const double neutrinoZ,
        const double jetX, const double jetY, const double jetZ, const double jetM,
        const double jetSF, const double metX, const double metY, HiggsSolverInfo * info
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

//    std::cout <<" A "<<chisq <<" "<< neutrinoX<<" "<<neutrinoY<<" "<< neutrinoZ <<" "<< jetSF  <<" "<< extraMetPar <<" "<< extraMetPerp<<" "<<" "<< jet.pt() <<" "<< jetM<<" " <<wlnu.mass() <<" "<< hww.mass()  <<std::endl;
    return chisq;

}

void HiggsSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par,
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

double HiggsSolver::hSolverMinimization(const ASTypes::CylLorentzVectorF& lep,
        const ASTypes::CylLorentzVectorF& jet, const ASTypes::CylLorentzVectorF& met,
        bool jetIsVirtual, const HWWParameters& pars, HiggsSolverInfo * info) {


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
    ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,
            std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
    return MomentumF(neutrino);
}


double HiggsSolver::posMETParErr     = 0.067;
double HiggsSolver::negMETParErr     = 0.16;
double HiggsSolver::metPerpErr       = 31;
double HiggsSolver::jetErr           = 0.11;
double HiggsSolver::onWlnuMeanJet    = 30;
double HiggsSolver::offWlnuMeanJet   = 80;
double HiggsSolver::onWlnuMeanWlnu   = 80;
double HiggsSolver::offWlnuMeanWlnu  = 41;
double HiggsSolver::offWlnuPosWlnuErr= 5;
double HiggsSolver::offWnluNegWlnuErr= 16;
double HiggsSolver::onWlnuWlnuErr    = 2.3;
double HiggsSolver::onWlnuHWWErr     = 8.3;
double HiggsSolver::offWlnuHWWErr    = 3;



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

void OneDimPDFWInterp::setInterp(double i){
    if(i < iLow) i = iLow;
    else if(i > iHigh) i = iHigh;
    for(int iB = 0; iB <= nB+1; ++iB){
        hCur->SetBinContent(iB,
                hM->GetBinContent(iB)*i
                +hB->GetBinContent(iB));
    }
}

double OneDimPDFWInterp::getProbability(const double x) const {
    return hCur->GetBinContent(hCur->FindFixBin(x));
}

double BASEPDF::extapDown(const double x, const double bound, const double probAtBound,
        const double integBound)  {
    return probAtBound*std::exp((x-bound)*probAtBound/integBound);
}

double BASEPDF::extapUp(const double x, const double bound, const double probAtBound,
        const double integBound)  {
    return probAtBound*std::exp((bound-x)*probAtBound/integBound);
}
double BASEPDF::interpolate(const double x, const double x1, const double x2,
        const double y1, const double y2) {
    return y1 + (x-x1)*((y2-y1)/(x2-x1));
}

void OneDimPDFWInterpAndExtrap::setup(TFile * inFile, const std::string& hLowName,
        const std::string& hHighName, const double lowValue, const double highValue, bool verbose ){
    OneDimPDFWInterp::setup(inFile,hLowName,hHighName,lowValue,highValue,verbose);
    minX = hCur->GetBinCenter(1);
    maxX = hCur->GetBinCenter(nB);
}


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
void OneDimPDFWExtrap::setup(TFile * inFile, const std::string& hName, bool verbose ){
    hCur = TObjectHelper::getObject<TH1>(inFile,hName,verbose);
    hCur->SetDirectory(0);

    nB = hCur->GetNbinsX();
    bW = hCur->GetBinWidth(1);
    minX = hCur->GetBinCenter(1);
    maxX = hCur->GetBinCenter(nB);
}


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

double TwoDimPDF::getProbability(const double x, const double y) const {
    if(std::isinf(x)||std::isinf(y)) return 0;
    if(std::isnan(x)||std::isnan(y)) return 0;

//    std::cout <<"!!! -> (" << nBX<<","<<minX<<","<<maxX<<") ("<< nBY<<","<<minY<<","<<maxY<<") ("
//            <<x <<","<<y<<") ";



    //Do the corners first
    if(x <= minX){
        if(y<= minY){

//            std::cout <<" 1: " << extapDown(x,minX,
//                    h->GetBinContent(1,1),h->GetBinContent(0,1))<<":"<<
//                            extapDown(y,minY,h->GetBinContent(1,1),h->GetBinContent(1,0))<<std::endl;
            return bW*extapDown(x,minX,
                    h->GetBinContent(1,1),h->GetBinContent(0,1))
            *  extapDown(y,minY,
                    h->GetBinContent(1,1),h->GetBinContent(1,0));

        }
        if(y>= maxY){

//            std::cout <<" 2: " << extapDown(x,minX,
//                    h->GetBinContent(1,nBY),h->GetBinContent(0,nBY))
//            <<":"<<  extapUp(y,maxY,
//                    h->GetBinContent(1,nBY),h->GetBinContent(1,nBY+1))<<std::endl;

            return bW*extapDown(x,minX,
                    h->GetBinContent(1,nBY),h->GetBinContent(0,nBY))
            *  extapUp(y,maxY,
                    h->GetBinContent(1,nBY),h->GetBinContent(1,nBY+1));


        }
    }
    if(x >= maxX){
        if(y<= minY){
//            std::cout <<" 3: " << extapUp(x,maxX,
//                    h->GetBinContent(nBX,1),h->GetBinContent(nBX+1,1))
//            <<":"<<  extapDown(y,minY,
//                    h->GetBinContent(nBX,1),h->GetBinContent(nBX,0))<<std::endl        ;
                            return bW*extapUp(x,maxX,
                    h->GetBinContent(nBX,1),h->GetBinContent(nBX+1,1))
            *  extapDown(y,minY,
                    h->GetBinContent(nBX,1),h->GetBinContent(nBX,0));
        }
        if(y>= maxY){
//            std::cout <<" 4: " << extapUp(x,maxX,
//                    h->GetBinContent(nBX,nBY),h->GetBinContent(nBX+1,nBY))
//            <<":"<<  extapUp(y,maxY,
//                    h->GetBinContent(nBX,nBY),h->GetBinContent(nBX,nBY+1))<<std::endl;
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
//            std::cout <<" 5: " << yEP<<" "<< yEI<<" "<< extapDown(x,minX,yEP,yEI)<<std::endl;
            return bW*extapDown(x,minX,yEP,yEI);
        }

//        std::cout <<" 6: " << yEP<<" "<< yEI<<" "<< extapUp(x,maxX,yEP,yEI)<<std::endl;
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
//            std::cout <<" 7: " << xEP<<" "<< xEP<<" "<< extapDown(y,minY,xEP,xEP)<<std::endl;
            return bW*extapDown(y,minY,xEP,xEP);
        }
//        std::cout <<" 8: " << xEP<<" "<< xEP<<" "<< extapUp(y,maxY,xEP,xEP)<<std::endl;
        return bW*extapUp(y,maxY,xEP,xEI);
    }

    //And everything else
//    std::cout <<" 9: " << h->Interpolate(x,y)<<std::endl;
    return bW*h->Interpolate(x,y);
}
void HiggsLiFunction::setObservables(const MomentumF& inL,
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

void HiggsLiFunction::setIterationStorage(const double * p){

// testting with the parameters being the neutrino momentum, not extra met
//        neutX = p[EMET_PERP];
//        neutY = p[EMET_PAR];
//
//    auto getPerp =[&](const double momx, const double momy)->double{
//        return momx*hwwPerpNormX+momy*hwwPerpNormY;};
//    auto getPar =[&](const double momx, const double momy)->double{
//        return momx*hwwParNormX+momy*hwwParNormY;};
//
//    neutPerp = getPerp(neutX,neutY);
//    neutPar  = getPar(neutX,neutY);
//    extraMetPerp = metPerp - neutPerp;
//    extraMetPar = metPar - neutPar;



//    wqqSF = p[WQQ_RES];
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



//testing the old chi2 method in this framework
double HiggsLiFunction::test(const double* p){
//    std::cout << neutE <<" -> ";
    setIterationStorage(p);
//    std::cout << jetM <<std::endl;

    const double metParErr =   extraMetPar >= 0 ? HiggsSolver::posMETParErr :
            HiggsSolver::negMETParErr;
    const double metParChiSq =  (extraMetPar/hwwMag)*(extraMetPar/hwwMag)/(metParErr*metParErr);


    const double metPerpError = HiggsSolver::metPerpErr;
    const double metPerpChiSq = extraMetPerp*extraMetPerp/(metPerpError*metPerpError);

    const double jetSFChiSq = (wqqSF-1.0)*(wqqSF-1.0)/(HiggsSolver::jetErr*HiggsSolver::jetErr);

    const double meanWlnuMass = jetM > 60 ? HiggsSolver::offWlnuMeanWlnu : HiggsSolver::onWlnuMeanWlnu;
    const double wlnuMassError = jetM > 60
            ? (wlnu.mass() > HiggsSolver::offWlnuMeanWlnu ? HiggsSolver::offWlnuPosWlnuErr : HiggsSolver::offWnluNegWlnuErr)
                    : HiggsSolver::onWlnuWlnuErr;
    const double wlnuMassChiSq = std::pow( (wlnu.mass() -meanWlnuMass)/wlnuMassError,2  );

    const double hWWMassError = jetM > 60 ? HiggsSolver::offWlnuHWWErr : HiggsSolver::onWlnuHWWErr;
    const double hWWMassChiSq = std::pow( (hww.mass() -125)/hWWMassError,2  );


    const double chisq= hWWMassChiSq +wlnuMassChiSq +metPerpChiSq +metParChiSq+jetSFChiSq ;

//    std::cout <<" B "<< chisq<<" "<<neutrino.px()<<" "<<neutrino.py()<<" "<<neutrino.pz()<<" "<<wqqSF <<" "<< extraMetPar/hwwMag <<" "<< extraMetPerp<<" "<< scaledQQJet.pt() <<" "<< jetM<<" " <<" "<< wlnu.mass() <<" "<< hww.mass()  <<" "<<std::endl;
//    std::cout <<" B "<<jetM<<" "<<neutrino.pt()<<" "<<neutE <<" "<< neutrino.E() <<" "<<lepton.pt() <<" "<<lepton.mass() <<" "<< hww.mass() <<" "<< wlnu.mass() <<" "<< chisq <<std::endl;
    return chisq;
}


double HiggsLiFunction::operator()(const double* p){

    setIterationStorage(p);


    LL = 0;
    LL += std::log(pdfs[HiggsLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[HiggsLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[HiggsLi::WQQ_RES]->getProbability(p[WQQ_RES]));
//    LL += std::log(pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hww.mass(),wlnu.mass()));
    LL += std::log(pdfs[HiggsLi::WLNU_MASS]->getProbability(wlnu.mass()));
    LL += std::log(pdfs[HiggsLi::HWW_MASS]->getProbability(hww.mass()));


    // Testing replacing certain terms with the chisq
    //    const double jetSFChiSq = (wqqSF-1.0)*(wqqSF-1.0)/(HiggsSolver::jetErr*HiggsSolver::jetErr);
    //    LL += -1*jetSFChiSq/2.0;
    //    LL += std::log(pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hww.mass(),wlnu.mass()));
    //    const double hWWMassError = jetM > 60 ? HiggsSolver::offWlnuHWWErr : HiggsSolver::onWlnuHWWErr;
    //    const double hWWMassChiSq = std::pow( (hww.mass() -125)/hWWMassError,2  );
    //    const double meanWlnuMass = jetM > 60 ? HiggsSolver::offWlnuMeanWlnu : HiggsSolver::onWlnuMeanWlnu;
    //    const double wlnuMassError = jetM > 60
    //            ? (wlnu.mass() > HiggsSolver::offWlnuMeanWlnu ? HiggsSolver::offWlnuPosWlnuErr : HiggsSolver::offWnluNegWlnuErr)
    //                    : HiggsSolver::onWlnuWlnuErr;
    //    const double wlnuMassChiSq = std::pow( (wlnu.mass() -meanWlnuMass)/wlnuMassError,2  );
    //
    //    LL += -1*hWWMassChiSq/2.0;
    //    LL += -1*wlnuMassChiSq/2.0;



//    std::cout <<" **** "<< -2.0*LL <<" "<<neutX<<" "<<neutY<<" "<<p[NEUT_Z]<<" "<<p[WQQ_RES]<<" "<< hww.mass() <<" "<< wlnu.mass()
//<<" ("<<p[EMET_PERP]<<","<<pdfs[HiggsLi::EMET_PERP]->getProbability(p[EMET_PERP])<<","<<std::log(pdfs[HiggsLi::EMET_PERP]->getProbability(p[EMET_PERP]))
//    <<") ("<<p[EMET_PAR]<<","<<pdfs[HiggsLi::EMET_PAR]->getProbability(p[EMET_PAR])<<","<<std::log(pdfs[HiggsLi::EMET_PAR]->getProbability(p[EMET_PAR]))
//    <<") ("<<pdfs[HiggsLi::WQQ_RES]->getProbability(p[WQQ_RES])<<","<<std::log(pdfs[HiggsLi::WQQ_RES]->getProbability(p[WQQ_RES]))
//    <<") ("<<mP<<","<<std::log(mP)<<")"<<std::endl;


//    std::cout <<" **** "<< -2.0*LL <<" "<<p[EMET_PERP]<<" "<<p[EMET_PAR]<<" "<<p[NEUT_Z]<<" "<<p[WQQ_RES]<<" "<< hww.mass() <<" "<< wlnu.mass()
//<<" ("<<pdfs[HiggsLi::EMET_PERP]->getProbability(p[EMET_PERP])<<","<<std::log(pdfs[HiggsLi::EMET_PERP]->getProbability(p[EMET_PERP]))
//    <<") ("<<pdfs[HiggsLi::EMET_PAR]->getProbability(p[EMET_PAR])<<","<<std::log(pdfs[HiggsLi::EMET_PAR]->getProbability(p[EMET_PAR]))
//    <<") ("<<pdfs[HiggsLi::WQQ_RES]->getProbability(p[WQQ_RES])<<","<<std::log(pdfs[HiggsLi::WQQ_RES]->getProbability(p[WQQ_RES]))
//    <<") ("<<mP<<","<<std::log(mP)<<")"<<std::endl;


    return -2.0*LL;
}


//Testing with the neutrino as a parameter
//double HiggsLiFunction::operator()(const double* p){
//    neutE = std::sqrt(p[NEUT_X]*p[NEUT_X]+p[NEUT_Y]*p[NEUT_Y]+p[NEUT_Z]*p[NEUT_Z]);
//
//
//    neutrino.SetPxPyPzE(p[NEUT_X],p[NEUT_Y],p[NEUT_Z],neutE);
//
//    neutPerp = neutrino.px()*hwwPerpNormX+neutrino.py()*hwwPerpNormY;
//    neutPar  = neutrino.px()*hwwParNormX+neutrino.py()*hwwParNormY;
//
//    wqqSF = 1.0/(p[QQJET_RES]+1.0);
//    jetE   = std::sqrt(
//            wqqSF*wqqSF*(qqJet.px()*qqJet.px()+qqJet.py()*qqJet.py()+qqJet.pz()*qqJet.pz())
//            +jetM*jetM);
//    scaledQQJet.SetPxPyPzE(wqqSF*qqJet.px(),wqqSF*qqJet.py(),wqqSF*qqJet.pz(),jetE);
//
//    wlnu = lepton + neutrino;
//    hww = wlnu + scaledQQJet;
//
//    extraMetPar = (metPar-neutPar)/hwwMag;
//    extraMetPerp = (metPerp-neutPerp);
//
//    double mP = pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hww.mass(),wlnu.mass());
//
//    LL = 0;
//    LL += std::log(pdfs[HiggsLi::EMET_PERP]->getProbability(extraMetPerp));
//    LL += std::log(pdfs[HiggsLi::EMET_PAR]->getProbability(extraMetPar));
//    LL += std::log(pdfs[HiggsLi::WQQ_RES]->getProbability(p[QQJET_RES]));
////    LL += std::log(pdfs[HiggsLi::HWW_WLNU_MASS]->getProbability(hww.mass(),wlnu.mass()));
//    LL += std::log(mP);
//
//    std::cout <<" **** "<< -2.0*LL <<" "<<p[NEUT_X]<<" "<<p[NEUT_Y]<<" "<<p[NEUT_Z]<<" "<<p[QQJET_RES]<<" "<< hww.mass() <<" "<< wlnu.mass()
//<<" ("<<extraMetPerp<<","<<pdfs[HiggsLi::EMET_PERP]->getProbability(extraMetPerp)<<","<<std::log(pdfs[HiggsLi::EMET_PERP]->getProbability(extraMetPerp))
//    <<") ("<<extraMetPar<<","<<pdfs[HiggsLi::EMET_PAR]->getProbability(extraMetPar)<<","<<std::log(pdfs[HiggsLi::EMET_PAR]->getProbability(extraMetPar))
//    <<") ("<<pdfs[HiggsLi::WQQ_RES]->getProbability(p[QQJET_RES])<<","<<std::log(pdfs[HiggsLi::WQQ_RES]->getProbability(p[QQJET_RES]))
//    <<") ("<<mP<<","<<std::log(mP)<<")"<<std::endl;
//
//
//    return -2.0*LL;
//}

HiggsLi::HiggsLi(const std::string& dataDir) : dataDir(dataDir),
        osqq_functor(&osqq_function,&HiggsLiFunction::operator (),4),
        vqq_functor(&vqq_function,&HiggsLiFunction::operator (),4),
        osqq_minimizer( ROOT::Minuit::kSimplex,4 ),
        vqq_minimizer( ROOT::Minuit::kSimplex,4 )
{

    auto setupMinimizer=[&](Minimizer& min, ROOT::Math::Functor& f) {
//        min.SetMaxFunctionCalls(1000000);
//        min.SetMaxIterations(100000);
//        min.SetTolerance(100);
        min.SetFunction(f);
//        min.SetPrintLevel(3);
        resetParameters(min);

    };







    setupMinimizer(vqq_minimizer,vqq_functor);
    setupMinimizer(osqq_minimizer,osqq_functor);

    osqq_pdfs.resize(NPDFS);
    osqq_pdfs[EMET_PERP]    .reset(new OneDimPDFWInterpAndExtrap());
    osqq_pdfs[EMET_PAR]     .reset(new OneDimPDFWInterpAndExtrap());
    osqq_pdfs[WQQ_RES]      .reset(new OneDimPDFWInterpAndExtrap());
    osqq_pdfs[WQQ_SDMASS]   .reset(new OneDimPDFWInterp());
//    osqq_pdfs[HWW_WLNU_MASS].reset(new TwoDimPDF());
    osqq_pdfs[WLNU_MASS].reset(new OneDimPDFWInterpAndExtrap());
    osqq_pdfs[HWW_MASS].reset(new OneDimPDFWInterpAndExtrap());

    vqq_pdfs.resize(NPDFS);
    vqq_pdfs[EMET_PERP].reset(new OneDimPDFWInterpAndExtrap());
    vqq_pdfs[EMET_PAR].reset(new OneDimPDFWInterpAndExtrap());
    vqq_pdfs[WQQ_RES].reset(new OneDimPDFWInterpAndExtrap());
    vqq_pdfs[WQQ_SDMASS].reset(new OneDimPDFWInterp());
//    vqq_pdfs[HWW_WLNU_MASS].reset(new TwoDimPDF());
    vqq_pdfs[WLNU_MASS].reset(new OneDimPDFWInterpAndExtrap());
    vqq_pdfs[HWW_MASS].reset(new OneDimPDFWInterpAndExtrap());

    osqq_function.pdfs.resize(NPDFS);
    osqq_function.pdfs[EMET_PERP]     = osqq_pdfs[EMET_PERP]    ;
    osqq_function.pdfs[EMET_PAR]      = osqq_pdfs[EMET_PAR]     ;
    osqq_function.pdfs[WQQ_RES]       = osqq_pdfs[WQQ_RES]      ;
    osqq_function.pdfs[WQQ_SDMASS]    = osqq_pdfs[WQQ_SDMASS]   ;
//    osqq_function.pdfs[HWW_WLNU_MASS] = osqq_pdfs[HWW_WLNU_MASS];
    osqq_function.pdfs[WLNU_MASS] = osqq_pdfs[WLNU_MASS];
    osqq_function.pdfs[HWW_MASS] = osqq_pdfs[HWW_MASS];
    osqq_function.jetM = 80;
    vqq_function.pdfs.resize(NPDFS);
    vqq_function.pdfs[EMET_PERP]     = vqq_pdfs[EMET_PERP]    ;
    vqq_function.pdfs[EMET_PAR]      = vqq_pdfs[EMET_PAR]     ;
    vqq_function.pdfs[WQQ_RES]       = vqq_pdfs[WQQ_RES]      ;
    vqq_function.pdfs[WQQ_SDMASS]    = vqq_pdfs[WQQ_SDMASS]   ;
//    vqq_function.pdfs[HWW_WLNU_MASS] = vqq_pdfs[HWW_WLNU_MASS];
    vqq_function.pdfs[WLNU_MASS] = vqq_pdfs[WLNU_MASS];
    vqq_function.pdfs[HWW_MASS] = vqq_pdfs[HWW_MASS];
    vqq_function.jetM = 31;

}

void HiggsLi::resetParameters(Minimizer& min,
        const double neutZ,const double bounds) {

//    min.SetLimitedVariable(HiggsLiFunction::EMET_PERP,"EMET_PERP",0,30,-1000,1000);
//    min.SetLimitedVariable(HiggsLiFunction::EMET_PAR,"EMET_PAR",0,0.1,-5,5);
//    min.SetLimitedVariable(HiggsLiFunction::NEUT_Z  ,"NEUT_Z",neutZ,500,neutZ-bounds,neutZ+bounds);
//    min.SetLimitedVariable(HiggsLiFunction::WQQ_RES ,"QQJET_RES",0,0.5,-0.9,5);


    min.SetVariable(HiggsLiFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(HiggsLiFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetVariable(HiggsLiFunction::NEUT_Z  ,"NEUT_Z",neutZ,500);
    min.SetLimitedVariable(HiggsLiFunction::WQQ_RES ,"QQJET_RES",0,0.5,-0.9,5);


//    min.SetLimitedVariable(HiggsLiFunction::EMET_PERP,"EMET_PERP" ,neutZ,500,-3000,3000); //nx
//    min.SetLimitedVariable(HiggsLiFunction::EMET_PAR,"EMET_PAR"  ,bounds,500,-3000,3000);//ny
//    min.SetLimitedVariable(HiggsLiFunction::NEUT_Z  ,"NEUT_Z",0,500,-3000,3000);
//    min.SetLimitedVariable(HiggsLiFunction::WQQ_RES ,"QQJET_RES",1,2    ,-5,5);


}

//void HiggsLi::resetParameters(ROOT::Minuit2::Minuit2Minimizer& min,
//        const double metX, const double metY, const double metZ) {
//
//    min.SetLimitedVariable(HiggsLiFunction::NEUT_X,"NEUT_X",metX,100,-3000,3000);
//    min.SetLimitedVariable(HiggsLiFunction::NEUT_Y,"NEUT_Y",metY,100,-3000,3000);
//    min.SetLimitedVariable(HiggsLiFunction::NEUT_Y,"NEUT_Z",metZ,100,-3000,3000);
//    min.SetLimitedVariable(HiggsLiFunction::QQJET_RES,"QQJET_RES",0,.2,-0.9,5);
//}
void HiggsLi::setup(std::string fileName, const double ptCorB_,const double ptCorM_,
        bool verbose){
    TFile * inFile = TObjectHelper::getFile(dataDir+fileName,"read",verbose);

    hwwPT = TObjectHelper::getObject<TH1>(inFile,"signal_avgHWWMag",verbose);
    const double lowM = hwwPT->GetBinContent(1);
    const double highM = hwwPT->GetBinContent(2);

    ptCorB = ptCorB_;
    ptCorM = ptCorM_;



    auto setup1D = [&] (int entry, const std::string& vName){
        ((OneDimPDFWInterp*)(vqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_low_vqq_")+vName,
                std::string("signal_high_vqq_")+vName,lowM,highM,verbose);
        ((OneDimPDFWInterp*)(osqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_low_osqq_")+vName,
                std::string("signal_high_osqq_")+vName,lowM,highM,verbose);
    };
    auto setup2D = [&] (int entry, const std::string& vName){
        ((TwoDimPDF*)(vqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_vqq_")+vName,verbose);
        ((TwoDimPDF*)(osqq_pdfs[entry].get()))-> setup(inFile,
                std::string("signal_osqq_")+vName,verbose);
    };

    setup1D(EMET_PERP,"extraMetPerp");
    setup1D(EMET_PAR,"extraMetParRelhwwMag");
    setup1D(WQQ_RES,"wqqPTRes");
    setup1D(WQQ_SDMASS,"qqSDMassCoarse");
    setup1D(WLNU_MASS,"Wlnu");
    setup1D(HWW_MASS,"hWW");
//    setup2D(HWW_WLNU_MASS,"hWW_v_Wlnu");

    delete inFile;
}

double HiggsLi::getCorrHWWPT(const double recoPT) const{
    double boundPT = recoPT;
    if(recoPT < 300) boundPT = 300;
    if(recoPT > 2000) boundPT = 2000;
    double corPT = (-1*ptCorB + std::sqrt(ptCorB*ptCorB+4.*ptCorM*boundPT))/(2.*ptCorM);
    return recoPT * corPT/boundPT; //apply as a correction due to the limted range;
}

void HiggsLi::minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
        double qqSDMass, HiggsSolverInfoDebug& out){
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

    auto doFit = [&](
            Minimizer& min, HiggsLiFunction& f,
            std::vector<std::shared_ptr<BASEPDF>>& pdfs, HiggsSolInfoDebug& out  ){

        f.setObservables(lepton,met,qqJet);
        resetParameters(min,0,3000 );
//        resetParameters(min,met.x(),met.y());

        out.minOut=min.Minimize();
        out.noSDLikli = min.MinValue();
        out.likeli = out.noSDLikli - 2.0*std::log(pdfs[WQQ_SDMASS]->getProbability(qqSDMass));

        const double *xs = min.X();
        f.setIterationStorage(xs);
        out.neutrino = f.neutrino;
        out.ptRes = xs[HiggsLiFunction::WQQ_RES];
        out.wqqjet = f.scaledQQJet;
        out.wlnu = f.wlnu;
        out.hWW = f.hww;

        out.emetperp = xs[HiggsLiFunction::EMET_PERP];
        out.emetpar = xs[HiggsLiFunction::EMET_PAR];

    };
//    if(qqSDMass < 60){
//        vqq_minimizer.SetFunction(vqq_functor);
//        doFit(vqq_minimizer,vqq_function,vqq_pdfs,out.vqq_sol);
//        out.min_sol = out.vqq_sol;
//    } else {
//        osqq_minimizer.SetFunction(osqq_functor);
//        doFit(osqq_minimizer,osqq_function,osqq_pdfs,out.osqq_sol);
//        out.min_sol = out.osqq_sol;
//    }
    osqq_minimizer.SetFunction(osqq_functor);
    doFit(osqq_minimizer,osqq_function,osqq_pdfs,out.osqq_sol);
    vqq_minimizer.SetFunction(vqq_functor);
    doFit(vqq_minimizer,vqq_function,vqq_pdfs,out.vqq_sol);
    out.min_sol = out.vqq_sol.likeli < out.osqq_sol.likeli ? out.vqq_sol : out.osqq_sol;
}


void BkgLiFunction::setObservables(const MomentumF& inL,
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

void BkgLiFunction::setIterationStorage(const double * p){

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


double BkgLiFunction::operator()(const double* p){
    setIterationStorage(p);
    LL = 0;
    LL += std::log(pdfs[BkgLi::EMET_PERP]->getProbability(p[EMET_PERP]));
    LL += std::log(pdfs[BkgLi::EMET_PAR]->getProbability(p[EMET_PAR]));
    LL += std::log(pdfs[BkgLi::WLNU_MASS]->getProbability(wlnu.mass()));
    LL += std::log(pdfs[BkgLi::HWW_MASS]->getProbability(hww.mass()));
    return -2.0*LL;
}



BkgLi::BkgLi(const std::string& dataDir) : dataDir(dataDir),
        functor(&function,&BkgLiFunction::operator (),4),
        minimizer( ROOT::Minuit::kSimplex,4 )
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

    function.pdfs.resize(NPDFS);
    function.pdfs[EMET_PERP]     = pdfs[EMET_PERP]    ;
    function.pdfs[EMET_PAR]      = pdfs[EMET_PAR]     ;
    function.pdfs[WLNU_MASS]     = pdfs[WLNU_MASS]      ;
    function.pdfs[HWW_MASS]      = pdfs[HWW_MASS]   ;
}

void BkgLi::resetParameters(Minimizer& min) {
    min.SetVariable(BkgLiFunction::EMET_PERP,"EMET_PERP",0,100);
    min.SetVariable(BkgLiFunction::EMET_PAR,"EMET_PAR",0,0.5);
    min.SetVariable(BkgLiFunction::NEUT_Z  ,"NEUT_Z",0,500);
    min.SetFixedVariable(3,"DUMMY",0);
}
void BkgLi::setup(std::string fileName, bool verbose){
    TFile * inFile = TObjectHelper::getFile(dataDir+fileName,"read",verbose);

//    auto setup1D = [&] (int entry, const std::string& vName){
//        ((OneDimPDFWExtrap*)(pdfs[entry].get()))-> setup(inFile,
//                std::string("ttbarPW_")+vName,verbose);
//    };
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


void BkgLi::minimize(const MomentumF& lepton, const MomentumF& met, const MomentumF& qqJet,
        HiggsSolverInfoDebug& out){

    const double oHWWMag = (lepton.p4()+met.p4()+qqJet.p4()).pt();
    auto interp = [&](PDFList p) {
      ((OneDimPDFWInterp*)(pdfs[p].get()))->setInterp(oHWWMag);

    };
    interp(EMET_PERP);
    interp(EMET_PAR);
    interp(WLNU_MASS);
    interp(HWW_MASS);

    auto doFit = [&](
            Minimizer& min, BkgLiFunction& f,HiggsSolInfoDebug& out  ){

        f.setObservables(lepton,met,qqJet);
        resetParameters(min);

        out.minOut=min.Minimize();

        out.likeli = min.MinValue();
        const double *xs = min.X();
        f.setIterationStorage(xs);
        out.neutrino = f.neutrino;
        out.wqqjet = f.qqJet;
        out.wlnu = f.wlnu;
        out.hWW = f.hww;
        out.emetperp = xs[BkgLiFunction::EMET_PERP];
        out.emetpar = xs[BkgLiFunction::EMET_PAR];

    };

    minimizer.SetFunction(functor);
    doFit(minimizer,function,out.min_sol);
}



}




