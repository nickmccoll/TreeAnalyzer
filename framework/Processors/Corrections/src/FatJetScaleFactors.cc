
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "DataFormats/interface/FatJet.h"

namespace TAna {
HbbFatJetScaleFactors::HbbFatJetScaleFactors(const std::string& dataDir, const std::string& sfFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    sdMassCorr_eta_0p0_0p8.reset(new  TObjectHelper::GraphEContainer(file,"sdMassCorr_eta_0p0_0p8",true,300,2400,verbose));
    sdMassCorr_eta_0p8_1p6.reset(new  TObjectHelper::GraphEContainer(file,"sdMassCorr_eta_0p8_1p6",true,300,2400,verbose));
    sdMassCorr_eta_1p6_2p4.reset(new  TObjectHelper::GraphEContainer(file,"sdMassCorr_eta_1p6_2p4",true,300,2400,verbose));
    jec_eta_0p0_0p6       .reset(new  TObjectHelper::TF1Container   (file,"jec_eta_0p0_0p6"       ,true,200,2000,verbose));
    jec_eta_0p6_1p4       .reset(new  TObjectHelper::TF1Container   (file,"jec_eta_0p6_1p4"       ,true,200,2000,verbose));
    jec_eta_1p4_1p8       .reset(new  TObjectHelper::TF1Container   (file,"jec_eta_1p4_1p8"       ,true,200,1000,verbose));
    jec_eta_1p8_2p0       .reset(new  TObjectHelper::TF1Container   (file,"jec_eta_1p8_2p0"       ,true,200,1000,verbose));
    jec_eta_2p0_2p4       .reset(new  TObjectHelper::TF1Container   (file,"jec_eta_2p0_2p4"       ,true,200,1000,verbose));
    delete file;
}
//--------------------------------------------------------------------------------------------------
float HbbFatJetScaleFactors::getMassScaleFactor(const float corrJetPT, const float jetAbsETA) const {
    if(jetAbsETA < 0.8)
        return sdMassCorr_eta_0p0_0p8->eval(corrJetPT);
    else if(jetAbsETA < 1.6)
        return sdMassCorr_eta_0p0_0p8->eval(corrJetPT);
    else
        return sdMassCorr_eta_0p8_1p6->eval(corrJetPT);
}
//--------------------------------------------------------------------------------------------------
float HbbFatJetScaleFactors::getJEC(const FatJet* fj) const {
    const float jetAbsETA = fj->absEta();
    const float jetPT = fj->pt();
    if(jetAbsETA < 0.6)
        return jec_eta_0p0_0p6->eval(jetPT);
    else if(jetAbsETA < 1.4)
        return jec_eta_0p6_1p4->eval(jetPT);
    else if(jetAbsETA < 1.8)
        return jec_eta_1p4_1p8->eval(jetPT);
    else if(jetAbsETA < 2.0)
        return jec_eta_1p8_2p0->eval(jetPT);
    else
        return jec_eta_2p0_2p4->eval(jetPT);
}
//--------------------------------------------------------------------------------------------------
float HbbFatJetScaleFactors::getCorrSDMass(const FatJet* fj) const {
    return fj->rawSdMom().mass() * getMassScaleFactor(fj->pt()*getJEC(fj), fj->absEta());
}
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
SoftDropMassScaleFactorsFromW::SoftDropMassScaleFactorsFromW(const std::string& dataDir, const std::string& sfFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    puppisd_corrGEN     .reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_gen",verbose) );
    puppisd_corrRECO_cen.reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_reco_0eta1v3",verbose) );
    puppisd_corrRECO_for.reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_reco_1v3eta2v5",verbose) );
    delete file;
}
//--------------------------------------------------------------------------------------------------
float SoftDropMassScaleFactorsFromW::getScaleFactor(const MomentumF* mom) const {
    return puppisd_corrGEN->eval( mom->pt() )
            *  ( mom->absEta() <= 1.3 ? puppisd_corrRECO_cen : puppisd_corrRECO_for  )->eval(mom->pt());
}
//--------------------------------------------------------------------------------------------------
float SoftDropMassScaleFactorsFromW::getCorrSDMass(const FatJet* fj) const {
    return fj->rawSdMom().mass() * getScaleFactor(fj);
}

}
