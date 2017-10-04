
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "DataFormats/interface/FatJet.h"

namespace TAna {
SoftDropMassScaleFactors::SoftDropMassScaleFactors(const std::string& dataDir, const std::string& sfFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    puppisd_corrGEN     .reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_gen",verbose) );
    puppisd_corrRECO_cen.reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_reco_0eta1v3",verbose) );
    puppisd_corrRECO_for.reset(new  TObjectHelper::TF1Container(file,"puppiJECcorr_reco_1v3eta2v5",verbose) );
    delete file;
}

float SoftDropMassScaleFactors::getScaleFactor(const MomentumF* mom) const {
    return puppisd_corrGEN->eval( mom->pt() )
            *  ( mom->absEta() <= 1.3 ? puppisd_corrRECO_cen : puppisd_corrRECO_for  )->eval(mom->pt());
}
float SoftDropMassScaleFactors::getCorrSDMass(const FatJet* fj) const {
    return fj->rawSdMom().mass() * getScaleFactor(fj);
}
}
