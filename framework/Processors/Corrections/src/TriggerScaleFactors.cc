
#include "Processors/Corrections/interface/TriggerScaleFactors.h"

namespace TAna {
TriggerScaleFactors::TriggerScaleFactors(const std::string& dataDir){
    dataDirectory = dataDir;
}

void TriggerScaleFactors::setParameters(const EventParameters& evtParam, bool verbose) {
	TFile *file = TObjectHelper::getFile(dataDirectory+evtParam.leptonCorrSFFile,"read",verbose);
    electronSFs_1l.reset(new  TObjectHelper::Hist1DContainer(file,"electronSF_lnuqq",verbose) );
    muonSFs_1l.reset(new  TObjectHelper::Hist1DContainer(file,"muonSF_lnuqq",verbose) );

    electronEffs_mcHT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"electronEffs_mcHT_lnulnu",verbose) );
    electronEffs_mcPT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"electronEffs_mcPT_lnulnu",verbose) );
    electronEffs_dataHT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"electronEffs_dataHT_lnulnu",verbose) );
    electronEffs_dataPT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"electronEffs_dataPT_lnulnu",verbose) );
    muonEffs_mcHT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"muonEffs_mcHT_lnulnu",verbose) );
    muonEffs_mcPT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"muonEffs_mcPT_lnulnu",verbose) );
    muonEffs_dataHT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"muonEffs_dataHT_lnulnu",verbose) );
    muonEffs_dataPT_2l.reset(new  TObjectHelper::Hist1DContainer(file,"muonEffs_dataPT_lnulnu",verbose) );

    delete file;
}

float TriggerScaleFactors::getElectronTriggerSF_1l(const float ht) const{
    return electronSFs_1l->getBinContentByValue(ht).val();
}
float TriggerScaleFactors::getMuonTriggerSF_1l(const float ht) const{
    return muonSFs_1l->getBinContentByValue(ht).val();
}
float TriggerScaleFactors::getSingleLeptonTriggerSF(const float ht, const bool leadingLepIsMuon) const {
    if(leadingLepIsMuon) return getMuonTriggerSF_1l(ht) ;
    return getElectronTriggerSF_1l(ht);
}

float TriggerScaleFactors::getMuonTriggerEff_2l(const float value, const bool getHT, const bool isData) const {
	if(getHT) {
		if(isData) return muonEffs_dataHT_2l->getBinContentByValue(value).val();
		else       return muonEffs_mcHT_2l->getBinContentByValue(value).val();
	} else {
		if(isData) return muonEffs_dataPT_2l->getBinContentByValue(value).val();
		else       return muonEffs_mcPT_2l->getBinContentByValue(value).val();
	}
}

float TriggerScaleFactors::getElectronTriggerEff_2l(const float value, const bool getHT, const bool isData) const {
	if(getHT) {
		if(isData) return electronEffs_dataHT_2l->getBinContentByValue(value).val();
		else       return electronEffs_mcHT_2l->getBinContentByValue(value).val();
	} else {
		if(isData) return electronEffs_dataPT_2l->getBinContentByValue(value).val();
		else       return electronEffs_mcPT_2l->getBinContentByValue(value).val();
	}
}

float TriggerScaleFactors::getDileptonTriggerSF(const float ht, const float pt2,
		const bool lep1IsMuon, const bool lep2IsMuon) const {

	float eff_ht_mc=0, eff_pt2_mc=0, eff_ht_data=0, eff_pt2_data=0;

	if(lep1IsMuon) {
		eff_ht_mc   = getMuonTriggerEff_2l(ht,true,false);
		eff_ht_data = getMuonTriggerEff_2l(ht,true,true);
	} else {
		eff_ht_mc   = getElectronTriggerEff_2l(ht,true,false);
		eff_ht_data = getElectronTriggerEff_2l(ht,true,true);
	}

	if(lep2IsMuon) {
		eff_pt2_mc   = getMuonTriggerEff_2l(pt2,false,false);
		eff_pt2_data = getMuonTriggerEff_2l(pt2,false,true);
	} else {
		eff_pt2_mc   = getElectronTriggerEff_2l(pt2,false,false);
		eff_pt2_data = getElectronTriggerEff_2l(pt2,false,true);
	}

	float mcEff = eff_ht_mc + eff_pt2_mc - eff_ht_mc*eff_pt2_mc;
	float dataEff = eff_ht_data + eff_pt2_data - eff_ht_data*eff_pt2_data;

	return (dataEff / mcEff);
}
}



