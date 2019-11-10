
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "TreeReaders/interface/EventReader.h"
#include "DataFormats/interface/GenParticle.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {
namespace EventWeights {
float calcNormalizedEventWeight(const EventReader& reader_event, const float cross,
        const float numE, const float lumi) {
    return (reader_event.weight > 0 ? 1.0 : -1.0) *lumi * cross *1000 /numE;
}
float getNormalizedEventWeight(const EventReader& reader_event, const float cross,
        const float numE, const EventParameters& evtParam, const float genMtt, const int nLepsTT) {
    if(cross < 0 ||numE < 0) return reader_event.weight;

    float sgn = (reader_event.weight > 0 ? 1.0 : -1.0);

    if (evtParam.doTTBarStitching && reader_event.process.val() == FillerConstants::TTBAR) {
    	if (reader_event.sampParam.val() == 1000) {
    		if (genMtt >= 1010) {
    			if      (nLepsTT==0) return evtParam.ttbarXSecSF_1000toInf_nLep0 * evtParam.lumi * sgn;
    			else if (nLepsTT==1) return evtParam.ttbarXSecSF_1000toInf_nLep1 * evtParam.lumi * sgn;
    			else if (nLepsTT==2) return evtParam.ttbarXSecSF_1000toInf_nLep2 * evtParam.lumi * sgn;
    		} else return 0;
    	} else if (reader_event.sampParam.val() == 700) {
    		if (genMtt >= 710 && genMtt <= 960) {
    			if      (nLepsTT==0) return evtParam.ttbarXSecSF_700to1000_nLep0 * evtParam.lumi * sgn;
    			else if (nLepsTT==1) return evtParam.ttbarXSecSF_700to1000_nLep1 * evtParam.lumi * sgn;
    			else if (nLepsTT==2) return evtParam.ttbarXSecSF_700to1000_nLep2 * evtParam.lumi * sgn;
    		} else return 0;
    	} else {
    		if (genMtt >= 1010) {
    			if      (nLepsTT==0) return evtParam.ttbarXSecSF_1000toInf_nLep0 * evtParam.lumi * sgn;
    			else if (nLepsTT==1) return evtParam.ttbarXSecSF_1000toInf_nLep1 * evtParam.lumi * sgn;
    			else if (nLepsTT==2) return evtParam.ttbarXSecSF_1000toInf_nLep2 * evtParam.lumi * sgn;
    		} else if (genMtt >= 710 && genMtt <= 960) {
    			if      (nLepsTT==0) return evtParam.ttbarXSecSF_700to1000_nLep0 * evtParam.lumi * sgn;
    			else if (nLepsTT==1) return evtParam.ttbarXSecSF_700to1000_nLep1 * evtParam.lumi * sgn;
    			else if (nLepsTT==2) return evtParam.ttbarXSecSF_700to1000_nLep2 * evtParam.lumi * sgn;
    		}
    	}
    }
    return calcNormalizedEventWeight(reader_event,cross,numE,evtParam.lumi);
}
}


PUScaleFactors::PUScaleFactors(const std::string& dataDir){
	dataDirectory = dataDir;
}

void PUScaleFactors::setParameters(const EventParameters& evtParam, bool verbose) {
    TFile * file = TObjectHelper::getFile(dataDirectory+evtParam.puCorrSFFile,"read",verbose);
    nominalSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_nom",verbose) );
    downSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_down",verbose) );
    upSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_up",verbose) );
    delete file;
}

float PUScaleFactors::getCorrection(const unsigned int trueNumInteractions,
        const CorrHelp::CORRTYPE corrType) const {
    switch(corrType) {
    case CorrHelp::NOMINAL :
        return  nominalSF->getBinContentByValue(trueNumInteractions).val();
    case CorrHelp::NONE :
        return  1.0;
    case CorrHelp::UP :
        return  upSF->getBinContentByValue(trueNumInteractions).val();
    case CorrHelp::DOWN :
        return  downSF->getBinContentByValue(trueNumInteractions).val();
    default:
        return 1.0;
    }
}

TopPTWeighting::TopPTWeighting(const std::string& dataDir, const std::string& sfFile,
        bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    weightConsts.reset(new  TObjectHelper::Hist1DContainer(file,"weightConsts",verbose) );
    a  = weightConsts->getBinContentByBinNumber(1).val();
    b  = weightConsts->getBinContentByBinNumber(2).val();
    nf = weightConsts->getBinContentByBinNumber(3).val();
    delete file;
}

float TopPTWeighting::getCorrection(const ASTypes::size8 process,
        const SMDecayEvent& decayEvent) const{
    if(process != FillerConstants::TTBAR) return 1;
    return std::exp(a +b*getAvgPT(decayEvent))*nf;
}

float TopPTWeighting::getCorrectionNoNorm(const ASTypes::size8 process,
        const SMDecayEvent& decayEvent) const{
    if(process != FillerConstants::TTBAR) return 1;
    return std::exp(a +b*getAvgPT(decayEvent));
}
float TopPTWeighting::getAvgPT(const SMDecayEvent& decayEvent) const {
    if(decayEvent.topDecays.size() != 2) return 0;
    return std::sqrt( std::min(decayEvent.topDecays[0].top->pt(),float(800.0)) *
            std::min(decayEvent.topDecays[1].top->pt(),float(800.0)));
}


}



