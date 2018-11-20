
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "TreeReaders/interface/EventReader.h"
#include "DataFormats/interface/GenParticle.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {
namespace EventWeights {
float calcNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi) {
    return (reader_event.weight > 0 ? 1.0 : -1.0) *lumi * cross *1000 /numE;
}
float getNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi) {
    if(cross < 0 ||numE < 0) return reader_event.weight;
    return calcNormalizedEventWeight(reader_event,cross,numE,lumi);
}
float get4bXSecLimit(size mass) {
    static const std::vector<std::pair<size,float>> sigMassesAnd4bXSecLimit = {
            std::pair<size,float>(600, 120),
            std::pair<size,float>(800, 50),
            std::pair<size,float>(1000,20),
            std::pair<size,float>(1200,10),
            std::pair<size,float>(1400,6),
            std::pair<size,float>(1600,5),
            std::pair<size,float>(1800,4),
            std::pair<size,float>(2000,3),
            std::pair<size,float>(2500,2),
            std::pair<size,float>(3000,2),
            std::pair<size,float>(3500,2),
            std::pair<size,float>(4000,2),
            std::pair<size,float>(4500,2)
    };
    static const size nSigMassesAnd4bXSecLimit= sigMassesAnd4bXSecLimit.size();
    for(unsigned int iB = 0; iB < nSigMassesAnd4bXSecLimit; ++iB){
        if(mass == sigMassesAnd4bXSecLimit[iB].first) return sigMassesAnd4bXSecLimit[iB].second;
        if(mass < sigMassesAnd4bXSecLimit[iB].first  ){
            if(iB == 0) return  sigMassesAnd4bXSecLimit[iB].second;
            return sigMassesAnd4bXSecLimit[iB-1].second +
                    float(mass)*(sigMassesAnd4bXSecLimit[iB].second - sigMassesAnd4bXSecLimit[iB-1].second)
                    /float(sigMassesAnd4bXSecLimit[iB].first - sigMassesAnd4bXSecLimit[iB-1].first);
        }
    }
    return sigMassesAnd4bXSecLimit[nSigMassesAnd4bXSecLimit-1].second;
}
}


PUScaleFactors::PUScaleFactors(const std::string& dataDir, const std::string& sfFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    nominalSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_nominal",verbose) );
    downSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_down",verbose) );
    upSF.reset(new  TObjectHelper::Hist1DContainer(file,"puSF_up",verbose) );
    delete file;
}
float PUScaleFactors::getCorrection(const unsigned int trueNumInteractions, const CorrHelp::CORRTYPE corrType) const {
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

TopPTWeighting::TopPTWeighting(const std::string& dataDir, const std::string& sfFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+sfFile,"read",verbose);
    weightConsts.reset(new  TObjectHelper::Hist1DContainer(file,"weightConsts",verbose) );
    a  = weightConsts->getBinContentByBinNumber(1).val();
    b  = weightConsts->getBinContentByBinNumber(2).val();
    nf = weightConsts->getBinContentByBinNumber(3).val();
    delete file;
}

float TopPTWeighting::getCorrection(const ASTypes::size8 process,const SMDecayEvent& decayEvent) const{
    if(process != FillerConstants::TTBAR) return 1;
    return std::exp(a +b*getAvgPT(decayEvent))*nf;
}

float TopPTWeighting::getCorrectionNoNorm(const ASTypes::size8 process,const SMDecayEvent& decayEvent) const{
    if(process != FillerConstants::TTBAR) return 1;
    return std::exp(a +b*getAvgPT(decayEvent));
}
float TopPTWeighting::getAvgPT(const SMDecayEvent& decayEvent) const {
    if(decayEvent.topDecays.size() != 2) return 0;
    return std::sqrt( std::min(decayEvent.topDecays[0].top->pt(),float(800.0)) *
            std::min(decayEvent.topDecays[1].top->pt(),float(800.0)));
}


}



