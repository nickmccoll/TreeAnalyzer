
#include "../interface/BTagScaleFactors.h"
#include "DataFormats/interface/Jet.h"
#include "DataFormats/interface/FatJet.h"
#include "Processors/Corrections/interface/BTagCalibrationStandalone.h"



namespace TAna {
using namespace CorrHelp;
using namespace BTagging;

//_____________________________________________________________________________
BTagScaleFactors::BTagScaleFactors(BTagging::CSVWP maxWP, std::string ht,std::string lt, const std::string& dataDir, const std::string& sfFile, const std::string& effFile, bool verbose) :
                        calib(new BTagCalibration("CSVv2", dataDir+sfFile))
{

    auto makeNewReader = [&] (const BTagEntry::OperatingPoint op){
        calibReaders.emplace_back(new BTagCalibrationReader(op, "central",{"up", "down"}));
        calibReaders.back()->load(*calib, BTagEntry::FLAV_B,ht);
        calibReaders.back()->load(*calib, BTagEntry::FLAV_C,ht);
        calibReaders.back()->load(*calib, BTagEntry::FLAV_UDSG,lt);

    };
    if(maxWP >= CSV_L) makeNewReader(BTagEntry::OP_LOOSE );
    if(maxWP >= CSV_M) makeNewReader(BTagEntry::OP_MEDIUM);
    if(maxWP >= CSV_T) makeNewReader(BTagEntry::OP_TIGHT );

    TFile * efile = TObjectHelper::getFile(dataDir+effFile,"read",verbose);
    std::vector<TString> taggers {"loose","med","tight"};
    std::vector<TString> flvs {"b","c","l"};
    for(unsigned int iF = FLV_B; iF <= FLV_L; ++iF){
        efficiencies.push_back(std::vector< std::unique_ptr<TObjectHelper::Hist2DContainer> >());
        for(unsigned int iT = CSV_L; iT <= maxWP; ++iT){ //the tagger is displaced by 1 for the inclusive entry
            efficiencies[iF].emplace_back(new  TObjectHelper::Hist2DContainer(efile,TString::Format("%s_%s",flvs[iF].Data(),taggers[iT-1].Data()).Data(),verbose));
        }
    }
}
//_____________________________________________________________________________
BTagScaleFactors::~BTagScaleFactors() {}
//_____________________________________________________________________________
float BTagScaleFactors::getJetEff(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp) const {
    return efficiencies[flv][wp -1]->getBinContentByValue(pt, eta < 0 ? -eta : eta).val();
}
//_____________________________________________________________________________
float BTagScaleFactors::getJetSF(const BTagging::FLAVOR flv, const float pt, const float eta, const BTagging::CSVWP wp,
        CorrHelp::CORRTYPE corrT) const{
    //wp[0] = inclusive
    return calibReaders[wp-1]->eval_auto_bounds(systNames[corrT],BTagEntry::JetFlavor(flv), eta,pt  );   // BTagEntry::JetFlavor set to have same order as BTagging::Flavor
}
//_____________________________________________________________________________
float BTagScaleFactors::getJetCorr(const BaseRecoJet* jet,  CorrHelp::CORRTYPE lightT,  CorrHelp::CORRTYPE heavyT) const {
    const float pt = jet->pt();
    const float eta = jet->eta();
    const auto flv = jetFlavor(*jet);
    const float csv = jet->csv();
    const CorrHelp::CORRTYPE corrT =  flv == FLV_L ? lightT : heavyT;
    if(corrT == NONE) return 1.0;

    float lE = 1;
    float hE = 1;
    float lSF = 1;
    float hSF = 1;

    assignVals(flv,pt,eta,csv,corrT,lE,hE,lSF,hSF);

    const float origEff = lE - hE;
    const float newEff = lSF*lE - hSF*hE;
    if(origEff == 0 || newEff == 0 || origEff > 1.0 || newEff > 1.0){
        TString errStr = TString::Format(
                "BTagScaleFactors::getJetCorr() -> Bad efficiencies: corr(pt,eta,flv):lowEff,highEff,lowSF,highSF,origEff,newEff -> %u(%f,%f,%u)%f,%f,%f,%f,%f,%f",
                corrT,pt,eta,flv,lE,hE,lSF,hSF,origEff,newEff);
        throw std::invalid_argument(errStr.Data());
    }

    return newEff/origEff;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
JetBTagScaleFactors::JetBTagScaleFactors(const std::string& dataDir, const std::string& sfFile, const std::string& effFile, bool verbose) :
                BTagScaleFactors(CSV_T,"comb","incl", dataDir, sfFile,effFile,verbose) {}
//_____________________________________________________________________________
JetBTagScaleFactors::~JetBTagScaleFactors() {}
//_____________________________________________________________________________
void JetBTagScaleFactors::assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv, const CorrHelp::CORRTYPE corrT,
        float& lE,float& hE,float& lSF,float& hSF) const {
    auto gE = [&](const BTagging::CSVWP wp)->float{return getJetEff(flv,pt,eta,wp);};
    auto gS = [&](const BTagging::CSVWP wp)->float{return getJetSF(flv,pt,eta,wp,corrT);};

    if(csv <  CSVWP_VALS[CSV_L]){
        lE = 1.0;
        hE = gE(CSV_L);
        lSF = 1.0;
        hSF = gS(CSV_L);
    } else if(csv <  CSVWP_VALS[CSV_M]){
        lE = gE(CSV_L);
        hE = gE(CSV_M);
        lSF = gS(CSV_L);
        hSF = gS(CSV_M);
    } else if(csv <  CSVWP_VALS[CSV_T]){
        lE = gE(CSV_M);
        hE = gE(CSV_T);
        lSF = gS(CSV_M);
        hSF = gS(CSV_T);
    } else {
        lE = gE(CSV_T);
        hE = 0;
        lSF = gS(CSV_T);
        hSF = 0;
    }
}
//_____________________________________________________________________________
float JetBTagScaleFactors::getSF(const std::vector<const Jet*>& jets, CorrHelp::CORRTYPE lightT,  CorrHelp::CORRTYPE heavyT) const{
    float SF = 1.0;
    for(const auto* j : jets) SF *= getJetCorr(j,lightT,heavyT);
    return SF;
}
//_____________________________________________________________________________
//_____________________________________________________________________________

SubJetBTagScaleFactors::SubJetBTagScaleFactors(const std::string& dataDir, const std::string& sfFile, const std::string& effFile, bool verbose) :
                BTagScaleFactors(CSV_M,"lt","incl",dataDir, sfFile,effFile,verbose) {}
//_____________________________________________________________________________
SubJetBTagScaleFactors::~SubJetBTagScaleFactors() {}
//_____________________________________________________________________________
void SubJetBTagScaleFactors::assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv, const CorrHelp::CORRTYPE corrT,
        float& lE,float& hE,float& lSF,float& hSF) const {
    auto gE = [&](const BTagging::CSVWP wp)->float{return getJetEff(flv,pt,eta,wp);};
    auto gS = [&](const BTagging::CSVWP wp)->float{return getJetSF(flv,pt,eta,wp,corrT);};

    if(csv <  CSVWP_VALS[CSV_L]){
        lE = 1.0;
        hE = gE(CSV_L);
        lSF = 1.0;
        hSF = gS(CSV_L);
    } else if(csv <  CSVWP_VALS[CSV_M]){
        lE = gE(CSV_L);
        hE = gE(CSV_M);
        lSF = gS(CSV_L);
        hSF = gS(CSV_M);
    } else {
        lE = gE(CSV_M);
        hE = 0;
        lSF = gS(CSV_M);
        hSF = 0;

    }
}
//_____________________________________________________________________________
float SubJetBTagScaleFactors::getSF(const std::vector<const FatJet*>& fatJets, CorrHelp::CORRTYPE lightT,  CorrHelp::CORRTYPE heavyT) const{
    std::vector<const BaseRecoJet*> subjets;
    for(const auto* fj : fatJets){
        if(fj){
            for(const auto& sj : fj->subJets()){
                if(sj.absEta() < 2.4 && sj.pt() >= 30) subjets.push_back(static_cast<const BaseRecoJet*>(&sj) );
            }
        }
    }

    float SF = 1.0;
    for(const auto* j : subjets) SF *= getJetCorr(j,lightT,heavyT);
    return SF;
}
}
