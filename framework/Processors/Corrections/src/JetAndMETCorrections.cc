
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"



namespace TAna {
using namespace CorrHelp;
//--------------------------------------------------------------------------------------------------
JERCorrector::JERCorrector (const std::string& dataDir,const std::string& resFile,std::shared_ptr<TRandom3> rndGen, const CorrHelp::CORRTYPE cT )
: cT(cT),resObj(dataDir+resFile),rndGen(rndGen){

}
//--------------------------------------------------------------------------------------------------
float JERCorrector::getSF(const CorrHelp::CORRTYPE cT, const float absETA) const {
    if(cT == CorrHelp::NONE) return 1;
    if(absETA < 0.522 ) return (cT == CorrHelp::NOMINAL ? 1.1595 : (cT== CorrHelp::DOWN ? 1.095  : 1.224  ));
    if(absETA < 0.783 ) return (cT == CorrHelp::NOMINAL ? 1.1948 : (cT== CorrHelp::DOWN ? 1.1296 : 1.26   ));
    if(absETA < 1.131 ) return (cT == CorrHelp::NOMINAL ? 1.1464 : (cT== CorrHelp::DOWN ? 1.0832 : 1.2096 ));
    if(absETA < 1.305 ) return (cT == CorrHelp::NOMINAL ? 1.1609 : (cT== CorrHelp::DOWN ? 1.0584 : 1.2634 ));
    if(absETA < 1.740 ) return (cT == CorrHelp::NOMINAL ? 1.1278 : (cT== CorrHelp::DOWN ? 1.0292 : 1.2264 ));
    if(absETA < 1.930 ) return (cT == CorrHelp::NOMINAL ? 1.1000 : (cT== CorrHelp::DOWN ? 0.9921 : 1.2079 ));
    if(absETA < 2.043 ) return (cT == CorrHelp::NOMINAL ? 1.1426 : (cT== CorrHelp::DOWN ? 1.0212 : 1.264  ));
    if(absETA < 2.322 ) return (cT == CorrHelp::NOMINAL ? 1.1512 : (cT== CorrHelp::DOWN ? 1.0372 : 1.2652 ));
    if(absETA < 2.500 ) return (cT == CorrHelp::NOMINAL ? 1.2963 : (cT== CorrHelp::DOWN ? 1.0592 : 1.5334 ));
    if(absETA < 2.853 ) return (cT == CorrHelp::NOMINAL ? 1.3418 : (cT== CorrHelp::DOWN ? 1.1327 : 1.5509 ));
    if(absETA < 2.964 ) return (cT == CorrHelp::NOMINAL ? 1.7788 : (cT== CorrHelp::DOWN ? 1.578  : 1.9796 ));
    if(absETA < 3.139 ) return (cT == CorrHelp::NOMINAL ? 1.1869 : (cT== CorrHelp::DOWN ? 1.0626 : 1.3112 ));
    else                return (cT == CorrHelp::NOMINAL ? 1.1922 : (cT== CorrHelp::DOWN ? 1.0434 : 1.341  ));

}
//--------------------------------------------------------------------------------------------------
void JERCorrector::processJets(JetReader& jetreader,Met& met,const GenJetCollection& genjets, const float rho){
    if(cT == CORRTYPE::NONE) return;
    JMEStand::JetParameters parameters;
    parameters.setRho(std::min(40.0f,rho));
    std::vector<const GenJet*> gjptrs; gjptrs.reserve(genjets.size());
//    std::cout << "START! "<< rho <<std::endl;
    for(auto& gj: genjets) gjptrs.push_back(&gj);
    double deltaMX= 0;
    double deltaMY= 0;

    for(auto& j : jetreader.jets){
//        std::cout << j.pt() <<" "<< j.eta() <<"->";
        if(j.pt() < 15) {
//            std::cout <<std::endl;
            continue;
        }
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);
        auto oRawFact = j.toRawFactor();
        auto ptSF = correctJet(&j,gjptrs,parameters);
        deltaMX += ((*jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
        deltaMY += ((*jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
//        std::cout << std::sqrt((*jetreader.metUnc_rawPx)[j.index()]*(*jetreader.metUnc_rawPx)[j.index()] + (*jetreader.metUnc_rawPy)[j.index()]*(*jetreader.metUnc_rawPy)[j.index()] )/oRawFact <<","<<
//                (std::sqrt((*jetreader.metUnc_rawPx)[j.index()]*(*jetreader.metUnc_rawPx)[j.index()] + (*jetreader.metUnc_rawPy)[j.index()]*(*jetreader.metUnc_rawPy)[j.index()] )/oRawFact)*(1.0-ptSF) <<std::endl;
//        std::cout << j.pt() <<" "<< j.eta() <<std::endl;
    }

    ASTypes::CylLorentzVectorF deltaV(std::sqrt(deltaMX*deltaMX+deltaMY*deltaMY),0.0f, (deltaMX == 0.0 && deltaMY == 0.0) ? 0 : std::atan2(deltaMY, deltaMX),0 );
    met.p4() += deltaV;
//    std::cout << "END! "<< rho <<std::endl;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void JERCorrector::processFatJets(FatJetCollection& jets,const GenFatJetCollection& genjets, const float rho) {
    if(cT == CORRTYPE::NONE) return;
    JMEStand::JetParameters parameters;
    parameters.setRho(std::min(40.0f,rho));
    std::vector<const GenJet*> gjptrs; gjptrs.reserve(genjets.size());
    for(auto& gj: genjets) gjptrs.push_back(&gj);
    for(auto& j : jets){
        if(j.pt() < 15) continue;
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);
        correctJet(&j,gjptrs,parameters);
    }
    std::sort(jets.begin(),jets.end(),PhysicsUtilities::greaterPT<FatJet>());
}

//--------------------------------------------------------------------------------------------------
float JERCorrector::correctJet(Jet* jet, const std::vector<const GenJet*> genjets, JMEStand::JetParameters& params){
    const float jres = resObj.evaluateFormula(*resObj.getRecord(params),params);
    const float resSF = getSF(cT,jet->absEta());

    double nearDR = 0;
    int idx = PhysicsUtilities::findNearestDRDeref(*jet,genjets,nearDR);
//    std::cout <<"(";
//    if(idx >= 0) std::cout <<genjets[idx]->pt()<<",";
//    else std::cout << "-1,";
//    std::cout << jres <<","<<resSF<<",";
    float ptSF = 1.0;
    if(idx>= 0 && nearDR < 0.4/2.0 && std::fabs(jet->pt()-genjets[idx]->pt()) <3.0*jres*jet->pt() ){
        ptSF = 1 + (resSF-1.0)*(jet->pt()-genjets[idx]->pt())/jet->pt();
//        std::cout <<(jet->pt()-genjets[idx]->pt())/jet->pt() <<","<<ptSF<<")";
    } else{
        ptSF = 1 + rndGen->Gaus(0.,jres)*std::sqrt(std::max(0.0f,resSF*resSF-1));
//        std::cout <<ptSF<<")";
    }
        if(ptSF>0){
            jet->setP4(ptSF*jet->pt(),jet->eta(),jet->phi(),ptSF*jet->mass());
            jet->setRawFactor(jet->toRawFactor()/ptSF);
        } else{
            jet->setP4(0.0f,jet->eta(),jet->phi(),0.0f);
        }

    return std::max(0.0f,ptSF);
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
JESUncShifter::JESUncShifter (const CorrHelp::CORRTYPE cT)
: cT(cT){}
//--------------------------------------------------------------------------------------------------
void JESUncShifter::processJets(JetReader& jetreader,Met& met){
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    double deltaMX= 0;
    double deltaMY= 0;
    for(auto& j : jetreader.jets){
        auto oRawFact = j.toRawFactor();
        auto ptSF = correctJet(&j);
        if(j.pt() >= 15) {
            deltaMX += ((*jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
            deltaMY += ((*jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
        }
    }
    ASTypes::CylLorentzVectorF deltaV(std::sqrt(deltaMX*deltaMX+deltaMY*deltaMY),0.0f, (deltaMX == 0.0 && deltaMY == 0.0) ? 0 : std::atan2(deltaMY, deltaMX),0 );
    met.p4() += deltaV;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void JESUncShifter::processFatJets(FatJetCollection& jets){
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    for(auto& j : jets){
        correctJet(&j);
        for(unsigned int iSJ= 0; iSJ < j.nSubJets(); ++iSJ)
            correctJet(&j.subJet(iSJ));

    }
    std::sort(jets.begin(),jets.end(),PhysicsUtilities::greaterPT<FatJet>());
}
//--------------------------------------------------------------------------------------------------
float JESUncShifter::correctJet(BaseRecoJet* jet) const {
    float ptSF = 1.0;
    if(cT == CORRTYPE::DOWN) ptSF = 1.0 - jet->jecUnc();
    else if(cT == CORRTYPE::UP) ptSF = 1.0 + jet->jecUnc();
    if(ptSF>0){
        jet->setP4(ptSF*jet->pt(),jet->eta(),jet->phi(),ptSF*jet->mass());
        jet->setRawFactor(jet->toRawFactor()/ptSF);
    } else{
        jet->setP4(0.0f,jet->eta(),jet->phi(),0.0f);
    }
    return std::max(0.0f,ptSF);
}
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
METUncShifter::METUncShifter (const CorrHelp::CORRTYPE cT)
: cT(cT){}
//--------------------------------------------------------------------------------------------------
void METUncShifter::process(Met& met, const EventReader& eventreader) const {
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    //calculate deltaMet so it can be applied pos/neg and to possibly already corrected met
    const Met met_uncUp(ASTypes::CylLorentzVectorF(eventreader.met_unclUp_pt.val(),0,eventreader.met_unclUp_phi.val(),0));
    const Met met_std(ASTypes::CylLorentzVectorF(eventreader.met_pt.val(),0,eventreader.met_phi.val(),0));
    auto deltaM = met_uncUp.p4()-met_std.p4();

    if(cT == CORRTYPE::UP) met.p4() += deltaM;
    else met.p4() -= deltaM;
}
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
}


