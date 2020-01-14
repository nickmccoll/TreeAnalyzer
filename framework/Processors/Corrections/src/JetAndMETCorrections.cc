
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"



namespace TAna {
using namespace CorrHelp;
//--------------------------------------------------------------------------------------------------
JERCorrector::JERCorrector (const std::string& dataDir,std::shared_ptr<TRandom3> rndGen, const CorrHelp::CORRTYPE cT )
: dataDir(dataDir), cT(cT),rndGen(rndGen){

}
//--------------------------------------------------------------------------------------------------
void JERCorrector::setParameters(const JetParameters& param) {

    ak8Puppi_resObj  .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK8Puppi_resFile));
    ak8Puppi_sfObj   .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK8Puppi_sfFile ));
    ak4CHS_resObj    .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK4CHS_resFile  ));
    ak4CHS_sfObj     .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK4CHS_sfFile   ));
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
        if(j.pt() > 7000) continue;
        if(j.pt() < 15) {
//            std::cout <<std::endl;
            continue;
        }
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);
        auto oRawFact = j.toRawFactor();

        const float jres = ak4CHS_resObj->evaluateFormula(*ak4CHS_resObj
                ->getRecord(parameters),parameters);

        const float jSF =  ak4CHS_sfObj->getRecord(parameters)
                 ->getParametersValues()[getSFCount(cT)];

        auto ptSF = correctJet(&j,gjptrs,jres,jSF,0.4);

        deltaMX += ((jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
        deltaMY += ((jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
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
        if(j.pt() > 7000) continue;
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);

        const float jres = ak8Puppi_resObj->evaluateFormula(*ak8Puppi_resObj
                ->getRecord(parameters),parameters);
        const float jSF =  ak8Puppi_sfObj->getRecord(parameters)
                ->getParametersValues()[getSFCount(cT)];

        correctJet(&j,gjptrs,jres,jSF,0.8);
    }
    std::sort(jets.begin(),jets.end(),PhysicsUtilities::greaterPT<FatJet>());
}

//--------------------------------------------------------------------------------------------------
int JERCorrector::getSFCount(const CORRTYPE c) {
    switch(c){
    case NOMINAL:
        return 0;
    case DOWN:
        return 1;
    case UP:
        return 2;
    default:
        throw std::invalid_argument(std::string("JERCorrector cannot be used with NONE"));
    }
};
//--------------------------------------------------------------------------------------------------
float JERCorrector::correctJet(Jet* jet, const std::vector<const GenJet*> genjets, const float jres,
        const float resSF, const float coneR){

    double nearDR = 0;
    int idx = PhysicsUtilities::findNearestDRDeref(*jet,genjets,nearDR);
//    std::cout <<"(";
//    if(idx >= 0) std::cout <<genjets[idx]->pt()<<",";
//    else std::cout << "-1,";
//    std::cout << jres <<","<<resSF<<",";
    float ptSF = 1.0;
    if(idx>= 0 && nearDR < coneR/2.0 && std::fabs(jet->pt()-genjets[idx]->pt()) <3.0*jres*jet->pt() ){
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
            deltaMX += ((jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
            deltaMY += ((jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
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
HEM1516TestCorrector::HEM1516TestCorrector (const CorrHelp::CORRTYPE cT) : cT(cT){}
//--------------------------------------------------------------------------------------------------
void HEM1516TestCorrector::processJets(JetReader& jetreader, Met& met) {
    if(cT == CORRTYPE::NONE) return;
    float deltaPX = 0, deltaPY = 0;
    for(auto& j : jetreader.jets) {
    	if(j.pt() < 15) continue;
    	if(!j.passTightID()) continue;

    	// instructions given in email on Oct 20 2019
    	if(j.phi() > -1.57 && j.phi() < -0.87) {
    		float sf = 1.0;
    		if(j.eta() > -2.5 && j.eta() < -1.3) sf = 0.8;
    		else if (j.eta() > -3.0 && j.eta() < -2.5) sf = 0.65;
    		else continue;

    		const ASTypes::CartLorentzVector newMom(sf*j.px(),sf*j.py(),sf*j.pz(),sf*j.E());

    		deltaPX += (1-sf)*j.px();
    		deltaPY += (1-sf)*j.py();
    		j.setP4(newMom);
    	}
    }

    ASTypes::CartLorentzVector deltaMet(deltaPX,deltaPY,0.0f,std::sqrt(deltaPX*deltaPX + deltaPY*deltaPY));
    met.p4() += deltaMet;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void HEM1516TestCorrector::processFatJets(FatJetCollection& fatjets) {
    if(cT == CORRTYPE::NONE) return;

    for(auto& fj : fatjets) {
    	// * should be scaling down all fatjets before considering subjets
    	if(!fj.nSubJets()) continue;

    	float origSdEnergy = fj.sdMom().E();
    	for(unsigned int idx=0; idx < fj.nSubJets(); ++idx) {
    		if(fj.subJet(idx).pt() < 15) continue;
    		// should I do anything to ensure there are no leptons in this subjet??

    		if(fj.subJet(idx).phi() > -1.57 && fj.subJet(idx).phi() < -0.87) {
    			float sf = 1.0;
    			if(fj.subJet(idx).eta() > -1.3) continue;
    			if(fj.subJet(idx).eta() > -2.5 && fj.subJet(idx).eta() < -1.3) sf = 0.8;
    			else if (fj.subJet(idx).eta() > -3.0 && fj.subJet(idx).eta() < -2.5) sf = 0.65;

        		ASTypes::CartLorentzVector newMom(sf*fj.subJet(idx).px(),sf*fj.subJet(idx).py(),
        				sf*fj.subJet(idx).pz(),sf*fj.subJet(idx).E());
        		fj.subJet(idx).setP4(newMom);
    		}
    	}
    	float scaleFactor = fj.sdMom().E() / origSdEnergy;
    	if(scaleFactor < 1.0) {
    		ASTypes::CartLorentzVector newMom(scaleFactor*fj.px(),scaleFactor*fj.py(),scaleFactor*fj.pz(),scaleFactor*fj.E());
    		fj.setP4(newMom);
    	}
    }
    std::sort(fatjets.begin(),fatjets.end(),PhysicsUtilities::greaterPT<FatJet>());
}

}


