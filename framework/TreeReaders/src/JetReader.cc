
#include "TreeReaders/interface/JetReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

JetReader::JetReader(std::string branchName, bool isRealData, bool fillGenJets,  bool fillRecoJets) : BaseReader("JetReader",branchName),
        realData(isRealData),fillGenJets(fillGenJets),fillRecoJets(fillRecoJets)
{};

JetReader::~JetReader(){
    delete pt          ;
    delete eta         ;
    delete phi         ;
    delete mass        ;
    delete chef        ;
    delete toRawFact   ;
    delete metUnc_rawPx;
    delete metUnc_rawPy;
    delete csv         ;
    delete id          ;
    delete hadronFlavor;
    delete partonFlavor;
    delete JECUnc;
    delete genIDX      ;
    delete gen_pt      ;
    delete gen_eta     ;
    delete gen_phi     ;
    delete gen_mass    ;

}

void JetReader::setup(TreeReaderWrapper * wrapper){
//    if(fillRecoJets){
//        wrapper->setBranchAddressPre(branchName,"pt"          ,&pt          ,true);
//        wrapper->setBranchAddressPre(branchName,"eta"         ,&eta         ,true);
//        wrapper->setBranchAddressPre(branchName,"phi"         ,&phi         ,true);
//        wrapper->setBranchAddressPre(branchName,"mass"        ,&mass        ,true);
//        wrapper->setBranchAddressPre(branchName,"chef"        ,&chef        ,true);
//        wrapper->setBranchAddressPre(branchName,"toRawFact"   ,&toRawFact        ,true);
//        wrapper->setBranchAddressPre(branchName,"metUnc_rawPx",&metUnc_rawPx        ,true);
//        wrapper->setBranchAddressPre(branchName,"metUnc_rawPy",&metUnc_rawPy        ,true);
//        wrapper->setBranchAddressPre(branchName,"csv"         ,&csv         ,true);
//        wrapper->setBranchAddressPre(branchName,"id"          ,&id          ,true);
//
//        if(!realData){
//            wrapper->setBranchAddressPre(branchName,"hadronFlavor",&hadronFlavor,true);
//            wrapper->setBranchAddressPre(branchName,"partonFlavor",&partonFlavor,true);
//            wrapper->setBranchAddressPre(branchName,"JECUnc"      ,&JECUnc,true);
//        }
//        if(fillGenJets && !realData) {
//            wrapper->setBranchAddressPre(branchName,"genIDX"      ,&genIDX      ,true);
//        }
//    }
//
//
//    if(fillGenJets && !realData) {
//        wrapper->setBranchAddressPre(branchName,"gen_pt"      ,&gen_pt      ,true);
//        wrapper->setBranchAddressPre(branchName,"gen_eta"     ,&gen_eta     ,true);
//        wrapper->setBranchAddressPre(branchName,"gen_phi"     ,&gen_phi     ,true);
//        wrapper->setBranchAddressPre(branchName,"gen_mass"    ,&gen_mass    ,true);
//    }


}

void JetReader::processVars() {
    jets.clear();
    genJets.clear();

    std::vector<GenJet*> genInd(gen_pt->size(),0);

    if(!realData && fillGenJets){
        for(unsigned int iO = 0; iO < gen_pt->size(); ++iO){
            genJets.emplace_back(ASTypes::CylLorentzVectorF(gen_pt->at(iO),gen_eta->at(iO),gen_phi->at(iO),gen_mass->at(iO)),iO);
        }
        std::sort(genJets.begin(), genJets.end(), PhysicsUtilities::greaterPT<GenJet>());
        for(unsigned int iO = 0; iO < genJets.size(); ++iO){
            genInd[genJets[iO].index()] = &genJets[iO];
        }
    }

    if(fillRecoJets){
        for(unsigned int iO = 0; iO < pt->size(); ++iO){
            ASTypes::int8 hadronFlv =  0;
            ASTypes::int8 partonFlv =  0;
            float         jecUnc    =  0;
           if(!realData){
               hadronFlv = hadronFlavor->at(iO);
               partonFlv = partonFlavor->at(iO);
               jecUnc    = JECUnc->at(iO);
           }
            GenJet *gj = 0;
            if(fillGenJets && !realData && genIDX->at(iO) != 255 ){
                gj = genInd[genIDX->at(iO)];
            }

            jets.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),mass->at(iO)),iO,
                    toRawFact->at(iO),csv->at(iO),id->at(iO), hadronFlv,partonFlv,jecUnc,gj);
        }
        std::sort(jets.begin(), jets.end(), PhysicsUtilities::greaterPT<Jet>());
    }
}


}
