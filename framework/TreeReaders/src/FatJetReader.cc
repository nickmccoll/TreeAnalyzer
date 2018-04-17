
#include "TreeReaders/interface/FatJetReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

FatJetReader::FatJetReader(std::string branchName, bool isRealData, bool fillGenFatJets) : BaseReader("FatJetReader",branchName),
        realData(isRealData),fillGenFatJets(fillGenFatJets)
{};

FatJetReader::~FatJetReader(){
    delete pt          ;
    delete eta         ;
    delete phi         ;
    delete mass        ;
    delete toRawFact   ;
    delete csv         ;
    delete id          ;
    delete hadronFlavor;
    delete partonFlavor;
    delete JECUnc      ;
    delete genIDX      ;
    delete gen_pt      ;
    delete gen_eta     ;
    delete gen_phi     ;
    delete gen_mass    ;

    delete bbt               ;
    delete tau1              ;
    delete tau2              ;
    delete tau3              ;
    delete sj1_pt            ;
    delete sj1_eta           ;
    delete sj1_phi           ;
    delete sj1_mass          ;
    delete sj1_toRawFact     ;
    delete sj1_csv           ;
    delete sj1_JECUnc        ;
    delete sj1_hadronFlavor  ;
    delete sj1_partonFlavor  ;
    delete sj2_pt            ;
    delete sj2_eta           ;
    delete sj2_phi           ;
    delete sj2_mass          ;
    delete sj2_toRawFact     ;
    delete sj2_csv           ;
    delete sj2_JECUnc        ;
    delete sj2_hadronFlavor  ;
    delete sj2_partonFlavor  ;

}

void FatJetReader::setup(TreeReadingWrapper * wrapper){
    wrapper->setBranchAddressPre(branchName,"pt"          ,&pt          ,true);
    wrapper->setBranchAddressPre(branchName,"eta"         ,&eta         ,true);
    wrapper->setBranchAddressPre(branchName,"phi"         ,&phi         ,true);
    wrapper->setBranchAddressPre(branchName,"mass"        ,&mass        ,true);
    wrapper->setBranchAddressPre(branchName,"toRawFact"   ,&toRawFact   ,true);
    wrapper->setBranchAddressPre(branchName,"csv"         ,&csv         ,true);
    wrapper->setBranchAddressPre(branchName,"id"          ,&id          ,true);

    if(!realData) {
        wrapper->setBranchAddressPre(branchName,"hadronFlavor",&hadronFlavor,true);
        wrapper->setBranchAddressPre(branchName,"partonFlavor",&partonFlavor,true);
        wrapper->setBranchAddressPre(branchName,"JECUnc"      ,&JECUnc,true);
    }

    if(fillGenFatJets && !realData) {
        wrapper->setBranchAddressPre(branchName,"genIDX"      ,&genIDX      ,true);
        wrapper->setBranchAddressPre(branchName,"gen_pt"      ,&gen_pt      ,true);
        wrapper->setBranchAddressPre(branchName,"gen_eta"     ,&gen_eta     ,true);
        wrapper->setBranchAddressPre(branchName,"gen_phi"     ,&gen_phi     ,true);
        wrapper->setBranchAddressPre(branchName,"gen_mass"    ,&gen_mass    ,true);
    }

    wrapper->setBranchAddressPre(branchName,"bbt"              , &bbt              ,true);
    wrapper->setBranchAddressPre(branchName,"tau1"             , &tau1             ,true);
    wrapper->setBranchAddressPre(branchName,"tau2"             , &tau2             ,true);
    wrapper->setBranchAddressPre(branchName,"tau3"             , &tau3             ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_pt"           , &sj1_pt           ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_eta"          , &sj1_eta          ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_phi"          , &sj1_phi          ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_mass"         , &sj1_mass         ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_toRawFact"    , &sj1_toRawFact       ,true);
    wrapper->setBranchAddressPre(branchName,"sj1_csv"          , &sj1_csv          ,true);
    if(!realData) {
        wrapper->setBranchAddressPre(branchName,"sj1_JECUnc"   , &sj1_JECUnc ,true);
        wrapper->setBranchAddressPre(branchName,"sj1_hadronFlavor" , &sj1_hadronFlavor ,true);
        wrapper->setBranchAddressPre(branchName,"sj1_partonFlavor" , &sj1_partonFlavor ,true);
    }
    wrapper->setBranchAddressPre(branchName,"sj2_pt"           , &sj2_pt           ,true);
    wrapper->setBranchAddressPre(branchName,"sj2_eta"          , &sj2_eta          ,true);
    wrapper->setBranchAddressPre(branchName,"sj2_phi"          , &sj2_phi          ,true);
    wrapper->setBranchAddressPre(branchName,"sj2_mass"         , &sj2_mass         ,true);
    wrapper->setBranchAddressPre(branchName,"sj2_toRawFact"    , &sj2_toRawFact       ,true);
    wrapper->setBranchAddressPre(branchName,"sj2_csv"          , &sj2_csv          ,true);
    if(!realData) {
        wrapper->setBranchAddressPre(branchName,"sj2JECUnc"        , &sj2_JECUnc ,true);
        wrapper->setBranchAddressPre(branchName,"sj2_hadronFlavor" , &sj2_hadronFlavor ,true);
        wrapper->setBranchAddressPre(branchName,"sj2_partonFlavor" , &sj2_partonFlavor ,true);
    }




}

void FatJetReader::processVars() {
    jets.clear();
    genJets.clear();

    std::vector<GenFatJet*> genInd(gen_pt->size(),0);

    if(!realData && fillGenFatJets){
        for(unsigned int iO = 0; iO < gen_pt->size(); ++iO){
            genJets.emplace_back(ASTypes::CylLorentzVectorF(gen_pt->at(iO),gen_eta->at(iO),gen_phi->at(iO),gen_mass->at(iO)),iO);
        }
        std::sort(genJets.begin(), genJets.end(), PhysicsUtilities::greaterPT<GenFatJet>());
        for(unsigned int iO = 0; iO < genJets.size(); ++iO){
            genInd[genJets[iO].index()] = &genJets[iO];
        }
    }

    for(unsigned int iO = 0; iO < pt->size(); ++iO){
         ASTypes::int8 hadronFlv =  0;
         ASTypes::int8 partonFlv =  0;
         float         jecUnc    =  0;
        if(!realData){
            hadronFlv = hadronFlavor->at(iO);
            partonFlv = partonFlavor->at(iO);
            jecUnc    = JECUnc->at(iO);
        }


        GenFatJet *gj = 0;
        if(fillGenFatJets && !realData && genIDX->at(iO) != 255 ){
            gj = genInd[genIDX->at(iO)];
        }

        jets.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),mass->at(iO)),iO,
                toRawFact->at(iO),csv->at(iO),id->at(iO), hadronFlv,partonFlv,jecUnc,gj);
        jets.back().addFatJetInfo(bbt->at(iO),tau1->at(iO),tau2->at(iO),tau3->at(iO));

        auto addSubJet =[&](const std::vector<float>* pt,const std::vector<float>* eta,const std::vector<float>* phi,const std::vector<float>* mass, int iO,
                const std::vector<float>* toRaw,const std::vector<float>* csv, std::vector<ASTypes::int8> * hadronFlavor,
                std::vector<ASTypes::int8> * partonFlavor,const std::vector<float>* JECUnc) {
            if(pt->at(iO) == 0) return;
            ASTypes::int8 hadronFlv =  0;
            ASTypes::int8 partonFlv =  0;
            float         jecUnc    =  0;
           if(!realData){
               hadronFlv = hadronFlavor->at(iO);
               partonFlv = partonFlavor->at(iO);
               jecUnc    = JECUnc->at(iO);
           }
            SubJet jet(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),mass->at(iO)),-1,
                    toRaw->at(iO), csv->at(iO),hadronFlv,partonFlv,jecUnc);
            jets.back().addSubJet(jet);
        };
        addSubJet(sj1_pt,sj1_eta,sj1_phi,sj1_mass, iO,sj1_toRawFact,sj1_csv,sj1_hadronFlavor,sj1_partonFlavor,sj1_JECUnc);
        addSubJet(sj2_pt,sj2_eta,sj2_phi,sj2_mass, iO,sj2_toRawFact,sj2_csv,sj2_hadronFlavor,sj2_partonFlavor,sj2_JECUnc);

    }
    std::sort(jets.begin(), jets.end(), PhysicsUtilities::greaterPT<FatJet>());

}


}
