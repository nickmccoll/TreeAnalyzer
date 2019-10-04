
#include "TreeReaders/interface/FatJetReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{
//--------------------------------------------------------------------------------------------------
FatJetReader::FatJetReader(std::string branchName, bool isRealData,
        bool fillGenFatJets,bool fillBTagging, bool fillWTagging, bool fillLSFInfo) :
                BaseReader("FatJetReader",branchName),
                realData(isRealData),fillGenFatJets(fillGenFatJets),fillBTagging(fillBTagging),
                fillWTagging(fillWTagging),fillLSFInfo(fillLSFInfo)
{};

FatJetReader::~FatJetReader(){}

//--------------------------------------------------------------------------------------------------
void FatJetReader::setup(TreeReaderWrapper * wrapper){
    wrapper->setBranch(branchName,"pt"       ,pt       ,true);
    wrapper->setBranch(branchName,"eta"      ,eta      ,true);
    wrapper->setBranch(branchName,"phi"      ,phi      ,true);
    wrapper->setBranch(branchName,"mass"     ,mass     ,true);
    wrapper->setBranch(branchName,"toRawFact",toRawFact,true);
    wrapper->setBranch(branchName,"id"       ,id       ,true);
    if(fillBTagging){
        wrapper->setBranch(branchName,"bbt",bbt,true);
        wrapper->setBranch(branchName,"deep_MDZHbb",deep_MDZHbb,true);
        wrapper->setBranch(branchName,"deep_MDHbb" ,deep_MDHbb,true);
        wrapper->setBranch(branchName,"deep_Hbb"   ,deep_Hbb,true);
    }
    if(fillWTagging){
        wrapper->setBranch(branchName,"deep_W"   ,deep_W,true);
    }
    wrapper->setBranch(branchName,"tau1"     ,tau1     ,true);
    wrapper->setBranch(branchName,"tau2"     ,tau2     ,true);
    wrapper->setBranch(branchName,"sdMass_z0p15",sdMass_z0p15     ,true);
    wrapper->setBranch(branchName,"sdMass_z0p05",sdMass_z0p05     ,true);
    if(fillLSFInfo){
        wrapper->setBranch(branchName,"ecfN2"      ,ecfN2      ,true);
        wrapper->setBranch(branchName,"ecfM2"      ,ecfM2      ,true);
        wrapper->setBranch(branchName,"ecfD2"      ,ecfD2      ,true);
        wrapper->setBranch(branchName,"ecfN3"      ,ecfN3      ,true);
        wrapper->setBranch(branchName,"ecfU3"      ,ecfU3      ,true);
        wrapper->setBranch(branchName,"ecfU2"      ,ecfU2      ,true);
        wrapper->setBranch(branchName,"tau3"       ,tau3       ,true);
        wrapper->setBranch(branchName,"tau4"       ,tau4       ,true);
        wrapper->setBranch(branchName,"lsf3"       ,lsf3       ,true);
        wrapper->setBranch(branchName,"dRLep"      ,dRLep      ,true);
        wrapper->setBranch(branchName,"lepInJetMVA",lepInJetMVA,true);
    }

    if(!realData){
        wrapper->setBranch(branchName,"hadronFlavor",hadronFlavor,true);
        wrapper->setBranch(branchName,"partonFlavor",partonFlavor,true);
        wrapper->setBranch(branchName,"JECUnc"      ,JECUnc      ,true);
        if(fillGenFatJets){
            wrapper->setBranch(branchName,"genIDX"  ,genIDX   ,true);
            wrapper->setBranch(branchName,"gen_pt"  ,gen_pt   ,true);
            wrapper->setBranch(branchName,"gen_eta" ,gen_eta  ,true);
            wrapper->setBranch(branchName,"gen_phi" ,gen_phi  ,true);
            wrapper->setBranch(branchName,"gen_mass",gen_mass ,true);
        }
    }

    wrapper->setBranch(branchName,"sjIDX1"      ,sjIDX1      ,true);
    wrapper->setBranch(branchName,"sjnum"       ,sjnum       ,true);
    wrapper->setBranch(branchName,"sj_pt"       ,sj_pt       ,true);
    wrapper->setBranch(branchName,"sj_eta"      ,sj_eta      ,true);
    wrapper->setBranch(branchName,"sj_phi"      ,sj_phi      ,true);
    wrapper->setBranch(branchName,"sj_mass"     ,sj_mass     ,true);
    wrapper->setBranch(branchName,"sj_toRawFact",sj_toRawFact,true);
    if(fillBTagging){
        wrapper->setBranch(branchName,"sj_deep_csv",sj_deep_csv,true);
    }
    if(!realData){
        wrapper->setBranch(branchName,"sj_hadronFlavor",sj_hadronFlavor,true);
        wrapper->setBranch(branchName,"sj_partonFlavor",sj_partonFlavor,true);
        wrapper->setBranch(branchName,"sj_JECUnc"      ,sj_JECUnc      ,true);
    }



}

//--------------------------------------------------------------------------------------------------
void FatJetReader::processVars() {
    jets.clear();
    genJets.clear();

    std::vector<GenFatJet*> genInd;
    if(!realData && fillGenFatJets){
        genInd.resize(gen_pt.size());
        for(unsigned int iO = 0; iO < gen_pt.size(); ++iO){
            genJets.emplace_back(
                    ASTypes::CylLorentzVectorF(gen_pt[iO],gen_eta[iO],gen_phi[iO],gen_mass[iO]),iO);
        }
        std::sort(genJets.begin(), genJets.end(), PhysicsUtilities::greaterPT<GenFatJet>());
        for(unsigned int iO = 0; iO < genJets.size(); ++iO){
            genInd[genJets[iO].index()] = &genJets[iO];
        }
    }

    for(unsigned int iO = 0; iO < pt.size(); ++iO){
        jets.emplace_back(ASTypes::CylLorentzVectorF(pt[iO],eta[iO],phi[iO],mass[iO]),iO,
                toRawFact[iO],id[iO]);
        if(fillBTagging)
            jets.back().addFJBtagging(bbt[iO],deep_MDZHbb[iO],deep_MDHbb[iO],deep_Hbb[iO]);
        if(fillWTagging)
            jets.back().addWTaging(deep_W[iO]);
        jets.back().addSubStructure(tau1[iO],tau2[iO],sdMass_z0p15[iO],sdMass_z0p05[iO]);


        if(!realData){
            GenJet *gj = fillGenFatJets && genIDX[iO] != 255 ?
                    genInd[genIDX[iO]] : 0;
            jets.back().addMCInfo(hadronFlavor[iO],partonFlavor[iO],JECUnc[iO],gj);
        }

        for(unsigned int iSJ = sjIDX1[iO]; iSJ < sjIDX1[iO]+sjnum[iO]; ++iSJ ){
            SubJet subJet(
                    ASTypes::CylLorentzVectorF(sj_pt[iSJ],sj_eta[iSJ],sj_phi[iSJ],sj_mass[iSJ])
                    ,iSJ,sj_toRawFact[iSJ]);
            if(fillBTagging)
                subJet.addBTagging(sj_deep_csv[iSJ]);
            if(!realData)
                subJet.addMCInfo(sj_hadronFlavor[iSJ],sj_partonFlavor[iSJ],sj_JECUnc[iSJ]);
            jets.back().addSubJet(subJet);
        }
    }
    std::sort(jets.begin(), jets.end(), PhysicsUtilities::greaterPT<FatJet>());

}


}
