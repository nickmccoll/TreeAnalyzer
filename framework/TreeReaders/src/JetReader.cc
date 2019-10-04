
#include "TreeReaders/interface/JetReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{
//--------------------------------------------------------------------------------------------------
JetReader::JetReader(std::string branchName, bool isRealData, bool fillGenJets,
        bool fillRecoJets,bool fillBTagging) : BaseReader("JetReader",branchName),
        realData(isRealData),fillGenJets(fillGenJets),fillRecoJets(fillRecoJets),
        fillBTagging(fillBTagging)
{};

JetReader::~JetReader(){}
//--------------------------------------------------------------------------------------------------
void JetReader::setup(TreeReaderWrapper * wrapper){
    if(fillRecoJets){
        wrapper->setBranch(branchName,"pt"          ,pt          ,true);
        wrapper->setBranch(branchName,"eta"         ,eta         ,true);
        wrapper->setBranch(branchName,"phi"         ,phi         ,true);
        wrapper->setBranch(branchName,"mass"        ,mass        ,true);
        wrapper->setBranch(branchName,"toRawFact"   ,toRawFact   ,true);
        wrapper->setBranch(branchName,"metUnc_rawPx",metUnc_rawPx,true);
        wrapper->setBranch(branchName,"metUnc_rawPy",metUnc_rawPy,true);
        if(fillBTagging){
            wrapper->setBranch(branchName,"csv"         ,csv         ,true);
            wrapper->setBranch(branchName,"deep_csv"    ,deep_csv    ,true);
            wrapper->setBranch(branchName,"deep_flavor"    ,deep_flavor    ,true);
        }
        wrapper->setBranch(branchName,"id"          ,id          ,true);
        if(!realData){
            wrapper->setBranch(branchName,"hadronFlavor",hadronFlavor,true);
            wrapper->setBranch(branchName,"partonFlavor",partonFlavor,true);
            wrapper->setBranch(branchName,"JECUnc"      ,JECUnc      ,true);
            if(fillGenJets)wrapper->setBranch(branchName,"genIDX"  ,genIDX  ,true);
        }
    }
    if(!realData&&fillGenJets){
        wrapper->setBranch(branchName,"gen_pt"  ,gen_pt  ,true);
        wrapper->setBranch(branchName,"gen_eta" ,gen_eta ,true);
        wrapper->setBranch(branchName,"gen_phi" ,gen_phi ,true);
        wrapper->setBranch(branchName,"gen_mass",gen_mass,true);

    }
}
//--------------------------------------------------------------------------------------------------
void JetReader::processVars() {
    jets.clear();
    genJets.clear();

    std::vector<GenJet*> genInd;

    if(!realData && fillGenJets){
        genInd.resize(gen_pt.size());
        for(unsigned int iO = 0; iO < gen_pt.size(); ++iO){
            genJets.emplace_back(
                    ASTypes::CylLorentzVectorF(gen_pt[iO],gen_eta[iO],gen_phi[iO],gen_mass[iO]),iO);
        }
        std::sort(genJets.begin(), genJets.end(), PhysicsUtilities::greaterPT<GenJet>());
        for(unsigned int iO = 0; iO < genJets.size(); ++iO){
            genInd[genJets[iO].index()] = &genJets[iO];
        }
    }

    if(fillRecoJets){
        for(unsigned int iO = 0; iO < pt.size(); ++iO){

            jets.emplace_back(ASTypes::CylLorentzVectorF(pt[iO],eta[iO],phi[iO],mass[iO]),iO,
                    toRawFact[iO],id[iO]);
            if(fillBTagging)
                jets.back().addBTagging(deep_csv[iO],csv[iO],deep_flavor[iO]);

            if(!realData){
                GenJet *gj = fillGenJets && genIDX[iO] != 255 ?
                        genInd[genIDX[iO]] : 0;
                jets.back().addMCInfo(hadronFlavor[iO],partonFlavor[iO],JECUnc[iO],gj);
            }
        }
        std::sort(jets.begin(), jets.end(), PhysicsUtilities::greaterPT<Jet>());
    }
}


}
