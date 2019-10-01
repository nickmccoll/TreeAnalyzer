
#include "TreeReaders/interface/ElectronReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{
//--------------------------------------------------------------------------------------------------
ElectronReader::ElectronReader(std::string branchName, bool storeIDVars) :
        BaseReader("ElectronReader",branchName),storeIDVars(storeIDVars)
{};

ElectronReader::~ElectronReader(){}
//--------------------------------------------------------------------------------------------------
void ElectronReader::setup(TreeReaderWrapper * wrapper){
    wrapper->setBranch(branchName,"pt"          ,pt         ,true);
    wrapper->setBranch(branchName,"eta"         ,eta        ,true);
    wrapper->setBranch(branchName,"phi"         ,phi        ,true);
    wrapper->setBranch(branchName,"q"           ,q          ,true);
    wrapper->setBranch(branchName,"scEta"       ,scEta      ,true);
    wrapper->setBranch(branchName,"scE"         ,scE        ,true);
    wrapper->setBranch(branchName,"uncorPt"     ,uncorPt    ,true);
    wrapper->setBranch(branchName,"id"          ,id         ,true);
    wrapper->setBranch(branchName,"d0"          ,d0         ,true);
    wrapper->setBranch(branchName,"dz"          ,dz         ,true);
    wrapper->setBranch(branchName,"sip3D"       ,sip3D      ,true);
    wrapper->setBranch(branchName,"mvaID"       ,mvaID      ,true);
    wrapper->setBranch(branchName,"miniIso"     ,miniIso    ,true);
    wrapper->setBranch(branchName,"eaRelIso"    ,eaRelIso   ,true);
    wrapper->setBranch(branchName,"dRnorm"      ,dRnorm     ,true);
    wrapper->setBranch(branchName,"lepAct_o_pt" ,lepAct_o_pt,true);
    wrapper->setBranch(branchName,"sc_act_o_pt" ,sc_act_o_pt,true);
    wrapper->setBranch(branchName,"sc_dr_act"   ,sc_dr_act  ,true);
    if(storeIDVars){
        wrapper->setBranch(branchName,"passMedCutBased"      ,passMedCutBased      ,false);
        wrapper->setBranch(branchName,"passTightCutBased"    ,passTightCutBased    ,false);
        wrapper->setBranch(branchName,"full5x5_sigmaIetaIeta",full5x5_sigmaIetaIeta,false);
        wrapper->setBranch(branchName,"abs_dEtaSeed"         ,abs_dEtaSeed         ,false);
        wrapper->setBranch(branchName,"abs_dPhiIn"           ,abs_dPhiIn           ,false);
        wrapper->setBranch(branchName,"HoE"                  ,HoE                  ,false);
        wrapper->setBranch(branchName,"HoE_BC"               ,HoE_BC               ,false);
        wrapper->setBranch(branchName,"abs_1oEm1op"          ,abs_1oEm1op          ,false);
        wrapper->setBranch(branchName,"missInnerHits"        ,missInnerHits        ,false);
        wrapper->setBranch(branchName,"passConvVeto"         ,passConvVeto         ,false);
        wrapper->setBranch(branchName,"seeding"              ,seeding              ,false);
        wrapper->setBranch(branchName,"nSaturatedXtals"      ,nSaturatedXtals      ,false);
        wrapper->setBranch(branchName,"e2x5OverE5x5"         ,e2x5OverE5x5         ,false);
        wrapper->setBranch(branchName,"e1x5OverE5x5"         ,e1x5OverE5x5         ,false);
        wrapper->setBranch(branchName,"isolEmHadDepth1"      ,isolEmHadDepth1      ,false);
    }
}


//--------------------------------------------------------------------------------------------------
void ElectronReader::processVars() {
    electrons.clear();
    for(unsigned int iO = 0; iO < pt.size(); ++iO){
        electrons.emplace_back(ASTypes::CylLorentzVectorF(pt[iO],eta[iO],phi[iO],0),iO,
                q[iO],d0[iO],dz[iO],sip3D[iO]);
        electrons.back().setIsos(miniIso[iO],eaRelIso[iO],0,0);
        electrons.back().setSysts(dRnorm[iO],lepAct_o_pt[iO]);
        electrons.back().addElectronInfo(scEta[iO],scE[iO],mvaID[iO],id[iO],
                 sc_act_o_pt[iO],sc_dr_act[iO] );
    }
    std::sort(electrons.begin(), electrons.end(), PhysicsUtilities::greaterPT<Electron>());
}



}
