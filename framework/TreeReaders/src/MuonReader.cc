
#include "TreeReaders/interface/MuonReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{
//--------------------------------------------------------------------------------------------------
MuonReader::MuonReader(std::string branchName,bool isRealData) :
        BaseReader("MuonReader",branchName),realData(isRealData)
{};

MuonReader::~MuonReader(){}
//--------------------------------------------------------------------------------------------------
void MuonReader::setup(TreeReaderWrapper * wrapper){
    wrapper->setBranch(branchName,"pt"         ,pt         ,true);
    wrapper->setBranch(branchName,"eta"        ,eta        ,true);
    wrapper->setBranch(branchName,"phi"        ,phi        ,true);
    wrapper->setBranch(branchName,"q"          ,q          ,true);
    wrapper->setBranch(branchName,"id"         ,id         ,true);
    wrapper->setBranch(branchName,"d0"         ,d0         ,true);
    wrapper->setBranch(branchName,"dz"         ,dz         ,true);
    wrapper->setBranch(branchName,"sip3D"      ,sip3D      ,true);
    wrapper->setBranch(branchName,"miniIso"    ,miniIso    ,true);
    wrapper->setBranch(branchName,"dBRelISO"   ,dBRelISO   ,true);
    wrapper->setBranch(branchName,"ptRel"      ,ptRel      ,true);
    wrapper->setBranch(branchName,"ptRatio"    ,ptRatio    ,true);
    wrapper->setBranch(branchName,"dRnorm"     ,dRnorm     ,true);
    wrapper->setBranch(branchName,"lepAct_o_pt",lepAct_o_pt,true);
    if(realData) wrapper->setBranch(branchName,"simType",simType,true);
}
//--------------------------------------------------------------------------------------------------
void MuonReader::processVars() {
    muons.clear();
    for(unsigned int iO = 0; iO < pt.size(); ++iO){
        muons.emplace_back(ASTypes::CylLorentzVectorF(pt[iO],eta[iO],phi[iO],0),iO,
                q[iO],d0[iO],dz[iO],sip3D[iO]);
        muons.back().setIsos(miniIso[iO],dBRelISO[iO],ptRel[iO],ptRatio[iO]);
        muons.back().setSysts(dRnorm[iO],lepAct_o_pt[iO]);
        muons.back().setMuonInfo(id[iO]);
    }
    std::sort(muons.begin(), muons.end(), PhysicsUtilities::greaterPT<Muon>());
}



}
