
#include "TreeReaders/interface/MuonReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

MuonReader::MuonReader(std::string branchName) : BaseReader("MuonReader",branchName)
{};

MuonReader::~MuonReader(){
    delete pt        ;
    delete eta       ;
    delete phi       ;
    delete q         ;
    delete d0        ;
    delete dz        ;
    delete sip3D     ;
    delete miniIso   ;
    delete dBRelISO  ;
    delete id        ;
    delete dRnorm    ;
    delete PtRatioLepAct;

}

void MuonReader::setup(TreeReadingWrapper * wrapper){
    wrapper->setBranchAddressPre(branchName,"pt"       ,&pt       );
    wrapper->setBranchAddressPre(branchName,"eta"      ,&eta      );
    wrapper->setBranchAddressPre(branchName,"phi"      ,&phi      );
    wrapper->setBranchAddressPre(branchName,"q"        ,&q        );
    wrapper->setBranchAddressPre(branchName,"d0"       ,&d0       );
    wrapper->setBranchAddressPre(branchName,"dz"       ,&dz       );
    wrapper->setBranchAddressPre(branchName,"sip3D"    ,&sip3D    );
    wrapper->setBranchAddressPre(branchName,"miniIso"  ,&miniIso  );
    wrapper->setBranchAddressPre(branchName,"dBRelISO" ,&dBRelISO );
    wrapper->setBranchAddressPre(branchName,"id"       ,&id       );

    wrapper->setBranchAddressPre(branchName,"dRnorm"       ,&dRnorm       );
    wrapper->setBranchAddressPre(branchName,"PtRatioLepAct"       ,&PtRatioLepAct       );

}

void MuonReader::processVars() {
    muons.clear();
    for(unsigned int iO = 0; iO < pt->size(); ++iO){
        muons.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),0),iO,
                q->at(iO),d0->at(iO),dz->at(iO),sip3D->at(iO),miniIso->at(iO),dRnorm->at(iO),PtRatioLepAct->at(iO));
        muons.back().addMuonInfo(dBRelISO->at(iO),id->at(iO));
    }
    std::sort(muons.begin(), muons.end(), PhysicsUtilities::greaterPT<Muon>());
}



}
