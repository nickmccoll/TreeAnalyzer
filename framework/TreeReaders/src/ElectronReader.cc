
#include "TreeReaders/interface/ElectronReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

ElectronReader::ElectronReader(std::string branchName) : BaseReader("ElectronReader",branchName)
{};

ElectronReader::~ElectronReader(){
    delete pt        ;
    delete eta       ;
    delete phi       ;
    delete q         ;
    delete scEta     ;
    delete d0        ;
    delete dz        ;
    delete sip3D     ;
    delete mvaID     ;
    delete mvaID_cat ;
    delete miniIso   ;
    delete eaRelISO  ;
    delete id        ;

}

void ElectronReader::setup(TreeReadingWrapper * wrapper){
    wrapper->setBranchAddressPre(branchName,"pt"       ,&pt       );
    wrapper->setBranchAddressPre(branchName,"eta"      ,&eta      );
    wrapper->setBranchAddressPre(branchName,"phi"      ,&phi      );
    wrapper->setBranchAddressPre(branchName,"q"        ,&q        );
    wrapper->setBranchAddressPre(branchName,"scEta"    ,&scEta    );
    wrapper->setBranchAddressPre(branchName,"d0"       ,&d0       );
    wrapper->setBranchAddressPre(branchName,"dz"       ,&dz       );
    wrapper->setBranchAddressPre(branchName,"sip3D"    ,&sip3D    );
    wrapper->setBranchAddressPre(branchName,"mvaID"    ,&mvaID    );
    wrapper->setBranchAddressPre(branchName,"mvaID_cat",&mvaID_cat);
    wrapper->setBranchAddressPre(branchName,"miniIso"  ,&miniIso  );
    wrapper->setBranchAddressPre(branchName,"eaRelISO" ,&eaRelISO );
    wrapper->setBranchAddressPre(branchName,"id"       ,&id       );
}

void ElectronReader::processVars() {
    electrons.clear();
    for(unsigned int iO = 0; iO < pt->size(); ++iO){
        electrons.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),0),iO,
                q->at(iO),d0->at(iO),dz->at(iO),sip3D->at(iO),miniIso->at(iO));
        electrons.back().addElectronInfo(scEta->at(iO),mvaID->at(iO),eaRelISO->at(iO),id->at(iO));
    }
    std::sort(electrons.begin(), electrons.end(), PhysicsUtilities::greaterPT<Electron>());
}



}
