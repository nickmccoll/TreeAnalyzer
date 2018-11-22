
#include "TreeReaders/interface/PhotonReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

PhotonReader::PhotonReader(std::string branchName) : BaseReader("PhotonReader",branchName)
{};

PhotonReader::~PhotonReader(){
    delete pt        ;
    delete eta       ;
    delete phi       ;
    delete hadOvEm   ;
}

void PhotonReader::setup(TreeReaderWrapper * wrapper){
//    wrapper->setBranchAddressPre(branchName,"pt"        ,&pt         );
//    wrapper->setBranchAddressPre(branchName,"eta"       ,&eta        );
//    wrapper->setBranchAddressPre(branchName,"phi"       ,&phi        );
//    wrapper->setBranchAddressPre(branchName,"hadOvEm"   ,&hadOvEm    );
//    wrapper->setBranchAddressPre(branchName,"hadTowOvEm",&hadTowOvEm    );
}

void PhotonReader::processVars() {
    photons.clear();
    for(unsigned int iO = 0; iO < pt->size(); ++iO){
        photons.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),0),iO);
    }
    std::sort(photons.begin(), photons.end(), PhysicsUtilities::greaterPT<Photon>());
}
}
