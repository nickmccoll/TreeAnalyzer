
#include "TreeReaders/interface/GenParticleReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

GenParticleReader::GenParticleReader(std::string branchName) : BaseReader("GenParticleReader",branchName){};
GenParticleReader::~GenParticleReader(){
    delete pt       ;
    delete eta      ;
    delete phi      ;
    delete mass     ;
    delete status   ;
    delete pdgid    ;
    delete nmoms    ;
    delete firstmom ;
    delete ndaus    ;
    delete firstdau ;
    delete assoc    ;
}


void GenParticleReader::setup(TreeReaderWrapper * wrapper){
//    wrapper->setBranchAddressPre(branchName,"pt"      ,&pt      ,true);
//    wrapper->setBranchAddressPre(branchName,"eta"     ,&eta     ,true);
//    wrapper->setBranchAddressPre(branchName,"phi"     ,&phi     ,true);
//    wrapper->setBranchAddressPre(branchName,"mass"    ,&mass    ,true);
//    wrapper->setBranchAddressPre(branchName,"status"  ,&status  ,true);
//    wrapper->setBranchAddressPre(branchName,"pdgid"   ,&pdgid   ,true);
//    wrapper->setBranchAddressPre(branchName,"nmoms"   ,&nmoms   ,true);
//    wrapper->setBranchAddressPre(branchName,"firstmom",&firstmom,true);
//    wrapper->setBranchAddressPre(branchName,"ndaus"   ,&ndaus   ,true);
//    wrapper->setBranchAddressPre(branchName,"firstdau",&firstdau,true);
//    wrapper->setBranchAddressPre(branchName,"assoc"   ,&assoc   ,true);
}

void GenParticleReader::processVars() {
    genParticles.clear();
    genParticles.reserve(pt->size());

    for(unsigned int iP = 0; iP < pt->size(); ++iP){
      genParticles.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iP),eta->at(iP),phi->at(iP),mass->at(iP)),&genParticles);
      genParticles.back().setStorage(status->at(iP),pdgid->at(iP),nmoms->at(iP),firstmom->at(iP),ndaus->at(iP),firstdau->at(iP),assoc);
    }
}


}
