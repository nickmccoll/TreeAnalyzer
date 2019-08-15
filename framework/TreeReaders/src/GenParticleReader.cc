
#include "TreeReaders/interface/GenParticleReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

GenParticleReader::GenParticleReader(std::string branchName) : BaseReader("GenParticleReader",branchName){};
GenParticleReader::~GenParticleReader(){}


void GenParticleReader::setup(TreeReaderWrapper * wrapper){
    wrapper->setBranch(branchName,"pt"      ,pt      ,true);
    wrapper->setBranch(branchName,"eta"     ,eta     ,true);
    wrapper->setBranch(branchName,"phi"     ,phi     ,true);
    wrapper->setBranch(branchName,"mass"    ,mass    ,true);
    wrapper->setBranch(branchName,"status"  ,status  ,true);
    wrapper->setBranch(branchName,"pdgid"   ,pdgid   ,true);
    wrapper->setBranch(branchName,"nmoms"   ,nmoms   ,true);
    wrapper->setBranch(branchName,"firstmom",firstmom,true);
    wrapper->setBranch(branchName,"assoc"   ,assoc   ,true);
    wrapper->setBranch(branchName,"mTTBar"  ,mTTBar  ,true);
}

void GenParticleReader::processVars() {
    genParticles.clear();
    genParticles.reserve(pt.size());

    for(unsigned int iP = 0; iP < pt.size(); ++iP){
      genParticles.emplace_back(ASTypes::CylLorentzVectorF(pt[iP],eta[iP],phi[iP],mass[iP]),status[iP],pdgid[iP],&genParticles);
      for(ASTypes::size16 iM = 0; iM < nmoms[iP]; ++iM  ) genParticles.back().addMother(assoc[ firstmom[iP]+iM ] );
    }

    for(unsigned int iP = 0; iP < genParticles.size(); ++iP){
        for(ASTypes::size16 iM = 0; iM < genParticles[iP].numberOfMothers(); ++iM){
            genParticles[genParticles[iP].motherIndex(iM)].addDaughter(iP);
        }
    }
}


}
