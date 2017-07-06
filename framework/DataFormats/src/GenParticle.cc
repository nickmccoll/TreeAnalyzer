
#include "DataFormats/interface/GenParticle.h"

namespace TAna {
void GenParticle::setStorage(const ASTypes::size8 instatus, const int inpdgId,
    const stor nMoms, const stor firstMom, const stor nDaus, const stor firstDau, const std::vector<stor> * assocList){
  status_       = instatus    ;
  pdgId_        = inpdgId     ;
  nMoms_        = nMoms       ;
  firstMom_     = firstMom    ;
  nDaus_        = nDaus       ;
  firstDau_     = firstDau    ;
  assocList_    = assocList   ;
}

const GenParticle * GenParticle::mother(const stor idx)const{
  if(assocList_ == 0) throw std::invalid_argument("GenParticle::mother() -> storage not loaded!");
  if(idx >= nMoms_ ) throw std::invalid_argument("GenParticle::mother() -> invalid mother index");
  return &genParticles_->at(assocList_->at(firstMom_ + idx));
}

CandidateRef<GenParticle > GenParticle::motherRef(const stor idx) const {
  return CandidateRef<GenParticle >(mother(idx),assocList_->at(firstMom_ + idx));
}

const GenParticle * GenParticle::daughter(const stor idx)const{
  if(assocList_ == 0) throw std::invalid_argument("GenParticle::daughter() -> storage not loaded!");
  if(idx >= nDaus_ ) throw std::invalid_argument("GenParticle::daughter() -> invalid mother index");
  return &genParticles_->at(assocList_->at(firstDau_ + idx));
}

CandidateRef<GenParticle> GenParticle::daughterRef(const stor idx) const {
  return CandidateRef<GenParticle>(daughter(idx),assocList_->at(firstDau_ + idx));
}

}
