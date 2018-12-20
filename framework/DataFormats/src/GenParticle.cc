
#include "DataFormats/interface/GenParticle.h"
#include <TTreeReaderArray.h>

namespace TAna {
GenParticle::stor GenParticle::motherIndex(const stor mNum)const{
    if(mNum >= numberOfMothers() ) throw std::invalid_argument("GenParticle::mother() -> invalid mother index");
    return momIdxs[mNum];
}

const GenParticle * GenParticle::mother(const stor mNum)const{ return &(*genParticles_)[motherIndex(mNum)];}

CandidateRef<GenParticle > GenParticle::motherRef(const stor mNum) const {
  return CandidateRef<GenParticle >(mother(mNum),momIdxs[mNum]);
}

GenParticle::stor GenParticle::daughterIndex(const stor dNum)const{
    if(dNum >= numberOfDaughters() ) throw std::invalid_argument("GenParticle::mother() -> invalid daughter index");
    return dauIdxs[dNum];
}

const GenParticle * GenParticle::daughter(const stor dNum)const{return &(*genParticles_)[daughterIndex(dNum)];}

CandidateRef<GenParticle> GenParticle::daughterRef(const stor dNum) const {
  return CandidateRef<GenParticle>(daughter(dNum),dauIdxs[dNum]);
}

}
