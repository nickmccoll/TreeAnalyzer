#ifndef TREEANALYZER_DATAFORMATS_GENPARTICLE_H
#define TREEANALYZER_DATAFORMATS_GENPARTICLE_H

#include "DataFormats/interface/Momentum.h"

namespace TAna {
template <class Type> class CandidateRef;

class GenParticle : public MomentumF
{
public :
  typedef ASTypes::size16 stor;

  GenParticle() : status_(-1), pdgId_(-1), nMoms_(0), firstMom_(0),nDaus_(0),firstDau_(0),assocList_(0),genParticles_(0) {}

  template <class InputCoordSystem>
  GenParticle(ROOT::Math::LorentzVector<InputCoordSystem> inMomentum,std::vector<GenParticle> * genParticles = 0) : MomentumF(inMomentum),
  status_(-1), pdgId_(-1), nMoms_(0), firstMom_(0),nDaus_(0),firstDau_(0),assocList_(0),genParticles_(genParticles) {}
  ~GenParticle(){}

  void setGenPrtPtr(std::vector<GenParticle> * genParticles) {genParticles_ = genParticles;}
  void setStorage(const ASTypes::size8 instatus, const int inpdgId,
      const stor nMoms, const stor firstMom, const stor nDaus, const stor firstDau, const std::vector<stor> * assocList);

  int status()      const { return status_;}
  int pdgId()       const { return pdgId_; }
  int absPdgId()    const { return std::abs(pdgId_); }

  stor numberOfMothers() const {return nMoms_;}
  const GenParticle * mother(const stor idx)const;
  CandidateRef<GenParticle > motherRef(const stor idx) const;

  stor numberOfDaughters() const {return nDaus_;}
  const GenParticle * daughter(const stor idx)const;
  CandidateRef<GenParticle> daughterRef(const stor idx) const;


  //Dummy functions used to add compatibility with reco::GenParticles
  int charge() const {return 0;}


protected :
  ASTypes::size8  status_  ;
  int    pdgId_   ;
  stor   nMoms_   ;
  stor   firstMom_;
  stor   nDaus_   ;
  stor   firstDau_;
  const std::vector<stor> * assocList_;
  const std::vector<GenParticle> * genParticles_;


};
typedef std::vector<GenParticle>     GenParticleCollection;


//dummy class to add compatibility with GenParticle
template <class Type>
class CandidateRef {
public:
  CandidateRef() : ptr(0), idx(0) {};
  CandidateRef(const Type * inptr, const unsigned int inidx) : ptr(inptr), idx(inidx) {}

  void set(const Type * inptr, const unsigned int inidx) {ptr = inptr; idx = inidx;}
  const Type* operator->() const {return  ptr;}
  const Type& operator* () const {return *ptr;}
  unsigned int key() const {return idx;}
  bool isNull() const {return ptr == 0;}

private:
  const Type * ptr;
  unsigned int idx;
};

typedef  CandidateRef<GenParticle> GenParticleRef;

}
#endif
