#ifndef TREEANALYZER_DATAFORMATS_GENPARTICLE_H
#define TREEANALYZER_DATAFORMATS_GENPARTICLE_H

#include "DataFormats/interface/Momentum.h"

template <typename T> class TTreeReaderArray;

namespace TAna {
template <class Type> class CandidateRef;

//--------------------------------------------------------------------------------------------------
// GenParticle
//--------------------------------------------------------------------------------------------------
class GenParticle : public MomentumF
{
public :
  typedef ASTypes::size16 stor;

  GenParticle() : status_(-1), pdgId_(-1),genParticles_(0) {}

  template <class InputCoordSystem>
  GenParticle(ROOT::Math::LorentzVector<InputCoordSystem> inMomentum, const ASTypes::size8 instatus=0, const int inpdgId = 0, std::vector<GenParticle> * genParticles = 0) : MomentumF(inMomentum),
  status_(instatus), pdgId_(inpdgId),genParticles_(genParticles) {}
  ~GenParticle(){}

  void setGenPrtPtr(std::vector<GenParticle> * genParticles) {genParticles_ = genParticles;}
  void addDaughter(const stor idx) {dauIdxs.push_back(idx); }
  void addMother(const stor idx) {momIdxs.push_back(idx); }

  int status()      const { return status_;}
  int pdgId()       const { return pdgId_; }
  int absPdgId()    const { return std::abs(pdgId_); }

  stor numberOfMothers() const {return momIdxs.size();}
  stor motherIndex(const stor mNum)const;
  const GenParticle * mother(const stor mNum)const;
  CandidateRef<GenParticle > motherRef(const stor mNum) const;

  stor numberOfDaughters() const {return dauIdxs.size();}
  stor daughterIndex(const stor dNum)const;
  const GenParticle * daughter(const stor dNum)const;
  CandidateRef<GenParticle> daughterRef(const stor dNum) const;


  //Dummy functions used to add compatibility with reco::GenParticles
  int charge() const {return 0;}


protected :
  ASTypes::size8  status_  ;
  int    pdgId_   ;
  const std::vector<GenParticle> * genParticles_;
  std::vector<stor> momIdxs;
  std::vector<stor> dauIdxs;




};
typedef std::vector<GenParticle>     GenParticleCollection;

//--------------------------------------------------------------------------------------------------
// CandidateRef...dummy class to add compatibility with CMSSW GenParticle
//--------------------------------------------------------------------------------------------------
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
