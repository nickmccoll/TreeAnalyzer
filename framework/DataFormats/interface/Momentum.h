//--------------------------------------------------------------------------------------------------
// 
// Momemntum.h
// 
// Momentum class, used as base class for custom objects.
// 
// 
//--------------------------------------------------------------------------------------------------

#ifndef TREEANALYZER_DATAFORMATS_MOMENTUM_H
#define TREEANALYZER_DATAFORMATS_MOMENTUM_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include<vector>

namespace TAna {

template <class CoordSystem>
class Momentum
{

public :
  Momentum() {}

  template <class InputCoordSystem>
  Momentum(ROOT::Math::LorentzVector<InputCoordSystem> inMomentum) : fMom(inMomentum) {}

  ~Momentum(){}

  // Functions to facilitate momentum operations
  float   pt()      const { return fMom.Pt();       }
  float   eta()     const { return fMom.Eta();      }
  float   absEta()  const { return std::fabs(fMom.Eta());}
  float   phi()     const { return fMom.Phi();      }
  float   mass()    const { return fMom.M();        }
  float   E()       const { return fMom.E();        }
  float   energy()  const { return fMom.E();        }
  float   Et()      const { return fMom.Et();       }
  float   mt()      const { return fMom.Mt();       }
  float   px()      const { return fMom.Px();       }
  float   py()      const { return fMom.Py();       }
  float   pz()      const { return fMom.Pz();       }
  float   p()       const { return fMom.P();        }
  float   y()       const { return fMom.Rapidity(); }
  float   theta()   const { return fMom.Theta();    }

  //Momentum getting and setting functions
  ROOT::Math::LorentzVector<CoordSystem>&        p4()         { return fMom; }
  const  ROOT::Math::LorentzVector<CoordSystem>& p4()   const { return fMom; }

  template< class Coords >
  void setP4(const ROOT::Math::LorentzVector<Coords> & v )    { fMom = v;    }
  template< class Coords >
  void setP4(const Coords pt, const Coords eta, const Coords phi, const Coords mass )
  { fMom = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Coords> >(pt,eta,phi,mass) ;    }

  //cout the momentum
  friend std::ostream& operator<<(std::ostream& os, const Momentum<CoordSystem>& m){
    os << "("<<m.pt()<<","<<m.eta()<<","<<m.phi()<<","<<m.mass()<<")";
    return os;
  }

private:
  ROOT::Math::LorentzVector<CoordSystem>	 fMom;

};

  typedef Momentum<ASTypes::CylLorentzCoord>  MomentumD;
  typedef Momentum<ASTypes::CylLorentzCoordF> MomentumF;

  typedef std::vector<MomentumD> MomentumDCollection;
  typedef std::vector<MomentumF> MomentumFCollection;

}

#endif
