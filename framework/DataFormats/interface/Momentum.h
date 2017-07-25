

#ifndef TREEANALYZER_DATAFORMATS_MOMENTUM_H
#define TREEANALYZER_DATAFORMATS_MOMENTUM_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include<vector>
using ASTypes::size;

namespace TAna {


//--------------------------------------------------------------------------------------------------
// Momentum: base class for object that has a momentum
//--------------------------------------------------------------------------------------------------

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
    float   rap()     const { return fMom.Rapidity(); }
    float   theta()   const { return fMom.Theta();    }
    float   x()       const { return fMom.x();        }
    float   y()       const { return fMom.y();        }
    float   z()       const { return fMom.z();        }

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

//--------------------------------------------------------------------------------------------------
// IndexedMomentum: Base class for object that has a momentum and an index (e.g. in a tree)
//--------------------------------------------------------------------------------------------------

template <class CoordSystem>
class IndexedMomentum : public Momentum<CoordSystem>
{
public :
    IndexedMomentum() : _idx(-1) {}

    template <class InputCoordSystem>
    IndexedMomentum(const ROOT::Math::LorentzVector<InputCoordSystem>& mom, const int idx = -1) : Momentum<CoordSystem>(mom), _idx(idx) {}

    ~IndexedMomentum(){}

    int   index()             const { return _idx;  }
    void  setIndex(const int& newIndex)   { _idx = newIndex;    }


protected :
    int   _idx;  //Index in Jet vector
};

typedef IndexedMomentum<ASTypes::CylLorentzCoord>  IndexedMomentumD;
typedef IndexedMomentum<ASTypes::CylLorentzCoordF> IndexedMomentumF;

typedef std::vector<IndexedMomentumD> IndexedMomentumDCollection;
typedef std::vector<IndexedMomentumF> IndexedMomentumFCollection;

}

#endif
