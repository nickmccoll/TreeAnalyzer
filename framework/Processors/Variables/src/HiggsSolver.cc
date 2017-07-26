#include "Processors/Variables/interface/HiggsSolver.h"

namespace TAna {

MomentumF HiggsSolver::getInvisible(const MomentumF& met, const MomentumF& vis, const double hMass){
    const double a = hMass*hMass - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
    const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
    const double B = -4*a* vis.z();
    const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
    const double delta = B*B -4*A*C;

    double metZ = 0;
    if(delta < 0) {
        metZ= -B/(2*A);
    } else {
        const double pos = (-B + std::sqrt(delta))/(2*A);
        const double neg = (-B - std::sqrt(delta))/(2*A);
        if(std::fabs(pos) < std::fabs(neg)) metZ = pos;
        else metZ = neg;
    }
    ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
    return MomentumF(neutrino);
}

}
