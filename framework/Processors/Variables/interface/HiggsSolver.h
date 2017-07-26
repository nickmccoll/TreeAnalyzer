
#ifndef PROCESSORS_VARIABLES_HIGGSSOLVER_H
#define PROCESSORS_VARIABLES_HIGGSSOLVER_H


#include "DataFormats/interface/Momentum.h"


namespace TAna {
namespace HiggsSolver {
MomentumF getInvisible(const MomentumF& met, const MomentumF& vis, const double hMass = 125);
}



}


#endif

