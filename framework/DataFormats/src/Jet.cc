
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {
bool Jet::passPUID()    const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_PU);}
bool Jet::passLooseID() const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_LOOSE);}
bool Jet::passTightID() const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_TIGHT);}

}
