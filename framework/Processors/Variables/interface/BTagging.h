
#ifndef PROCESSORS_VARIABLES_BTAGGING_H
#define PROCESSORS_VARIABLES_BTAGGING_H




namespace TAna {

namespace BTagging{
enum  CSVWP { CSV_INCL, CSV_L, CSV_M, CSV_T};
const float CSVWP_VALS[] = {-100,0.5426,0.8484,0.9535};

template<typename Jet>
bool isLooseCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_L];}

template<typename Jet>
bool isMediumCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_M];}

template<typename Jet>
bool isTightCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_T];}

template<typename Jet>
bool csvTagged(const Jet& jet, CSVWP wp) { return jet.csv() >=  CSVWP_VALS[wp];}


enum  BBTWP { BBT_INCL, BBT_L, BBT_M1, BBT_M2, BBT_T};
const float BBT_VALS[] = {-100,0.3,0.6,0.8,0.9};
}


}

#endif

