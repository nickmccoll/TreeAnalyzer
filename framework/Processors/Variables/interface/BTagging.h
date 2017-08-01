
#ifndef PROCESSORS_VARIABLES_BTAGGING_H
#define PROCESSORS_VARIABLES_BTAGGING_H




namespace TAna {

namespace BTagging{
enum  CSVWP { CSV_INCL, CSV_L, CSV_M, CSV_T};
const float CSVWP_VALS[] = {-100,0.5426,0.8484,0.9535};
enum  BBTWP { BBT_INCL, BBT_L, BBT_M1, BBT_M2, BBT_T};
const float BBT_VALS[] = {-100,0.3,0.6,0.8,0.9};
}


}

#endif

