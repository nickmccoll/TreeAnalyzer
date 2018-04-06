#include "../predTools/CutConstants.h"
#include <vector>
#include "TFile.h"
#include "HistoPlotting/include/Plotter.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;



void plotCRBkgTests(int step = 0){
    hhFilename +="_CR";
    std::string filename = hhFilename;
//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000},true,true);
//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{30,210},false,true);

//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_LMT_I_full","mu_LMT_I_full","emu_L_I_full","emu_M_I_full","emu_T_I_full","emu_LMT_LP_full","emu_LMT_HP_full"},{700,4000},true,true);

//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_LMT_I_full","mu_LMT_I_full","emu_L_I_full","emu_M_I_full","emu_T_I_full","emu_LMT_LP_full","emu_LMT_HP_full"},{30,80,100,160,210},false,true);
//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"emu_LMT_I_full","emu_L_I_full","emu_M_I_full","emu_T_I_full"},{30,80,100,160,210},false,true);
    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"mu_LMT_I_full","mu_L_I_full"},{30,80,100,160,210},false,true);
//    test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_LMT_I_full","mu_LMT_I_full"},{30,80,100,160,210},false,true);
}
