#if !defined(__CINT__) || defined(__MAKECINT__)
#include "DataCardMaker.h"
using namespace CutConstants;
using namespace ASTypes;

void go(const std::string& signalName, const std::string& filename) {
    const std::string inputDir = "../inputs/";
    const std::string fPF =inputDir+filename;
    const std::string category = "std";
    std::string cmd = "combineCards.py ";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l == lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]||p == purSels[PUR_LMT]) continue;
        if(h != hadSels[HAD_FULL] ) continue;

        auto card = DataCardMaker(l,p,"13TeV",1,category);
        const std::string cat =  l+"_"+p +"_13TeV";
        cmd += std::string(" ")+ cat +"=datacard_"+category+"_"+cat +".txt";

        //Make search variables
        card.addVar(MOD_MJ,100,0,1000,false);
        card.addVar(MOD_MR,1000,0,10000,false);
        card.addVar(MOD_MS,2000,true);

        //Systematic vars
        card.addVar("CMS_scale_prunedj",0,-.1,1,false);
        card.addVar("CMS_res_prunedj",0,-.5,0.5,false);
        card.addVar("CMS_scale_j",0,-.1,1,false);
        card.addVar("CMS_scale_MET",0,-.1,1,false);
        card.addVar("CMS_res_j",0,-.5,0.5,false);
        card.addVar("CMS_res_MET",0,-.5,0.5,false);

        if(signalName=="radHH"){
            card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, fPF + "_"+signalName +"_"+l+"_"+p+"_"+h +"_2D_fit.json",
                    {{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},{{"CMS_scale_j",1},{"CMS_scale_MET",1}},{{"CMS_res_j",1},{"CMS_res_MET",1}}, p == purSels[PUR_L],MOD_MS);
            card.addParametricYieldWithUncertainty(signalName,0,fPF + "_"+signalName +"_"+l+"_"+p+"_"+h +"_yield.json",1,"CMS_tau21_PtDependence","log("+MOD_MS+"/600)",0.041,MOD_MS);
        } else throw std::invalid_argument("makeCard::go() -> Bad parsing");

        //qg bkg
        auto qgSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_"+n+"_"+l+"_"+p ));};
        card.addHistoShapeFromFile(bkgSels[BKG_QG],{MOD_MJ,MOD_MR}, fPF + "_"+bkgSels[BKG_QG] +"_"+l+"_"+p+"_"+h +"_2D_template.root","histo",{qgSyst("PTX"),qgSyst("OPTX"),qgSyst("PTY"),qgSyst("OPTY")});
        card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+l+"_"+p+"_"+h+"_"+hhMCS);
        //lost t/w bkg
        auto lostTWSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_"+n+"_"+l+"_"+p ));};
        card.addHistoShapeFromFile(bkgSels[BKG_LOSTTW],{MOD_MJ,MOD_MR}, fPF + "_"+bkgSels[BKG_LOSTTW] +"_"+l+"_"+p+"_"+h +"_2D_template.root","histo",{lostTWSyst("PTX"),lostTWSyst("OPTX"),lostTWSyst("PTY"),lostTWSyst("OPTY")});
        card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+l+"_"+p+"_"+h+"_"+hhMCS);

        //mW
        card.add1DBKGParametricShape(bkgSels[BKG_MW],MOD_MJ,fPF+"_"+bkgSels[BKG_MW]+"_"+lepSels[LEP_EMU]+"_"+purSels[PUR_LMT]+"_"+hadSels[HAD_NONE]+"_fit.json",{{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},MOD_MR,MOD_MJ);
        auto lostMWSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_MW]+"_"+n+"_"+l+"_"+p ));};
        card.addHistoShapeFromFile(bkgSels[BKG_MW],{MOD_MR}, fPF + "_"+bkgSels[BKG_MW] +"_"+l+"_"+p+"_"+h +"_template.root","histo",{lostMWSyst("PT"),lostMWSyst("OPT")},false,0,MOD_MR);
        card.conditionalProduct(bkgSels[BKG_MW],bkgSels[BKG_MW] + "_"+MOD_MJ,MOD_MR,bkgSels[BKG_MW] + "_"+MOD_MR);
        card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+l+"_"+p+"_"+h+"_"+hhMCS);

        //mT
        card.add1DBKGParametricShape(bkgSels[BKG_MT],MOD_MJ,fPF+"_"+bkgSels[BKG_MT]+"_"+lepSels[LEP_EMU]+"_"+p+"_"+hadSels[HAD_NONE]+"_fit.json",{{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},MOD_MR,MOD_MJ);
        auto lostMTSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_MT]+"_"+n+"_"+l+"_"+p ));};
        card.addHistoShapeFromFile(bkgSels[BKG_MT],{MOD_MR}, fPF + "_"+bkgSels[BKG_MT] +"_"+l+"_"+p+"_"+h +"_template.root","histo",{lostMTSyst("PT"),lostMTSyst("OPT")},false,0,MOD_MR);
        card.conditionalProduct(bkgSels[BKG_MT],bkgSels[BKG_MT] + "_"+MOD_MJ,MOD_MR,bkgSels[BKG_MT]+ "_"+MOD_MR);
        card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+l+"_"+p+"_"+h+"_"+hhMCS);

        //Data
//        card.importBinnedData(filename + "_data.root","data_"+l+"_"+p+"_"+h,{MOD_MJ,MOD_MR});
        card.importBinnedData(fPF + "_pd.root","data_"+l+"_"+p+"_"+h+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

        //Systematics
        //luminosity
        card.addSystematic("CMS_lumi","lnN",{StrFlt("radHH",1.026)});
        //kPDF uncertainty for the signal
        card.addSystematic("CMS_pdf","lnN",{StrFlt("radHH",1.01)});
        //
        //lepton efficiency
        card.addSystematic("CMS_eff_"+l,"lnN",{StrFlt("radHH",1.1)});
        //
        //W+jets cross section in acceptance-dominated by pruned mass
        card.addSystematic("CMS_"+filename + bkgSels[BKG_QG]+"_norm_"+l+"_"+p,"lnN",{StrFlt(bkgSels[BKG_QG],1.5)});
        card.addSystematic("CMS_"+filename + bkgSels[BKG_LOSTTW]+"_norm_"+l+"_"+p,"lnN",{StrFlt(bkgSels[BKG_LOSTTW],1.5)});
        card.addSystematic("CMS_"+filename + bkgSels[BKG_MW]+"_norm_"+l+"_"+p,"lnN",{StrFlt(bkgSels[BKG_MW],1.5)});
        card.addSystematic("CMS_"+filename + bkgSels[BKG_MT]+"_norm_"+l+"_"+p,"lnN",{StrFlt(bkgSels[BKG_MT],1.5)});
        //
        //tau21
        card.addSystematic("CMS_"+filename+"_tau21_eff","lnN",{StrFlt("radHH",1+0.14)});
        //Btag
        card.addSystematic("CMS_btag_fake","lnN",{StrFlt("radHH",1+0.02)});
        if(p == purSels[PUR_L])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1-0.03)});
        if(p == purSels[PUR_M])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1+0.03)});
        if(p == purSels[PUR_T])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1+0.06)});

        //pruned mass scale
        card.addParamSystematic("CMS_scale_j",0.0,0.02);
        card.addParamSystematic("CMS_res_j",0.0,0.05);
        card.addParamSystematic("CMS_scale_prunedj",0.0,0.0094);
        card.addParamSystematic("CMS_res_prunedj",0.0,0.2);
        card.addParamSystematic("CMS_scale_MET",0.0,0.02);
        card.addParamSystematic("CMS_res_MET",0.0,0.01);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_PTX_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_OPTX_"+l+"_"+p,0.0,0.6);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_PTY_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_OPTY_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_PTX_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_OPTX_"+l+"_"+p,0.0,0.6);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_PTY_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_OPTY_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_MW]+"_PT_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_MW]+"_OPT_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_MT]+"_PT_"+l+"_"+p,0.0,0.333);
        card.addParamSystematic("CMS_"+filename+"_"+bkgSels[BKG_MT]+"_OPT_"+l+"_"+p,0.0,0.333);

        card.makeCard();
    }
    std::cout << cmd<<std::endl;

}
#endif

void makeCard(std::string signalString = "radHH"){
    go(signalString,hhFilename);
}
