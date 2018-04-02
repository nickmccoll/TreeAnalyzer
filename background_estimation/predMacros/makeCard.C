#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
using namespace CutConstants;
using namespace ASTypes;

void go(const std::string& signalName, const std::string& filename) {
    const std::string inputDir = "../inputs/";
    const std::string fPF =inputDir+filename;
    const std::string category = "std";
    std::string cmd = "combineCards.py ";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_LMT]|| b==btagCats[BTAG_I] ) continue;
        if(p == purCats[PURE_I] ) continue;
        if(h != hadCuts[HAD_FULL] ) continue;

        auto card = DataCardMaker(l,b+"_"+p +"_"+h ,"13TeV",1,category);

        const std::string cat = l +"_"+b+"_"+p +"_"+h;
        cmd += std::string(" ")+ cat +"_13TeV=datacard_"+category+"_"+cat +"_13TeV.txt";

        auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& p, const std::string& h,  const std::string& pf) -> std::string
                {return fPF + "_"+proc +"_"+l +"_"+b+"_"+p +"_"+h +"_"+ pf; };
        auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };


        auto qgSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_QG]+"_"+n+"_"+cat ));};
        auto lostTWSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_LOSTTW]+"_"+n+"_"+cat ));};
        auto mwSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_MW]+"_"+n+"_"+cat));};
        auto mtSyst =[&](const std::string& n)->StrStr{return StrStr(n,std::string("CMS_"+filename+"_"+bkgSels[BKG_MT]+"_"+n+"_"+cat ));};


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
            card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, inputName(signalName,"2D_fit.json"),
                    {{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},{{"CMS_scale_j",1},{"CMS_scale_MET",1}},{{"CMS_res_j",1},{"CMS_res_MET",1}}, b == btagCats[BTAG_L],MOD_MS);
//            card.add2DSignalParametricShapeNoCond(signalName,MOD_MJ,MOD_MR, inputName(signalName,"2D_fit.json"),
//                    {{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},{{"CMS_scale_j",1},{"CMS_scale_MET",1}},{{"CMS_res_j",1},{"CMS_res_MET",1}}, b == btagCats[BTAG_L],MOD_MS);
            card.addParametricYieldWithUncertainty(signalName,0,inputName(signalName,"yield.json"),1,"CMS_tau21_PtDependence","log("+MOD_MS+"/600)",0.041,MOD_MS);
        } else throw std::invalid_argument("makeCard::go() -> Bad parsing");

        //qg bkg
        card.addHistoShapeFromFile(bkgSels[BKG_QG],{MOD_MJ,MOD_MR}, inputName(bkgSels[BKG_QG],"2D_template.root"),"histo",{qgSyst("PTX"),qgSyst("OPTX"),qgSyst("PTY"),qgSyst("OPTY")});
        card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS);
        //lost t/w bkg
        card.addHistoShapeFromFile(bkgSels[BKG_LOSTTW],{MOD_MJ,MOD_MR},inputName(bkgSels[BKG_LOSTTW],"2D_template.root"),"histo",{lostTWSyst("PTX"),lostTWSyst("OPTX"),lostTWSyst("PTY"),lostTWSyst("OPTY")});
        card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS);

        //mW
        card.add1DBKGParametricShape(bkgSels[BKG_MW],MOD_MJ,fullInputName(bkgSels[BKG_MW],lepCats[LEP_EMU],btagCats[BTAG_LMT],purCats[PURE_I],hadCuts[HAD_NONE],"fit.json"),{{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MW],{MOD_MR}, inputName(bkgSels[BKG_MW],"template.root"),"histo",{mwSyst("PT"),mwSyst("OPT")},false,0,MOD_MR);
        card.conditionalProduct(bkgSels[BKG_MW],bkgSels[BKG_MW] + "_"+MOD_MJ,MOD_MR,bkgSels[BKG_MW] + "_"+MOD_MR);
        card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS);

        //mT
        card.add1DBKGParametricShape(bkgSels[BKG_MT],MOD_MJ,fullInputName(bkgSels[BKG_MT],lepCats[LEP_EMU],b,purCats[PURE_I],hadCuts[HAD_NONE],"fit.json"),{{"CMS_scale_prunedj",1}},{{"CMS_res_prunedj",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MT],{MOD_MR}, inputName(bkgSels[BKG_MT],"template.root"),"histo",{mtSyst("PT"),mtSyst("OPT")},false,0,MOD_MR);
        card.conditionalProduct(bkgSels[BKG_MT],bkgSels[BKG_MT] + "_"+MOD_MJ,MOD_MR,bkgSels[BKG_MT]+ "_"+MOD_MR);
        card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS);

        //Data
//        card.importBinnedData(filename + "_data.root","data_"+l+"_"+p+"_"+h,{MOD_MJ,MOD_MR});
        card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

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
        card.addSystematic("CMS_"+filename +"_"+ bkgSels[BKG_QG]+"_norm_"+cat,"lnN",{StrFlt(bkgSels[BKG_QG],1.5)});
        card.addSystematic("CMS_"+filename +"_"+ bkgSels[BKG_LOSTTW]+"_norm_"+cat,"lnN",{StrFlt(bkgSels[BKG_LOSTTW],1.5)});
        card.addSystematic("CMS_"+filename +"_"+ bkgSels[BKG_MW]+"_norm_"+cat,"lnN",{StrFlt(bkgSels[BKG_MW],1.5)});
        card.addSystematic("CMS_"+filename +"_"+ bkgSels[BKG_MT]+"_norm_"+cat,"lnN",{StrFlt(bkgSels[BKG_MT],1.5)});
        //
        //tau21
        if(p == purCats[PURE_HP])
            card.addSystematic("CMS_"+filename+"_tau21_eff","lnN",{StrFlt("radHH",1+0.14)});
        if(p == purCats[PURE_LP])
            card.addSystematic("CMS_"+filename+"_tau21_eff","lnN",{StrFlt("radHH",1-0.33)});
        //Btag
        card.addSystematic("CMS_btag_fake","lnN",{StrFlt("radHH",1+0.02)});
        if(p == btagCats[BTAG_L])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1-0.03)});
        if(p == btagCats[BTAG_M])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1+0.03)});
        if(p == btagCats[BTAG_T])
            card.addSystematic("CMS_btag_eff" ,"lnN",{StrFlt("radHH",1+0.06)});

        //pruned mass scale
        card.addParamSystematic("CMS_scale_j",0.0,0.02);
        card.addParamSystematic("CMS_res_j",0.0,0.05);
        card.addParamSystematic("CMS_scale_prunedj",0.0,0.0094);
        card.addParamSystematic("CMS_res_prunedj",0.0,0.2);
        card.addParamSystematic("CMS_scale_MET",0.0,0.02);
        card.addParamSystematic("CMS_res_MET",0.0,0.01);
        card.addParamSystematic(qgSyst("PTX") .second,0.0,0.333);
        card.addParamSystematic(qgSyst("OPTX").second,0.0,0.6);
        card.addParamSystematic(qgSyst("PTY") .second,0.0,0.333);
        card.addParamSystematic(qgSyst("OPTY").second,0.0,0.333);
        card.addParamSystematic(lostTWSyst("PTX") .second,0.0,0.333);
        card.addParamSystematic(lostTWSyst("OPTX").second,0.0,0.6);
        card.addParamSystematic(lostTWSyst("PTY") .second,0.0,0.333);
        card.addParamSystematic(lostTWSyst("OPTY").second,0.0,0.333);
        card.addParamSystematic(mwSyst("PT") .second,0.0,0.333);
        card.addParamSystematic(mwSyst("OPT") .second,0.0,0.333);
        card.addParamSystematic(mtSyst("PT") .second,0.0,0.333);
        card.addParamSystematic(mtSyst("OPT") .second,0.0,0.333);

        card.makeCard();
    }
    std::cout << cmd <<" > combinedCard.txt"<<std::endl;

}
#endif

void makeCard(std::string signalString = "radHH"){
    go(signalString,hhFilename);
}
