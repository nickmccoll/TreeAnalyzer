#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
using namespace CutConstants;
using namespace ASTypes;


std::string getSystName (const std::string& prefix, const std::string& proc,const std::string& name, const std::string& sel){
    std::string sn = prefix;
    if(proc.size() )sn += "_"+proc;
    if(name.size() )sn += "_"+name;
    if(sel.size() )sn += "_"+sel;
    return sn;
};

void go(const int insig, const std::string& filename, const std::string& mainDir,  REGION reg, bool simpleSignal) {
    const std::string sigInputDir =  mainDir + (simpleSignal ? "/signalInputsNoCond/" : "/signalInputs/");
    const std::string sfPF =sigInputDir + filename;
    const std::string signalName = signals[insig];
    std::string inputDir = mainDir;

    std::string fPF;
    switch(reg){
    case REG_SR:
        fPF = mainDir + "/bkgInputs/" + filename;
        break;
    case REG_TOPCR:
        fPF = mainDir + "/bkgInputsTopCR/" + filename + "_TopCR";
        break;
    case REG_QGCR:
        fPF = mainDir + "/bkgInputsQGCR/" + filename + "_QGCR";
        break;
    }

    const std::string category = "std";
    std::string cmd = "combineCards.py ";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(p == purCats[PURE_I] ) continue;
        if(h != hadCuts[HAD_FULL] ) continue;


        auto card = DataCardMaker(l,b+"_"+p +"_"+h ,"13TeV",1,category);
//        std::vector<double> newYBins =
//        {700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
//                1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,
//                1500,2000,4000
//        };
//        card.rebinY(newYBins);
//        card.rebinY(13,700,2000);

        const std::string cat = l +"_"+b+"_"+p +"_"+h;
        cmd += std::string(" ")+ category +"_"+cat +"_13TeV=datacard_"+category+"_"+cat +"_13TeV.txt";

        auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& p, const std::string& h,  const std::string& pf) -> std::string
                {return fPF + "_"+proc +"_"+l +"_"+b+"_"+p +"_"+h +"_"+ pf; };
        auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
        auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

        auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
            return getSystName("" + filename,proc,name, sel == "-1" ? cat : sel  );
        };

        //Make search variables
        card.addVar(MOD_MJ,100,0,1000,false);
        card.addVar(MOD_MR,1000,0,10000,false);
        card.addVar(MOD_MS,2000,true);

        //---------------------------------------------------------------------------------------------------
        //Get rates and contributions for backgrounds
        //---------------------------------------------------------------------------------------------------
        card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS,true);
        double rate_lostTW = card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS,true);
        double rate_mw =     card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS,true);
        double rate_mt =     card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS,true);

        //---------------------------------------------------------------------------------------------------
        //Add Systematics first since the param systs need to have the variables added to the workspace
        //---------------------------------------------------------------------------------------------------
        //luminosity
        card.addSystematic("yield","lnN",{{signalName,1.0354}});//lumi = 2.5 jer = 1,  jes = 1 met = 0.5, pdf= 2

        //lepton efficiency
        if(l==lepCats[LEP_E])
            card.addSystematic("eff_"+l,"lnN",{{signalName,1.065}});
        else
            card.addSystematic("eff_"+l,"lnN",{{signalName,1.0576}});
        //
        //tau21
        card.addParamSystematic("tau21_PtDependence",0.0,0.041);
        if(p == purCats[PURE_HP])
            card.addSystematic(systName("","tau21_eff",""),"lnN",{{signalName,1+0.14}});
        if(p == purCats[PURE_LP])
            card.addSystematic(systName("","tau21_eff",""),"lnN",{{signalName,1-0.33}});
        //Btag
        card.addSystematic("btag_fake","lnN",{{signalName,1+0.01}});
        card.addParamSystematic("btag_eff",0.0,0.1);

        //pruned mass scale
        card.addParamSystematic("hh_scale",0.0,0.0122); // jes 1 jer 0.5 met 0.5
        card.addParamSystematic("hh_res",0.0,0.045); // jes 2 jer 4 met 0.5
        card.addParamSystematic("hbb_scale",0.0,0.0094);
        card.addParamSystematic("hbb_res",0.0,0.2);
        //KDE shape systematics
        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"PTX",b) ,0.0,0.5);
        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"OPTX",b),0.0,1.0);
        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"PTY")   ,0.0,1.0);
        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"OPTY")  ,0.0,1.0);
//        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"PT2Y") ,0.0,1.0);
        //top HH resolution and scale
        card.addParamSystematic(systName("top","res"  ) ,0.0,0.20);
        card.addParamSystematic(systName("top","scale") ,0.0,0.25);
        card.addParamSystematic(systName("top","mt_rel_scale",b) ,0.0,0.25);
        card.addParamSystematic(systName("top","lostmw_rel_scale",b) ,0.0,0.25);
        //top Hbb resolution and scale
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"PTX" ,b),0.0,0.3);
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"OPTX",b),0.0,0.6);
        //Normalization
        card.addSystematic(systName(bkgSels[BKG_QG],"norm")  ,"lnN",{{bkgSels[BKG_QG],1.5}});
        card.addSystematic(systName("top","norm")            ,"lnN",{{bkgSels[BKG_LOSTTW],1.25},{bkgSels[BKG_MW],1.25},{bkgSels[BKG_MT],1.25}});
        if(reg == REG_QGCR){
            card.addSystematic(systName("top","tFrac",b)         ,"lnN",{{bkgSels[BKG_MT],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_mt/rate_mw}});
            card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_lostTW/rate_mw}});
        } else {
            card.addSystematic(systName("top","wFrac",b)         ,"lnN",{{bkgSels[BKG_MW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_mw/rate_mt}});
            card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_lostTW/rate_mt}});
        }



        //---------------------------------------------------------------------------------------------------
        //Signal
        //---------------------------------------------------------------------------------------------------
//        if(signalName==signalName){
            //Conditional template
            if(!simpleSignal){
                card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                        {{"hbb_scale",1}},{{"hbb_res",1}},{{"hh_scale",1}},{{"hh_res",1}}, b == btagCats[BTAG_L],MOD_MS);
            }else {
                //Non conditional template
                card.add2DSignalParametricShapeNoCond(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                        {{"hbb_scale",1}},{{"hbb_res",1}},{{"hh_scale",1}},{{"hh_res",1}}, b == btagCats[BTAG_L],MOD_MS);
            }
            std::string tau21Form = "(1.0+tau21_PtDependence*"+ (p == purCats[PURE_HP] ? "log("+MOD_MS+"/1000)" : "((0.054/0.041)*(-log("+MOD_MS+"/1000)))")+")";
            std::string brealForm = "(1.0+btag_eff*";
            if(b == btagCats[BTAG_L]) brealForm+= "(0.07-3.34*10^(-4)*"+MOD_MS+"+5.95*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_M]) brealForm+="(-1.75+0.27*log("+MOD_MS+")))";
            if(b == btagCats[BTAG_T]) brealForm+="(-1.21+0.27*log("+MOD_MS+")))";
            std::string uncForm = tau21Form+"*"+brealForm;
            card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),1,uncForm,{"tau21_PtDependence","btag_eff"}
                            ,MOD_MS);
//        } else throw std::invalid_argument("makeCard::go() -> Bad parsing");

        //---------------------------------------------------------------------------------------------------
        //QG
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts qgKDESysts;
        qgKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_QG],"PTX",b),"1"  }});
        qgKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_QG],"OPTX",b),"1"  }});
        qgKDESysts.addSyst("PTY",{{systName(bkgSels[BKG_QG],"PTY"),"1"  }});
        qgKDESysts.addSyst("OPTY",{{systName(bkgSels[BKG_QG],"OPTY"),"1"  }});
//        qgKDESysts.addSyst("PT2Y",{{systName(bkgSels[BKG_QG],"PT2Y"),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_QG],{MOD_MJ,MOD_MR}, inputName(bkgSels[BKG_QG],"2D_template.root"),"histo",qgKDESysts);

        //---------------------------------------------------------------------------------------------------
        //Lost t/W
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts twKDESysts;
        twKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_LOSTTW],"PTX",b),"1"  }});
        twKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_LOSTTW],"OPTX",b),"1"  }});
        twKDESysts.addSyst("PTY",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b),"1"}});
        twKDESysts.addSyst("OPTY",{{systName("top","res"  ),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_LOSTTW],{MOD_MJ,MOD_MR},inputName(bkgSels[BKG_LOSTTW],"2D_template.root"),"histo",twKDESysts);

        //---------------------------------------------------------------------------------------------------
        //mW
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mwKDESysts;
        mwKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b),"1"}});
        mwKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MW],MOD_MJ,inputName(bkgSels[BKG_MW],"MJJ_SFFit.json"),{{"hbb_scale",1}},{{"hbb_res",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MW],{MOD_MR}, inputName(bkgSels[BKG_MW],"MVV_template.root"),"histo",mwKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MW],bkgSels[BKG_MW] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MW] + "_"+MOD_MR);

        //---------------------------------------------------------------------------------------------------
        //mt
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mtKDESysts;
        mtKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","mt_rel_scale",b),"1"}});
        mtKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MT],MOD_MJ,inputName(bkgSels[BKG_MT],"MJJ_SFFit.json"),{{"hbb_scale",1}},{{"hbb_res",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MT],{MOD_MR}, inputName(bkgSels[BKG_MT],"MVV_template.root"),"histo",mtKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MT],bkgSels[BKG_MT] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MT]+ "_"+MOD_MR);

        //---------------------------------------------------------------------------------------------------
        //Data
        //---------------------------------------------------------------------------------------------------
        card.importBinnedData(fPF + "_data_distributions.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

//        card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

        card.makeCard();
    }
    std::ofstream outFile("comp.sh",std::ios::out|std::ios::trunc);
    outFile << cmd <<" > combinedCard.txt";
    outFile.close();

    std::cout << cmd <<" > combinedCard.txt"<<std::endl;

}
#endif

void makeCard(int inreg = REG_SR, int insig = RADION,    bool condSignal= true){
    std::cout <<" <<<<< "<< inreg <<" "<< condSignal <<" "<<signals[insig]<<std::endl;
    REGION reg = REGION(inreg);
    if(reg == REG_QGCR) btagCats = qgBtagCats;
    std::string mainDir = "../../";
    go(insig,hhFilename,mainDir,reg,!condSignal);
}
