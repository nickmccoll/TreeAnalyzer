#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
using namespace CutConstants;
using namespace ASTypes;


std::string getSystName (const std::string& prefix, const std::string& proc,const std::string& name, const std::string& sel){
    std::string sn = prefix;

    if(proc.size() ){
        if(sn.size()) sn+="_";
        sn += proc;
    }
    if(name.size() ){
        if(sn.size()) sn+="_";
        sn += name;
    }
    if(sel.size() ){
        if(sn.size()) sn+="_";
        sn += sel;
    }
    return sn;
};

void go(const int insig, const std::string& filename, const std::string& mainDir, REGION reg, int channel){
    const std::string sigInputDir =  mainDir +  "/signalInputs/";
    const std::string sfPF = sigInputDir + filename;
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
    case REG_NONTOPCR:
        fPF = mainDir + "/bkgInputsNonTopCR/" + filename + "_NonTopCR";
        break;
    }

    const std::string category = "std";
    std::string cmd = "combineCards.py ";

    if(channel == 0 || channel == 1) {
        for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats) for(const auto& h :hadCuts){
            if(l == lepCats[LEP_EMU] ) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(p == purCats[PURE_I] ) continue;
            if(h != hadCuts[HAD_FULL] ) continue;


            auto card = DataCardMaker(l,b+"_"+p +"_"+h ,"13TeV",1,category);

            const std::string cat = l +"_"+b+"_"+p +"_"+h;
            cmd += std::string(" ")+ category +"_"+cat +"_13TeV=datacard_"+category+"_"+cat +"_13TeV.txt";

            auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& p, const std::string& h, const std::string& pf) -> std::string
                    {return fPF + "_"+proc +"_"+l +"_"+b+"_"+p +"_"+h +"_"+ pf; };
            auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
            auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

            auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
                return getSystName("",proc,name, sel == "-1" ? l +"_"+b+"_"+p : sel  );
            };

            //Make search variables
            card.addVar(MOD_MJ,100,0,1000,false);
            card.addVar(MOD_MR,1000,0,10000,false);
            card.addVar(MOD_MS,2000,true);

            //---------------------------------------------------------------------------------------------------
            //Get rates and contributions for backgrounds
            //---------------------------------------------------------------------------------------------------
            card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_lostTW = card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_mw =     card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_mt =     card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS,1.0,true);

            //---------------------------------------------------------------------------------------------------
            //Add Systematics first since the param systs need to have the variables added to the workspace
            //---------------------------------------------------------------------------------------------------
            //luminosity
            card.addSystematic("yield","lnN",{{signalName,1.0391}});//lumi = 2.5 pdf= 2, PU = 0.5, btagfake=1

            //lepton efficiency
            if(l==lepCats[LEP_E])
                card.addSystematic("eff_"+l,"lnN",{{signalName,1.0602}}); //2% trigger / 5.5% for reco  / 1.4% ID / ISO 0.2%
            else
                card.addSystematic("eff_"+l,"lnN",{{signalName,1.0566}}); //2% trigger / ID 1%  /  ISO 5.2%
            //
            //tau21
            card.addParamSystematic("tau21_PtDependence",0.0,0.041);
            if(p == purCats[PURE_HP])
                card.addSystematic(systName("","tau21_eff",""),"lnN",{{signalName,1+0.14}});
            if(p == purCats[PURE_LP])
                card.addSystematic(systName("","tau21_eff",""),"lnN",{{signalName,1-0.33}});
            //Btag
            card.addParamSystematic("btag_eff",0.0,0.1);

            //pruned mass scale

    //        card.addParamSystematic("hh_scale",0.0,0.0122); // jes 1 jer 0.5 met 0.5
    //        card.addParamSystematic("hh_res",0.0,0.045); // jes 2 jer 4 met 0.5
            card.addParamSystematic("unclust",0.0,0.01);
            card.addParamSystematic("jes",0.0,0.01);
            card.addParamSystematic("jer",0.0,0.01);

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

            if(reg == REG_NONTOPCR){
                card.addSystematic(systName("top","tFrac",b)         ,"lnN",{{bkgSels[BKG_MT],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_mt/rate_mw}});
                card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_lostTW/rate_mw}});
            } else {
                card.addSystematic(systName("top","wFrac",b)         ,"lnN",{{bkgSels[BKG_MW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_mw/rate_mt}});
                card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_lostTW/rate_mt}});
            }

            //---------------------------------------------------------------------------------------------------
            //Signal
            //---------------------------------------------------------------------------------------------------
            card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                                 {{"hbb_scale",1}},{{"hbb_res",1}},{{"unclust",0.5},{"jes",1},{"jer",0.5}},{{"unclust",0.5},{"jes",2},{"jer",5}}, b == btagCats[BTAG_L],MOD_MS);
            std::string tau21Form = "(1.0+tau21_PtDependence*"+ (p == purCats[PURE_HP] ? "log("+MOD_MS+"/1000)" : "((0.054/0.041)*(-log("+MOD_MS+"/1000)))")+")";
            std::string brealForm = "(1.0+btag_eff*";
            if(b == btagCats[BTAG_L]) brealForm+= "(0.22-4.7*10^(-4)*"+MOD_MS+"+9.4*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_M]) brealForm+="(-0.27+3.3*10^(-4)*"+MOD_MS+"-3.7*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_T]) brealForm+="(0.22+3.8*10^(-4)*"+MOD_MS+"-4.9*10^(-8)*"+MOD_MS+"^(2)))";
            std::string jetForm = "(1.0+unclust)*(1.0+jer)*(1.0+0.5*jes)";
            std::string uncForm = tau21Form+"*"+brealForm+"*"+jetForm;
            double tau21Corr = p == purCats[PURE_HP]  ? 1.03 : 0.95;
            card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
                    tau21Corr,uncForm,{"tau21_PtDependence","btag_eff","unclust","jer","jes"},MOD_MS);
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
//            card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

            card.makeCard();
        }
    }
    if (channel == 0 || channel == 2) {
        for(const auto& l :dilepCats) for(const auto& b :btagCats) for(const auto& s :selCuts){
            if(l == dilepCats[LEP_INCL] ) continue;
            if(b == btagCats[BTAG_LMT]) continue;
            if(s != selCuts[SEL_FULL] ) continue;

            auto card = DataCardMaker(l,b+"_"+s,"13TeV",1,category);

            const std::string cat = l +"_"+b +"_"+s;
            cmd += std::string(" ")+ category +"_"+cat +"_13TeV=datacard_"+category+"_"+cat +"_13TeV.txt";

            auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& s, const std::string& pf) -> std::string
                    {return fPF + "_"+proc +"_"+l +"_"+b +"_"+s +"_"+ pf; };
            auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
            auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

            auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
                return getSystName("",proc,name, sel == "-1" ? l +"_"+b : sel  );
            };

            //Make search variables
            card.addVar(MOD_MJ,100,0,1000,false);
            card.addVar(MOD_MR,1000,0,10000,false);
            card.addVar(MOD_MS,2000,true);

            //---------------------------------------------------------------------------------------------------
            //Get rates and contributions for backgrounds
            //---------------------------------------------------------------------------------------------------
            card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_lostTW = card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_mw =     card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS,1.0,true);
            double rate_mt =     card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS,1.0,true);

            //---------------------------------------------------------------------------------------------------
            //Add Systematics first since the param systs need to have the variables added to the workspace
            //---------------------------------------------------------------------------------------------------
            //luminosity
            card.addSystematic("yield","lnN",{{signalName,1.0391}});//lumi = 2.5 pdf= 2, PU = 0.5, btagfake=1

            //lepton efficiency
            if(l==lepCats[LEP_E])
                card.addSystematic("eff_"+l,"lnN",{{signalName,1.0602}}); //2% trigger / 5.5% for reco  / 1.4% ID / ISO 0.2%
            else
                card.addSystematic("eff_"+l,"lnN",{{signalName,1.0566}}); //2% trigger / ID 1%  /  ISO 5.2%

            //Btag
            card.addParamSystematic("btag_eff",0.0,0.1);

            //pruned mass scale

    //        card.addParamSystematic("hh_scale",0.0,0.0122); // jes 1 jer 0.5 met 0.5
    //        card.addParamSystematic("hh_res",0.0,0.045); // jes 2 jer 4 met 0.5
            card.addParamSystematic("unclust",0.0,0.01);
            card.addParamSystematic("jes",0.0,0.01);
            card.addParamSystematic("jer",0.0,0.01);

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

            if(reg == REG_NONTOPCR){
                card.addSystematic(systName("top","tFrac",b)         ,"lnN",{{bkgSels[BKG_MT],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_mt/rate_mw}});
                card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MW],1.0-0.25*rate_lostTW/rate_mw}});
            } else {
                card.addSystematic(systName("top","wFrac",b)         ,"lnN",{{bkgSels[BKG_MW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_mw/rate_mt}});
                card.addSystematic(systName("top","lostFrac",b)      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.25},{bkgSels[BKG_MT],1.0-0.25*rate_mt/rate_lostTW}}); // changed to mt / tw
            }
            //---------------------------------------------------------------------------------------------------
            //Signal
            //---------------------------------------------------------------------------------------------------
            card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                                 {{"hbb_scale",1}},{{"hbb_res",1}},{{"unclust",0.5},{"jes",1},{"jer",0.5}},{{"unclust",0.5},{"jes",2},{"jer",5}}, b == btagCats[BTAG_L],MOD_MS);
            std::string brealForm = "(1.0+btag_eff*";
            if(b == btagCats[BTAG_L]) brealForm+= "(0.22-4.7*10^(-4)*"+MOD_MS+"+9.4*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_M]) brealForm+="(-0.27+3.3*10^(-4)*"+MOD_MS+"-3.7*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_T]) brealForm+="(0.22+3.8*10^(-4)*"+MOD_MS+"-4.9*10^(-8)*"+MOD_MS+"^(2)))";
            std::string jetForm = "(1.0+unclust)*(1.0+jer)*(1.0+0.5*jes)";
            std::string uncForm = brealForm+"*"+jetForm;
            card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
                    1.0,uncForm,{"btag_eff","unclust","jer","jes"},MOD_MS);

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
//            card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

            card.makeCard();
        }
    }

    std::ofstream outFile("comp.sh",std::ios::out|std::ios::trunc);

    outFile << cmd <<" > combinedCard.txt";
    outFile.close();

    std::cout << cmd <<" > combinedCard.txt"<<std::endl;

}
#endif

void makeCard(int inreg = REG_SR, int insig = RADION, int channel = 1){
	if (channel != 0 && channel != 1 && channel != 2) {
		std::cout<<"channel needs to be either 1 (single lep), 2 (dilep), or 0 (both)"<<std::endl;
		return;
	}
    std::cout <<" <<<<< "<< inreg <<" " <<" "<<signals[insig]<<std::endl;
    REGION reg = REGION(inreg);
    if(reg == REG_NONTOPCR) btagCats = qgBtagCats;
    std::string mainDir = "../";
    go(insig,hhFilename,mainDir,reg,channel);
}
