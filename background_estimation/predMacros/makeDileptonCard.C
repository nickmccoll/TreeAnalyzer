#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
#include "../predTools/CutConstants.h"
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
    printf("debug0\n");

    const std::string category = "std";
    std::string cmd = "combineCards.py ";
    for(const auto& l :dilepCats) for(const auto& b :btagCats) for(const auto& s :selCuts){
        if(l == dilepCats[LEP_INCL] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(s != selCuts[SEL_FULL] ) continue;


        auto card = DataCardMaker(l,b +"_"+s ,"13TeV",1,category);
//        std::vector<double> newYBins =
//        {700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
//                1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,
//                1500,2000,4000
//        };
//        card.rebinY(newYBins);
//        card.rebinY(13,700,2000);

        const std::string cat = l +"_"+b +"_"+s;
        cmd += std::string(" ")+ category +"_"+cat +"_13TeV=datacard_"+category+"_"+cat +"_13TeV.txt";

        auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& h,  const std::string& pf) -> std::string
                {return fPF + "_"+proc +"_"+l +"_"+b +"_"+h +"_"+ pf; };
        auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
        auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

        auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
            return getSystName("",proc,name, sel == "-1" ? l +"_"+b : sel  );
        };

        //Make search variables
        card.addVar(MOD_MJ,100,0,1000,false);
        card.addVar(MOD_MR,1000,0,10000,false);
        card.addVar(MOD_MS,2000,true);
        printf("debug0.1\n");

        //---------------------------------------------------------------------------------------------------
        //Get rates and contributions for backgrounds
        //---------------------------------------------------------------------------------------------------
        double rate_nontop = card.addFixedYieldFromFile(llBkgSels[BKG_NONTOP],1,fPF+"_"+llBkgSels[BKG_NONTOP]+"_distributions.root",llBkgSels[BKG_NONTOP]+"_"+cat+"_"+hhMCS,true);
        double rate_top = card.addFixedYieldFromFile(llBkgSels[BKG_TOP],1,fPF+"_"+llBkgSels[BKG_TOP]+"_distributions.root",llBkgSels[BKG_TOP]+"_"+cat+"_"+hhMCS,true);

        //---------------------------------------------------------------------------------------------------
        //Add Systematics first since the param systs need to have the variables added to the workspace
        //---------------------------------------------------------------------------------------------------
        //luminosity
        card.addSystematic("yield","lnN",{{signalName,1.0391}});//lumi = 2.5 pdf= 2, PU = 0.5, btagfake=1

        //lepton efficiency
        if(l==dilepCats[LEP_OF])
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
        card.addParamSystematic(systName(llBkgSels[BKG_NONTOP]    ,"PTX",b) ,0.0,0.5);
        card.addParamSystematic(systName(llBkgSels[BKG_NONTOP]    ,"OPTX",b),0.0,1.0);
        card.addParamSystematic(systName(llBkgSels[BKG_NONTOP]    ,"PTY")   ,0.0,1.0);
		card.addParamSystematic(systName(llBkgSels[BKG_NONTOP]    ,"OPTY")  ,0.0,1.0);
//        card.addParamSystematic(systName(bkgSels[BKG_QG]    ,"PT2Y") ,0.0,1.0);

        //Normalization
        card.addSystematic(systName(llBkgSels[BKG_NONTOP],"norm")  ,"lnN",{{llBkgSels[BKG_NONTOP],1.5}});
        card.addSystematic(systName("top","norm")            ,"lnN",{{llBkgSels[BKG_TOP],1.25}});


        printf("debug0.2\n");


        //---------------------------------------------------------------------------------------------------
        //Signal
        //---------------------------------------------------------------------------------------------------
//        if(signalName==signalName){
            //Conditional template
            if(!simpleSignal){
//                card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
//                        {{"hbb_scale",1}},{{"hbb_res",1}},{{"unlc",1}},{{"hh_res",1}}, b == btagCats[BTAG_L],MOD_MS);
                card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                        {{"hbb_scale",1}},{{"hbb_res",1}},{{"unclust",0.5},{"jes",1},{"jer",0.5}},{{"unclust",0.5},{"jes",2},{"jer",5}}, b == btagCats[BTAG_L],MOD_MS);
            }else {
                //Non conditional template
                card.add2DSignalParametricShapeNoCond(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
                        {{"hbb_scale",1}},{{"hbb_res",1}},{{"unclust",0.5},{"jes",1},{"jer",0.5}},{{"unclust",0.5},{"jes",2},{"jer",5}}, b == btagCats[BTAG_L],MOD_MS);
            }
            std::string brealForm = "(1.0+btag_eff*";
            if(b == btagCats[BTAG_L]) brealForm+= "(0.22-4.7*10^(-4)*"+MOD_MS+"+9.4*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_M]) brealForm+="(-0.27+3.3*10^(-4)*"+MOD_MS+"-3.7*10^(-8)*"+MOD_MS+"^(2)))";
            if(b == btagCats[BTAG_T]) brealForm+="(0.22+3.8*10^(-4)*"+MOD_MS+"-4.9*10^(-8)*"+MOD_MS+"^(2)))";
            std::string jetForm = "(1.0+unclust)*(1.0+jer)*(1.0+0.5*jes)";
            std::string uncForm = brealForm+"*"+jetForm;
            card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),1.0,uncForm,{"btag_eff","unclust","jer","jes"}
                            ,MOD_MS);
//        } else throw std::invalid_argument("makeCard::go() -> Bad parsing");
            printf("debug0.3\n");

        //---------------------------------------------------------------------------------------------------
        //NONTOP
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts qgKDESysts;
        qgKDESysts.addSyst("PTX",{{systName(llBkgSels[BKG_NONTOP],"PTX",b),"1"  }});
        qgKDESysts.addSyst("OPTX",{{systName(llBkgSels[BKG_NONTOP],"OPTX",b),"1"  }});
        qgKDESysts.addSyst("PTY",{{systName(llBkgSels[BKG_NONTOP],"PTY"),"1"  }});
        qgKDESysts.addSyst("OPTY",{{systName(llBkgSels[BKG_NONTOP],"OPTY"),"1"  }});
//        qgKDESysts.addSyst("PT2Y",{{systName(bkgSels[BKG_QG],"PT2Y"),"1"  }});
        card.addHistoShapeFromFile(llBkgSels[BKG_NONTOP],{MOD_MJ,MOD_MR}, inputName(llBkgSels[BKG_NONTOP],"2D_template.root"),"histo",qgKDESysts);
        printf("debug0.4\n");

        //---------------------------------------------------------------------------------------------------
        //TOP
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts twKDESysts;
        printf("debug0.41\n");

//        twKDESysts.addSyst("PTX",{{systName(llBkgSels[BKG_REALB],"PTX",b),"1"  }});
//        printf("debug0.42\n");
//
//        twKDESysts.addSyst("OPTX",{{systName(llBkgSels[BKG_REALB],"OPTX",b),"1"  }});
//        printf("debug0.43\n");

//        twKDESysts.addSyst("PTY",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b),"1"}});
//        printf("debug0.44\n");
//
//        twKDESysts.addSyst("OPTY",{{systName("top","res"  ),"1"  }});
        printf("debug0.45\n");

        card.addHistoShapeFromFile(llBkgSels[BKG_TOP],{MOD_MJ,MOD_MR},inputName(llBkgSels[BKG_TOP],"2D_template.root"),"histo",twKDESysts);
        printf("debug0.5\n");

        //---------------------------------------------------------------------------------------------------
        //Data
        //---------------------------------------------------------------------------------------------------
        card.importBinnedData(fPF + "_data_distributions.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

//        card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});
        card.makeCard();
        printf("debug1\n");

    }
    std::ofstream outFile("comp.sh",std::ios::out|std::ios::trunc);
    outFile << cmd <<" > combinedCard.txt";
    outFile.close();

    std::cout << cmd <<" > combinedCard.txt"<<std::endl;

}
#endif

void makeDileptonCard(int inreg = REG_SR, int insig = RADION,    bool condSignal= true){
    std::cout <<" <<<<< "<< inreg <<" "<< condSignal <<" "<<signals[insig]<<std::endl;
    REGION reg = REGION(inreg);
    if(reg == REG_QGCR) btagCats = qgBtagCats;
    std::string mainDir = "/Users/brentstone/Dropbox/Physics/HHbbWW/BEtrees/Dilepton17/";
    go(insig,hhFilename,mainDir,reg,!condSignal);
}
