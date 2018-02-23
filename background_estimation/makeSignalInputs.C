
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "makeJSON.C"
#include "makePlots.C"
#include "InputsHelper.h"
#include "FunctionFitter.C"

void makeSignalFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    int nameIDX =  inputFile.find("XXX", 0);
    std::vector<PlotVar> vars;
    std::vector<PlotSel> sels;
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        sels.emplace_back(l +"_"+p +"_"+h,
                l.cut +"&&"+p.cut+"&&"+h.cut);
    }
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    }

    for(const auto& sM : signalMassBins){
        std::vector<PlotSamp> samps = { {name +"_m"+ASTypes::int2Str(sM),"1.0"}};
        std::string outFileName=filename+"_"+name+"_m"+ASTypes::int2Str(sM) + (doIncl ? "_inclM_distributions.root" : "_distributions.root");
        std::string inputName = inputFile;
        inputName.replace(nameIDX,3,ASTypes::int2Str(sM));
        MakePlots a(inputName,outFileName,samps,sels,vars,cut,nomW.cut);
    }
    std::string compiledFile =  filename+"_"+name + (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    std::string allFiles = filename+"_"+name+"_m*" + (doIncl ? "_inclM_distributions.root" : "_distributions.root");

    gSystem->Exec((std::string("hadd -f ")+ compiledFile + " " + allFiles).c_str());
    gSystem->Exec((std::string("rm ") + allFiles).c_str());
}
void makeSignal1DShapes(const std::string& name, const std::string& filename, const std::string& catName, const std::string& fitName, bool fitMJJ, CJSON* prevJSON, TFile* iF, bool doExpo){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;
    auto setup1DFit = [&](const TH1* hbbH, double HHMass){
        std::string pF = fitMJJ ? "SMJJ" : "SMVV";
        auto vN=[&](std::string var)->std::string{return var+pF;};
        fitters.emplace_back(new CBFunctionFitter(hbbH,doExpo,pF,{ fitMJJ ? "MJJ" : "MVV"}));

        auto fitter = &* fitters.back();
        if(fitMJJ){
            fitter->setVar(vN("mean")     ,125,90,180);
            fitter->setVar(vN("sigma")       ,10,5,20);
        } else {
            fitter->setVar(vN("mean")  ,HHMass   ,HHMass -200,HHMass+200);
            fitter->setVar(vN("sigma") ,HHMass*0.05,HHMass*0.025,HHMass*0.2);
        }

        fitter->setVar(vN("n")   ,  5  ,1,6);
        fitter->setVar(vN("n2")  ,5,3,6);
        fitter->setConst(vN("n")  ,1);
        fitter->setConst(vN("n2")  ,1);

        if(prevJSON){
            fitter->setVar(vN("alpha")     ,prevJSON->evalFunc(vN("alpha")  ,HHMass) ,0.1,2);
            fitter->setConst(vN("alpha"),1);
            fitter->setVar(vN("alpha2")  ,prevJSON->evalFunc(vN("alpha2")  ,HHMass),0.1,3);
            fitter->setConst(vN("alpha2")  ,1);
        } else {
            fitter->setVar(vN("alpha")     ,0.9 ,0.1,2);
            fitter->setVar(vN("alpha2")  ,1.5,0.1,3);
        }
        if(doExpo){
            fitter->setVar(vN("slope")  ,-1,-10,0);
            fitter->setVar(vN("fE")  ,0.1,0,0.75);
        }
        if(fitMJJ){
            fitter->w->var("MJJ")->setRange("fit",30,210);
            fitter->w->var("MJJ")->setRange("coef",30,210);
        } else {
            fitter->w->var("MVV")->setRange("fit",HHMass*.75,HHMass*1.50);
            fitter->w->var("MVV")->setRange("coef",HHMass*.75,HHMass*1.50);
        }

        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
    };


    for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
        std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
        std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
        auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
        if(hbb_hh_H ==0) continue;
        auto hbb_H = fitMJJ ? projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS) :  projY(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS,115,135);
        setup1DFit(&*hbb_H,signalMassBins[iS]);
        plotter.addFit(&*fitters.back(),signalMassBins[iS],ptName);
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");

}

void makeSignalMJJShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit1stIt";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_LTMB] ) continue;
        bool doExpo = p == purSels[PUR_L];
        const std::string catName = l+"_"+p+"_"+h;
        makeSignal1DShapes(name,filename,catName,fitName,true,0,iF,doExpo);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += p == purSels[PUR_L] ? " -minX 550 -maxX 4550 " :" -minX 550 -maxX 4550 ";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0 ";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+  ( doExpo  ?  jsonArgsExpo : jsonArgsStd ));
    }
}
void makeSignalMJJShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_LTMB] ) continue;
        bool doExpo = p == purSels[PUR_L];
        const std::string catName = l+"_"+p+"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_MJJ_fit1stIt.json");
        oldJSON.fillFunctions("MH");
        makeSignal1DShapes(name,filename,catName,fitName,true,&oldJSON,iF,doExpo);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 550 -maxX 4550";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0 ";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4 ";
        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+  ( doExpo  ?  jsonArgsExpo : jsonArgsStd ));
        newJSON.replaceEntry(std::string("alphaS")+"MJJ", oldJSON.getP(std::string("alphaS")+"MJJ") );
        newJSON.replaceEntry(std::string("alpha2S")+"MJJ", oldJSON.getP(std::string("alpha2S")+"MJJ") );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}
void makeSignalMVVShapes1D(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MVV_fit";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l == lepSels[LEP_EMU] ) continue;
//        if(!(p == purSels[PUR_LMT] || p == purSels[PUR_M]) ) continue;
        if(!(h == hadSels[HAD_LTMB] || h == hadSels[HAD_FULL]) ) continue;
        const std::string catName = l+"_"+p+"_"+h;
        makeSignal1DShapes(name,filename,catName,fitName,false,0,iF,false);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += p == purSels[PUR_L] ? " -minX 550 -maxX 4550 " :" -minX 550 -maxX 4550 ";
        std::string jsonArgsStd = " -g meanSMVV:pol1,sigmaSMVV:pol1,alphaSMVV:pol1,alpha2SMVV:laur3,nSMVV:pol0,n2SMVV:pol0 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+ jsonArgsStd);
    }
}

void makeSignal2DShapes(const std::string& name, const std::string& filename, std::string& catName, std::string& fitName, CJSON* mjjJSON,CJSON* mvvJSON, bool fixMVVinTerms, bool doExpo, TFile* iF){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;

    auto setup2DFit = [&](const TH2* hbbH, const double HHMass){
        std::string pF = "S";
        auto pnX =[&] (std::string v) ->std::string{return v+pF+"MJJ";};
        auto pnY =[&] (std::string v) ->std::string{return v+pF+"MVV";};

        fitters.emplace_back(new CBFunctionFitter2D(hbbH,doExpo,pF,{"MJJ","MVV"}));
        auto fitter = &* fitters.back();
        fitter->w->var("MJJ")->setRange("fit",30,210);
        fitter->w->var("MVV")->setRange("fit",HHMass*.75,HHMass*1.25);
        fitter->w->var("MVV")->setRange("coef",minInclHHMass,maxInclHHMass);
        fitter->w->var("MJJ")->setRange("coef",30,210);

        if(mjjJSON){
            fitter->  setVar(pnX("mean")  ,mjjJSON->evalFunc(pnX("mean")  ,HHMass),90,180);
            fitter->  setVar(pnX("sigma") ,mjjJSON->evalFunc(pnX("sigma") ,HHMass),5,20);
            fitter->  setVar(pnX("alpha") ,mjjJSON->evalFunc(pnX("alpha") ,HHMass),0.1,3);
            fitter->  setVar(pnX("alpha2"),mjjJSON->evalFunc(pnX("alpha2"),HHMass),0.1,3);
            fitter->  setVar(pnX("n")     ,mjjJSON->evalFunc(pnX("n")     ,HHMass),1,6);
            fitter->  setVar(pnX("n2")    ,mjjJSON->evalFunc(pnX("n2")    ,HHMass),3,6);
            fitter->setConst(pnX("mean")  ,1);
            fitter->setConst(pnX("sigma") ,1);
            fitter->setConst(pnX("alpha") ,1);
            fitter->setConst(pnX("alpha2"),1);
            fitter->setConst(pnX("n")     ,1);
            fitter->setConst(pnX("n2")    ,1);
            if(doExpo){
                fitter->setVar(pnX("slope") ,mjjJSON->evalFunc(pnX("slope") ,HHMass),-10,0);
                fitter->setVar(pnX("fE")    ,mjjJSON->evalFunc(pnX("fE")    ,HHMass),0,0.75);
                fitter->setConst(pnX("slope")     ,1);
                fitter->setConst(pnX("fE")    ,1);
            }
        }else {
            fitter->  setVar(pnX("mean")  ,125,90,180);
            fitter->  setVar(pnX("sigma") ,10,5,20);
            fitter->  setVar(pnX("alpha") ,0.9 ,0.1,2);
            fitter->  setVar(pnX("alpha2"),1.5,0.1,3);
            fitter->  setVar(pnX("n")     ,5  ,1,6);
            fitter->  setVar(pnX("n2")    ,5,3,6);
            fitter->setConst(pnX("n")     ,1);
            fitter->setConst(pnX("n2")   ,1);
            if(doExpo){
                fitter->setVar(pnX("slope"),-1,-10,0);
                fitter->setVar(pnX("fE")   ,0.1,0,0.75);
            }
        }

        fitter->setVar(pnY("mean")  ,HHMass   ,HHMass -200,HHMass+200);
        fitter->setVar(pnY("sigma") ,HHMass*0.05,HHMass*0.025,HHMass*0.2);
        fitter->setVar(pnY("alpha") ,1.5,0.1,3);
        fitter->setVar(pnY("alpha2"),1.5,0.1,3);
        fitter->setVar(pnY("n")     ,5,1,6);
        fitter->setVar(pnY("n2")    ,5,1,6);
        fitter->setVar(pnY("maxS")  ,2.5 ,0,5);
        fitter->setVar(pnY("mean_p1")  ,.044,.01,.10);
        fitter->setVar(pnY("sigma_p1") ,0,0,1);
        fitter->setConst(pnY("n")     ,1);
        fitter->setConst(pnY("n2")    ,1);
        fitter->setConst(pnY("maxS")  ,1);
//        if(doExpo){
//            fitter->setVar(pnY("meanE")  ,HHMass   ,HHMass -200,HHMass+200);
//            fitter->setVar(pnY("sigmaE") ,HHMass*0.05,HHMass*0.025,HHMass*0.2);
//        }

        if(mvvJSON){
            fitter->setVar(pnY("alpha") ,mvvJSON->evalFunc(pnY("alpha") ,HHMass),0.1,3);
            fitter->setVar(pnY("alpha2"),mvvJSON->evalFunc(pnY("alpha2"),HHMass),0.1,3);
            fitter->setConst(pnY("alpha")     ,1);
            fitter->setConst(pnY("alpha2")    ,1);
            if(fixMVVinTerms){
                fitter->setVar(pnY("mean_p1")  ,mvvJSON->evalFunc(pnY("mean_p1") ,HHMass),.01,.10);
                fitter->setVar(pnY("sigma_p1") ,mvvJSON->evalFunc(pnY("sigma_p1") ,HHMass),0,1);
                fitter->setConst(pnY("mean_p1")     ,1);
                fitter->setConst(pnY("sigma_p1")    ,1);
                fitter->setVar(pnY("mean")  ,mvvJSON->evalFunc(pnY("mean"),HHMass),HHMass -200,HHMass+200);
                fitter->setVar(pnY("sigma") ,mvvJSON->evalFunc(pnY("sigma"),HHMass),HHMass*0.025,HHMass*0.2);
                fitter->setConst(pnY("mean")     ,1);
                fitter->setConst(pnY("sigma")    ,1);
            }
        }

        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
//        fitter->w->writeToFile((filename +"_"+name+"_"+catName+"_"+ASTypes::flt2Str(HHMass)+"_temp.root").c_str(),true);
        std::cout<< "---------------------------END!-------------------------------"<<std::endl;

    };

    for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
        std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
        std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
        auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
        if(hbb_hh_H ==0) continue;
        setup2DFit(&*hbb_hh_H,signalMassBins[iS]);
        plotter.addFit2D(&*fitters.back(),signalMassBins[iS],ptName);
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");
}


void makeSignal2DShapesFirstIteration(const std::string& name, const std::string& filename){
    std::string fitName = "2D_fit1stIt";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l == lepSels[LEP_EMU] ) continue;
                if(p !=  purSels[PUR_L] ) continue;
        //        if(p !=  purSels[PUR_LMT] ) continue;
//        if(!(p ==  purSels[PUR_LMT] ||p ==  purSels[PUR_L]) ) continue;
        if(h !=  hadSels[HAD_LTMB] ) continue;

        bool doExpo = p ==  purSels[PUR_L];

        std::string catName = l+"_"+p+"_"+h;
        std::string mjjcatName = lepSels[LEP_EMU] +"_"+p+"_"+hadSels[HAD_LTMB];
        std::string mvvcatName = l +"_"+purSels[PUR_LMT]+"_"+hadSels[HAD_LTMB];
        CJSON mjjJSON(     filename+"_"+name+"_"+mjjcatName+"_MJJ_fit.json");
        mjjJSON.fillFunctions("MH");
        CJSON mvvJSON(     filename+"_"+name+"_"+mvvcatName+"_MVV_fit.json");
        mvvJSON.fillFunctions("MH");
        makeSignal2DShapes(name,filename,catName,fitName,&mjjJSON,&mvvJSON,false,doExpo,iF);
        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 550 -maxX 4550 -v ";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4,meanESMVV:pol1,sigmaESMVV:pol1";
        std::string jsonArgsVV0 = "meanSMVV:pol1,sigmaSMVV:pol1,alphaSMVV:pol1,alpha2SMVV:laur3,nSMVV:pol0,n2SMVV:pol0";
        std::string jsonArgsVV1 = "maxSSMVV:pol0,mean_p1SMVV:pol2,sigma_p1SMVV:pol1";
        std::string jsonArgs = argsP1 + (doExpo ? jsonArgsExpo : jsonArgsStd ) + ","+jsonArgsVV0+","+jsonArgsVV1;
        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        newJSON.replaceEntries(mjjJSON);
        newJSON.replaceEntry("alphaSMVV", mvvJSON.getP("alphaSMVV") );
        newJSON.replaceEntry("alpha2SMVV", mvvJSON.getP("alpha2SMVV") );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}

void makeSignal2DShapesSecondIteration(const std::string& name, const std::string& filename){
    std::string fitName = "2D_fit";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l == lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I] || p == purSels[PUR_LMT]) continue;
        if(h != hadSels[HAD_FULL] ) continue;
        bool doExpo = p == purSels[PUR_L];

        std::string catName = l+"_"+p+"_"+h;
        std::string mjjCatName = lepSels[LEP_EMU]+"_"+p+"_"+hadSels[HAD_LTMB];
        std::string mvvCatName = l+"_"+purSels[PUR_LMT]+"_"+hadSels[HAD_LTMB];

        CJSON mjjJSON(     filename+"_"+name+"_"+mjjCatName+"_MJJ_fit.json");
        mjjJSON.fillFunctions("MH");
        CJSON mvvJSON(     filename+"_"+name+"_"+mvvCatName+"_2D_fit1stIt.json");
        mvvJSON.fillFunctions("MH");
        makeSignal2DShapes(name,filename,catName,fitName,&mjjJSON,&mvvJSON,true,doExpo,iF);


        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 550 -maxX 4550";
        std::string jsonArgsStd = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:laur4,nSMJJ:pol0,n2SMJJ:pol0";
        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4";
//        std::string jsonArgsExpo = " -g meanSMJJ:laur4,sigmaSMJJ:laur4,alphaSMJJ:laur4,alpha2SMJJ:pol1,nSMJJ:pol0,n2SMJJ:pol0,slopeSMJJ:laur4,fESMJJ:pol4,meanESMVV:pol1,sigmaESMVV:pol1";
        std::string jsonArgsVV0 = "meanSMVV:pol1,sigmaSMVV:pol1,alphaSMVV:pol1,alpha2SMVV:laur3,nSMVV:pol0,n2SMVV:pol0";
        std::string jsonArgsVV1 = "maxSSMVV:pol0,mean_p1SMVV:pol2,sigma_p1SMVV:pol1";
        std::string jsonArgs = argsP1 + (doExpo ? jsonArgsExpo : jsonArgsStd ) + ","+jsonArgsVV0+","+jsonArgsVV1;


        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        newJSON.replaceEntries(mjjJSON);
        newJSON.replaceEntry("alphaSMVV", mvvJSON.getP("alphaSMVV") );
        newJSON.replaceEntry("alpha2SMVV", mvvJSON.getP("alpha2SMVV") );
        newJSON.replaceEntry("mean_p1SMVV", mvvJSON.getP("mean_p1SMVV") );
        newJSON.replaceEntry("sigma_p1SMVV", mvvJSON.getP("sigma_p1SMVV") );
        newJSON.replaceEntry("meanSMVV"   , mvvJSON.getP("meanSMVV") );
        newJSON.replaceEntry("sigmaSMVV"  , mvvJSON.getP("sigmaSMVV") );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}





void go(std::string treeDir) {
    std::string filename = hhFilename;
    std::string signalTrees = treeDir + "/out_radion_hh_bbinc_mXXX_0.root";
    std::string name = radionSig;

    //    makeSignalFittingDistributions(name,filename,signalTrees,hhRange.cut+"&&"+hbbRange.cut,false);
    //    makeSignalMJJShapes1stIt(name,filename);
    //    makeSignalMJJShapes2ndIt(name,filename);
//    makeSignalMVVShapes1D(name,filename);
//    makeSignal2DShapesFirstIteration(name,filename);
    makeSignal2DShapesSecondIteration(name,filename);


    //        makeSignalFittingDistributions(name,filename,signalTrees,hhInclRange.cut+"&&"+hbbInclRange.cut,true);
    //            makeSignal2DShapes(name,filename);

    //    makeSignalMJJShapes1stIt(name,filename);
    //        makeSignalMJJShapes(name,filename);
    //            makeSignal1DMVVShapes(name,filename);
    //        makeSignal2DShapesFirstIteration(name,filename);



}
#endif

void makeSignalInputs(std::string treeDir = "trees/"){
    go(treeDir);
}
