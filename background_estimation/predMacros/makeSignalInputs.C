
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "../predTools/makeJSON.C"
#include "../predTools/makePlots.C"
#include "../predTools/InputsHelper.h"
#include "../predTools/FunctionFitter.C"

void makeSignalFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    int nameIDX =  inputFile.find("XXX", 0);
    std::vector<PlotVar> vars;
    std::vector<PlotSel> sels;
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        sels.emplace_back(l +"_"+b+"_"+p +"_"+h,
                l.cut +"&&"+b.cut+"&&"+p.cut+"&&"+h.cut);
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
void makeSignalYields(const std::string& name, const std::string& filename, const double BR= 2*0.5824*(.2137+.002619)){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        FunctionParameterPlotter plotter;
        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        TGraphErrors* yieldGraph = new TGraphErrors;
        int n = 0;
        for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
            std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
            auto hh_H = TObjectHelper::getObject<TH1>(iF,hName+"_"+hhMCS,false,false);
            if(hh_H ==0) {
                std::cout << hName+"_"+hhMCS<<std::endl;
                continue;
            }
            double error = 0;
            double integral = hh_H->IntegralAndError(1,hh_H->GetNbinsX(),error);
            yieldGraph->SetPoint(n,signalMassBins[iS],integral*BR);
            yieldGraph->SetPointError(n,0.0,error*BR);
            n++;
        }
        plotter.addFit(yieldGraph,"yield");
        plotter.write(filename+"_"+name+"_"+catName+"_yield.root");
        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_yield.root ";
        argsP1 += " -minX 700 -maxX 3600 ";
        argsP1 += " -g yield:laur5 ";
        argsP1 += " -var "+MOD_MS+" ";
        MakeJSON(filename+"_"+name+"_"+catName+"_yield.json",argsP1);
        delete yieldGraph;
    }

    iF->Close();
    delete iF;
}


std::string vnMJ(const std::string& var){return var+MOD_MJ;};
std::string vnMR(const std::string& var){return var+MOD_MR;};
const std::string fitMJJStd  =vnMJ("mean")+":laur4,"+vnMJ("sigma")+":laur4,"+vnMJ("alpha")+":laur4,"+vnMJ("alpha2")+":laur4,"+vnMJ("n")+":pol0,"+vnMJ("n2")+":pol0";
const std::string fitMJJExpo =vnMJ("mean")+":laur4,"+vnMJ("sigma")+":laur4,"+vnMJ("alpha")+":laur4,"+vnMJ("alpha2")+":pol1," +vnMJ("n")+":pol0,"+vnMJ("n2")+":pol0,"+vnMJ("slope")+":laur4,"+vnMJ("fE")+":pol4";
const std::string fitMVV     =vnMR("mean")+":pol1," +vnMR("sigma")+":pol1," +vnMR("alpha")+":laur3," +vnMR("alpha2")+":laur3,"+vnMR("n")+":pol0,"+vnMR("n2")+":pol0";
const std::string fitCond    =vnMR("maxS")+":pol0,"+vnMR("mean_p1")+":pol2,"+vnMR("sigma_p1")+":laur3";
//const std::string fitCond    =vnMR("mean_p1")+":pol2";
const std::string fitCondMVV =vnMJ("maxS")+":pol0,"+vnMJ("mean_p1")+":pol2,"+vnMJ("sigma_p1")+":pol1";


void makeSignal1DShapes(const std::string& name, const std::string& filename, const std::string& catName, const std::string& fitName, bool fitMJJ, CJSON* prevJSON, TFile* iF, bool doExpo, bool restrictX = true){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;
    const std::string modStr = fitMJJ ? MOD_MJ : MOD_MR;
    auto setup1DFit = [&](const TH1* hbbH, double HHMass){
        auto vN=[&](std::string var)->std::string{return fitMJJ ? vnMJ(var):vnMR(var);};
        fitters.emplace_back(new CBFunctionFitter(hbbH,doExpo,modStr,{modStr}));
        auto fitter = &* fitters.back();

        if(fitMJJ){
            fitter->w->var(modStr.c_str())->setRange("fit",30,210);
            fitter->w->var(modStr.c_str())->setRange("coef",30,210);
        } else {
            //            fitter->w->var(modStr.c_str())->setRange("fit",std::max(HHMass*.75,minHHMass),std::min(HHMass*1.5,maxHHMass));
            fitter->w->var(modStr.c_str())->setRange("fit",minHHMass,maxHHMass);
            fitter->w->var(modStr.c_str())->setRange("coef",minHHMass,maxHHMass);
        }
        //        fitter->w->pdf((std::string("model")+modStr).c_str())->fixAddCoefRange("coef",true);


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


        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::SumCoefRange("coef"),RooFit::NumCPU(8)});
    };


    for(unsigned int iS = 0; iS < signalMassBins.size(); ++iS){
        std::string hName   = name+"_m"+ASTypes::int2Str(signalMassBins[iS])+"_"+catName;
        std::string ptName  = std::string("m") + ASTypes::flt2Str(signalMassBins[iS]);
        auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
        if(hbb_hh_H ==0) continue;
        auto hbb_H = fitMJJ ? projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS) :  (restrictX ? projY(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS,115,135) : projY(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS)) ;
        setup1DFit(&*hbb_H,signalMassBins[iS]);
        plotter.addFit(&*fitters.back(),signalMassBins[iS],ptName);
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");

}


void makeSignalMJJShapes1stIt(const std::string& name, const std::string& filename){
    //    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root"); //Changing to cond on MR
    const std::string fitName = "MJJ_fit1stIt";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l != lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_I] ) continue;
        if(p != purCats[PURE_I]) continue;
        if(h != hadCuts[HAD_LTMB]) continue;
        bool doExpo = b == btagCats[BTAG_L];
        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        makeSignal1DShapes(name,filename,catName,fitName,true,0,iF,doExpo);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+" -var "+MOD_MS +" ";
        argsP1 += doExpo ? " -minX 700 -maxX 3800 " :" -minX 700 -maxX 3800 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" -g "+  ( doExpo  ?  fitMJJExpo : fitMJJStd ));
    }
}
void makeSignalMJJShapes2ndIt(const std::string& name, const std::string& filename){
    //    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root"); //Changing to cond on MR
    const std::string fitName = "MJJ_fit";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l != lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_I] ) continue;
        if(p != purCats[PURE_I]) continue;
        if(h != hadCuts[HAD_LTMB]) continue;
        bool doExpo = b == btagCats[BTAG_L];

        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_MJJ_fit1stIt.json");
        oldJSON.fillFunctions(MOD_MS);
        makeSignal1DShapes(name,filename,catName,fitName,true,&oldJSON,iF,doExpo);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 700 -maxX 3800 "+" -var "+MOD_MS +" ";
        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" -g "+  ( doExpo  ?  fitMJJExpo : fitMJJStd ));
        newJSON.replaceEntry(vnMJ("alpha"), oldJSON.getP(vnMJ("alpha")) );
        newJSON.replaceEntry(vnMJ("alpha2"), oldJSON.getP(vnMJ("alpha2")) );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}

void makeSignalMVVShapes1D(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    const std::string fitName = "MVV_fit1stIt";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b != btagCats[BTAG_LMT] ) continue;
        if(p != purCats[PURE_I] ) continue;
        if(h != hadCuts[HAD_LTMB]) continue;
        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        makeSignal1DShapes(name,filename,catName,fitName,false,0,iF,false,true);
        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+" -var "+MOD_MS +" ";
        argsP1 += " -minX 700 -maxX 3800 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" -g "+ fitMVV);
    }
}

void makeSignalMVVShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    const std::string fitName = "MVV_fit1stIt";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b != btagCats[BTAG_LMT] ) continue;
//        if(p != purCats[PURE_I] ) continue;
        if(!(h == hadCuts[HAD_LTMB] || h == hadCuts[HAD_FULL])  ) continue;
        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        makeSignal1DShapes(name,filename,catName,fitName,false,0,iF,false,false);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+" -var "+MOD_MS +" ";
        argsP1 += " -minX 700 -maxX 3800 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" -g "+ fitMVV);
    }
}
void makeSignalMVVShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    const std::string fitName = "MVV_fit";
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b != btagCats[BTAG_LMT] ) continue;
        //        if(p != purCats[PURE_I] ) continue;
        if(!(h == hadCuts[HAD_LTMB] || h == hadCuts[HAD_FULL])  ) continue;
        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_MVV_fit1stIt.json");
        oldJSON.fillFunctions(MOD_MS);
        makeSignal1DShapes(name,filename,catName,fitName,false,&oldJSON,iF,false,false);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+" -var "+MOD_MS +" ";
        argsP1 +=  " -minX 700 -maxX 3800 ";
        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" -g "+fitMVV);
        newJSON.replaceEntry(vnMJ("alpha"), oldJSON.getP(vnMJ("alpha")) );
        newJSON.replaceEntry(vnMJ("alpha2"), oldJSON.getP(vnMJ("alpha2")) );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");


    }
}

void makeSignal2DShapes(const std::string& name, const std::string& filename,const std::string& catName, std::string& fitName, CJSON* mjjJSON,CJSON* mvvJSON, bool fixMVVinTerms, bool doExpo, TFile* iF, const bool cond = true){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;

    auto setup2DFit = [&](const TH2* hbbH, const double HHMass){
        auto pnX =[&] (std::string v) ->std::string{return vnMJ(v);};
        auto pnY =[&] (std::string v) ->std::string{return vnMR(v);};
        if(cond)
            fitters.emplace_back(new CBFunctionFitter2D(hbbH,doExpo,"",{MOD_MJ,MOD_MR}));
        else
            fitters.emplace_back(new CBFunctionFitter2DNoCond(hbbH,doExpo,"",{MOD_MJ,MOD_MR}));
        auto fitter = &* fitters.back();
        fitter->w->var(MOD_MR.c_str())->setRange("coef",minHHMass,maxHHMass);
        fitter->w->var(MOD_MJ.c_str())->setRange("coef",30,210);
        fitter->w->var(MOD_MJ.c_str())->setRange("fit",30,210);
        fitter->w->var(MOD_MR.c_str())->setRange("fit",std::max(HHMass*.75,minHHMass),std::min(HHMass*1.25,maxHHMass));

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
        fitter->setConst(pnY("n")     ,1);
        fitter->setConst(pnY("n2")    ,1);
        if(cond){
            fitter->setVar(pnY("maxS")  ,2.5 ,0,5);
            fitter->setVar(pnY("mean_p1")  ,.044,.01,.10);
            fitter->setVar(pnY("sigma_p1") ,0,0,1);
            fitter->setConst(pnY("maxS")  ,1);
        }


        if(mvvJSON){
            fitter->setVar(pnY("alpha") ,mvvJSON->evalFunc(pnY("alpha") ,HHMass),0.1,3);
            fitter->setVar(pnY("alpha2"),mvvJSON->evalFunc(pnY("alpha2"),HHMass),0.1,3);
            fitter->setVar(pnY("mean")  ,mvvJSON->evalFunc(pnY("mean"),HHMass),HHMass -200,HHMass+200);
            fitter->setVar(pnY("sigma") ,mvvJSON->evalFunc(pnY("sigma"),HHMass),HHMass*0.025,HHMass*0.2);
            fitter->setConst(pnY("alpha")     ,1);
            fitter->setConst(pnY("alpha2")    ,1);
            if(fixMVVinTerms){
                if(cond){
                    fitter->setVar(pnY("mean_p1")  ,mvvJSON->evalFunc(pnY("mean_p1") ,HHMass),.01,.10);
                    fitter->setVar(pnY("sigma_p1") ,mvvJSON->evalFunc(pnY("sigma_p1") ,HHMass),0,1);
                    fitter->setConst(pnY("mean_p1")     ,1);
                    fitter->setConst(pnY("sigma_p1")    ,1);
                }
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

void makeSignal2DShapesCondMVV(const std::string& name, const std::string& filename, std::string& catName, std::string& fitName, CJSON* mjjJSON,CJSON* mvvJSON, bool fixCorrTerms, bool doExpo, TFile* iF){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;

    auto setup2DFit = [&](const TH2* hbbH, const double HHMass){
        auto pnX =[&] (std::string v) ->std::string{return vnMJ(v);};
        auto pnY =[&] (std::string v) ->std::string{return vnMR(v);};

        fitters.emplace_back(new CBFunctionFitter2DCondMVV(hbbH,doExpo,"",{MOD_MJ,MOD_MR}));
        auto fitter = &* fitters.back();
        fitter->w->var(MOD_MJ.c_str())->setRange("fit",30,210);
        fitter->w->var(MOD_MR.c_str())->setRange("fit",std::max(HHMass*.75,minHHMass),std::min(HHMass*1.5,maxHHMass));
        fitter->w->var(MOD_MR.c_str())->setRange("coef",minHHMass,maxHHMass);
        fitter->w->var(MOD_MJ.c_str())->setRange("coef",30,210);

        if(mvvJSON){
            fitter->setVar(pnY("mean")  ,mvvJSON->evalFunc(pnY("mean"),HHMass),HHMass -200,HHMass+200);
            fitter->setVar(pnY("sigma") ,mvvJSON->evalFunc(pnY("sigma"),HHMass),HHMass*0.025,HHMass*0.2);
            fitter->setVar(pnY("alpha") ,mvvJSON->evalFunc(pnY("alpha") ,HHMass),0.1,3);
            fitter->setVar(pnY("alpha2"),mvvJSON->evalFunc(pnY("alpha2"),HHMass),0.1,3);
            fitter->setVar(pnY("n")     ,mvvJSON->evalFunc(pnY("n"),HHMass),1,6);
            fitter->setVar(pnY("n2")    ,mvvJSON->evalFunc(pnY("n2"),HHMass),1,6);
            fitter->setConst(pnY("mean")     ,1);
            fitter->setConst(pnY("sigma")   ,1);
            fitter->setConst(pnY("alpha")     ,1);
            fitter->setConst(pnY("alpha2")   ,1);
        } else {
            fitter->setVar(pnY("mean")  ,HHMass   ,HHMass -200,HHMass+200);
            fitter->setVar(pnY("sigma") ,HHMass*0.05,HHMass*0.025,HHMass*0.2);
            fitter->setVar(pnY("alpha") ,1.5,0.1,3);
            fitter->setVar(pnY("alpha2"),1.5,0.1,3);
            fitter->setVar(pnY("n")     ,5,1,6);
            fitter->setVar(pnY("n2")    ,5,1,6);
        }
        fitter->setConst(pnY("n")     ,1);
        fitter->setConst(pnY("n2")   ,1);


        fitter->  setVar(pnX("mean")  ,125,90,180);
        fitter->  setVar(pnX("sigma") ,10,5,20);
        fitter->  setVar(pnX("alpha") ,0.9 ,0.1,2);
        fitter->  setVar(pnX("alpha2"),1.5,0.1,3);
        fitter->  setVar(pnX("n")     ,5  ,1,6);
        fitter->  setVar(pnX("n2")    ,5,3,6);
        fitter->setVar(pnX("maxS")  ,2.5 ,0,5);
        fitter->setVar(pnX("mean_p1")  ,.044,.01,.10);
        //        fitter->setVar(pnX("sigma_p1")  ,0,0,1);
        if(doExpo){
            fitter->setVar(pnX("slope"),-1,-10,0);
            fitter->setVar(pnX("fE")   ,0.1,0,0.75);
        }

        if(mjjJSON){
            fitter->  setVar(pnX("alpha") ,mjjJSON->evalFunc(pnX("alpha") ,HHMass),0.1,2);
            fitter->  setVar(pnX("alpha2"),mjjJSON->evalFunc(pnX("alpha2"),HHMass),0.1,3);
            fitter->  setVar(pnX("n")     ,mjjJSON->evalFunc(pnX("n")     ,HHMass),1,6);
            fitter->  setVar(pnX("n2")    ,mjjJSON->evalFunc(pnX("n2")    ,HHMass),3,6);
            fitter->setConst(pnX("alpha")     ,1);
            fitter->setConst(pnX("alpha2")     ,1);
            if(doExpo){
                fitter->setVar(pnX("slope"),mjjJSON->evalFunc(pnX("slope")     ,HHMass),-10,0);
                fitter->setVar(pnX("fE")   ,mjjJSON->evalFunc(pnX("fE")     ,HHMass),0,0.75);
                fitter->setConst(pnX("slope")     ,1);
                fitter->setConst(pnX("fE")     ,1);
            }
            if(fixCorrTerms){
                fitter->  setVar(pnX("mean_p1") ,mjjJSON->evalFunc(pnX("mean_p1") ,HHMass),.01,.10);
                fitter->  setVar(pnX("sigma_p1") ,mjjJSON->evalFunc(pnX("sigma_p1") ,HHMass),0,1);
                fitter->  setVar(pnX("maxS"),mjjJSON->evalFunc(pnX("maxS"),HHMass),0,5);
                fitter->setConst(pnX("mean_p1")     ,1);
                fitter->setConst(pnX("sigma_p1")     ,1);
            }
        }
        fitter->setConst(pnX("maxS")     ,1);
        fitter->setConst(pnX("n")     ,1);
        fitter->setConst(pnX("n2")   ,1);

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

//void makeSignal2DShapesCondMVVFirstIteration(const std::string& name, const std::string& filename){
//    std::string fitName = "2D_fit1stIt";
//    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
//    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
//        if(p !=  purSels[PUR_LMT] ) continue;
//        if(h !=  hadSels[HAD_LTMB] ) continue;
//
//        bool doExpo = p ==  purSels[PUR_L];
//
//        std::string catName = l+"_"+p+"_"+h;
//        std::string mjjcatName = lepSels[LEP_EMU] +"_"+p+"_"+hadSels[HAD_LTMB];
//        std::string mvvcatName = l +"_"+purSels[PUR_LMT]+"_"+hadSels[HAD_LTMB];
//        CJSON mjjJSON(     filename+"_"+name+"_"+mjjcatName+"_MJJ_fit1stIt.json");
//        mjjJSON.fillFunctions(MOD_MS);
//        CJSON mvvJSON(     filename+"_"+name+"_"+mvvcatName+"_MVV_fit.json");
//        mvvJSON.fillFunctions(MOD_MS);
//        makeSignal2DShapesCondMVV(name,filename,catName,fitName,&mjjJSON,&mvvJSON,false,doExpo,iF);
//        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 700 -maxX 3800 -var "+MOD_MS+" ";
//        std::string jsonArgs = argsP1 + " -g " +(doExpo ? fitMJJExpo : fitMJJStd ) + ","+fitMVV+","+fitCondMVV;
//        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
//        newJSON.replaceEntries(mvvJSON);
////        newJSON.replaceEntry(vnMJ("alpha"), mjjJSON.getP(vnMJ("alpha")) );
////        newJSON.replaceEntry(vnMJ("alpha2"), mjjJSON.getP(vnMJ("alpha2")) );
//        if(doExpo){
//            newJSON.replaceEntry(vnMJ("fE"), mjjJSON.getP(vnMJ("fE")) );
//            newJSON.replaceEntry(vnMJ("slope"), mjjJSON.getP(vnMJ("slope")) );
//        }
//        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
//    }
//}


//void makeSignal2DShapesFirstIteration(const std::string& name, const std::string& filename){
//    std::string fitName = "2D_fit1stIt";
//    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
//    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
//        if(l == lepSels[LEP_EMU] ) continue;
//        if(p !=  purSels[PUR_LMT] ) continue;
//        if(h !=  hadSels[HAD_LTMB] ) continue;
//
//        bool doExpo = p ==  purSels[PUR_L];
//
//        std::string catName = l+"_"+p+"_"+h;
//        std::string mjjcatName = lepSels[LEP_EMU] +"_"+p+"_"+hadSels[HAD_LTMB];
//        std::string mvvcatName = l +"_"+purSels[PUR_LMT]+"_"+hadSels[HAD_LTMB];
//        CJSON mjjJSON(     filename+"_"+name+"_"+mjjcatName+"_MJJ_fit.json");
//        mjjJSON.fillFunctions(MOD_MS);
//        CJSON mvvJSON(     filename+"_"+name+"_"+mvvcatName+"_MVV_fit1stIt.json");
//        mvvJSON.fillFunctions(MOD_MS);
//        makeSignal2DShapes(name,filename,catName,fitName,&mjjJSON,&mvvJSON,false,doExpo,iF);
//        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 700 -maxX 3800 -var "+MOD_MS+" ";
//        std::string jsonArgs = argsP1 + " -g " +(doExpo ? fitMJJExpo : fitMJJStd ) + ","+fitMVV+","+fitCond;
//        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
//        newJSON.replaceEntries(mjjJSON);
//        newJSON.replaceEntry(vnMR("alpha"), mvvJSON.getP(vnMR("alpha")) );
//        newJSON.replaceEntry(vnMR("alpha2"), mvvJSON.getP(vnMR("alpha2")) );
//        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
//    }
//}

void makeSignal2DShapesSecondIteration(const std::string& name, const std::string& filename){
    std::string fitName = "2D_fit";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_I] || b == btagCats[BTAG_LMT] ) continue;
        if(p == purCats[PURE_I]) continue;
        if(h != hadCuts[HAD_FULL] ) continue;
        bool doExpo = b == btagCats[BTAG_L];

        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        std::string mjjCatName = lepCats[LEP_EMU]+"_"+b+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB];
        std::string mvvCatName = l+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB];

        CJSON mjjJSON(     filename+"_"+name+"_"+mjjCatName+"_MJJ_fit.json");
        mjjJSON.fillFunctions(MOD_MS);
        CJSON mvvJSON(     filename+"_"+name+"_"+mvvCatName+"_MVV_fit1stIt.json");
        mvvJSON.fillFunctions(MOD_MS);
        makeSignal2DShapes(name,filename,catName,fitName,&mjjJSON,&mvvJSON,false,doExpo,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 700 -maxX 3800 -var "+MOD_MS+" ";
        std::string jsonArgs = argsP1 +" -g "+ (doExpo ? fitMJJExpo : fitMJJStd ) + ","+fitMVV+","+fitCond;


        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        newJSON.replaceEntries(mjjJSON);
        newJSON.replaceEntry(vnMR("alpha"   ), mvvJSON.getP(vnMR("alpha"   )));
        newJSON.replaceEntry(vnMR("alpha2"  ), mvvJSON.getP(vnMR("alpha2"  )));
        //        newJSON.replaceEntry(vnMR("mean_p1" ), mvvJSON.getP(vnMR("mean_p1" )));
        //        newJSON.replaceEntry(vnMR("sigma_p1"), mvvJSON.getP(vnMR("sigma_p1")));
        //        newJSON.replaceEntry(vnMR("mean"    ), mvvJSON.getP(vnMR("mean"    )));
        //        newJSON.replaceEntry(vnMR("sigma"   ), mvvJSON.getP(vnMR("sigma"   )));
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}


void combine2DShapesNoCond(const std::string& name, const std::string& filename){
    std::string fitName = "2D_fit";
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_distributions.root");
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats)  for(const auto& h :hadCuts){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_I] || b == btagCats[BTAG_LMT] ) continue;
        if(p == purCats[PURE_I]) continue;
        if(h != hadCuts[HAD_FULL] ) continue;
        bool doExpo = b == btagCats[BTAG_L];

        const std::string catName = l +"_"+b+"_"+p +"_"+h;
        std::string mjjCatName = lepCats[LEP_EMU]+"_"+b+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LTMB];
        std::string mvvCatName = l+"_"+btagCats[BTAG_LMT]+"_"+p+"_"+hadCuts[HAD_FULL];

        CJSON mjjJSON(     filename+"_"+name+"_"+mjjCatName+"_MJJ_fit.json");
        mjjJSON.fillFunctions(MOD_MS);
        CJSON mvvJSON(     filename+"_"+name+"_"+mvvCatName+"_MVV_fit1stIt.json");
        mvvJSON.fillFunctions(MOD_MS);
        makeSignal2DShapes(name,filename,catName,fitName,&mjjJSON,&mvvJSON,false,doExpo,iF,false);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root "+ " -minX 700 -maxX 3800 -var "+MOD_MS+" ";
        std::string jsonArgs = argsP1 +" -g "+ (doExpo ? fitMJJExpo : fitMJJStd ) + ","+fitMVV;


        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",jsonArgs);
        newJSON.replaceEntries(mjjJSON);
        newJSON.replaceEntries(mvvJSON);
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}





void go(int step,std::string treeDir) {
    std::string filename = hhFilename;
    std::string signalTrees = treeDir + "/out_radion_hh_bbinc_mXXX_0.root";
    std::string name = radionSig;
    if(step == 0){
        makeSignalFittingDistributions(name,filename,signalTrees,hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeSignalFittingDistributions(name,filename,signalTrees,hhRange.cut+"&&"+hbbRange.cut,false);
        makeSignalYields(name,filename);
    }

    if(step == 1){
        makeSignalMJJShapes1stIt(name,filename);
        makeSignalMJJShapes2ndIt(name,filename);
        makeSignalMVVShapes1D(name,filename);
        makeSignal2DShapesSecondIteration(name,filename);
    }

    if(step == 2){
        makeSignalMJJShapes1stIt(name,filename);
        makeSignalMJJShapes2ndIt(name,filename);
        makeSignalMVVShapes1stIt(name,filename);
//        makeSignalMVVShapes2ndIt(name,filename);
        combine2DShapesNoCond(name,filename);
    }




    //old stuff
    //    makeSignal2DShapesFirstIteration(name,filename);
    //    makeSignalMVVShapes1stIt(name,filename);
    //    makeSignalMVVShapes2ndIt(name,filename);
    //    makeSignalMJJShapes1stIt(name,filename);
    //    makeSignal2DShapesCondMVVFirstIteration(name,filename);
}
#endif

void makeSignalInputs(int step = 0,std::string treeDir = "../trees/"){
    go(step,treeDir);
}
