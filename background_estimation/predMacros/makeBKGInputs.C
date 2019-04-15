
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/make2DDetectorParam.C"
#include "../predTools/make1DTemplateWithScaledKernels.C"
#include "../predTools/make1DTemplateWithAdaKern.C"
#include "../predTools/make2DTemplateWithAdaKern.C"
#include "../predTools/fit2DTemplate.C"
#include "../predTools/fit1DTemplate.C"
#include "../predTools/makeJSON.C"
#include "../predTools/mergeHistosToPDF2D.C"
#include "../predTools/makePlots.C"
#include "../predTools/InputsHelper.h"
#include "../predTools/cutHistos1D.C"
#include "../predTools/FunctionFitter.C"
#include "../predTools/fitBKGRatio.C"
#include "../predTools/fit2DFitKDETemplate.C"
#include "TRandom3.h"


std::string getTTBarSF(const std::string& filename){
    CJSON json(filename+"_ttbarSF.json");
    std::string qToW = json.getP(0).second;
    auto replace = [&](const std::string& vn, const std::string tf1n){
        std:size_t index = 0;
        while (true) {
            index = qToW.find(vn, index);
            if (index == std::string::npos) break;
            qToW.replace(index, vn.size(), tf1n);
            index += 1;
        }
    };
    replace(MOD_MS,hhMCS);
    return std::string("(")+processes[TTBAR].cut+"?("+qToW+"):1.0)";
}

std::string getQCDSF(const std::string& name, const std::string& filename,
        const LEPCats l,const PURCats p,const HADCuts h, double sf=1.0 ){
    CJSON json(filename+"_"+name+"_"
            +lepCats[l]+"_"+inclBtagCat+"_"+purCats[p]+"_"+hadCuts[h]
                                                                   +"_QCDSF.json");
    std::string qToW = json.getP(0).second;
    auto replace = [&](const std::string& vn, const std::string tf1n){
        std:size_t index = 0;
        while (true) {
            index = qToW.find(vn, index);
            if (index == std::string::npos) break;
            qToW.replace(index, vn.size(), tf1n);
            index += 1;
        }
    };
    replace(MOD_MS,hhMCS);
    return std::string("(")+processes[WJETS].cut+"?(1+("+flt2Str(sf)+")*(" +qToW+")):1.0)";
}

std::string getMVVExpoMinMaxStr(std::string sample){
    float eMin, eMax;
    if(strFind(sample,bkgSels[BKG_QG])){
        eMin = 2500;
        eMax = 4500;
    } else if(strFind(sample,bkgSels[BKG_LOSTTW])){
        eMin = 2000;
        eMax = 3000;
    } else if(strFind(sample,bkgSels[BKG_MW])){
        eMin = 1500;
        eMax = 2500;
    } else{
        eMin = 2000 ;
        eMax = 3500 ;
    }
    return std::string(" -emin ") + flt2Str(eMin) + " -emax "+flt2Str(eMax) + " ";
}

const float hh_scaleUnc = 0.0005; // coef of 1 means that we get a 100% scale at 2 TeV;
const float hh_resUnc   = 1400; // coef of 1 means that we get a 200% scale at 700 GeV

const float hbb_scaleUnc = 0.005; // coef of 1 means that we get a 100% scale at 200 GeV;
const float hbb_resUnc   = 30;   // coef of 1 means that we get a 100% scale at 30 GeV;

//Originals////
//JJ ->   -hs 0.00714 -hr 45
//VV ->   -hs 0.0003 -hr 1200


void makeBackgroundShapesMJJAdaKernel(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& baseSel="1.0",
        bool addQCDSF = false,float khxs = 1,float khxc = 5){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU )) continue;
        if(strFind(name,bkgSels[BKG_QG])||strFind(name,bkgSels[BKG_LOSTTW])){
            if(!ci.is(BTAG_LMT)) continue;
        } else {
            if(ci.is(BTAG_LMT)) continue;
        }
        if(!ci.is(PURE_I  )) continue;
        if(!ci.is(HAD_LTMB)) continue;

        std::string tempFile=filename+"_"+name+"_"+ci.name()+"_MJJ_incl_template.root";
        std::string cut =  std::string("(")+baseSel+"&&"+ci.cut()+")";
        if(addQCDSF) cut += "*"+getQCDSF(name,filename,ci.l,ci.p,ci.h);
        std::string args = std::string("-v -n histo ")+" -x "+hbbMCS.cut+" -g hbbGenMass " +
                " -xb "+getHbbBinningString(true)+ " -s "+cut+" -w "+nomW.cut
                + " -khs "+ flt2Str(khxs) +" -khc "+ flt2Str(khxc);
        args += " -ks 1.5 -kr 1.5 -hs " +flt2Str(hbb_scaleUnc)+ " -hr " + flt2Str(hbb_resUnc)+ " ";
        args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT -t nsSo ";
        make1DTemplateWithAdaKern(inputFile,tempFile, args);
    }
}

void makeBackgroundShapesMVVAdaKernel(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& baseSel="1.0",bool addQCDSF = false,
        float khxs = 1,float khxc = 5){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    CatIterator ci;
    while(ci.getBin()){
        if(strFind(name,bkgSels[BKG_QG])){
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(HAD_LTMB)) continue;
        }
        if(strFind(name,bkgSels[BKG_LOSTTW])){
            if(!ci.is(LEP_EMU )) continue;
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(HAD_LTMB)) continue;
        }
        if(strFind(name,bkgSels[BKG_MT])||strFind(name,bkgSels[BKG_MW])){
            if(!ci.is(LEP_EMU )) continue;
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(PURE_I  )) continue;
            if(!ci.is(HAD_LT  )) continue;
        }
        std::string tempFile=filename+"_"+name+"_"+ci.name()+"_MVV_incl_template.root";
        std::string cut =  std::string("(")+baseSel+"&&"+ci.cut()+")";
        std::string weight = std::string("(") + nomW.cut+")";
        if(addQCDSF){
            weight += "*"+getQCDSF(name,filename,ci.l,ci.p,ci.h);
        }
        std::string args = std::string("-v -n histo ")+" -x "+hhMCS.cut+" -g hbbGenMass " +
                " -xb "+getHHBinningString(true)+" -s "+cut+" -w "+weight + " -khs "+ flt2Str(khxs)
                +" -khc "+ flt2Str(khxc);
        args += " -ks 1.5 -kr 1.5 -hs "+flt2Str(hh_scaleUnc) +" -hr "+flt2Str(hh_resUnc)+" ";
        args += std::string(" -doS ") + getMVVExpoMinMaxStr(name);
        args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT -t nsSo ";
        make1DTemplateWithAdaKern(inputFile,tempFile, args);
    }

}

void makeBackgroundShapes2DConditional(const std::string name, const std::string filename,
        const std::string inputFile, const std::string baseSel="1.0", bool addQCDSF = false,
        bool xIsCond = false, float khxs = 1,float khxc = 5,float khys = 1,float khyc = 5) {

    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU )) continue;
        if(!ci.is(BTAG_LMT)) continue;
        if(!ci.is(PURE_I  )) continue;
        if(!ci.is(HAD_LTMB  )) continue;

        const std::string tempFile=filename+"_"+name+"_"+ci.name()+"_2D_cond_template.root";
        std::string cut =  std::string("(")+baseSel+"&&"+ci.cut()+")";
        if(addQCDSF) cut += "*"+getQCDSF(name,filename,ci.l,ci.p,ci.h);

        std::string args = std::string("-v -n histo ") + " -vx "+ hbbMCS.cut+ " -vy "+hhMCS.cut
                + " -xb "+getHbbBinningString(true)+" -yb "+getHHBinningString(true)
                +  " -s "+ cut +" -w "+nomW.cut+" ";
        args+=std::string(" -khxs ")+ flt2Str(khxs) +" -khxc "+ flt2Str(khxc)
                        +" -khys "+ flt2Str(khys) +" -khyc "+ flt2Str(khyc) + " ";
        args += "-eopt y  "+ getMVVExpoMinMaxStr(name);
        if(xIsCond){
            args+=" -hs "+ flt2Str(hbb_scaleUnc) + " -hr " + flt2Str(hbb_resUnc) + " ";
            args+=" -xIsCond ";
        } else {
            args+=" -hs "+flt2Str(hh_scaleUnc) +" -hr "+flt2Str(hh_resUnc)+" ";
        }

        make2DTemplateWithAdaKern(inputFile,tempFile, args);
    }
}


void mergeBackgroundShapes(const std::string& name, const std::string& filename, bool xIsCond){

    CatIterator ci;
    while(ci.getBin()){
        if(strFind(name,bkgSels[BKG_QG])){
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(HAD_LTMB)) continue;
        } else {
            if(!ci.is(LEP_EMU)) continue;
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(HAD_LTMB)) continue;
        }
        std::string inFile1D=
                filename+"_"+name+"_"+ ci.name()+"_"+(xIsCond? "MVV":"MJJ") +"_incl_template.root";
        std::string inFile2D=filename+"_"+name+"_"
                +lepCats[LEP_EMU] +"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I] +"_"+hadCuts[HAD_LTMB]
                +"_2D_cond_template.root";
        std::string mergedTemp=filename+"_"+name+"_"+ci.name()+"_2D_merged_template.root";
        std::string args = std::string("-v -n histo ");
        if(xIsCond){ args += " -xIsCond ";}
        args += " -in1D "+  inFile1D + " -in2D "+ inFile2D
                + " -sX Scale:ScaleX,Res:ResX,PT:PTX,OPT:OPTX -sY PT:PTY,OPT:OPTY,PT2:PT2Y "
                + " -xb "+getHbbBinningString(false)+" -yb "+getHHBinningString(false);
        mergeHistosToPDF2D(mergedTemp, args);
    }
}

void cutMVVTemplate(const std::string& name, const std::string& filename){
    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU))continue;
        if(!ci.is(BTAG_LMT))continue;
        if(!ci.is(PURE_I))continue;
        if(!ci.is(HAD_LT))continue;
        std::string inFileX=filename+"_"+name+"_"+ci.name()+"_MVV_incl_template.root";
        std::string outFileX=filename+"_"+name+"_"+ci.name()+"_MVV_cut_template.root";
        std::string args = std::string("-v -n histo ") + " -i "+  inFileX
                + " -s Scale:Scale,Res:Res,PT:PT,OPT:OPT "+ " -xb "+getHHBinningString(false);
        cutHistos1D(outFileX, args);
    }
}




void getQCDScaleFactor(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& cut="1.0"){
    std::vector<PlotVar> vars;
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,inclHHMassBins);
    std::vector<PlotSamp> samps = {
            {name+"_"+processes[TTBAR],processes[TTBAR].cut},
            {name+"_"+processes[WJETS],processes[WJETS].cut},
            {name+"_"+processes[QCD],processes[QCD].cut},
            {name+"_"+processes[OTHER],processes[OTHER].cut}
    };

    //have to use the old iteration to get the inclusive category
    std::vector<CutStr > srPCrBtagCats = btagCats;
    srPCrBtagCats.push_back(inclBtagCat);

    std::vector<PlotSel> sels;
    for(const auto& l :lepCats) for(const auto& b :srPCrBtagCats)
        for(const auto& p :purCats)  for(const auto& h :hadCuts){
            sels.emplace_back(l +"_"+b+"_"+p +"_"+h,
                    l.cut +"&&"+b.cut+"&&"+p.cut+"&&"+h.cut);
        }

    std::string distName=filename+"_"+name+ "_QCDSF_distributions.root";
    MakePlots a(inputFile,distName,samps,sels,vars,cut,nomW.cut);

    for(const auto& l :lepCats) for(const auto& b :srPCrBtagCats)
        for(const auto& p :purCats)  for(const auto& h :hadCuts){
            const std::string catName = l +"_"+b+"_"+p +"_"+h;
            std::string args = "-nF "+distName +" -nH " + catName + "_"+hhMCS.cut
                    + " -nN "+ name+"_"+processes[QCD] + " -nD "+ name+"_"+processes[WJETS]
                    + " -b 0,500,600,700,800,900,1000,1100,1250,1500,2500,3100,5050 -lb "
                    +" -g mix -fVar "+ MOD_MS + " -fMin "  +flt2Str(minHHMass)
                    + " -fMax "  +flt2Str(maxHHMass);
            fitBKGRatio((filename+"_"+name+"_" +catName +"_QCDSF.json").c_str(),args);
        }
}



void makeFittingDistributions(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& cut="1.0",
        bool doIncl = true, bool addQCDSF = false,
        const std::vector<std::pair<std::string,std::string>>& samples={}){
    std::vector<PlotVar> vars;
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,
                hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,
                hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,inclHHMassBins );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,
                hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,
                hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,hhMassBins );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,
                hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,hhMassBins );
    }
    std::vector<PlotSamp> samps = { {name,"1.0"}};
    if(samples.size()){
        samps.clear();
        for(const auto& s: samples)
            samps.emplace_back(s.first,s.second);
    }
    std::vector<PlotSel> sels;
    CatIterator ci;
    while(ci.getBin()){
        if(addQCDSF){
            sels.emplace_back(ci.name(),
                    ci.cut()+"*"+getQCDSF(name,filename,ci.l,ci.p,ci.h));
            sels.emplace_back( std::string("noQCDSF_") +  ci.name(),ci.cut());
        } else {
            sels.emplace_back(ci.name(),ci.cut());
        }
    }
    std::string outFileName=filename+"_"+name+
            (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,nomW.cut);
}

void fitBackgroundShapes2DConditional(std::string name, const std::string& filename,
        bool xIsCond, const std::string& fittedName =""){
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    CatIterator ci;
    while(ci.getBin()){
        std::string hName = name+"_"+ci.name()+"_"+hbbMCS+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+ci.name()+"_2D_template.root";

        std::string tempFile= filename+"_"+name+"_"
                +lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[ci.p]+"_"+hadCuts[HAD_LTMB]
                +"_2D_merged_template.root";
        if(strFind(name,bkgSels[BKG_QG])){
            tempFile=filename+"_"+name+"_"
                    +lepCats[ci.l]+"_"+btagCats[BTAG_LMT]+"_"+purCats[ci.p]+"_"+hadCuts[HAD_LTMB]
                    +"_2D_merged_template.root";
        }
        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PTX,OPTX,PTY,OPTY "
                +" -sA PT2Y " + " -fH " + distFileName + " -nH "+ hName;
        if(xIsCond) args+=" -xCy ";
        fit2DTemplate(outName,args);

    }
}
void fitBackgroundShapesMVV(std::string name, const std::string& filename,
        const std::string& fittedName =""){
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    CatIterator ci;
    while(ci.getBin()){
        std::string hName = name+"_"+ci.name()+"_"+hhMCS;
        std::string inName = filename + "_"+name+"_"
                +lepCats[LEP_EMU]+"_"+btagCats[BTAG_LMT]+"_"+purCats[PURE_I]+"_"+hadCuts[HAD_LT]
                +"_MVV_cut_template.root";
        std::string outName = filename + "_"+name+"_"+ci.name()+"_MVV_template.root";
        std::string args = std::string("-v ") + "-fT "+ inName+" -nT histo -s PT,OPT "
                + " -fH " + distFileName + " -nH "+ hName;
        fit1DTemplate(outName,args);
    }
}

void makeBKG1DShapes(const std::string& name, const std::string& filename,
        const std::string& catName, const std::string& fitName,
        bool isW, CJSON* prevJSON, TFile* iF){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;
    auto setup1DFit = [&](const TH1* hbbH, double HHMass)->FunctionFitter*{
        auto vN=[&](std::string var)->std::string{return var;};
        fitters.emplace_back(new CBFunctionFitter(hbbH,{},false,"",{MOD_MJ}));

        FunctionFitter* fitter = fitters.back().get();
        if(isW){
            fitter->setVar(vN("mean")     ,90,80,100);
            fitter->setVar(vN("sigma")     ,8,5,15);
            fitter->setVar(vN("alpha")     ,1.18,0.1,10);
            fitter->setVar(vN("alpha2")  ,0.9,0.1,5);
            fitter->setVar(vN("n")  ,5,1,6);
            fitter->setVar(vN("n2")  , 2,1,6);
            fitter->setConst(vN("n")  ,1);
            fitter->setConst(vN("n2")  ,1);
            fitter->w->var(MOD_MJ.c_str())->setRange("fit",30,160);
        } else {
            fitter->setVar(vN("mean")     ,180,120,195);
            fitter->setVar(vN("sigma")     ,15.8,14,20);
            fitter->setVar(vN("alpha")     ,1.5 ,0.1,2);
            fitter->setVar(vN("alpha2")  ,1.5,0.5,3);
            fitter->setVar(vN("n")  ,5,1,6);
            fitter->setVar(vN("n2")  , 5,1,6);
            fitter->setConst(vN("n")  ,1);
            fitter->setConst(vN("n2")  ,1);
            fitter->w->var(MOD_MJ.c_str())->setRange("fit",70,230);
        }
        if(prevJSON){
            fitter->setVar(vN("alpha")     ,prevJSON->evalFunc(vN("alpha")  ,HHMass) ,0.1,10);
            fitter->setConst(vN("alpha"),1);
            fitter->setVar(vN("alpha2")  ,prevJSON->evalFunc(vN("alpha2")  ,HHMass),0.1,10);
            fitter->setConst(vN("alpha2")  ,1);
        }
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::NumCPU(8)});
        return fitter;
    };


    std::string hName = name+"_"+catName;
    auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
    if(hbb_hh_H==0) return;
    auto hhH  = projY(&*hbb_hh_H,hName+"_"+hhMCS);

    for(unsigned int iP = 0; iP < resPTBins.size() -1; ++iP){
        std::string ptName  = flt2Str(resPTBins[iP]) +"to"+flt2Str(resPTBins[iP+1]);
        auto hbbH = projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS,resPTBins[iP],resPTBins[iP+1]);
        double mean = getMean(&*hhH,resPTBins[iP],resPTBins[iP+1] );
        plotter.addFit(setup1DFit(&*hbbH,mean),mean,ptName);
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");
}


void makeResWMJJShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit1stIt";

    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU)) continue;
        if(!ci.is(PURE_I)) continue;
        if(!ci.is(HAD_NONE)) continue;
        makeBKG1DShapes(name,filename,ci.name(),fitName,true,0,iF);

        std::string argsP1 = std::string("-i ")+filename+"_"+name+"_"+ci.name()+"_"+fitName+".root"
                +" -var "+MOD_MR+" ";
        argsP1 += " -minX 600 -maxX 3000 ";
        std::string jsonArgsStd =
                " -g mean:laur2,sigma:laur2,alpha:laur4,alpha2:laur3,n:pol0,n2:pol0 ";
        MakeJSON(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json",argsP1+" "+  jsonArgsStd );
    }
}

void makeResWMJJShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit";
    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU)) continue;
        if(!ci.is(PURE_I)) continue;
        if(!ci.is(HAD_NONE)) continue;
        CJSON oldJSON(     filename+"_"+name+"_"+ci.name()+"_MJJ_fit1stIt.json");
        oldJSON.fillFunctions(MOD_MR);
        makeBKG1DShapes(name,filename,ci.name(),fitName,true,&oldJSON,iF);

        std::string argsP1 = std::string("-i ")+filename+"_"+name+"_"+ci.name()+"_"+fitName+".root"
                +" -var "+MOD_MR+" ";
        argsP1 += " -minX 600 -maxX 3000 ";
        std::string jsonArgsStd =
                " -g mean:laur2,sigma:laur2,alpha:laur4,alpha2:laur3,n:pol0,n2:pol0 ";

        CJSON newJSON = getJSON(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json",
                argsP1+" "+jsonArgsStd);
        newJSON.replaceEntry("alpha", oldJSON.getP("alpha") );
        newJSON.replaceEntry("alpha2", oldJSON.getP("alpha2") );
        newJSON.write(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json");
    }
}


void makeResTopMJJShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit1stIt";

    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU)) continue;
        if(!ci.is(PURE_I)) continue;
        if(!ci.is(HAD_NONE)) continue;
        makeBKG1DShapes(name,filename,ci.name(),fitName,false,0,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+ci.name()+"_"+fitName+".root"
                +" -var "+MOD_MR+" ";
        argsP1 += " -minX 700 -maxX 3500 ";
        std::string jsonArgsStd =
                " -g mean:laur3,sigma:laur2,alpha:laur4,alpha2:laur3,n:pol0,n2:pol0 ";
        MakeJSON(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json",argsP1+" "+  jsonArgsStd );


    }
}

void makeResTopMJJShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "MJJ_fit";
    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_EMU )) continue;
        if(!ci.is(PURE_I  )) continue;
        if(!ci.is(HAD_NONE)) continue;
        CJSON oldJSON(     filename+"_"+name+"_"+ci.name()+"_MJJ_fit1stIt.json");
        oldJSON.fillFunctions(MOD_MR);
        makeBKG1DShapes(name,filename,ci.name(),fitName,false,&oldJSON,iF);

        std::string argsP1 = std::string("-i ")+filename+"_"+name+"_"+ci.name()+"_"+fitName+".root"
                +" -var "+MOD_MR+" ";
        argsP1 += " -minX 700 -maxX 3000 ";
        std::string jsonArgsStd=
                " -g mean:laur3,sigma:laur2,alpha:laur4,alpha2:laur3,n:pol0,n2:pol0 ";

        CJSON newJSON = getJSON(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json"
                ,argsP1+" "+jsonArgsStd);
        newJSON.replaceEntry("alpha", oldJSON.getP("alpha") );
        newJSON.replaceEntry("alpha2", oldJSON.getP("alpha2") );
        newJSON.write(filename+"_"+name+"_"+ci.name()+"_"+fitName+".json");
    }
}


void fitMJJSF(const std::string& name, const std::string& filename){
    std::string fittingFileName = filename+"_"+name+"_distributions.root";
    CatIterator ci;
    while(ci.getBin()){
        std::string jsonFile = filename+"_"+name+"_"+lepCats[LEP_EMU]+"_";
        jsonFile += strFind(name,bkgSels[BKG_MW]) ?btagCats[BTAG_LMT]:btagCats[ci.b];
        jsonFile+=std::string("_")+purCats[PURE_I]+"_"+hadCuts[HAD_NONE] +"_MJJ_fit.json";

        std::string kdeFile = filename + "_"+name+"_"+ci.name()+"_MVV_template.root";
        std::string outFile = filename + "_"+name+"_"+ci.name()+"_MJJ_SFFit.json";

        std::string args = std::string(" -v -fT ")+kdeFile+" -nT histo "+ " -varX "+MOD_MJ
                +" -varY "+MOD_MR+ " -fJS "+jsonFile + " -sS 0.05 -sR 0.2 "
                +" -fH "+ fittingFileName+ " -nH "+ name+"_"+ci.name()+"_"+hbbMCS+"_"+hhMCS+" ";
        if(strFind(name,bkgSels[BKG_MT])) args+= " -sA1 0.2 -mTop ";

        fit2DFitKDETemplate(outFile,args);
    }
}

void convertFuncFitTo2DTemplate(const std::string& name, const std::string& filename,
        const std::string& funcParamPostfix=""){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    CatIterator ci;
    while(ci.getBin()){
        std::string jsonFile = filename + "_"+name+"_"+ci.name() +"_MJJ_SFFit.json";

        std::unique_ptr<TAxis> xAx(new TAxis(nHbbMassBins*10,minHbbMass,maxHbbMass));

        CBFunctionFitter xFit(0,{xAx.get()},false,funcParamPostfix,{MOD_MJ});

        CJSON json(jsonFile);
        json.fillFunctions(MOD_MR);
        auto pn = [&](const std::string& v) ->std::string{return v + funcParamPostfix;};
        auto * iF =  TObjectHelper::getFile(filename + "_"+name+"_"+ci.name()+"_MVV_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH1>(iF,"histo",false,false);
        if(hh_H==0) continue;
        oF->cd();
        TH2 * h_2D = new TH2F((name+"_"+ci.name()).c_str(),
                (std::string(";")+hbbMCS.title+";"+hhMCS.title).c_str(),
                nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,&hhMassBins[0]);
        for(int iY = 1; iY <= hh_H->GetNbinsX(); ++iY){
            double yNorm = hh_H->GetBinContent(iY);
            double hhV = hh_H->GetBinCenter(iY);

            xFit.setVar(pn("mean"  ),json.evalFunc(pn("mean"  ),hhV));
            xFit.setVar(pn("sigma" ),json.evalFunc(pn("sigma" ),hhV));
            xFit.setVar(pn("alpha" ),json.evalFunc(pn("alpha" ),hhV));
            xFit.setVar(pn("alpha2"),json.evalFunc(pn("alpha2"),hhV));
            xFit.setVar(pn("n"     ),json.evalFunc(pn("n"     ),hhV));
            xFit.setVar(pn("n2"    ),json.evalFunc(pn("n2"    ),hhV));
            auto xHist = xFit.pdf1D(name+"_"+ci.name()+ASTypes::int2Str(iY)+"_hbb");

            for(int iX = 1; iX <= h_2D->GetNbinsX(); ++iX){
                double xNorm = xHist->Integral( 10*iX -9,10*iX  );
                h_2D->SetBinContent(iX,iY,xNorm*yNorm);
            }
            delete xHist;
        }
        h_2D->Write();
        iF->Close();
    }
    oF->Close();
}

void compile2DTemplatesForDebug(const std::string& name, const std::string& filename){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    CatIterator ci;
    while(ci.getBin()){
        auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_"+ci.name()+"_2D_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH2>(iF,"histo",false,false);
        if(hh_H==0) continue;
        oF->cd();
        hh_H->SetXTitle(hbbMCS.title.c_str());
        hh_H->SetYTitle(hhMCS.title.c_str());
        hh_H->Write((name+"_"+ci.name()).c_str());
        iF->Close();
    }
    oF->Close();
}

void makePseudoData(const std::string& name, const std::string& filename, const double sf = 1){
    TRandom3 * rand = new TRandom3(1234);
    std::vector<TFile*> yieldFiles;
    std::vector<TFile*> tempFiles;
    for(unsigned int bkg = BKG_QG; bkg <= BKG_MT; ++bkg){
        yieldFiles.push_back(TObjectHelper::getFile(
                filename+"_"+bkgSels[bkg]+"_distributions.root"));
        tempFiles.push_back(TObjectHelper::getFile(
                filename+"_"+bkgSels[bkg]+"_2D_template_debug.root"));
    }
    TFile *oF = new TFile((filename + "_"+name+".root").c_str(),"recreate");
    CatIterator ci;
    while(ci.getBin()){
        if(ci.is(LEP_EMU)) continue;
        if(ci.is(PURE_I)) continue;
        if(!ci.is(HAD_FULL)) continue;
        double totProb = 0;
        TH2* totH = 0;
        for(unsigned int bkg = BKG_QG; bkg <= BKG_MT; ++bkg){
            auto y_H = TObjectHelper::getObject<TH2>(
                    yieldFiles[bkg],bkgSels[bkg]+"_"+ci.name()+"_"+hbbMCS+"_"+hhMCS);
            auto t_H = TObjectHelper::getObject<TH2>(tempFiles[bkg],bkgSels[bkg]+"_"+ci.name());
            std::cout <<bkgSels[bkg]+"_"+ci.name() +" -> "<< y_H->Integral()<<" "<<t_H->Integral();
            if(totH) std::cout<<" "<< totH->Integral();
            std::cout <<"\n";
            t_H->Scale(y_H->Integral()/t_H->Integral());
            if(totH == 0) totH = (TH2*)t_H->Clone();
            else totH->Add(&*t_H);
        }

        for(int iX = 1; iX <= totH->GetNbinsX(); ++iX)
            for(int iY = 1; iY <= totH->GetNbinsY(); ++iY){
                double mean = totH->GetBinContent(iX,iY);
                int val = rand->Poisson(sf*mean);
                totH->SetBinContent(iX,iY,val);
                totH->SetBinError(iX,iY,std::sqrt(float(val)));
            }
        totH->Write((std::string("data_")+ci.name()+"_"+hbbMCS+"_"+hhMCS).c_str() );
    }

    oF->Close();
    for(auto* f: yieldFiles)f->Close();
    for(auto* f: tempFiles)f->Close();
    delete rand;
}

void makeDataDistributions(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    std::vector<PlotVar> vars;
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,
                nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,
                hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,hhMassBins );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,
                nHbbMassBins,minHbbMass,maxHbbMass,
                hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,hhMassBins );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,
                nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,hhMassBins );
    }
    std::vector<PlotSamp> samps = { {name,"1.0"}};
    std::vector<PlotSel> sels;
    CatIterator ci;
    while(ci.getBin()){
        if(!ci.is(HAD_FULL)) continue;
        sels.emplace_back(ci.name(),ci.cut());
    }
    std::string outFileName=filename+"_"+name+ (
            doIncl ? "_inclM_distributions.root" : "_distributions.root");
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,"1.0");
}

void go(int modelToDo, std::string treeDir) {
    std::string filename = hhFilename;
    if(modelToDo == -1){
        makeDataDistributions("data",hhFilename,treeDir+"../betrees_data.root",
                hhRange.cut+"&&"+hbbRange.cut,false);
        return;
    }
    if(modelToDo ==5){
        return;
    }

    //Turn on TTBar scaling
    nomW.cut = nomW.cut +"*"+getTTBarSF("../supportInputs/HHlnujj");
    std::string treeAreaIncl = treeDir + "../betrees_mc.root";

    if(modelToDo == BKG_QG)
    {
        std::string name = bkgSels[BKG_QG];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = bkgSels[BKG_QG].cut + "&&"+ aQCD.cut;
        getQCDScaleFactor(name,filename,  treeAreaIncl, bkgSels[BKG_QG].cut+"&&"+hbbRange.cut);
        makeFittingDistributions(name,filename,treeArea,
                hhInclRange.cut+"&&"+hbbInclRange.cut,true,true,
                {{bkgSels[BKG_QG],genSel},{bkgSels[BKG_QG]+"_wQCD",bkgSels[BKG_QG].cut}}
        );
        makeFittingDistributions(name,filename,treeArea,
                hhRange.cut+"&&"+hbbRange.cut,false,true,
                {{bkgSels[BKG_QG],genSel},{bkgSels[BKG_QG]+"_wQCD",bkgSels[BKG_QG].cut}}
        );

        makeBackgroundShapes2DConditional(name,filename,treeArea,
                genSel,true,true,0.4,1,0.75,3);//P(hbb|hh)
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut,true,1,3);
        mergeBackgroundShapes(name,filename,true);
        fitBackgroundShapes2DConditional(name,filename,true);

        compile2DTemplatesForDebug(name,filename);
    }

    if(modelToDo == BKG_LOSTTW){
        std::string name = bkgSels[BKG_LOSTTW];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = bkgSels[BKG_LOSTTW].cut;
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);


        makeBackgroundShapes2DConditional(name,filename,treeArea,
                genSel,false,true,0.5,3,0.75,5);//P(hbb|hh) 0.5,2,0.75,8 old values
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut,false,1,4);
        mergeBackgroundShapes(name,filename,true);
        fitBackgroundShapes2DConditional(name,filename,true);

        compile2DTemplatesForDebug(name,filename);
    }

    if(modelToDo == BKG_MW){
        std::string name = bkgSels[BKG_MW];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = bkgSels[BKG_MW].cut;
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);

        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        fitBackgroundShapesMVV(name,filename);

        makeResWMJJShapes1stIt(name,filename);
        makeResWMJJShapes2ndIt(name,filename);
        fitMJJSF(name,filename);
        convertFuncFitTo2DTemplate(name,filename);
    }

    if(modelToDo == BKG_MT){
        std::string name = bkgSels[BKG_MT];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = bkgSels[BKG_MT].cut;
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);

        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        fitBackgroundShapesMVV(name,filename);

        makeResTopMJJShapes1stIt(name,filename);
        makeResTopMJJShapes2ndIt(name,filename);
        fitMJJSF(name,filename);
        convertFuncFitTo2DTemplate(name,filename);
    }

    //Make pseudo data
    if(int(modelToDo) == 4){//4
        makePseudoData("pd",filename,1);
    }

}
#endif

void makeBKGInputs(int bkgToDo = BKG_QG, std::string treeDir = "../trees/bkgCompLMT/"){
    go(bkgToDo,treeDir);
}
