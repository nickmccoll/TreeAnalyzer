
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "make2DDetectorParam.C"
#include "make1DTemplateWithScaledKernels.C"
#include "make1DTemplateWithAdaKern.C"
#include "make2DTemplateWithAdaKern.C"
#include "fit2DTemplate.C"
#include "fit1DTemplate.C"
#include "makeJSON.C"
#include "mergeHistosToPDF2D.C"
#include "makePlots.C"
#include "InputsHelper.h"
#include "cutHistos1D.C"
#include "FunctionFitter.C"

void makeDetectorParam(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0"){
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    std::string xArgs = std::string("-v -n xHisto -x ")+ hbbMCS.cut +"/hbbGenMass "+ " -xb 500,0,5 "+" -s "+cut +" -w "+nomW.cut+ " -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000";
    std::string yArgs = std::string("-v -n yHisto -x ")+ hhMCS.cut +"/genhhMass "+ " -xb 500,0,5 "+" -s "+cut +" -w "+nomW.cut+ " -y hbbGenPT -ylb -yb 200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000,5000";
    make2DDetectorParam(inputFile,std::string("x_")+resFile,xArgs);
    make2DDetectorParam(inputFile,std::string("y_")+resFile,yArgs);
    gSystem->Exec(TString::Format("hadd -f %s *_%s",resFile.c_str(),resFile.c_str()));
    gSystem->Exec(TString::Format("rm *_%s",resFile.c_str()));
}

void makeBackgroundShapesMJJAdaKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0",float khxs = 1,float khxc = 5){
    std::string tempFile=filename+"_"+name+"_incl_template.root";
    std::string resFile=filename+"_"+name+"_detectorResponse.root";
    std::string args = std::string("-v -n histo ")+" -x "+hbbMCS.cut+" -g hbbGenMass -xb "+hbbInclBinning.cut+ " -s "+cut+" -w "+nomW.cut+ " -khs "+ flt2Str(khxs) +" -khc "+ flt2Str(khxc);
    args += " -kss -ks 1.5 -kr 1.5 -hs 0.00714 -hr 45 ";
    args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT ";
    make1DTemplateWithAdaKern(inputFile,tempFile, args);
}
void makeBackgroundShapesMVVAdaKernel(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0",float khxs = 1,float khxc = 5){
    std::string tempFile=filename+"_"+name+"_incl_template.root";
    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    float eMin = (name.find("mt") != std::string::npos) ?  2000:1500;
    float eMax = (name.find("mt") != std::string::npos) ?  4500:2500;

    std::string args = std::string("-v -n histo ")+" -x "+hhMCS.cut+" -g hbbGenMass -xb "+hhInclBinning.cut+" -s "+cut+" -w "+nomW.cut+ " -khs "+ flt2Str(khxs) +" -khc "+ flt2Str(khxc);
    args += " -kss -ks 1.5 -kr 1.5 -hs 0.0003 -hr 1200 ";
    args += std::string(" -doS -emin ") + flt2Str(eMin) + " -emax "+flt2Str(eMax) + " ";
    args += std::string(" -vsf ")+resFile+ " -vsh scalexHisto -vsv hbbGenPT ";
    make1DTemplateWithAdaKern(inputFile,tempFile, args);
}

void makeBackgroundShapesMVVConditional(const std::string name, const std::string filename,  const std::string inputFile, const std::string cut="1.0", float khxs = 1,float khxc = 5,float khys = 1,float khyc = 5) {
    std::string tempFile=filename+"_"+name+"_incl_COND2D_template.root";
    float eMin = (name.find("lost") != std::string::npos) ? 1500 : 2000;
    float eMax = (name.find("lost") != std::string::npos) ? 2500 : 4500;
    std::string args = std::string("-v -n histo ") + " -vx "+ hhMCS.cut+ " -vy "+hbbMCS.cut+ " -xb "+hhInclBinning.cut+" -yb "+hbbInclBinning.cut+ " -ycb 50,60,80,100,120,140,160,180,200,220,250 "+ "-s "+ cut +" -w "+nomW.cut;
    args+=std::string(" -khxs ")+ flt2Str(khxs) +" -khxc "+ flt2Str(khxc) +" -khys "+ flt2Str(khys) +" -khyc "+ flt2Str(khyc) + " -hss ";
    args+=std::string(" -hs 0.0003 -hr 1200 ") + " -emin "+ flt2Str(eMin)+ " -emax "+ flt2Str(eMax);
    make2DTemplateWithAdaKern(inputFile,tempFile, args);
}


void mergeBackgroundShapes(const std::string& name, const std::string& filename){
    std::string inFileX=filename+"_"+name+"_incl_template.root";
    std::string inFileY=filename+"_"+name+"_incl_COND2D_template.root";
    std::string rootFile=filename+"_"+name+"_2D_template.root";
    std::string args = std::string("-v -n histo ") + " -inX "+  inFileX + " -inY "+ inFileY + " -sX Scale:ScaleX,Res:ResX,PT:PTX,OPT:OPTX -sY PT:PTY,OPT:OPTY "+ " -xb "+hbbBinning.cut+" -yb "+hhBinning.cut;
    mergeHistosToPDF2D(rootFile, args);
}

void cutMVVTemplate(const std::string& name, const std::string& filename){
    std::string inFileX=filename+"_"+name+"_incl_template.root";
    std::string rootFile=filename+"_"+name+"_template.root";
    std::string args = std::string("-v -n histo ") + " -i "+  inFileX + " -s Scale:Scale,Res:Res,PT:PT,OPT:OPT "+ " -xb "+hhBinning.cut;
    cutHistos1D(rootFile, args);
}



void makeFittingDistributions(const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut="1.0", bool doIncl = true){
    std::vector<PlotVar> vars;
    if(doIncl){
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nInclHbbMassBins,minInclHbbMass,maxInclHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nInclHHMassBins,minInclHHMass,maxInclHHMass );
    } else {
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
        vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
        vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    }
    std::vector<PlotSamp> samps = { {name,"1.0"}};
    std::vector<PlotSel> sels;
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        sels.emplace_back(l +"_"+p +"_"+h,
                l.cut +"&&"+p.cut+"&&"+h.cut);
    }
    std::string outFileName=filename+"_"+name+ (doIncl ? "_inclM_distributions.root" : "_distributions.root");
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,nomW.cut);
}

void fitBackgroundShapesMVVConditional(std::string name, const std::string& filename, const std::string& fittedName =""){
    std::string tempFile=filename+"_"+name+"_2D_template.root";
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string hName = name+"_"+l+"_"+p+"_"+h+"_"+hbbMCS+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+l+"_"+p+"_"+h+"_2D_template.root";
        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PTX,OPTX,PTY,OPTY " + " -fH " + distFileName + " -nH "+ hName;
        fit2DTemplate(outName,args);
    }
}
void fitBackgroundShapesMVV(std::string name, const std::string& filename, const std::string& fittedName =""){
    std::string tempFile=filename+"_"+name+"_template.root";
    //Different name in case we want to fit on a different selection with some template
    if(fittedName.size())name = fittedName;
    std::string distFileName=filename+"_"+name+"_distributions.root";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string hName = name+"_"+l+"_"+p+"_"+h+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+l+"_"+p+"_"+h+"_template.root";
        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PT,OPT " + " -fH " + distFileName + " -nH "+ hName;
        fit1DTemplate(outName,args);
    }
}

void makeBKG1DShapes(const std::string& name, const std::string& filename, const std::string& catName, const std::string& fitName,  bool isW, CJSON* prevJSON, TFile* iF){
    FunctionParameterPlotter plotter;
    std::vector<std::unique_ptr<FunctionFitter>> fitters;
    auto setup1DFit = [&](const TH1* hbbH, double HHMass){
        std::string pF = isW ? "W" : "T";
        auto vN=[&](std::string var)->std::string{return var+pF;};
        fitters.emplace_back(new CBFunctionFitter(hbbH,false,pF,{ "MJJ"}));

        auto fitter = &* fitters.back();
        if(isW){
            fitter->setVar(vN("mean")     ,90,80,100);
            fitter->setVar(vN("sigma")     ,8,5,15);
            fitter->setVar(vN("alpha")     ,1.18,0.1,10);
            fitter->setVar(vN("alpha2")  ,0.9,0.1,5);
            fitter->setVar(vN("n")  ,5,1,6);
            fitter->setVar(vN("n2")  , 2,1,6);
            fitter->setConst(vN("n")  ,1);
            fitter->setConst(vN("n2")  ,1);
            fitter->w->var("MJJ")->setRange("fit",30,160);
        } else {
            fitter->setVar(vN("mean")     ,180,120,195);
            fitter->setVar(vN("sigma")     ,15.8,14,20);
            fitter->setVar(vN("alpha")     ,1 ,0.1,2);
            fitter->setVar(vN("alpha2")  ,1.5,0.5,3);
            fitter->setVar(vN("n")  ,5,1,6);
            fitter->setVar(vN("n2")  , 5,1,6);
            fitter->setConst(vN("n")  ,1);
            fitter->setConst(vN("n2")  ,1);
            fitter->w->var("MJJ")->setRange("fit",70,230);
        }
        if(prevJSON){
            fitter->setVar(vN("alpha")     ,prevJSON->evalFunc(vN("alpha")  ,HHMass) ,0.1,10);
            fitter->setConst(vN("alpha"),1);
            fitter->setVar(vN("alpha2")  ,prevJSON->evalFunc(vN("alpha2")  ,HHMass),0.1,10);
            fitter->setConst(vN("alpha2")  ,1);
        }
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::NumCPU(8)});
        fitter->fit({RooFit::SumW2Error(1),RooFit::Range("fit"),RooFit::NumCPU(8)});
    };


    std::string hName = name+"_"+catName;
    auto hbb_hh_H = TObjectHelper::getObject<TH2>(iF,hName+"_"+hbbMCS+"_"+hhMCS,false,false);
    if(hbb_hh_H==0) return;
    auto hhH  = projY(&*hbb_hh_H,hName+"_"+hhMCS);

    for(unsigned int iP = 0; iP < resPTBins.size() -1; ++iP){
        std::string ptName  = flt2Str(resPTBins[iP]) +"to"+flt2Str(resPTBins[iP+1]);
        auto hbbH = projX(&*hbb_hh_H,hName+"_"+ptName +"_"+hbbMCS,resPTBins[iP],resPTBins[iP+1]);
        double mean = getMean(&*hhH,resPTBins[iP],resPTBins[iP+1] );
        setup1DFit(&*hbbH,mean);
        plotter.addFit(&*fitters[iP],mean,ptName);
    }
    plotter.write(filename+"_"+name+"_"+catName+"_"+fitName+".root");
}


void makeResWMJJShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "fit1stIt";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_NONE] ) continue;
        const std::string catName = l+"_"+p+"_"+h;
        makeBKG1DShapes(name,filename,catName,fitName,true,0,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += " -minX 500 -maxX 3000 ";
        std::string jsonArgsStd = " -g meanW:laur2,sigmaW:laur2,alphaW:laur4,alpha2W:laur3,nW:pol0,n2W:pol0 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+  jsonArgsStd );
    }
}

void makeResWMJJShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "fit";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_NONE] ) continue;
        const std::string catName = l+"_"+p+"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_fit1stIt.json");
        oldJSON.fillFunctions("MH");
        makeBKG1DShapes(name,filename,catName,fitName,true,&oldJSON,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += " -minX 500 -maxX 3000 ";
        std::string jsonArgsStd = " -g meanW:laur2,sigmaW:laur2,alphaW:laur4,alpha2W:laur3,nW:pol0,n2W:pol0 ";

        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+jsonArgsStd);
        newJSON.replaceEntry("alphaW", oldJSON.getP("alphaW") );
        newJSON.replaceEntry("alpha2W", oldJSON.getP("alpha2W") );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}


void makeResTopMJJShapes1stIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "fit1stIt";

    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_NONE] ) continue;
        const std::string catName = l+"_"+p+"_"+h;
        makeBKG1DShapes(name,filename,catName,fitName,false,0,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += " -minX 500 -maxX 3500 ";
        std::string jsonArgsStd = " -g meanT:laur3,sigmaT:laur2,alphaT:laur4,alpha2T:laur3,nT:pol0,n2T:pol0 ";
        MakeJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+  jsonArgsStd );
    }
}

void makeResTopMJJShapes2ndIt(const std::string& name, const std::string& filename){
    auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_inclM_distributions.root");
    const std::string fitName = "fit";
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        if(l != lepSels[LEP_EMU] ) continue;
        if(p == purSels[PUR_I]) continue;
        if(h != hadSels[HAD_NONE] ) continue;
        const std::string catName = l+"_"+p+"_"+h;
        CJSON oldJSON(     filename+"_"+name+"_"+catName+"_fit1stIt.json");
        oldJSON.fillFunctions("MH");
        makeBKG1DShapes(name,filename,catName,fitName,false,&oldJSON,iF);

        std::string argsP1 = std::string("-i ")+ filename+"_"+name+"_"+catName+"_"+fitName+".root ";
        argsP1 += " -minX 500 -maxX 3000 ";
        std::string jsonArgsStd= " -g meanT:laur3,sigmaT:laur2,alphaT:laur4,alpha2T:laur3,nT:pol0,n2T:pol0 ";

        CJSON newJSON = getJSON(filename+"_"+name+"_"+catName+"_"+fitName+".json",argsP1+" "+jsonArgsStd);
        newJSON.replaceEntry("alphaT", oldJSON.getP("alphaT") );
        newJSON.replaceEntry("alpha2T", oldJSON.getP("alpha2T") );
        newJSON.write(filename+"_"+name+"_"+catName+"_"+fitName+".json");
    }
}

void convertFuncFitTo2DTemplate(const std::string& name, const std::string& filename,const std::string& funcParamPostfix){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        std::string jsonFile = filename+"_"+name+"_"+lepSels[LEP_EMU]+"_";
        jsonFile += name.find("w")!= std::string::npos ?purSels[PUR_LMT]:p;
        jsonFile+=std::string("_")+hadSels[HAD_NONE] +"_fit.json";

        CBFunctionFitter xFit(0,false,funcParamPostfix,{"MJJ"});
        xFit.w->var("MJJ")->setMin(minHbbMass);
        xFit.w->var("MJJ")->setMax(maxHbbMass);
        xFit.w->var("MJJ")->setBins(nHbbMassBins*10);
        xFit.w->var("MJJ")->setVal((minHbbMass+maxHbbMass)/2.);

        CJSON json(jsonFile);
        json.fillFunctions("MH");
        auto pn = [&](const std::string& v) ->std::string{return v + funcParamPostfix;};
        auto * iF =  TObjectHelper::getFile(filename + "_"+name+"_"+l+"_"+p+"_"+h+"_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH1>(iF,"histo",false,false);
        if(hh_H==0) continue;
        oF->cd();
        TH2 * h_2D = new TH2F((name+"_"+l+"_"+p+"_"+h).c_str(),(std::string(";")+hbbMCS.title+";"+hhMCS.title).c_str(),nHbbMassBins,minHbbMass,maxHbbMass,nHHMassBins,minHHMass,maxHHMass);
        for(int iY = 1; iY <= hh_H->GetNbinsX(); ++iY){
            double yNorm = hh_H->GetBinContent(iY);
            double hhV = hh_H->GetBinCenter(iY);

            xFit.setVar(pn("mean"  ),json.evalFunc(pn("mean"  ),hhV));
            xFit.setVar(pn("sigma" ),json.evalFunc(pn("sigma" ),hhV));
            xFit.setVar(pn("alpha" ),json.evalFunc(pn("alpha" ),hhV));
            xFit.setVar(pn("alpha2"),json.evalFunc(pn("alpha2"),hhV));
            xFit.setVar(pn("n"     ),json.evalFunc(pn("n"     ),hhV));
            xFit.setVar(pn("n2"    ),json.evalFunc(pn("n2"    ),hhV));
            auto xHist = xFit.pdf1D(name+"_"+l+"_"+p+"_"+h+ASTypes::int2Str(iY)+"_hbb");
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
    for(const auto& l :lepSels) for(const auto& p :purSels) for(const auto& h :hadSels){
        auto * iF =  TObjectHelper::getFile(filename+"_"+name+"_"+l+"_" +p +"_"+h+"_2D_template.root");
        if(iF==0) continue;
        auto hh_H = TObjectHelper::getObject<TH2>(iF,"histo__x_y",false,false);
        if(hh_H==0) continue;
        oF->cd();
        hh_H->SetXTitle(hbbMCS.title.c_str());
        hh_H->SetYTitle(hhMCS.title.c_str());
        hh_H->Write((name+"_"+l+"_"+p+"_"+h).c_str());
        iF->Close();
    }
    oF->Close();
}

void go(BKGModels modelToDo, std::string treeDir) {
    std::string filename = hhFilename;
    std::string treeArea = treeDir + "/betrees_bkg.root";
    std::string signalTrees = treeDir + "/out_radion_hh_bbinc_mXXX_0.root";

    if(modelToDo == BKG_QG)
    {
        std::string name = bkgSels[BKG_QG];
        std::string genSel = bkgSels[BKG_QG].cut + "&&"+ aQCD.cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+ hadSels[HAD_LTMB].cut;
        makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        //        //MVV
        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,3,0.5,1);//x = hh
        //                makeBackgroundShapesMVVConditional(name+"_xs_0p75_xc_2_ys_0p75_yc_2,filename,treeArea,baseSel,0.75,2,0.75,2);//old
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hhInclRange.cut,true);
        //
        //        //MJJ
        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel+"&&"+hhRange.cut);
        mergeBackgroundShapes(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVVConditional(name,filename);

        //        compile2DTemplatesForDebug(name,filename);
    }

    if(modelToDo == BKG_LOSTTW){
        std::string name = bkgSels[BKG_LOSTTW];
        std::string genSel = bkgSels[BKG_LOSTTW].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LTMB].cut;
        makeDetectorParam(name,filename,treeArea, genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        //MVV
        makeBackgroundShapesMVVConditional(name,filename,treeArea,baseSel,0.75,8,0.5,2);//x = hh
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hhInclRange.cut,true);
        //MJJ
        makeBackgroundShapesMJJAdaKernel(name,filename,treeArea,baseSel+"&&"+hhRange.cut);
        mergeBackgroundShapes(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVVConditional(name,filename);

        compile2DTemplatesForDebug(name,filename);
    }

    if(modelToDo == BKG_MW){
        std::string name = bkgSels[BKG_MW];
        std::string genSel = bkgSels[BKG_MW].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LB].cut;
        //MVV
        makeDetectorParam(name,filename,treeArea,genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,baseSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVV(name,filename);

        //MJJ
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeResWMJJShapes1stIt(name,filename);
        makeResWMJJShapes2ndIt(name,filename);
        convertFuncFitTo2DTemplate(name,filename,"W");
    }

    if(modelToDo == BKG_MT){
        std::string name = bkgSels[BKG_MT];
        std::string genSel = bkgSels[BKG_MT].cut;
        std::string baseSel = genSel + "&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&" +hadSels[HAD_LB].cut;
        //MVV
        makeDetectorParam(name,filename,treeArea,genSel + "&&"+ hhRange.cut+"&&"+hbbRange.cut+"&&"+lepSels[LEP_EMU].cut+"&&"+purSels[PUR_LMT].cut+"&&"+hadSels[HAD_NONE].cut);
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,baseSel+"&&"+hbbRange.cut);
        cutMVVTemplate(name,filename);
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);
        fitBackgroundShapesMVV(name,filename);
        //
        //        //MJJ
        makeFittingDistributions(name,filename,treeArea,genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeResTopMJJShapes1stIt(name,filename);
        makeResTopMJJShapes2ndIt(name,filename);
        convertFuncFitTo2DTemplate(name,filename,"T");
    }

}
#endif

void makeBKGInputs(int bkgToDo = BKG_QG, std::string treeDir = "trees/"){
    go(static_cast<BKGModels>(bkgToDo),treeDir);
}
