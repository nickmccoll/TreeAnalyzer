
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
#include "../predTools/DileptonCutConstants.h"

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

std::string getMVVExpoMinMaxStr(std::string sample){
    float eMin, eMax;
    if(strFind(sample,llBkgSels[BKG_MISB])){
        eMin = 2000;
        eMax = 3000;
    } else if(strFind(sample,llBkgSels[BKG_REALB])){
        eMin = 2000;
        eMax = 3000;
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


void makeBackgroundShapesMVVAdaKernel(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& baseSel="1.0", float khxs = 1,float khxc = 5){

    std::string resFile=filename+"_"+name+"_detectorResponse.root";

    DilepCatIterator ci;
    while(ci.getBin()){
        if(strFind(name,llBkgSels[BKG_MISB])){
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(SEL_RPhiB)) continue;
        }
        if(strFind(name,llBkgSels[BKG_REALB])){
            if(!ci.is(LEP_INCL )) continue;
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(SEL_RPhiB)) continue;
        }

        std::string tempFile=filename+"_"+name+"_"+ci.name()+"_MVV_incl_template.root";
        std::string cut =  std::string("(")+baseSel+"&&"+ci.cut()+")";
        std::string weight = std::string("(") + nomW.cut+")";

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
        const std::string inputFile, const std::string baseSel="1.0", bool xIsCond = false,
        float khxs = 1,float khxc = 5,float khys = 1,float khyc = 5) {

    DilepCatIterator ci;
    while(ci.getBin()){
        if(!ci.is(LEP_INCL )) continue;
        if(!ci.is(BTAG_LMT)) continue;
        if(!ci.is(SEL_RPhiB  )) continue;

        const std::string tempFile=filename+"_"+name+"_"+ci.name()+"_2D_cond_template.root";
        std::string cut =  std::string("(")+baseSel+"&&"+ci.cut()+")";

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

    DilepCatIterator ci;
    while(ci.getBin()){
        if(strFind(name,bkgSels[BKG_MISB])){
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(SEL_RPhiB)) continue;
        } else {
            if(!ci.is(LEP_INCL)) continue;
            if(!ci.is(BTAG_LMT)) continue;
            if(!ci.is(SEL_RPhiB)) continue;
        }
        std::string inFile1D=
                filename+"_"+name+"_"+ ci.name()+"_"+(xIsCond? "MVV":"MJJ") +"_incl_template.root";
        std::string inFile2D=filename+"_"+name+"_"
                +dilepCats[LEP_INCL] +"_"+btagCats[BTAG_LMT]+"_"+selCuts[SEL_RPhiB]
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



void makeFittingDistributions(const std::string& name, const std::string& filename,
        const std::string inputFile, const std::string& cut="1.0", bool doIncl = true,
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
    DilepCatIterator ci;
    while(ci.getBin()){

        sels.emplace_back(ci.name(),ci.cut());

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

    DilepCatIterator ci;
    while(ci.getBin()){
        std::string hName = name+"_"+ci.name()+"_"+hbbMCS+"_"+hhMCS;
        std::string outName = filename + "_"+name+"_"+ci.name()+"_2D_template.root";

        std::string tempFile= filename+"_"+name+"_"
                +dilepCats[LEP_INCL]+"_"+btagCats[BTAG_LMT]+"_"+selCuts[SEL_RPhiB]
                +"_2D_merged_template.root";
        if(strFind(name,bkgSels[BKG_MISB])){
            tempFile=filename+"_"+name+"_"
                    +dilepCats[ci.l]+"_"+btagCats[BTAG_LMT]+"_"+selCuts[SEL_RPhiB]
                    +"_2D_merged_template.root";
        }

        std::string args = std::string("-v ") + "-fT "+ tempFile+" -nT histo -s PTX,OPTX,PTY,OPTY "
                +" -sA PT2Y " + " -fH " + distFileName + " -nH "+ hName;
        if(xIsCond) args+=" -xCy ";

        fit2DTemplate(outName,args);

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


void compile2DTemplatesForDebug(const std::string& name, const std::string& filename){
    TFile *oF = new TFile((filename + "_"+name+"_2D_template_debug.root").c_str(),"recreate");
    DilepCatIterator ci;
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
    DilepCatIterator ci;
    while(ci.getBin()){
//        if(!ci.is(SEL_FULL)) continue;
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
    if(modelToDo == 5){
        return;
    }

    //Turn on TTBar scaling
    nomW.cut = nomW.cut+"*" +getTTBarSF("../supportInputs/HHlnujj") ;
    std::string treeAreaIncl = treeDir + "../betrees_mc.root";

    if(modelToDo == BKG_MISB)
    {
        std::string name = llBkgSels[BKG_MISB];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = llBkgSels[BKG_MISB].cut + "&&"+ aQCD.cut;
        makeFittingDistributions(name,filename,treeArea,
                hhInclRange.cut+"&&"+hbbInclRange.cut,true,
                {{llBkgSels[BKG_MISB],genSel}}
        );

        makeFittingDistributions(name,filename,treeArea,
                hhRange.cut+"&&"+hbbRange.cut,false,
                {{llBkgSels[BKG_MISB],genSel}}
        );

        makeBackgroundShapes2DConditional(name,filename,treeArea,
                genSel,true,0.4,1,0.75,3);//P(hbb|hh)

        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut,1,3);

        mergeBackgroundShapes(name,filename,true);

        fitBackgroundShapes2DConditional(name,filename,true);

        compile2DTemplatesForDebug(name,filename);

    }

    if(modelToDo == BKG_REALB){
        std::string name = llBkgSels[BKG_REALB];
        std::string treeArea = treeDir + "/betrees_" +name+".root";
        std::string genSel = llBkgSels[BKG_REALB].cut;
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhInclRange.cut+"&&"+hbbInclRange.cut,true);
        makeFittingDistributions(name,filename,treeArea,
                genSel+ "&&"+ hhRange.cut+"&&"+hbbRange.cut,false);


        makeBackgroundShapes2DConditional(name,filename,treeArea,
                genSel,true,0.5,3,0.75,5);//P(hbb|hh) 0.5,2,0.75,8 old values
        makeBackgroundShapesMVVAdaKernel(name,filename,treeArea,genSel+"&&"+hbbRange.cut,1,4);
        mergeBackgroundShapes(name,filename,true);
        fitBackgroundShapes2DConditional(name,filename,true);

        compile2DTemplatesForDebug(name,filename);
    }

    //Make pseudo data
    if(int(modelToDo) == 4){//4
        makePseudoData("pd",filename,1);
    }

}
#endif

void makeBKGInputs_2l(int bkgToDo = BKG_MISB, std::string treeDir = "../bkgCompLMT/"){
    go(bkgToDo,treeDir);
}
