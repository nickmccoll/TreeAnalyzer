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


