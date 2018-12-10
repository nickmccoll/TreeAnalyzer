#include "plotTestHelper.h"
#include "TH1.h"
#include "TH2.h"

std::vector<TObject*> writeables;


void testHHPDFFits(std::string name, std::string filename, const std::vector<std::string> sels) {

    TFile * fd = new TFile((filename+"_"+name+"_distributions.root").c_str(),"read");

    for(const auto& s : sels){
        TH1 * hd = 0;
        fd->GetObject(  (name+"_"+s+"_hhMass").c_str(),hd);
        if(hd==0) continue;
        std::vector<TH1*> hs;
        std::vector<TString> hNs;


        TFile * ff = new TFile((filename+"_"+name+"_"+s+"_MVV_template.root").c_str(),"read");
        TH1* h = 0;
        ff->GetObject("histo",h);hs.push_back(h);hNs.push_back("SR Template: "+ s);
        ff->GetObject("originalPDF",h);hs.push_back(h);hNs.push_back("Baseline KDE before MC fit");

        for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(hd->Integral()/hs[iH]->Integral());

        Plotter * p = new Plotter();
        p->addHist(hd,"SR MC");
        for(unsigned int iH = 0; iH < hs.size(); ++iH){
            TH1 * h1D = hs[iH];
            for(int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
            p->addHist(h1D,hNs[iH],-1,1,4,20,1,false,true,false,"E");
        }
        p->setMinMax(.01,hs[0]->Integral());
        p->setUnderflow(false);
        p->setOverflow(false);
        p->rebin(8);
        p->setXTitle(hhMCS.title);
        p->setYTitle("N. of events");
        p->addText(getCategoryLabel(s).c_str(),0.15,0.88,0.03);

        p->setBotMinMax(0,2);
        // auto * c = p->drawSplitRatio(0,"stack",false,false,s);
        // c->GetPad(1)->SetLogy();
        // c->GetPad(1)->Update();
        auto * c = p->draw(false,(name + "_"+s+"_HHFIT").c_str());
        c->SetLogy();
        c->Update();
        writeables.push_back(c);

    }

}



std::vector<TObject*> testBKG1DFits(std::string name, std::string filename, std::string varName, std::string fitName, const std::vector<std::string>& sels){
    std::vector<std::string> canNames;
    for(unsigned int iM = 0; iM+1 < resPTBins.size(); ++iM){
        canNames.push_back(std::string("can_") +flt2Str(resPTBins[iM]) + "to"+flt2Str(resPTBins[iM+1]));
    }
    return test1DFits(name,filename,varName,fitName,sels,canNames);
}


class Dummy {
public:
    Dummy(const std::string& outName = "") : outName(outName) {};
    ~Dummy() {
        if(outName.size()){
            TFile * f = new TFile((outName+".root").c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
                w->Print((outName +"_"+w->GetTitle() +".pdf").c_str());
            }
            f->Close();
        }
    }
    std::string outName;
};






void plotResBkgTests(int step = 0, bool doMT = true, int inreg = REG_SR,  std::string outName = ""){
    REGION reg = REGION(inreg);

    std:: string inName =  "bkgInputs" ;
    auto srList = getSRList(reg);
    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
    }
    else if(reg == REG_QGCR){
        inName =  "bkgInputsQGCR";
        hhFilename +="_QGCR";
    }
    std::string filename = inName +"/"+hhFilename;

    CutStr mod = bkgSels [doMT ? BKG_MT : BKG_MW];
    if(outName.size()){
        outName += std::string("/") + mod;
        if(reg == REG_TOPCR) outName +=  "_TopCR";
        else if(reg == REG_QGCR) outName +=  "_QGCR";
    }

    std::vector<std::string> mtMJJBinning = {"emu_LMT_I_none"};
    if(reg != REG_QGCR && doMT) mtMJJBinning ={"emu_L_I_none","emu_M_I_none","emu_T_I_none"};

    switch(step){
    case 0:
        if(outName.size()) outName += "_MVV_temp";
        writeables = test1DKern(mod,filename,"MVV",{"emu_LMT_I_lt"});
        break;
    case 1:
        if(outName.size()) outName += "_MVV_fit";
        testHHPDFFits(mod,filename,srList);
        break;
    case 2:
        if(outName.size()) outName += "_MJJ_fit1stIt";
        writeables = testBKG1DFits(mod,filename,"","MJJ_fit1stIt",mtMJJBinning);
        break;
    case 3:
        if(outName.size()) outName += "_MJJ_fit";
        writeables = testBKG1DFits(mod,filename,"","MJJ_fit",mtMJJBinning);
        break;
    case 4:
        if(outName.size()) outName += "_MJJ_SFFit";
        writeables = test2DFits(mod,filename,srList,{700,4000},true,2,"MJJ_SFFit.json.root");
        break;
    case 5:
        if(outName.size()) outName += "_2DComp";
        writeables = test2DModel({mod},filename,srList,{700,4000});
        break;
    case 6:
        if(outName.size()) outName += "_all_2DComp";
        writeables = test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },
              filename,{"emu_LMT_I_full"},{700,4000});
        writeables = test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },
              filename,srList,{30,210},false);

    }

    Dummy d(outName);
}
