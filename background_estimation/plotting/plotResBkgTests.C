#include "plotTestHelper.h"
#include "TH1.h"
#include "TH2.h"

std::vector<TObject*> writeables;


void testHHKern(std::string name, std::string filename) {
    bool withRatio = true;
    TFile *f = new TFile((filename + "_"+name +"_incl_template.root").c_str(),"read");
    std::vector<TH1*> hs;
    std::vector<std::string> hNs;
    TH1* h = 0;
    f->GetObject("histo_data",h);hs.push_back(h);hNs.push_back("MC");
    f->GetObject("histo",h);     hs.push_back(h);hNs.push_back("KDE");
        f->GetObject("histo_KDE",h);     hs.push_back(h);hNs.push_back("KDE w/o expo. tail smoothing");

    int binL = hs[0]->FindFixBin(minHHMass);
    int binH = hs[0]->FindFixBin(maxHHMass);
    for(unsigned int iH = 1; iH < hs.size(); ++iH) hs[iH]->Scale(hs[0]->Integral(binL,binH)/hs[iH]->Integral(binL,binH));

    Plotter * p = new Plotter();
    Plotter * pf = new Plotter();
    for(unsigned int iH = 0; iH < hs.size(); ++iH){
        if(iH == 0){
            p->addHist(hs[iH],hNs[iH].c_str());
            pf->addHist(hs[iH],hNs[iH].c_str());
        }
        else  {
            p->addHistLine(hs[iH],hNs[iH].c_str());
            pf->addHistLine(hs[iH],hNs[iH].c_str());
        }
    }
    auto setupPlotter = [&](Plotter * p, std::string name, double rebin =-1){
        p->setMinMax(.0001,(rebin < 0 ? 1.0 : rebin) * hs[0]->Integral()/4);
        p->setUnderflow(false);
        p->setOverflow(false);
        p->setBotMinMax(0,2);
        p->setYTitle("N. of events");
        p->setXTitle(hhMCS.title.c_str());
        if(rebin > 0) p->rebin(rebin);
        if(withRatio){
            auto * c = p->drawSplitRatio(1,"stack",false,false,name);
            c->GetPad(1)->SetLogy();
            c->GetPad(1)->Update();
            writeables.push_back(c);
        } else {
            auto * c = p->draw(false,name);
            c->SetLogy();
            c->Update();
            writeables.push_back(c);
        }
    };

    setupPlotter(p,name + "_HHKDE_C",5);
    setupPlotter(pf,name + "_HHKDE_F");
}

void testHHPDFFits(std::string name, std::string filename) {
    std::vector<std::string> sels = {"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"};

    TFile * fd = new TFile((filename+"_"+name+"_distributions.root").c_str(),"read");
    TFile * fo = new TFile((filename+"_"+name+"_template.root").c_str(),"read");
    TH1 * hof = 0;
    fo->GetObject("histo",hof);

    for(const auto& s : sels){
        TH1 * hd = 0;
        fd->GetObject(  (name+"_"+s+"_hhMass").c_str(),hd);
        if(hd==0) continue;
        std::vector<TH1*> hs;
        std::vector<TString> hNs;
        hs.push_back((TH1*)hof->Clone());
        hNs.push_back("Baseline KDE before MC fit");

        TFile * ff = new TFile((filename+"_"+name+"_"+s+"_template.root").c_str(),"read");
        TH1* h = 0;
        ff->GetObject("histo",h);hs.push_back(h);hNs.push_back("Nominal SR PDF");

        for(unsigned int iH = 0; iH < hs.size(); ++iH) hs[iH]->Scale(hd->Integral()/hs[iH]->Integral());

        Plotter * p = new Plotter();
        p->addHist(hd,"SR MC");
        for(unsigned int iH = 0; iH < hs.size(); ++iH){
            TH1 * h1D = hs[iH];
            for(int iX = 1; iX <= h1D->GetNbinsX(); ++iX)h1D->SetBinError(iX,0);
            p->addHist(h1D,hNs[iH],-1,1,4,20,1,false,true,false,"E");
        }
        p->setMinMax(.0001,hs[0]->Integral());
        p->setUnderflow(false);
        p->setOverflow(false);
        p->rebin(8);
        p->setXTitle("HH mass [GeV]");
        p->setYTitle("N. of events");
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
            TFile * f = new TFile(outName.c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
            }
            f->Close();
        }
    }
    std::string outName;
};






void plotResBkgTests(int step = 0,bool doMT = true, bool doSR = true, std::string outName = ""){
    //    hhFilename +="_CR";
    std:: string inName = doSR ? "bkgInputs" : "bkgInputsCR";
    std::string filename = inName +"/"+hhFilename;


    CutStr mod = bkgSels [doMT ? BKG_MT : BKG_MW];
    if(outName.size()){
        if(step >= 5) outName += std::string("/All")  +  (doSR ? "" : "CR");
        else outName += std::string("/") + mod +  (doSR ? "" : "CR");
    }

    switch(step){
    case 0:
        if(outName.size()) outName += "_MVV_temp.root";
        testHHKern(mod,filename);
        break;
    case 1:
        if(outName.size()) outName += "_MVV_fit.root";
        testHHPDFFits(mod,filename);
        break;
    case 2:
        if(outName.size()) outName += "_MJJ_fit1stIt.root";
        if(doMT) writeables = testBKG1DFits(mod,filename,"","fit1stIt",{"emu_L_I_none","emu_M_I_none","emu_T_I_none"});
        else  writeables = testBKG1DFits(mod,filename,"","fit1stIt",{"emu_LMT_I_none"});
        break;
    case 3:
        if(outName.size()) outName += "_MJJ_fit.root";
        if(doMT) writeables = testBKG1DFits(mod,filename,"","fit",{"emu_L_I_none","emu_M_I_none","emu_T_I_none"});
        else  writeables = testBKG1DFits(mod,filename,"","fit",{"emu_LMT_I_none"});
        break;
    case 4:
        if(outName.size()) outName += "_2DComp.root";
        writeables = test2DModel({mod},filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});
        break;
    case 5:
        if(outName.size()) outName += "_2DComp.root";
        writeables = test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },
              filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});

    }

    Dummy d(outName);


//
//    if(step == 0)testHHKern(bkgSels[BKG_MW],filename);
//    if(step == 1)testHHPDFFits(bkgSels[BKG_MW],filename);
//    if(step == 2)testBKG1DFits(bkgSels[BKG_MW],filename,"","fit1stIt",{"emu_LMT_I_none"});
//    if(step == 3)testBKG1DFits(bkgSels[BKG_MW],filename,"","fit",{"emu_LMT_I_none"});
//    if(step == 4)test2DModel({bkgSels[BKG_MW] },filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});
//
//    if(step == 5)testHHKern(bkgSels[BKG_MT],filename);
//    if(step == 6)testHHPDFFits(bkgSels[BKG_MT],filename);
//    if(step == 7)testBKG1DFits(bkgSels[BKG_MT],filename,"","fit1stIt",{"emu_L_I_none","emu_M_I_none","emu_T_I_none"});
//    if(step == 8)testBKG1DFits(bkgSels[BKG_MT],filename,"","fit",{"emu_L_I_none","emu_M_I_none","emu_T_I_none"});
//    if(step == 9)test2DModel({bkgSels[BKG_MT] },filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});
//    if(step == 10)test2DModel({bkgSels[BKG_QG],bkgSels[BKG_LOSTTW],bkgSels[BKG_MW],bkgSels[BKG_MT] },filename,{"e_L_LP_full","mu_L_LP_full","e_M_LP_full","mu_M_LP_full","e_T_LP_full","mu_T_LP_full","e_L_HP_full","mu_L_HP_full","e_M_HP_full","mu_M_HP_full","e_T_HP_full","mu_T_HP_full"},{700,4000});
}
