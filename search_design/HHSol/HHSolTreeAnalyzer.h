
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "DataFormats/interface/Momentum.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "TPRegexp.h"
using namespace TAna;
using namespace ASTypes;

class HHSolTreeAnalyzer : public BaseTreeAnalyzer {
public:
    HHSolTreeAnalyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) :
    BaseTreeAnalyzer(fileName,treeName,treeInt,randSeed) {

        TPRegexp r1(".*radion.*");
        if(r1.Match(fileName)) isRadionSMP=true;


    }

    virtual void loadVariables() override {
        setBranch("","process"       ,process        ,true);
        setBranch("","sampParam"     ,sampParam      ,true);
        setBranch("","dhType"        ,dhType         ,true);
        setBranch("","hbbCat"        ,hbbCat         ,true);
        setBranch("","hbbDeepAK8"    ,hbbDeepAK8     ,true);
        setBranch("","isMuon"        ,isMuon         ,true);
        setBranch("","weight"        ,weight_        ,true);
        setBranch("","nAK4Btags"     ,nAK4Btags      ,true);
        setBranch("","hh_orig"       ,hh_orig        ,true);
        setBranch("","hh_chi2"       ,hh_chi2        ,true);
        setBranch("","md"            ,md             ,true);
        setBranch("","chi2"          ,chi2           ,true);
        setBranch("","wqqDR"         ,wqqDR          ,false);
        setBranch("","qqJet_pt"      ,qqJet_pt       ,true);
        setBranch("","qqJet_eta"     ,qqJet_eta      ,true);
        setBranch("","qqJet_phi"     ,qqJet_phi      ,true);
        setBranch("","qqJet_mass"    ,qqJet_mass     ,true);
        setBranch("","qqJet_SDmass"  ,qqJet_SDmass   ,true);
        setBranch("","qqJet_t2ot1"   ,qqJet_t2ot1    ,true);
        setBranch("","lep_pt"        ,lep_pt         ,true);
        setBranch("","lep_eta"       ,lep_eta        ,true);
        setBranch("","lep_phi"       ,lep_phi        ,true);
        setBranch("","met_pt"        ,met_pt         ,true);
        setBranch("","met_phi"       ,met_phi        ,true);
        setBranch("","bbJet_pt"      ,bbJet_pt       ,true);
        setBranch("","bbJet_eta"     ,bbJet_eta      ,true);
        setBranch("","bbJet_phi"     ,bbJet_phi      ,true);
        setBranch("","bbJet_mass"    ,bbJet_mass     ,true);
        setBranch("","bbJet_SDmass"  ,bbJet_SDmass   ,true);
        setBranch("","true_neut_pt"  ,true_neut_pt   ,true);
        setBranch("","true_neut_eta" ,true_neut_eta  ,true);
        setBranch("","true_neut_phi" ,true_neut_phi  ,true);
        setBranch("","true_neut_mass",true_neut_mass ,true);
        setBranch("","true_lep_pt"   ,true_lep_pt    ,true);
        setBranch("","true_lep_eta"  ,true_lep_eta   ,true);
        setBranch("","true_lep_phi"  ,true_lep_phi   ,true);
        setBranch("","true_jet_pt"   ,true_jet_pt    ,true);
        setBranch("","true_jet_eta"  ,true_jet_eta   ,true);
        setBranch("","true_jet_phi"  ,true_jet_phi   ,true);
        setBranch("","true_jet_mass" ,true_jet_mass  ,true);
    }

    virtual bool runEvent() override {


        if (*process == FillerConstants::SIGNAL)
            smpName = TString::Format("%s_m%i",isRadionSMP ? "radion" : "graviton",int(*sampParam));
        else smpName = FillerConstants::MCProcessNames[*process];

        qqJet.setP4(*qqJet_pt,*qqJet_eta,*qqJet_phi,*qqJet_mass);
        lepton.setP4(*lep_pt,*lep_eta,*lep_phi,float(0));
        met.setP4(*met_pt,float(0),*met_phi,float(0));
        bbJet.setP4(*bbJet_pt,*bbJet_eta,*bbJet_phi,*bbJet_mass);
        true_neut.setP4(*true_neut_pt,*true_neut_eta,*true_neut_phi,*true_neut_mass);
        true_lep.setP4(*true_lep_pt,*true_lep_eta,*true_lep_phi,float(0));
        true_jet.setP4(*true_jet_pt,*true_jet_eta,*true_jet_phi,*true_jet_mass);
        true_W.p4() = true_lep.p4() + true_neut.p4();
        true_H.p4() = true_W.p4() + true_jet.p4();
        true_hwwMag = true_H.pt();

        auto getPerp =[&](const MomentumF& mom)->double{
            return mom.px()*hwwPerpNormX+mom.py()*hwwPerpNormY;
        };
        auto getPar =[&](const MomentumF& mom)->double{
            return mom.px()*hwwParNormX+mom.py()*hwwParNormY;
        };

        hwwParX      = qqJet.px() + lepton.px() + met.px();
        hwwParY      = qqJet.py() + lepton.py() + met.py();
        hwwMag       = std::sqrt(hwwParX*hwwParX+hwwParY*hwwParY);
        hwwParNormX  = hwwParX/hwwMag;
        hwwParNormY  = hwwParY/hwwMag;
        hwwPerpNormX = -1*hwwParNormY;
        hwwPerpNormY = hwwParNormX;
        metPerp      = getPerp(met);
        metPar       = getPar(met);
        neutPerp     = getPerp(true_neut);
        neutPar      = getPar(true_neut);
        extraMetPerp = metPerp-neutPerp;
        extraMetPar  = metPar-neutPar;

        isVirtualWqq = true_jet.mass() < true_W.mass();
        weight = *weight_;

        if(*hbbDeepAK8 < 0 ) sanitizedDeepAK8Tag = 0;
        else if(*hbbDeepAK8 > 1) sanitizedDeepAK8Tag = 1;
        else sanitizedDeepAK8Tag = *hbbDeepAK8;

        return true;
    }


    bool isSignal() const {return *process == FillerConstants::SIGNAL;}


    //Reader info
    rd_size8 process        ;
    rd_int   sampParam      ;
    rd_size8 dhType         ;
    rd_size8 hbbCat         ;
    rd_float hbbDeepAK8     ;
    rd_size8 isMuon         ;
    rd_float weight_        ;
    rd_size8 nAK4Btags      ;
    rd_float hh_orig        ;
    rd_float hh_chi2        ;
    rd_float md             ;
    rd_float chi2           ;
    rd_float wqqDR          ;
    rd_float qqJet_pt       ;
    rd_float qqJet_eta      ;
    rd_float qqJet_phi      ;
    rd_float qqJet_mass     ;
    rd_float qqJet_SDmass   ;
    rd_float qqJet_t2ot1    ;
    rd_float lep_pt         ;
    rd_float lep_eta        ;
    rd_float lep_phi        ;
    rd_float met_pt         ;
    rd_float met_phi        ;
    rd_float bbJet_pt       ;
    rd_float bbJet_eta      ;
    rd_float bbJet_phi      ;
    rd_float bbJet_mass     ;
    rd_float bbJet_SDmass   ;
    rd_float true_neut_pt   ;
    rd_float true_neut_eta  ;
    rd_float true_neut_phi  ;
    rd_float true_neut_mass ;
    rd_float true_lep_pt    ;
    rd_float true_lep_eta   ;
    rd_float true_lep_phi   ;
    rd_float true_jet_pt    ;
    rd_float true_jet_eta   ;
    rd_float true_jet_phi   ;
    rd_float true_jet_mass  ;

    //Momenta
    MomentumF qqJet;
    MomentumF lepton;
    MomentumF met;
    MomentumF bbJet;
    MomentumF true_neut;
    MomentumF true_lep;
    MomentumF true_jet;
    MomentumF true_W;
    MomentumF true_H;
    float true_hwwMag      =0;

    float weight      =0;
    float hwwParX     =0;
    float hwwParY     =0;
    float hwwMag      =0;
    float hwwParNormX =0;
    float hwwParNormY =0;
    float hwwPerpNormX=0;
    float hwwPerpNormY=0;
    float metPerp     =0;
    float metPar      =0;
    float neutPerp    =0;
    float neutPar     =0;
    float extraMetPerp=0;
    float extraMetPar =0;

    bool isVirtualWqq = false;

    TString smpName;
    bool isRadionSMP = false;

    double sanitizedDeepAK8Tag = 0;



};

#endif
