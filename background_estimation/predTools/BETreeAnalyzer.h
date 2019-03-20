
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "TPRegexp.h"

using namespace TAna;
using namespace ASTypes;

class BETreeAnalyzer : public BaseTreeAnalyzer {
public:
    BETreeAnalyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
        auto match = r1.MatchS(fileName);
        const Int_t nrSubStr = match->GetLast()+1;
        if(nrSubStr>1){
            signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
        }
    }

    void loadVariables() override {
        if(isRealData()){
            setBranchAddress("","dataset"  ,&dataset,true);
            setBranchAddress("","dataRun"  ,&dataRun,true);
        } else {
            setBranchAddress("","process",&process,true);
            setBranchAddress("","dhType" ,&dhType ,true);
            setBranchAddress("","xsec"   ,&xsec   ,true);
            setBranchAddress("","trig_N" ,&trig_N ,true);
            setBranchAddress("","pu_N"   ,&pu_N   ,true);
            setBranchAddress("","lep_N"  ,&lep_N  ,true);
            setBranchAddress("","btag_N" ,&btag_N ,true);
        }
        setBranchAddress("","ht"        ,&ht       ,true);
        setBranchAddress("","isMuon"    ,&isMuon   ,true);
        setBranchAddress("","hbbMass"   ,&hbbMass  ,true);
        setBranchAddress("","hbbPT"     ,&hbbPT    ,true);
        setBranchAddress("","hbbNSJs"   ,&hbbNSJs  ,true);
        setBranchAddress("","hbbCSVCat" ,&hbbCSVCat,true);
        setBranchAddress("","hbbTau2o1" ,&hbbTau2o1,true);
        setBranchAddress("","hhMass"    ,&hhMass   ,true);
        setBranchAddress("","wlnuDR"    ,&wlnuDR   ,true);
        setBranchAddress("","wwDM"      ,&wwDM     ,true);
        setBranchAddress("","hwwPT"     ,&hwwPT    ,true);
        setBranchAddress("","wjjCSVCat" ,&wjjCSVCat,true);
        setBranchAddress("","wjjTau2o1" ,&wjjTau2o1,true);
        setBranchAddress("","wjjMass"   ,&wjjMass  ,true);
        setBranchAddress("","wjjPT"     ,&wjjPT    ,true);
        setBranchAddress("","wjjNSJs"   ,&wjjNSJs  ,true);
        setBranchAddress("","wlnuPT"    ,&wlnuPT   ,true);
        setBranchAddress("","nAK4Btags" ,&nAK4Btags,true);
        setBranchAddress("","minBtagMT" ,&minBtagMT,true);

        if(!isRealData()){
            setBranchAddress("","hbbGenPT"   ,&hbbGenPT  ,true);
            setBranchAddress("","hbbGenMass" ,&hbbGenMass,true);
            setBranchAddress("","hbbWQuark"  ,&hbbWQuark ,true);
            setBranchAddress("","hbbWEQuark" ,&hbbWEQuark,true);
            setBranchAddress("","genhhMass"  ,&genhhMass ,true);
            setBranchAddress("","genhhMass2" ,&genhhMass2,true);

        }

    }

    virtual bool runEvent() override {
        if(isRealData()) smpName = "data";
        else if (process == FillerConstants::SIGNAL) smpName = TString::Format("m%i",signal_mass);
        else smpName = FillerConstants::MCProcessNames[process];
        weight = isRealData() ? 1.0 : xsec*trig_N*pu_N*lep_N*btag_N;
        return true;
    }
    bool isSignal() const {return process == FillerConstants::SIGNAL;}


    //Shortcuts
    float           weight     =0; //std weight after all corrections
    int             signal_mass=0;
    TString         smpName  = "";



    //Event information and weights
    size16 process    = 0;
    size8  dhType     = 0;
    size16 dataset    = 0;
    size16 dataRun    = 0;
    float  xsec       = 0;
    float  trig_N     = 0;
    float  pu_N       = 0;
    float  lep_N      = 0;
    float  btag_N     = 0;
    float  ht         = 0;
    size8  isMuon     = 0;
    float  hbbMass    = 0;
    float  hbbPT      = 0;
    size8  hbbNSJs    =0;
    size8  hbbCSVCat  = 0;
    float  hbbTau2o1  = 0;
    float  hhMass     = 0;
    float  wlnuDR     = 0;
    float  wwDM       = 0;
    float  hwwPT      = 0;
    size8  wjjCSVCat  = 0;
    float  wjjTau2o1  = 0;
    float  wjjMass    = 0;
    float  wjjPT      = 0;
    size8  wjjNSJs    =0;
    float  wlnuPT     = 0;
    size8  nAK4Btags  = 0;
    float  minBtagMT  = 0;
    float  hbbGenPT    =0;
    float  hbbGenMass  =0;
    size8  hbbWQuark   =0;
    size8  hbbWEQuark   =0;
    float  genhhMass   =0;
    float  genhhMass2   =0;
};

#endif
