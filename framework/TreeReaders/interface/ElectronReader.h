#ifndef TREEANALYZER_TREEREADERS_ElectronREADER_H
#define TREEANALYZER_TREEREADERS_ElectronREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Electron.h"

namespace TAna{
//--------------------------------------------------------------------------------------------------
class ElectronReader: public BaseReader {
public:
    ElectronReader(std::string branchName, bool storeIDVars = false);
	virtual ~ElectronReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings
	bool storeIDVars;
public:
	//branches from the tree
	ra_float pt                    ;
	ra_float eta                   ;
	ra_float phi                   ;
	ra_int8 q                      ;
	ra_float scEta                 ;
	ra_float scE                   ;
	ra_float uncorPt               ;
	ra_size16 id                  ;
	ra_float d0                    ;
	ra_float dz                    ;
	ra_float sip3D                ;
	ra_float mvaID                 ;
	ra_float miniIso               ;
	ra_float eaRelIso              ;
	ra_float dRnorm                ;
	ra_float lepAct_o_pt           ;
	ra_float sc_act_o_pt           ;
	ra_float sc_dr_act            ;
	ra_float tthMVA                ;
	ra_size16 passMedCutBased      ;
	ra_size16 passTightCutBased    ;
	ra_float  full5x5_sigmaIetaIeta;
	ra_float  abs_dEtaSeed         ;
	ra_float  abs_dPhiIn           ;
	ra_float  HoE                  ;
	ra_float  HoE_BC               ;
	ra_float  abs_1oEm1op          ;
	ra_size8  missInnerHits        ;
	ra_size8  passConvVeto         ;
	ra_size8  seeding              ;
	ra_size8  nSaturatedXtals      ;
	ra_float  e2x5OverE5x5         ;
	ra_float  e1x5OverE5x5         ;
	ra_float  isolEmHadDepth1      ;

	//objects created in process
    ElectronCollection electrons;

};
}

#endif
