#ifndef TREEANALYZER_TREEREADERS_MuonREADER_H
#define TREEANALYZER_TREEREADERS_MuonREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Muon.h"

namespace TAna{
//--------------------------------------------------------------------------------------------------
class MuonReader: public BaseReader {
public:
    MuonReader(std::string branchName,bool isRealData);
	virtual ~MuonReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
    //settings
    bool realData   ;

public:
	//branches from the tree
     ra_float pt           ;
     ra_float eta          ;
     ra_float phi          ;
     ra_size  id           ;
     ra_int8  q            ;
     ra_float d0           ;
     ra_float dz           ;
     ra_float sip3D        ;
     ra_float miniIso      ;
     ra_float dBRelISO     ;
     ra_float ptRel        ;
     ra_float ptRatio      ;
     ra_float dRnorm       ;
     ra_float lepAct_o_pt  ;
     ra_int8  simType      ;


	//objects created in process
    MuonCollection muons;

};
}

#endif
