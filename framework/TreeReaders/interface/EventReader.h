#ifndef TREEANALYZER_TREEREADERS_EVENTREADER_H
#define TREEANALYZER_TREEREADERS_EVENTREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Met.h"

using ASTypes::size8 ;
using ASTypes::size16;
using ASTypes::size64;

namespace TAna{
class EventReader: public BaseReader {
public:
    EventReader(std::string branchName, bool isRealData);
	virtual ~EventReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

	//settings
	const bool realData;

public:
	//branches from the tree
    rd_size     run           ;
    rd_size     lumi          ;
    rd_size64   event         ;
    rd_size8    goodVtx       ;
    rd_size16   npv           ;
    rd_float    rho           ;
    rd_float    met_pt        ;
    rd_float    met_phi       ;
    rd_float    met_sig       ;
    rd_float    met_unclUp_pt ;
    rd_float    met_unclUp_phi;
    rd_float    met_raw_pt    ;
    rd_float    met_raw_phi   ;
    rd_float    nTruePUInts   ;
    rd_float    genWeight     ;
    rd_size8    process       ;
    rd_size8    dataEra       ;
    rd_size8    dataset       ;
    rd_size8    dataRun       ;
    ra_float    genWeights    ;

	rd_size   metFilters         ;
	rd_size64 triggerAccepts     ;
	rd_size64 triggerPrescales   ;


	//objects created in process
	Met met;
	Met rawMet;
	float weight =1;



};
}

#endif
