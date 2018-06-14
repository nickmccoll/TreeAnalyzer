#ifndef TREEANALYZER_TREEREADERS_MuonREADER_H
#define TREEANALYZER_TREEREADERS_MuonREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Muon.h"

namespace TAna{
class MuonReader: public BaseReader {
public:
    MuonReader(std::string branchName);
	virtual ~MuonReader();
	virtual void setup(TreeReadingWrapper * wrapper);
	virtual void processVars();

private:
	//settings

public:
	//branches from the tree
     std::vector<float>          * pt       = new std::vector<float> ;
     std::vector<float>          * eta      = new std::vector<float> ;
     std::vector<float>          * phi      = new std::vector<float> ;
     std::vector<ASTypes::int8  >* q        = new std::vector<ASTypes::int8  >;
     std::vector<float >         * d0       = new std::vector<float >;
     std::vector<float >         * dz       = new std::vector<float >;
     std::vector<float >         * sip3D    = new std::vector<float >;
     std::vector<float >         * miniIso  = new std::vector<float >;
     std::vector<float >         * dBRelISO = new std::vector<float >;
     std::vector<ASTypes::size16>* id       = new std::vector<ASTypes::size16>;

     std::vector<float>			 * dRnorm   = new std::vector<float> ;
     std::vector<float>			 * PtRatioLepAct   = new std::vector<float> ;

	//objects created in process
    MuonCollection muons;

};
}

#endif
