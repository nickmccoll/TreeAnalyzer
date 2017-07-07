#ifndef TREEANALYZER_TREEREADERS_ElectronREADER_H
#define TREEANALYZER_TREEREADERS_ElectronREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Electron.h"

namespace TAna{
class ElectronReader: public BaseReader {
public:
    ElectronReader(std::string branchName);
	virtual ~ElectronReader();
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
     std::vector<float >         * scEta    = new std::vector<float >;
     std::vector<float >         * d0       = new std::vector<float >;
     std::vector<float >         * dz       = new std::vector<float >;
     std::vector<float >         * sip3D    = new std::vector<float >;
     std::vector<float >         * mvaID    = new std::vector<float >;
     std::vector<ASTypes::size8 >* mvaID_cat= new std::vector<ASTypes::size8 >;
     std::vector<float >         * miniIso  = new std::vector<float >;
     std::vector<float >         * eaRelISO = new std::vector<float >;
     std::vector<ASTypes::size16>* id       = new std::vector<ASTypes::size16>;

	//objects created in process
    ElectronCollection electrons;

};
}

#endif
