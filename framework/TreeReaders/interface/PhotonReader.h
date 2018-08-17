#ifndef TREEANALYZER_TREEREADERS_PhotonReader_H
#define TREEANALYZER_TREEREADERS_PhotonReader_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Photon.h"

namespace TAna{
class PhotonReader: public BaseReader {
public:
    PhotonReader(std::string branchName);
	virtual ~PhotonReader();
	virtual void setup(TreeReadingWrapper * wrapper);
	virtual void processVars();

private:
	//settings
public:
	//branches from the tree
     std::vector<float>          * pt          = new std::vector<float> ;
     std::vector<float>          * eta         = new std::vector<float> ;
     std::vector<float>          * phi         = new std::vector<float> ;
     std::vector<float >         * hadOvEm     = new std::vector<float >;
     std::vector<float >         * hadTowOvEm  = new std::vector<float >;

	//objects created in process
    PhotonCollection photons;
};
}

#endif
