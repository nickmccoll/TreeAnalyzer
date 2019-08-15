#ifndef TREEANALYZER_TREEREADERS_GenParticleREADER_H
#define TREEANALYZER_TREEREADERS_GenParticleREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/GenParticle.h"
namespace TAna{
class GenParticleReader: public BaseReader {
public:

    GenParticleReader(std::string branchName);
	virtual ~GenParticleReader();

	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings

public:
	//branches from the tree

    ra_float  pt       ;
    ra_float  eta      ;
    ra_float  phi      ;
    ra_float  mass     ;
    ra_size8  status   ;
    ra_int    pdgid    ;
    ra_size16 nmoms    ;
    ra_size16 firstmom ;
    ra_size16 assoc    ;
    rd_float  mTTBar   ;


	//objects created in process
    GenParticleCollection genParticles;


};
}

#endif
