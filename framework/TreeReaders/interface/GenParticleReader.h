#ifndef TREEANALYZER_TREEREADERS_GenParticleREADER_H
#define TREEANALYZER_TREEREADERS_GenParticleREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/GenParticle.h"
namespace TAna{
class GenParticleReader: public BaseReader {
public:
    typedef ASTypes::size16 stor;

    GenParticleReader(std::string branchName);
	virtual ~GenParticleReader();

	virtual void setup(TreeReadingWrapper * wrapper);
	virtual void processVars();

private:
	//settings

public:
	//branches from the tree

    std::vector<float >*          pt       = new std::vector<float >;
    std::vector<float >*          eta      = new std::vector<float >;
    std::vector<float >*          phi      = new std::vector<float >;
    std::vector<float >*          mass     = new std::vector<float >;
    std::vector<ASTypes::size8 >* status   = new std::vector<ASTypes::size8 >;
    std::vector<int   >*          pdgid    = new std::vector<int   >;
    std::vector<stor  >*          nmoms    = new std::vector<stor  >;
    std::vector<stor  >*          firstmom = new std::vector<stor  >;
    std::vector<stor  >*          ndaus    = new std::vector<stor  >;
    std::vector<stor  >*          firstdau = new std::vector<stor  >;
    std::vector<stor  >*          assoc    = new std::vector<stor  >;


	//objects created in process
    GenParticleCollection genParticles;


};
}

#endif
