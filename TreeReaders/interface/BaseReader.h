#ifndef TREEANALYZER_TREEREADERS_BASEREADER_H
#define TREEANALYZER_TREEREADERS_BASEREADER_H
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"

class TreeReadingWrapper;
namespace TAna{
class BaseReader {
public:
	BaseReader(std::string branchName) : branchName(branchName) {};
	virtual ~BaseReader() {}
	virtual void setup(TreeReadingWrapper * wrapper) = 0; //set active branches
	virtual void process() = 0; //build objects
	bool isProcessed() const {return processStatus;}

	std::string branchName ="";
	bool processStatus     =false;

};
}

#endif
