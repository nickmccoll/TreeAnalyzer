#ifndef TREEANALYZER_TREEREADERS_BASEREADER_H
#define TREEANALYZER_TREEREADERS_BASEREADER_H
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
using ASTypes::size  ;
class TreeReadingWrapper;
namespace TAna{
class BaseReader {
public:
	BaseReader(std::string readerName,std::string branchName) : readerName(readerName), branchName(branchName) {};
	virtual ~BaseReader() {}
	//external functions
	virtual void initialize(TreeReadingWrapper * wrapper) final;

	virtual void setup(TreeReadingWrapper * wrapper) = 0; //set active branches
	virtual void processVars() = 0; //build objects

	const std::string readerName ="";
	const std::string branchName ="";

};
}

#endif
