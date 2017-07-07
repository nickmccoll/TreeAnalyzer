#ifndef TREEANALYZER_TREEREADERS_TemplateREADER_H
#define TREEANALYZER_TREEREADERS_TemplateREADER_H
#include "TreeReaders/interface/BaseReader.h"

namespace TAna{
class TemplateReader: public BaseReader {
public:
    TemplateReader(std::string branchName);
	virtual ~TemplateReader();
	virtual void setup(TreeReadingWrapper * wrapper);
	virtual void processVars();

private:
	//settings

public:
	//branches from the tree


	//objects created in process



};
}

#endif
