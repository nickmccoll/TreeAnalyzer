
#include "TreeReaders/interface/TemplateReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

TemplateReader::TemplateReader(std::string branchName) : BaseReader("TemplateReader",branchName){};
TemplateReader::~TemplateReader() {}

void TemplateReader::setup(TreeReadingWrapper * wrapper){
    double a;
    wrapper->setBranchAddressPre(branchName,"run"               , &a              , false);

}

void TemplateReader::processVars() {
}


}
