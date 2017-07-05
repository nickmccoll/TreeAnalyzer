
#include "TreeReaders/interface/TemplateReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

TemplateReader::TemplateReader(std::string branchName) : BaseReader("TemplateReader",branchName){};

void TemplateReader::setup(TreeReadingWrapper * wrapper){
    wrapper->getEntries();

//    wrapper->setBranchAddressPre(branchName,"run"               , &run              , false);

}

void TemplateReader::processVars() {
}


}
