
#ifndef TREEANALYZER_TREEANALYZER_BASETREECOPIER_H
#define TREEANALYZER_TREEANALYZER_BASETREECOPIER_H

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
namespace TAna {

//--------------------------------------------------------------------------------------------------
// EventAnalyzers to define how copiers are processed
//--------------------------------------------------------------------------------------------------

//For every entry copy a new entry if fillEvent returns true
class CopierEventAnalyzer : public BaseEventAnalyzer {
public:
    CopierEventAnalyzer() {};
    virtual ~CopierEventAnalyzer() {};
    virtual void analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency = 10000, int numEvents = -1, int startEvent = 0);
};

//Only fill the new tree when the user calls fillOutTree();
class ManualCopierEventAnalyzer : public BaseEventAnalyzer {
public:
    ManualCopierEventAnalyzer() {};
    virtual ~ManualCopierEventAnalyzer() {};
    virtual void analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency = 10000, int numEvents = -1, int startEvent = 0);
};

}


#endif
