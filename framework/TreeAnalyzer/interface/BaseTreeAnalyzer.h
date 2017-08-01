#ifndef TREEANALYZER_TREEANALYZER_BASETREEANALYZER_H
#define TREEANALYZER_TREEANALYZER_BASETREEANALYZER_H

#include "TRandom3.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/TreeInterface/interface/TreeWriter.h"
#include "TreeReaders/interface/BaseReader.h"
#include <memory>

using ASTypes::size;

namespace TAna {
//Function for processing tree type, expandable for more
enum TreeType {TREE_DATA, TREE_MC, TREE_OTHER};
TreeType getTreeType(const int inputInt);
std::string getTreeTypeName(const int inputInt);


//Event processor, tells the TreeAnalyzer how to
//process events. This one only reads the tree
class BaseTreeAnalyzer;
class BaseEventAnalyzer {
public:
	BaseEventAnalyzer() {};
	virtual ~BaseEventAnalyzer() {};
	virtual void analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency = 10000, int numEvents = -1, int startEvent = -1);
};



class BaseTreeAnalyzer {
public:
	BaseTreeAnalyzer(std::string fileName, std::string treeName, int inputTreeType, size randomSeed = 0);
	virtual ~BaseTreeAnalyzer();

    // Function that user calls to setup the output tree
	// MUST be called before "analyze" if you are making a new tree
	// can override if you are doing something special
	// Otherwise use a copy option:
	// ALL -> Copy all branches from the input tree
	// LOADED -> Copy only the loaded branches
	// NONE -> don't copy any
	enum TreeCopyingOptions {COPY_ERROR, COPY_ALL, COPY_LOADED, COPY_NONE};
    virtual void initializeTreeCopy(std::string outFileName, TreeCopyingOptions copyOptions, std::string directory = "treeMaker");

	// Function user calls to start the run, initializes and
	// starts the EventAnalyzer
	virtual void analyze(int reportFrequency = 10000, int numEvents = -1, int startEvent = -1);

    // Sets up the event analyzer, overload if you want a
    // different one
    virtual BaseEventAnalyzer * setupEventAnalyzer() {return new BaseEventAnalyzer();}


	//--------------------------------------------------------------------------------------------------
	// Functions the user must define in every job
	//--------------------------------------------------------------------------------------------------
	//Place all variable loading, include TreeReaders here
	virtual void loadVariables() = 0;
	//Place event analysis code here, the return value
	//has meaning in some context (e.g. tree copying)
	//This is handled by the EventAnalyzer
	virtual bool runEvent() = 0;

    //--------------------------------------------------------------------------------------------------
    // Additional functions the user must define when you run a tree making job
    //--------------------------------------------------------------------------------------------------
	//Book extra variables in the output tree
	virtual void bookOutputVariables();

	//--------------------------------------------------------------------------------------------------
	// Variable loading
	//--------------------------------------------------------------------------------------------------
	void load(std::shared_ptr<BaseReader> inReader);
	template<typename varType>
	void setBranchAddress(const std::string  branchName,const std::string varName, varType **var, bool require = false, bool verbose = true)
	{ tree.setBranchAddressPre(branchName,varName,var,require,verbose); }
	template<typename varType>
	void setBranchAddress(const std::string  branchName, const std::string varName, varType *var, bool require = false, bool verbose = true)
	{ tree.setBranchAddressPre(branchName,varName,var,require,verbose); }


	//--------------------------------------------------------------------------------------------------
	// Helper functions for event processing
	//--------------------------------------------------------------------------------------------------
	void setupReaders();
	void processReaders();
	bool nextEvent(const int reportFrequency)  {return tree.readEvent(eventNumber, reportFrequency);}
	void setEventNumber(const int newEventNumber) {eventNumber = newEventNumber;}
	void setSampleInfo(float inXSec, float inNumE);
	void setLumi(float inLumi);
	float xsec() const {return _xsec;}
	float nSampEvt() const {return _numSampleEvents;}
	float lumi() const {return _lumi;}

    //--------------------------------------------------------------------------------------------------
    // Helper functions for tree writing
    //--------------------------------------------------------------------------------------------------
	virtual void setupOutTree();
	void bookOutTree()   {outTree->book();}
    void resetOutData()  {outTree->reset();}
    void fillOutTree()   {outTree->fillTree(); }

	//--------------------------------------------------------------------------------------------------
	// Information
	//--------------------------------------------------------------------------------------------------
    std::shared_ptr<TRandom3> getRndGen()  { return randGen;}
	int  getEventNumber() const { return eventNumber;  }
	int  getEntries()     const { return tree.getEntries(); }
	bool isRealData()     const {return treeType == TREE_DATA;}


protected:
    const TreeType           treeType;
	TreeReadingWrapper       tree;
	int                      eventNumber; //current event number
	std::shared_ptr<TRandom3> randGen;

private:
	//Standard variables for normalization
	float _xsec            = -1;
	float _numSampleEvents = -1;
	float _lumi            =  1;
protected:

	std::string        outTreeName = "";
	std::string        outTreeDirectory = "";

	std::unique_ptr<TreeWriter>        outTree;
    TreeCopyingOptions outTreeCopyOpt = COPY_ERROR;
	std::vector<std::shared_ptr<BaseReader> > readers;
};

}
#endif
