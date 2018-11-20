#include "../interface/BaseTreeAnalyzer.h"
#include<iostream>
using namespace std;

namespace TAna {
//--------------------------------------------------------------------------------------------------
TreeType getTreeType(const int inputInt) {
    switch( inputInt){
    case TREE_DATA:
        return TREE_DATA;
    case TREE_MC:
        return TREE_MC;
    default:
        return TREE_OTHER;
    }
}
std::string getTreeTypeName(const int inputInt) {
    switch( inputInt){
    case TREE_DATA:
        return "DATA";
    case TREE_MC:
        return "MC";
    default:
        return "OTHER";
    }
}

void BaseEventAnalyzer::analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency, int numEvents, int startEvent){
    cout << " ++  Running over " << (numEvents < 0 ? "all" : TString::Format("at most %i",numEvents).Data()) << " events";
    if(startEvent >= 0 ) cout << ", starting with event: "<< startEvent;
    cout <<endl;
    ana->loadVariables();
    ana->setupReaders();
    ana->setEventRange(startEvent,numEvents);

    while(ana->nextEvent(reportFrequency)){
        ana->processReaders();
        ana->runEvent();
    }

}

//--------------------------------------------------------------------------------------------------
BaseTreeAnalyzer::BaseTreeAnalyzer(std::string fileName, std::string treeName, int inputTreeType, size randomSeed) :
        treeType(getTreeType(inputTreeType)),tree(fileName,treeName), eventNumber(0), randGen (std::make_shared<TRandom3>(randomSeed))
{
    //if the enviornment variable is set, use it
    char* datap = getenv ("TREEANALYZER_DATA");
    if(datap != NULL) dataDirectory = datap;

    std::cout << " \033[1;34m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m"  << std::endl;
    std::cout << " ++  Setting up BaseTreeAnalyzer"<<std::endl;
    std::cout << " ++  Will be analyzing a "<<getTreeTypeName(inputTreeType) <<" tree"<<std::endl;
    std::cout << " ++  Data dir: "<<dataDirectory<<std::endl;
    std::cout << " \033[1;34m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m"  << std::endl;

}

//--------------------------------------------------------------------------------------------------
BaseTreeAnalyzer::~BaseTreeAnalyzer() {
    if(outTree){outTree->write();}
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::analyze(int reportFrequency, int numEvents, int startEvent) {
    auto * evtAna = setupEventAnalyzer();
    evtAna->analyzeEvent(this,reportFrequency,numEvents,startEvent);
    delete evtAna;
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::load(std::shared_ptr<BaseReader> reader) {
    readers.push_back(reader);
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::setupReaders() {for(auto& r : readers) r->initialize(&tree);}
//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::processReaders() {for(auto& r : readers) r->processVars();}
//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::initializeTreeCopy(std::string outFileName, TreeCopyingOptions copyOptions, std::string directory) {
    outTreeName = outFileName;
    outTreeCopyOpt = copyOptions;
    outTreeDirectory =directory;
}
void BaseTreeAnalyzer::setupOutTree(){
    cout << " ++  Creating file " <<outTreeName;

    if(outTreeCopyOpt == COPY_ERROR)
        throw std::invalid_argument("BaseTreeAnalyzer::setupOutTree() -> Must call TreeCopyingOptions() before you call analyze()");

    TFile * outFile = new TFile(outTreeName.c_str(),"RECREATE");
    outFile->cd();

    TDirectory *cdtof = 0;
    if(outTreeDirectory != ""){
        cdtof = outFile->mkdir(outTreeDirectory.c_str());
        cdtof->cd();
    }

    TTree * newTree = 0;
    switch(outTreeCopyOpt){
    case COPY_ALL :
        tree.getTree()->SetBranchStatus("*",1);
        newTree = tree.getTree()->CloneTree(0);
        std::cout <<", and will copy all branches from the input tree.\n";
        break;
    case COPY_LOADED :
        newTree = tree.getTree()->CloneTree(0);
        std::cout <<", and will copy only loaded branches from the input tree.\n";
        break;
    default:
        TString name = tree.getTree()->GetName();
        if(name.Contains("/")){ //strip the directory
            const int start = name.Last('/') +1;
            if(start != name.Length()){
                name = name(start, name.Length() - start);
            }
        }
        newTree = new TTree(name,tree.getTree()->GetTitle());
        std::cout <<", and will copy no branches from the input tree.\n";
    }
    outTree.reset(new TreeWriter(outFile, newTree,cdtof));
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::bookOutputVariables() {
    throw std::invalid_argument("BaseTreeAnalyzer::bookOutputVariables() -> Must define if you are running an EventAnalyzer that makes a new tree!");
};
//--------------------------------------------------------------------------------------------------

void BaseTreeAnalyzer::setSampleInfo(float inXSec, float inNumE)
{_xsec =inXSec;_numSampleEvents=inNumE;
    if(inXSec > 0 ){
        std::cout << " ++  Sample cross-section is set to: "<< inXSec << " [pb]"<<std::endl;
        std::cout << " ++  Number of events produced in this sample: "<< inNumE <<std::endl;
    }

}
//--------------------------------------------------------------------------------------------------

void BaseTreeAnalyzer::setLumi(float inLumi) {
    _lumi=inLumi;
    std::cout << " ++  Luminosity is set to: "<< inLumi <<" [fb-1]" <<std::endl;
}


}
