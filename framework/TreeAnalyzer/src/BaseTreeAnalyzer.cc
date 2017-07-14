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
    if(startEvent >= 0 ){
        ana->setEventNumber(startEvent);
        if(numEvents >= 0 ) numEvents += startEvent;
    }
    while(ana->nextEvent(reportFrequency)){
        if(numEvents >= 0 && ana->getEventNumber() >= numEvents+1) return;
        ana->processReaders();
        ana->runEvent();
        ana->setEventNumber(ana->getEventNumber() +1);
    }
}

//--------------------------------------------------------------------------------------------------
BaseTreeAnalyzer::BaseTreeAnalyzer(std::string fileName, std::string treeName, int inputTreeType, size randomSeed) :
        treeType(getTreeType(inputTreeType)),tree(fileName,treeName), eventNumber(0), randGen (new TRandom3(randomSeed))
{
    std::cout << " \033[1;34m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m"  << std::endl;
    std::cout << " ++  Setting up BaseTreeAnalyzer"<<std::endl;
    std::cout << " ++  Will be analyzing a "<<getTreeTypeName(inputTreeType) <<" tree"<<std::endl;
    std::cout << " \033[1;34m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m"  << std::endl;

}

//--------------------------------------------------------------------------------------------------
BaseTreeAnalyzer::~BaseTreeAnalyzer() {
    for(auto * r : readers) delete r;
    delete randGen;
    if(outTree){outTree->write();}
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::analyze(int reportFrequency, int numEvents, int startEvent) {
    auto * evtAna = setupEventAnalyzer();
    evtAna->analyzeEvent(this,reportFrequency,numEvents,startEvent);
    delete evtAna;
}

//--------------------------------------------------------------------------------------------------
BaseReader* BaseTreeAnalyzer::load(BaseReader * reader) {
    readers.push_back(reader);
    return readers.back();
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::setupReaders() {for(auto * r : readers) r->initialize(&tree);}
//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::processReaders() {for(auto * r : readers) r->processVars();}
//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::initializeTreeCopy(std::string outFileName, TreeCopyingOptions copyOptions) {
    outTreeName = outFileName;
    outTreeCopyOpt = copyOptions;
}
void BaseTreeAnalyzer::setupOutTree(){
    cout << " ++  Creating file " <<outTreeName;

    if(outTreeCopyOpt == COPY_ERROR)
        throw std::invalid_argument("BaseTreeAnalyzer::setupOutTree() -> Must call TreeCopyingOptions() before you call analyze()");

    TFile * outFile = new TFile(outTreeName.c_str(),"RECREATE");
    outFile->cd();

    TTree * newTree = 0;
    switch(outTreeCopyOpt){
    case COPY_ALL :
        tree.getTree()->SetBranchStatus("*",1);
        newTree = tree.getTree()->CloneTree(0);
        std::cout <<", and will copy all branches from the input tree.";
        break;
    case COPY_LOADED :
        newTree = tree.getTree()->CloneTree(0);
        std::cout <<", and will copy only loaded branches from the input tree.";
        break;
    default:
        newTree = new TTree(tree.getTree()->GetName(),tree.getTree()->GetTitle());
        std::cout <<", and will copy no branches from the input tree.";
    }
    outTree = new TreeWriter(outFile, newTree);
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::bookOutputVariables() {
    throw std::invalid_argument("BaseTreeAnalyzer::bookOutputVariables() -> Must define if you are running an EventAnalyzer that makes a new tree!");
};

}
