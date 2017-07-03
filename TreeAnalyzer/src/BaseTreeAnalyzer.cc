#include "../interface/BaseTreeAnalyzer.h"
#include<iostream>
using namespace std;

namespace TAna {

//--------------------------------------------------------------------------------------------------
void BaseEventAnalyzer::analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency, int numEvents, int startEvent){
    cout << "Running over " << (numEvents < 0 ? "all" : TString::Format("at most %i",numEvents).Data()) << " events";
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
        ana->runEvent();
        ana->setEventNumber(ana->getEventNumber() +1);
    }
}

//--------------------------------------------------------------------------------------------------
BaseTreeAnalyzer::BaseTreeAnalyzer(std::string fileName, std::string treeName, size randomSeed) :
		        tree(fileName,treeName), eventNumber(0), randGen (new TRandom3(randomSeed))
{}

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
void BaseTreeAnalyzer::setupReaders() {for(auto * r : readers) r->setup(&tree);}

void BaseTreeAnalyzer::initializeTreeCopy(std::string outFileName, TreeCopyingOptions copyOptions) {
    outTreeName = outFileName;
    outTreeCopyOpt = copyOptions;
}
void BaseTreeAnalyzer::setupOutTree(){
    if(outTreeCopyOpt == COPY_ERROR)
        throw std::invalid_argument("BaseTreeAnalyzer::setupOutTree() -> Must call TreeCopyingOptions() before you call analyze()");

    TFile * outFile = new TFile(outTreeName.c_str(),"RECREATE");
    outFile->cd();

    TTree * newTree = 0;
    switch(outTreeCopyOpt){
    case COPY_ALL :
        tree.getTree()->SetBranchStatus("*",1);
        newTree = tree.getTree()->CloneTree(0);
        break;
    case COPY_LOADED :
        newTree = tree.getTree()->CloneTree(0);
        break;
    default:
        newTree = new TTree(tree.getTree()->GetName(),tree.getTree()->GetTitle());
    }
    outTree = new TreeWriter(outFile, newTree);
}

//--------------------------------------------------------------------------------------------------
void BaseTreeAnalyzer::bookOutputVariables() {
    throw std::invalid_argument("BaseTreeAnalyzer::bookOutputVariables() -> Must define if you are running an EventAnalyzer that makes a new tree!");
};

}
