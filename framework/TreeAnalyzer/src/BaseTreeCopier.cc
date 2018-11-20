#include "../interface/BaseTreeCopier.h"


using namespace TAna;

//--------------------------------------------------------------------------------------------------
void CopierEventAnalyzer::analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency, int numEvents, int startEvent){

  std::cout << " ++  Running over " << (numEvents < 0 ? "all" : TString::Format("at most %i",numEvents).Data()) << " events";
  if(startEvent >= 0 ) std::cout << ", starting with event: "<< startEvent;
  std::cout <<std::endl;

  ana->loadVariables();
  ana->setupReaders();

  ana->setupOutTree();
  ana->bookOutputVariables();
  ana->bookOutTree();

  ana->setEventRange(startEvent,numEvents);
  while(ana->nextEvent(reportFrequency)){
      ana->processReaders();
      ana->resetOutData();
      if(ana->runEvent())
          ana->fillOutTree();
  }
}

//--------------------------------------------------------------------------------------------------
void ManualCopierEventAnalyzer::analyzeEvent(BaseTreeAnalyzer * ana, int reportFrequency, int numEvents, int startEvent){

  std::cout << " ++  Running over " << (numEvents < 0 ? "all" : TString::Format("at most %i",numEvents).Data()) << " events";
  if(startEvent >= 0 ) std::cout << ", starting with event: "<< startEvent;
  std::cout <<std::endl;

  ana->loadVariables();
  ana->setupReaders();

  ana->setupOutTree();
  ana->bookOutputVariables();
  ana->bookOutTree();

  ana->setEventRange(startEvent,numEvents);
  while(ana->nextEvent(reportFrequency)){
      ana->processReaders();
      ana->resetOutData();
      ana->runEvent();
  }
}
