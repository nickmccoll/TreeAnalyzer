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

  if(startEvent >= 0 ){
      ana->setEventNumber(startEvent);
      if(numEvents >= 0 ) numEvents += startEvent;
  }

  while(ana->nextEvent(reportFrequency)){
      if(numEvents >= 0 && ana->getEventNumber() >= numEvents+1) return;
      ana->processReaders();
      ana->resetOutData();
      if(!ana->runEvent()) continue;
      ana->fillOutTree();
      ana->setEventNumber(ana->getEventNumber() +1);
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

  if(startEvent >= 0 ){
      ana->setEventNumber(startEvent);
      if(numEvents >= 0 ) numEvents += startEvent;
  }

  while(ana->nextEvent(reportFrequency)){
      if(numEvents >= 0 && ana->getEventNumber() >= numEvents+1) return;
      ana->processReaders();
      ana->resetOutData();
      ana->runEvent();
      ana->setEventNumber(ana->getEventNumber() +1);
  }
}
