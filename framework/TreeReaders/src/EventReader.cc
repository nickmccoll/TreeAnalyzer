
#include "TreeReaders/interface/EventReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

EventReader::EventReader(std::string branchName, bool isRealData) : BaseReader("EventReader",branchName),realData(isRealData) {};
EventReader::~EventReader() {}

void EventReader::setup(TreeReaderWrapper * wrapper){

    wrapper->setBranch(branchName,"run"               , run              , false);
    wrapper->setBranch(branchName,"lumi"              , lumi             , false);
    wrapper->setBranch(branchName,"event"             , event            , false);
    wrapper->setBranch(branchName,"goodVtx"           , goodVtx          , false);
    wrapper->setBranch(branchName,"npv"               , npv              , false);
    wrapper->setBranch(branchName,"rho"               , rho              , false);
    wrapper->setBranch(branchName,"met_pt"            , met_pt           , true );
    wrapper->setBranch(branchName,"met_phi"           , met_phi          , true );
    wrapper->setBranch(branchName,"met_sig"           , met_sig          , false);
    wrapper->setBranch(branchName,"met_unclUp_pt"     , met_unclUp_pt    , false);
    wrapper->setBranch(branchName,"met_unclUp_phi"    , met_unclUp_phi   , false);
    wrapper->setBranch(branchName,"met_raw_pt"        , met_raw_pt       , true);
    wrapper->setBranch(branchName,"met_raw_phi"       , met_raw_phi      , true);
    wrapper->setBranch(branchName,"dataEra"           , dataEra         , false);
    if(!realData){
        wrapper->setBranch(branchName,"nTruePUInts"       ,nTruePUInts  , false);
        wrapper->setBranch(branchName,"genWeight"         ,genWeight    , false);
        wrapper->setBranch(branchName,"process"           ,process      , false);
        wrapper->setBranch(branchName,"genWeights"        ,genWeights   , false);
    } else {
        wrapper->setBranch(branchName,"dataset"           , dataset      , false);
        wrapper->setBranch(branchName,"dataRun"           , dataRun      , false);
    }

    wrapper->setBranch(branchName,"metFilters"         , metFilters      , false);
    wrapper->setBranch(branchName,"triggerAccepts"    , triggerAccepts   , false);
    wrapper->setBranch(branchName,"triggerPrescales"  , triggerPrescales , false);


}

void EventReader::processVars() {
    met.setP4(*met_pt,float(0),*met_phi,float(0));
    rawMet.setP4(*met_raw_pt,float(0),*met_raw_phi,float(0));
    if(realData){
        weight = 1;
    } else {
        weight = *genWeight;
    }
}


}
