#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_EVENTSELECTION_INTERFACE_EVENTSELECTION_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_EVENTSELECTION_INTERFACE_EVENTSELECTION_H_
#include <vector>
#include "Configuration/interface/ReaderConstants.h"
namespace TAna {
class EventReader;
class Lepton;
namespace EventSelection {

bool passMuonTriggerSuite2016(const EventReader& reader_event);
bool passElectronTriggerSuite2016(const EventReader& reader_event);
bool passTriggerSuite2016(const EventReader& reader_event);

bool passMuonTriggerSuite2017(const EventReader& reader_event);
bool passElectronTriggerSuite2017(const EventReader& reader_event);
bool passTriggerSuite2017(const EventReader& reader_event);

bool passMuonTriggerSuite2018(const EventReader& reader_event);
bool passElectronTriggerSuite2018(const EventReader& reader_event);
bool passTriggerSuite2018(const EventReader& reader_event);

bool alwaysTrue(const EventReader& reader_event);

bool passTriggerPreselection(const EventParameters& params, const EventReader& reader_event,const float ht, const std::vector<const Lepton    *>& selectedLeptons);
bool passEventFilters(const EventParameters& params, const EventReader& reader_event);




}
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
