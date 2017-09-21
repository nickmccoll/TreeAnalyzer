#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_EVENTSELECTION_INTERFACE_EVENTSELECTION_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_EVENTSELECTION_INTERFACE_EVENTSELECTION_H_
#include <vector>
namespace TAna {
class EventReader;
class Lepton;
namespace EventSelection {

bool passEventFilters(const EventReader& reader_event);
bool passTriggerSuite(const EventReader& reader_event);
bool passMuonTriggerSuite(const EventReader& reader_event);
bool passElectronTriggerSuite(const EventReader& reader_event);
bool passTriggerPreselection(const EventReader& reader_event,const float ht, const std::vector<const Lepton    *>& selectedLeptons);


}
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
