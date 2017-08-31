#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_
namespace TAna {
class EventReader;
namespace EventWeights {
    //function to calculate weight normalized to some lumi
    float calcNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi);
    //Returns event weight normalized to some lumi
    //If this is a data event, set cross section to -1
    //Will first check to see if normWeight is loaded
    //If it is not (was added in post processing), will calculate it on the fly
    //Otherwise will output lumi*normWeight
    float getNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi = 1);

    float get4bXSecLimit(unsigned int  mass);

}
}



#endif /* FRAMEWORK_PROCESSORS_INTERFACE_EVENTWEIGHTS_H_ */
