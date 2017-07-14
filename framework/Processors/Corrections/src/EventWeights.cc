
#include "Processors/Corrections/interface/EventWeights.h"
#include "TreeReaders/interface/EventReader.h"
namespace TAna {
namespace EventWeights {
    float calcNormalizedEventWeight(const EventReader * reader_event, const float cross, const float numE, const float lumi) {
        return (reader_event->weight > 0 ? 1.0 : -1.0) *lumi * cross *1000 /numE;
    }
    float getNormalizedEventWeight(const EventReader * reader_event, const float cross, const float numE, const float lumi) {
        if(reader_event->normWeightLoaded) return lumi*reader_event->nomrmWeight;
        if(cross < 0 ||numE < 0) return reader_event->weight;
        return calcNormalizedEventWeight(reader_event,cross,numE,lumi);
    }
}

}



