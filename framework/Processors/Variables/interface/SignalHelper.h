
#ifndef PROCESSORS_VARIABLES_SIGNALHELPER_H
#define PROCESSORS_VARIABLES_SIGNALHELPER_H
#include <vector>
#include "TString.h"

#include "DataFormats/interface/Momentum.h"
#include "Configuration/interface/ReaderConstants.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"

namespace TAna {
class Lepton;
class ElectronReader;
class MuonReader;
class FatJet;
class GenParticle;

class SignalHelper {
private:
//	const FatJet* recoHbb=0;
    std::shared_ptr<ElectronReader   > elReader ;
    std::shared_ptr<MuonReader       > muReader ;
    DiHiggsEvent diHiggsEvt;

public:
	double minElRecoPt = 0;
	double maxElRecoEta = 99;
	double minMuRecoPt = 0;
	double maxMuRecoEta = 99;
	const GenParticle* genlep1=0;
	const GenParticle* genlep2=0;
	const Lepton* recolep1=0;
	const Lepton* recolep2=0;
	DiHiggsEvent::DECAYTYPE type = DiHiggsEvent::BAD;
	std::vector<const Muon*> candMuons;
	std::vector<const Electron*> candElectrons;

	SignalHelper(DiHiggsEvent dhEvt, std::shared_ptr<MuonReader> reader_muon,std::shared_ptr<ElectronReader> reader_electron);
	void setRecoLeptons(double matchDR);
	Lepton* getMatchedLepton(const GenParticle* genLep,double maxDR);
	const GenParticle* getGenLepFromTau(const GenParticle* tau);
	bool hasMatchedSingleLep();
	bool hasMatchedDileps();

};

}
#endif

