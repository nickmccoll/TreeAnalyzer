#ifndef TREEANALYZER_TREEREADERS_JetREADER_H
#define TREEANALYZER_TREEREADERS_JetREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Jet.h"

namespace TAna{
//--------------------------------------------------------------------------------------------------
// JetReader
//--------------------------------------------------------------------------------------------------
class JetReader: public BaseReader {
public:
    JetReader(std::string branchName, bool isRealData,
            bool fillGenJets = true,  bool fillRecoJets = true);
	virtual ~JetReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings
	bool realData   ;
	bool fillGenJets;
	bool fillRecoJets;
public:
	//branches from the tree
    ra_float pt                 ;
    ra_float eta                ;
    ra_float phi                ;
    ra_float mass               ;
    ra_float toRawFact          ;
    ra_float metUnc_rawPx       ;
    ra_float metUnc_rawPy       ;
    ra_float csv                ;
    ra_float deep_csv           ;
    ra_size8 id                 ;
    ra_int8  hadronFlavor       ;
    ra_int8  partonFlavor       ;
    ra_float JECUnc             ;
    ra_size8 genIDX             ;
    ra_float gen_pt             ;
    ra_float gen_eta            ;
    ra_float gen_phi            ;
    ra_float gen_mass           ;

	//objects created in process
    JetCollection jets;
    GenJetCollection genJets;


};
}

#endif
