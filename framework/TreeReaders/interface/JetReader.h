#ifndef TREEANALYZER_TREEREADERS_JetREADER_H
#define TREEANALYZER_TREEREADERS_JetREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Jet.h"

namespace TAna{
class JetReader: public BaseReader {
public:
    JetReader(std::string branchName, bool isRealData, bool fillGenJets = true,  bool fillRecoJets = true);
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
    std::vector<float>          * pt             = new std::vector<float>;
    std::vector<float>          * eta            = new std::vector<float>;
    std::vector<float>          * phi            = new std::vector<float>;
    std::vector<float>          * mass           = new std::vector<float>;
    std::vector<float>          * toRawFact      = new std::vector<float> ;
    std::vector<float>          * chef           = new std::vector<float> ;
    std::vector<float>          * metUnc_rawPx   = new std::vector<float> ;
    std::vector<float>          * metUnc_rawPy   = new std::vector<float> ;
    std::vector<float>          * csv            = new std::vector<float>;
    std::vector<ASTypes::size8> * id             = new std::vector<ASTypes::size8>;
    std::vector<ASTypes::int8>  * hadronFlavor   = new std::vector<ASTypes::int8> ;
    std::vector<ASTypes::int8>  * partonFlavor   = new std::vector<ASTypes::int8> ;
    std::vector<float>          * JECUnc         = new std::vector<float>         ;
    std::vector<ASTypes::size8> * genIDX         = new std::vector<ASTypes::size8>;
    std::vector<float>          * gen_pt         = new std::vector<float>;
    std::vector<float>          * gen_eta        = new std::vector<float>;
    std::vector<float>          * gen_phi        = new std::vector<float>;
    std::vector<float>          * gen_mass       = new std::vector<float>;

	//objects created in process
    JetCollection jets;
    GenJetCollection genJets;


};
}

#endif
