#ifndef TREEANALYZER_TREEREADERS_FatJetREADER_H
#define TREEANALYZER_TREEREADERS_FatJetREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/FatJet.h"

namespace TAna{
class FatJetReader: public BaseReader {
public:
    FatJetReader(std::string branchName, bool isRealData, bool fillGenFatJets = true);
	virtual ~FatJetReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings
	bool realData   ;
	bool fillGenFatJets;
public:
	//branches from the tree
    std::vector<float>          * pt             = new std::vector<float>;
    std::vector<float>          * eta            = new std::vector<float>;
    std::vector<float>          * phi            = new std::vector<float>;
    std::vector<float>          * mass           = new std::vector<float>;
    std::vector<float>          * toRawFact      = new std::vector<float>;
    std::vector<float>          * csv            = new std::vector<float>;
    std::vector<ASTypes::size8> * id             = new std::vector<ASTypes::size8>;
    std::vector<ASTypes::int8>  * hadronFlavor   = new std::vector<ASTypes::int8> ;
    std::vector<ASTypes::int8>  * partonFlavor   = new std::vector<ASTypes::int8> ;
    std::vector<float>          * JECUnc         = new std::vector<float> ;
    std::vector<ASTypes::size8> * genIDX         = new std::vector<ASTypes::size8>;
    std::vector<float>          * gen_pt         = new std::vector<float>;
    std::vector<float>          * gen_eta        = new std::vector<float>;
    std::vector<float>          * gen_phi        = new std::vector<float>;
    std::vector<float>          * gen_mass       = new std::vector<float>;

    std::vector<float>          * bbt               = new std::vector<float>;
    std::vector<float>          * tau1              = new std::vector<float>;
    std::vector<float>          * tau2              = new std::vector<float>;
    std::vector<float>          * tau3              = new std::vector<float>;
    std::vector<float>          * sj1_pt            = new std::vector<float>;
    std::vector<float>          * sj1_eta           = new std::vector<float>;
    std::vector<float>          * sj1_phi           = new std::vector<float>;
    std::vector<float>          * sj1_mass          = new std::vector<float>;
    std::vector<float>          * sj1_toRawFact     = new std::vector<float>;
    std::vector<float>          * sj1_csv           = new std::vector<float>;
    std::vector<float>          * sj1_JECUnc        = new std::vector<float>;
    std::vector<ASTypes::int8>  * sj1_hadronFlavor  = new std::vector<ASTypes::int8>;
    std::vector<ASTypes::int8>  * sj1_partonFlavor  = new std::vector<ASTypes::int8>;
    std::vector<float>          * sj2_pt            = new std::vector<float>;
    std::vector<float>          * sj2_eta           = new std::vector<float>;
    std::vector<float>          * sj2_phi           = new std::vector<float>;
    std::vector<float>          * sj2_mass          = new std::vector<float>;
    std::vector<float>          * sj2_toRawFact     = new std::vector<float>;
    std::vector<float>          * sj2_csv           = new std::vector<float>;
    std::vector<float>          * sj2_JECUnc        = new std::vector<float>;
    std::vector<ASTypes::int8>  * sj2_hadronFlavor  = new std::vector<ASTypes::int8>;
    std::vector<ASTypes::int8>  * sj2_partonFlavor  = new std::vector<ASTypes::int8>;


	//objects created in process
    FatJetCollection jets;
    GenFatJetCollection genJets;


};
}

#endif
