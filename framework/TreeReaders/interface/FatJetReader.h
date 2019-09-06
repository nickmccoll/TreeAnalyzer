#ifndef TREEANALYZER_TREEREADERS_FatJetREADER_H
#define TREEANALYZER_TREEREADERS_FatJetREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/FatJet.h"

namespace TAna{
//--------------------------------------------------------------------------------------------------
// FatJetReader
//--------------------------------------------------------------------------------------------------
class FatJetReader: public BaseReader {
public:
    FatJetReader(std::string branchName, bool isRealData,
            bool fillGenFatJets = true, bool fillBTagging = true, bool fillLSFInfo = false);
	virtual ~FatJetReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings
	bool realData       ;
	bool fillGenFatJets ;
	bool fillBTagging   ;
	bool fillLSFInfo    ;
public:
	//branches from the tree
    ra_float pt                 ;
    ra_float eta                ;
    ra_float phi                ;
    ra_float mass               ;
    ra_float toRawFact          ;
    ra_size8 id                 ;
    ra_float bbt                ;
    ra_float tau1               ;
    ra_float tau2               ;

    ra_float ecfN2              ;
    ra_float ecfM2              ;
    ra_float ecfD2              ;
    ra_float ecfN3              ;
    ra_float ecfU3              ;
    ra_float ecfU2              ;
    ra_float tau3               ;
    ra_float tau4               ;
    ra_float lsf3               ;
    ra_float dRLep              ;
    ra_float lepInJetMVA        ;

    ra_int8  hadronFlavor       ;
    ra_int8  partonFlavor       ;
    ra_float JECUnc             ;
    ra_size8 genIDX             ;
    ra_float gen_pt             ;
    ra_float gen_eta            ;
    ra_float gen_phi            ;
    ra_float gen_mass           ;

    ra_size8 sjIDX1             ;
    ra_size8 sjnum              ;
    ra_float sj_pt              ;
    ra_float sj_eta             ;
    ra_float sj_phi             ;
    ra_float sj_mass            ;
    ra_float sj_toRawFact       ;
    ra_float sj_csv             ;
    ra_float sj_deep_csv        ;
    ra_int8  sj_hadronFlavor    ;
    ra_int8  sj_partonFlavor    ;
    ra_float sj_JECUnc          ;


	//objects created in process
    FatJetCollection jets;
    GenFatJetCollection genJets;


};
}

#endif
