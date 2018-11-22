#ifndef TREEANALYZER_TREEREADERS_ElectronREADER_H
#define TREEANALYZER_TREEREADERS_ElectronREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Electron.h"

namespace TAna{
class ElectronReader: public BaseReader {
public:
    ElectronReader(std::string branchName, bool fillSCs = false, bool fillReco = false);
	virtual ~ElectronReader();
	virtual void setup(TreeReaderWrapper * wrapper);
	virtual void processVars();

private:
	//settings
	bool fillSCs;
	bool fillReco;
public:
	//branches from the tree
     std::vector<float>          * pt          = new std::vector<float> ;
     std::vector<float>          * eta         = new std::vector<float> ;
     std::vector<float>          * phi         = new std::vector<float> ;
     std::vector<ASTypes::int8  >* q           = new std::vector<ASTypes::int8  >;
     std::vector<float >         * scEta       = new std::vector<float >;
     std::vector<float >         * d0          = new std::vector<float >;
     std::vector<float >         * dz          = new std::vector<float >;
     std::vector<float >         * sip3D       = new std::vector<float >;
     std::vector<float >         * mvaID       = new std::vector<float >;
     std::vector<ASTypes::size8 >* mvaID_cat   = new std::vector<ASTypes::size8 >;
     std::vector<float >         * miniIso     = new std::vector<float >;
     std::vector<float >         * eaRelISO    = new std::vector<float >;
     std::vector<ASTypes::size16>* id          = new std::vector<ASTypes::size16>;
     std::vector<float>			 * dRnorm      = new std::vector<float> ;
     std::vector<float>			 * lepAct_o_pt = new std::vector<float> ;
     std::vector<float>          * sc_act_o_pt = new std::vector<float> ;
     std::vector<float>          * sc_dr_act   = new std::vector<float> ;
     std::vector<ASTypes::size8> * reco_flag   = new std::vector<ASTypes::size8> ;

     std::vector<float>          * sccol_et    = new std::vector<float> ;
     std::vector<float>          * sccol_eta   = new std::vector<float> ;
     std::vector<float>          * sccol_phi   = new std::vector<float> ;

	//objects created in process
    ElectronCollection electrons;
    MomentumFCollection superclusters;

};
}

#endif
