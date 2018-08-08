
#include "TreeReaders/interface/ElectronReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


namespace TAna{

ElectronReader::ElectronReader(std::string branchName, bool fillSCs) : BaseReader("ElectronReader",branchName),
        fillSCs(fillSCs)
{};

ElectronReader::~ElectronReader(){
    delete pt        ;
    delete eta       ;
    delete phi       ;
    delete q         ;
    delete scEta     ;
    delete d0        ;
    delete dz        ;
    delete sip3D     ;
    delete mvaID     ;
    delete mvaID_cat ;
    delete miniIso   ;
    delete eaRelISO  ;
    delete id        ;
    delete dRnorm    ;
    delete lepAct_o_pt;
    delete sc_act_o_pt;
    delete sc_dr_act;
    delete sccol_et  ;
    delete sccol_eta ;
    delete sccol_phi ;
}

void ElectronReader::setup(TreeReadingWrapper * wrapper){
    wrapper->setBranchAddressPre(branchName,"pt"         ,&pt         );
    wrapper->setBranchAddressPre(branchName,"eta"        ,&eta        );
    wrapper->setBranchAddressPre(branchName,"phi"        ,&phi        );
    wrapper->setBranchAddressPre(branchName,"q"          ,&q          );
    wrapper->setBranchAddressPre(branchName,"scEta"      ,&scEta      );
    wrapper->setBranchAddressPre(branchName,"d0"         ,&d0         );
    wrapper->setBranchAddressPre(branchName,"dz"         ,&dz         );
    wrapper->setBranchAddressPre(branchName,"sip3D"      ,&sip3D      );
    wrapper->setBranchAddressPre(branchName,"mvaID"      ,&mvaID      );
    wrapper->setBranchAddressPre(branchName,"mvaID_cat"  ,&mvaID_cat  );
    wrapper->setBranchAddressPre(branchName,"miniIso"    ,&miniIso    );
    wrapper->setBranchAddressPre(branchName,"eaRelISO"   ,&eaRelISO   );
    wrapper->setBranchAddressPre(branchName,"id"         ,&id         );
    wrapper->setBranchAddressPre(branchName,"dRnorm"     ,&dRnorm     ,false);
    wrapper->setBranchAddressPre(branchName,"lepAct_o_pt",&lepAct_o_pt,false);
    wrapper->setBranchAddressPre(branchName,"sc_act_o_pt",&sc_act_o_pt,false);
    wrapper->setBranchAddressPre(branchName,"sc_dr_act"  ,&sc_dr_act  ,false);
    if(fillSCs){
        wrapper->setBranchAddressPre(branchName,"sccol_et"   ,&sccol_et   ,false);
        wrapper->setBranchAddressPre(branchName,"sccol_eta"  ,&sccol_eta  ,false);
        wrapper->setBranchAddressPre(branchName,"sccol_phi"  ,&sccol_phi  ,false);
    }

}

void ElectronReader::processVars() {
    electrons.clear();
    for(unsigned int iO = 0; iO < pt->size(); ++iO){
        electrons.emplace_back(ASTypes::CylLorentzVectorF(pt->at(iO),eta->at(iO),phi->at(iO),0),iO,
                q->at(iO),d0->at(iO),dz->at(iO),sip3D->at(iO),miniIso->at(iO),
                dRnorm->size() ? dRnorm->at(iO) : 0,lepAct_o_pt->size() ? lepAct_o_pt->at(iO) : 0);
        electrons.back().addElectronInfo(scEta->at(iO),mvaID->at(iO),mvaID_cat->at(iO),eaRelISO->at(iO),id->at(iO),
                sc_act_o_pt->size() ? sc_act_o_pt->at(iO) : 0,sc_dr_act->size() ? sc_dr_act->at(iO) : 0 );
    }
    std::sort(electrons.begin(), electrons.end(), PhysicsUtilities::greaterPT<Electron>());

    if(fillSCs){
        superclusters.clear();
        for(unsigned int iO = 0; iO < sccol_et->size(); ++iO){
            superclusters.emplace_back(ASTypes::CylLorentzVectorF(sccol_et->at(iO),sccol_eta->at(iO),sccol_phi->at(iO),0));
        }
        std::sort(superclusters.begin(), superclusters.end(), PhysicsUtilities::greaterPT<MomentumF>());
    }
}



}
