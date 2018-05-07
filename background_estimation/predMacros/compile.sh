#!/bin/bash
# ./compile.sh jobs/4_30_beTreees/ compiled/
inputdir=$1
outputdir=$2

hadd ${outputdir}/wjets.root                    ${inputdir}/out_wjets-lnu_*.root
hadd ${outputdir}/ttbar-2l-madgraph.root        ${inputdir}/out_ttbar-2l-*.root
hadd ${outputdir}/wz-lnu2q-amcatnlo.root        ${inputdir}/out_wz-lnu2q*.root
hadd ${outputdir}/tbarW-5f-l-powheg.root        ${inputdir}/out_tbarW-*.root 
hadd ${outputdir}/tW-5f-l-powheg.root           ${inputdir}/out_tW-*.root
hadd ${outputdir}/ww-lnuqq-powheg.root          ${inputdir}/out_ww-lnuqq-*.root
hadd ${outputdir}/t-tbar1l-madgraph.root        ${inputdir}/out_t-tbar1l-*.root
hadd ${outputdir}/dyll.root                     ${inputdir}/out_dyll-*.root
hadd ${outputdir}/tbar-t1l-madgraph.root        ${inputdir}/out_tbar-t1l-*.root
hadd ${outputdir}/ttW-amcatnlo.root             ${inputdir}/out_ttW-*.root
hadd ${outputdir}/ttZ-amcatnlo.root             ${inputdir}/out_ttZ-*.root
hadd ${outputdir}/qcd.root                      ${inputdir}/out_qcd_ht*.root
hadd ${outputdir}/data_mu-2016h-pr-v3.root      ${inputdir}/out_data_e-2016b-23sep16*.root   	
hadd ${outputdir}/data_mu-2016h-pr-v2.root      ${inputdir}/out_data_e-2016c-23sep16*.root   	
hadd ${outputdir}/data_mu-2016g-23sep16.root    ${inputdir}/out_data_e-2016d-23sep16*.root   	
hadd ${outputdir}/data_mu-2016f-23sep16.root    ${inputdir}/out_data_e-2016e-23sep16*.root   	
hadd ${outputdir}/data_mu-2016e-23sep16.root    ${inputdir}/out_data_e-2016f-23sep16*.root   	
hadd ${outputdir}/data_mu-2016d-23sep16.root    ${inputdir}/out_data_e-2016g-23sep16*.root   	
hadd ${outputdir}/data_mu-2016c-23sep16.root    ${inputdir}/out_data_e-2016h-pr-v2*.root     	
hadd ${outputdir}/data_mu-2016b-23sep16.root    ${inputdir}/out_data_e-2016h-pr-v3*.root     	
hadd ${outputdir}/data_met-2016h-pr-v3.root     ${inputdir}/out_data_mu-2016b-23sep16*.root   
hadd ${outputdir}/data_met-2016h-pr-v2.root     ${inputdir}/out_data_mu-2016c-23sep16*.root   
hadd ${outputdir}/data_met-2016g-23sep16.root   ${inputdir}/out_data_mu-2016d-23sep16*.root   
hadd ${outputdir}/data_met-2016f-23sep16.root   ${inputdir}/out_data_mu-2016e-23sep16*.root   
hadd ${outputdir}/data_met-2016e-23sep16.root   ${inputdir}/out_data_mu-2016f-23sep16*.root   
hadd ${outputdir}/data_met-2016d-23sep16.root   ${inputdir}/out_data_mu-2016g-23sep16*.root   
hadd ${outputdir}/data_met-2016c-23sep16.root   ${inputdir}/out_data_mu-2016h-pr-v2*.root     
hadd ${outputdir}/data_met-2016b-23sep16.root   ${inputdir}/out_data_mu-2016h-pr-v3*.root     
hadd ${outputdir}/data_jetht-2016h-pr-v3.root   ${inputdir}/out_data_jetht-2016b-23sep16*.root
hadd ${outputdir}/data_jetht-2016h-pr-v2.root   ${inputdir}/out_data_jetht-2016c-23sep16*.root
hadd ${outputdir}/data_jetht-2016g-23sep16.root ${inputdir}/out_data_jetht-2016d-23sep16*.root
hadd ${outputdir}/data_jetht-2016f-23sep16.root ${inputdir}/out_data_jetht-2016e-23sep16*.root
hadd ${outputdir}/data_jetht-2016e-23sep16.root ${inputdir}/out_data_jetht-2016f-23sep16*.root
hadd ${outputdir}/data_jetht-2016d-23sep16.root ${inputdir}/out_data_jetht-2016g-23sep16*.root
hadd ${outputdir}/data_jetht-2016c-23sep16.root ${inputdir}/out_data_jetht-2016h-pr-v2*.root  
hadd ${outputdir}/data_jetht-2016b-23sep16.root ${inputdir}/out_data_jetht-2016h-pr-v3*.root  
hadd ${outputdir}/data_e-2016h-pr-v3.root       ${inputdir}/out_data_met-2016b-23sep16*.root	
hadd ${outputdir}/data_e-2016h-pr-v2.root       ${inputdir}/out_data_met-2016c-23sep16*.root	
hadd ${outputdir}/data_e-2016g-23sep16.root     ${inputdir}/out_data_met-2016d-23sep16*.root	
hadd ${outputdir}/data_e-2016f-23sep16.root     ${inputdir}/out_data_met-2016e-23sep16*.root	
hadd ${outputdir}/data_e-2016e-23sep16.root     ${inputdir}/out_data_met-2016f-23sep16*.root	
hadd ${outputdir}/data_e-2016d-23sep16.root     ${inputdir}/out_data_met-2016g-23sep16*.root	
hadd ${outputdir}/data_e-2016c-23sep16.root     ${inputdir}/out_data_met-2016h-pr-v2*.root  	
hadd ${outputdir}/data_e-2016b-23sep16.root     ${inputdir}/out_data_met-2016h-pr-v3*.root  	
cp ${inputdir}/out_radion*.root  ${outputdir}/