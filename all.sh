#!/bin/bash
#####################
#run prepare.sh before using this script
#USEAGE:
# $1 temp_data_file
# $2 param_file_for_fitting tempareture
# $3 cluster_base_name e.g. A1795
# $4 global.cfg 
# $5 "old" or "new" input other than "old" is equiv to "new"
# $1-$5 parameters are required
# 
######################
SCRIP_DIR=`echo $0 | awk -F"/all.sh" '{print $1}'`
COSMO_DIR="${SCRIP_DIR}/../mass_profile_all/mass_profile"
if [ "$5" = "old"  ]; then
    TEMP_DATA_FILE=$1
    PARAM_FILE=$2
    NH=`cat $4 | grep ^nh | awk '{print $2}'`
    ABUND=`cat $4 | grep ^abund | awk '{print $2}'`
    SBP_CFG=`cat $4 | grep ^sbp_cfg | awk '{print $2}'`
    SBP_DATA_FILE=`cat $4 | grep ^radius_sbp_file | awk '{print $2}'`
    E_HIGH=7.0
    E_LOW=0.7
    SBP_TYPE="CNT"
    FLAG="OLD"
    ABUND_DUMP_FILE=${ABUND}
    NH_DUMP_FILE=${NH}
else
    GLOBAL_CFG=$4
    ABUND_DATA_FILE=`cat ${GLOBAL_CFG} | grep ^abund_data_file | awk '{print $2}'`
    E_HIGH=`cat ${GLOBAL_CFG} | grep ^sbp_eninfo | awk '{print $3}'`
    E_LOW=`cat ${GLOBAL_CFG} | grep ^sbp_eninfo | awk '{print $2}'`
    SBP_TYPE=`cat ${GLOBAL_CFG} | grep ^sbp_eninfo | awk '{print $4}'`
    TEMP_DATA_FILE=`cat ${GLOBAL_CFG} | grep ^temp_data_file | awk '{print $2}'`
    PARAM_FILE=`cat ${GLOBAL_CFG} | grep ^param_file | awk '{print $2}'`
    NH_DATA_FILE=`cat ${GLOBAL_CFG} | grep ^nh_data_file | awk '{print $2}'`
    SBP_DATA_FILE=`cat ${GLOBAL_CFG} | grep ^sbp_data_file | awk '{print $2}'`
    FLAG="NEW"
    ABUND_DUMP_FILE="A_dump_tmp.dat"
    NH_DUMP_FILE="NH_dump_tmp.dat"
fi
R500=`cat ${PARAM_FILE} | grep ^R500 | awk '{print $2}'`
R200=`cat ${PARAM_FILE} | grep ^R200 | awk '{print $2}'`
REDSHIFT=`cat ${PARAM_FILE} | grep '^z\s' | awk '{print $2}'`
CM_PER_PIXEL=`${COSMO_DIR}/calc_distance ${REDSHIFT} | grep ^cm_per_pixel | awk '{print $2}'`
KPC_PER_PIXEL=`${COSMO_DIR}/calc_distance ${REDSHIFT} | grep ^kpc_per_pixel | awk '{print $2}'`
if [ -e "param_tmp.txt" ]; then 
    rm param_tmp.txt
fi
${SCRIP_DIR}/../density/fit_zzh_model.py ${TEMP_DATA_FILE} ${PARAM_FILE} param_tmp.txt data_temp.txt

if [ -e "T_dump.dat" ] ; then
    mv T_dump.dat T_dump.dat_bak
fi

if [ -e "T_dump_tmp.dat" ] ; then 
    rm T_dump_tmp.dat
fi
if [ -e "${ABUND_DUMP_FILE}" ] ; then
    rm ${ABUND_DUMP_FILE}
fi
if [ -e "${NH_DUMP_FILE}" ] ; then
    rm ${NH_DUMP_FILE}
fi
if [ "${FLAG}" = "NEW"  ]; then
    ${SCRIP_DIR}/../density/calc_T_profile.py param_tmp.txt T_dump_tmp.dat ${ABUND_DATA_FILE} ${ABUND_DUMP_FILE} ${NH_DATA_FILE} ${NH_DUMP_FILE}
else 
    ${SCRIP_DIR}/../density/calc_T_profile.py param_tmp.txt T_dump_tmp.dat
fi

${SCRIP_DIR}/coolfunc_calc.sh  T_dump_tmp.dat ${ABUND_DUMP_FILE} ${NH_DUMP_FILE} ${REDSHIFT} cfunc_for_density_fit_cnt.txt ${E_LOW} ${E_HIGH} cfunc_for_density_fit_erg.txt
${SCRIP_DIR}/coolfunc_calc.sh  T_dump_tmp.dat ${ABUND_DUMP_FILE} ${NH_DUMP_FILE} ${REDSHIFT} cfunc_for_chandra_density_fit_cnt.txt 0.7 7.0 cfunc_for_chandra_density_fit_erg.txt
if [ "${FLAG}" = "OLD" ]; then
    cp cfunc_for_densit_fit_cnt.txt cfunc_for_density_fit.txt
fi
#${SCRIP_DIR}/calc_T_profile.py param_tmp.txt T_dump.dat
if [ -e "param_sbp.txt" ] ; then 
    mv param_sbp.txt param_sbp.txt_bak
fi
if [ -e "data_sbp.txt" ]; then
    mv data_sbp.txt data_sbp.txt_bak
fi
#${SCRIP_DIR}/../density/fit_sbp.py ${SBP_DATA_FILE} ${SBP_CFG} cfunc_for_density_fit.txt param_sbp.txt data_sbp.txt

if [ -e den_dump.dat ]; then
    mv den_dump.dat den_dump.dat_bak
fi
#cat param_sbp.txt | awk -F"]" '{print $1}' | awk -F"[" '{print $2}' >tmp
#mv tmp param_sbp.txt
#${SCRIP_DIR}/calc_den_profile.py param_sbp.txt den_dump.dat

#${SCRIP_DIR}/calc_entropy.py T_dump.dat den_dump.dat entropy.txt
#cat entropy.txt | awk -F"[" '{print $2}' | awk -F"]" '{print $1}' >tmp
#mv tmp entropy.txt

cp $2 param_entropy.txt
NEC_0=0.05
#NEC_0=`cat param_sbp.txt | awk '{print $1+$4+$7}' `
echo "nec ${NEC_0} 0.05" >> param_entropy.txt
echo "nth_A 0.452 0.1" >>param_entropy.txt
echo "nth_B 1.401 0.7" >>param_entropy.txt
echo "nth_gamma 1.628 0.8" >>param_entropy.txt
echo "s 3000 5000 2000 10000" >>param_entropy.txt
echo "tau 1.3 0.5 1 4" >>param_entropy.txt
echo "fg 0.13 0.03 0.01 0.25" >>param_entropy.txt
###end here
if [ ${FLAG} = "NEW" ]; then
    if [ ${SBP_TYPE} = "CNT" ]; then
        CFUNC_FILE="cfunc_for_density_fit_cnt.txt"
    else 
        CFUNC_FILE="cfunc_for_density_fit_erg.txt"
    fi
else
    CFUNC_FILE="cfunc_for_density_fit.txt"
fi
echo "${TEMP_DATA_FILE} ${SBP_DATA_FILE} ${CFUNC_FILE} ${CM_PER_PIXEL}"
${SCRIP_DIR}/entropy_mcmc.py ${TEMP_DATA_FILE} param_entropy.txt ${SBP_DATA_FILE} ${CFUNC_FILE} ${CM_PER_PIXEL} cfunc_for_chandra_density_fit_cnt.txt

mv temperature.pdf ${3}_temp.pdf
mv sbp.pdf ${3}_sbp.pdf
mv entropy.pdf ${3}_entropy.pdf
mv result.csv ${3}_result.csv
mv array_plt.json ${3}_plt.json
mv mass.pdf ${3}_mass.pdf
mv csbp.pdf ${3}_csbp.pdf
