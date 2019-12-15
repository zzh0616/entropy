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
SCRIP_DIR=`echo $0 | awk -F"/all_temp.sh" '{print $1}'`
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

${SCRIP_DIR}/coolfunc_calc.sh  T_dump_tmp.dat ${ABUND_DUMP_FILE} ${NH_DUMP_FILE} ${REDSHIFT} cfunc_for_density_fit_cnt.txt ${E_LOW} ${E_HIGH} cfunc_for_density_fit_erg.txt
${SCRIP_DIR}/coolfunc_calc.sh  T_dump_tmp.dat ${ABUND_DUMP_FILE} ${NH_DUMP_FILE} ${REDSHIFT} cfunc_for_chandra_density_fit_cnt.txt 0.7 7.0 cfunc_for_chandra_density_fit_erg.txt
if [ "${FLAG}" = "OLD" ]; then
    cp cfunc_for_densit_fit_cnt.txt cfunc_for_density_fit.txt
fi

