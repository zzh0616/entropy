#shell program to calculate lx & energy loss 
#using calc_cf_cooling.py
#should be run in the individual fitting dir
# $1 name of the cluster
# $2 "read" or "calc" the sum_array_info
SCRIP_DIR=`echo $0 | awk -F"/calc_lx" '{print $1}'`
Z=`cat param_zzh_for_py.txt | grep ^z | awk '{print $2}'`
DL_CM=`${SCRIP_DIR}/../mass_profile_all/mass_profile/calc_distance ${Z} | grep ^d_l_cm | awk '{print $2}'`
DM_CM=`echo ${DL_CM} ${Z} | awk '{print $1/(1+$2)}'`
ABUND=`cat global.cfg | grep ^abund | awk '{print $2}'`
NH=`cat global.cfg | grep ^nh | awk '{print $2}'`
NAME=$1
FLAG_CALC=$2
${SCRIP_DIR}/calc_Tdump.py ${NAME} "main"
${SCRIP_DIR}/coolfunc_for_lxcalc.sh  T_dump_tmp.dat A_dump_tmp.dat NH_dump_tmp.dat ${Z} cfunc_for_lxcalc.txt

#${SCRIP_DIR}/feedback.py ${NAME} cfunc_for_lxcalc.txt ${DM_CM} ${FLAG_CALC}


