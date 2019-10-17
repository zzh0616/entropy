#shell program to calculate lx & energy loss 
#using calc_cf_cooling.py
#should be run in the individual fitting dir
# $1 param result file(csv) for each cluster
# $2 saved data_array_file(json) from mcmc sampling
SCRIP_DIR=`echo $0 | awk -F"/calc_lx" '{print $1}'`
Z=`cat param_zzh_for_py.txt | grep ^z | awk '{print $2}'`
DL_CM=`${SCRIP_DIR}/../mass_profile_all/mass_profile/calc_distance ${Z} | grep ^d_l_cm | awk '{print $2}'`
DM_CM=`echo ${DL_CM} ${Z} | awk '{print $1/(1+$2)}'`
ABUND=`cat global.cfg | grep ^abund | awk '{print $2}'`
NH=`cat global.cfg | grep ^nh | awk '{print $2}'`

#${SCRIP_DIR}/coolfunc_for_lxcalc.sh  T_dump_tmp.dat ${ABUND} ${NH} ${Z} cfunc_for_lxcalc.txt >/dev/null

${SCRIP_DIR}/feedback.py $2 $1 cfunc_for_lxcalc.txt ${DM_CM}


