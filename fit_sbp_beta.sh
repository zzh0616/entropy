NAME=$1
NEW_PARAM_FILE=$2
CENTER_ONLY=$3
SBP_DATA_FILE="${NAME}_sbp.txt"
JSON_FILE="${NAME}_plt.json"
z=`cat param_zzh_for_py.txt | grep '^z\s' | awk '{print $2}'`
SBP_TYPE=`cat global.cfg | grep ^sbp_eninfo | awk '{print $4}'`
if [ ${SBP_TYPE} = "CNT" ] ; then
    CFUNC_FILE="cfunc_for_density_fit_cnt.txt"
else
    CFUNC_FILE="cfunc_for_density_fit_erg.txt"
fi

if [ ${NEW_PARAM_FILE} = "T" ]; then

if [ -e "param_beta.txt" ] ; then
    mv param_beta.txt param_beta.txt_bak
fi
BKG_UP=`cat ${SBP_DATA_FILE} | sort -g | tail -n1 | awk '{print $3*0.9}'`
BKG_C=`cat ${NAME}_result.csv | grep ^sbp_c | awk -F"," '{print $2}'`
echo "n01 0.03 0.000 1.001
beta1 0.7 0.1 2.0
rc1 30 1 1000
n02 0.003 0.0 1.0
rc2 200 10 3000
beta2 0.7 0.3 2.0
bkg ${BKG_C} 0 ${BKG_UP}" > param_beta.txt

fi

~/scripts/entropy/fit_sbp_beta.py ${SBP_DATA_FILE} param_beta.txt ${CFUNC_FILE} cfunc_for_chandra_density_fit_cnt.txt ${JSON_FILE} ${z}
