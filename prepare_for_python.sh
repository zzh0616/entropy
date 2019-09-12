# parpare for entropy fitting
# should be execute in the mass dir
# $1 is the name of the dir below the entropy dir
PROFILE_FILE=tcl_temp_profile.txt
DIR=$1
SCRIP_DIR=/mnt/data/mass_sample/scrips
SBP_CFG=`cat global.cfg | grep ^sbp_cfg | awk '{print $2}'`
Z=`cat ${SBP_CFG} | grep ^z | awk '{print $2}'`
KPC_PERAL_PIXEL=`${SCRIP_DIR}/mass_profile_all/mass_profile/calc_distance $Z | grep kpc_per_pixel | awk '{print $2}'`
NAME=`${SCRIP_DIR}/lyt/analyze_path.sh $(pwd) | grep ^NAME | head -n1 | awk '{print $2}'`
BASE_DIR=/mnt/data/mass_sample/entropy
mkdir /mnt/data/mass_sample/entropy/${DIR}/${NAME} 2>/dev/null
${SCRIP_DIR}/entropy/gen_zzh_param_for_python.sh `pwd` ${DIR}/${NAME}
cd ../spc/profile
cat ${PROFILE_FILE} |  awk -v "a=${KPC_PERAL_PIXEL}" '{print $1*a,$2*a,$3,$4}' > /mnt/data/mass_sample/entropy/${DIR}/${NAME}/${NAME}.txt 
cd ../../mass

MASS_CFG="global.cfg"
SBP_CFG=`cat ${MASS_CFG} | grep ^sbp_cfg | awk '{print $2}'`
SBP_DAT=`cat ${MASS_CFG} | grep ^radius_sbp_file | awk '{print $2}'`
CM_PER_PIXEL=`cat ${SBP_CFG} | grep ^cm_per_pixel | awk '{print $2}'`
KPC_PER_PIXEL=`echo ${CM_PER_PIXEL} | awk '{print $1/3.0857e21}'`
cp ${MASS_CFG} ${BASE_DIR}/${DIR}/${NAME}
cp ${SBP_CFG} ${BASE_DIR}/${DIR}/${NAME}
cp ${SBP_DAT} ${BASE_DIR}/${DIR}/${NAME}

R200m=`${SCRIP_DIR}/entropy/calc_r200m.sh `
cd ${BASE_DIR}/${DIR}/${NAME}
sed -i '/^#/'d ${SBP_DAT}
cat ${SBP_DAT} | awk -v a=${KPC_PER_PIXEL} '{print $1*a,$2*a,$3,$4}' >tmp
mv tmp ${SBP_DAT}
echo ${R200m} >r200m.txt
