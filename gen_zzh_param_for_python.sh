#####
#useage:
#${SCRIP_DIR}/gen_zzh_param.sh the_mass_fitting_dir the_temperature_fitting_dir
#generate param_zzh.txt file
########
CFG_FILE=param_zzh_for_py.txt
SCRIP_DIR=`echo $0 | awk -F"/gen_zzh_param.sh" '{print $1}'`
cd $1
FIT_BASE_DIR="/home/astro/astro_data/mass_sample/entropy"
R500=`cat final_result.txt | grep ^r500 | awk '{print $2}'` #kpc
RS=`cat nfw_param_center.* | grep ^rs | awk '{print $2}' | head -n1` #kpc
RHO_ORI=`cat nfw_param_center.* | grep ^rho0 | awk '{print sqrt($2*$2)}' | head -n1` #Msun/kpc^3
G=6.673e-11 
KPC=3.08562e19 #m
MSUN=2e30 #kg
#T_AVE=`cat ../spc/profile/norm_profile_note.txt | grep ^T_AVE | awk '{print $3}'` #kev
#if [ -z ${T_AVE} ] ; then
#    read -p "T_AVE not found, please input->" -a T_AVE
#fi
PI=3.1415926
MP=1.02e-27 #kg
KEV=1.602e-16 #J
#echo  ${PI} ${G} ${MP} ${RHO} ${R500} ${MSUN} ${KPC} ${T_AVE} ${KEV} 
N1=`echo ${PI} ${G} ${MP} ${KEV} ${KPC} ${MSUN} | awk '{print 8*$1*$2*0.61*$3*$6/$4/3/$5}'` #kev*kpc/Msun
RS_ERR=`echo ${RS} | awk '{print $1*0.2}'` #kpc
RHO=`echo ${RHO_ORI} ${MSUN} ${KPC}| awk '{print $1}'` #Msun/kpc^3
RHO_ERR=`echo ${RHO} | awk '{print $1*0.2}'`
N2=1.0
N2_ERR=0.2
R200=`cat final_result.txt | grep ^r200 | awk '{print $2}'`
M200=`cat final_result.txt | grep ^m200 | awk '{print $2}'`

T200=`echo ${R200} ${M200} ${G} ${KPC} ${MSUN} ${KEV}  ${MP}| awk '{print $3*$2*0.61*$7/2/$1*$5/$4/$6}'`
SBP_CFG=`cat global.cfg | grep ^sbp_cfg | awk '{print $2}'`
Z=`cat ${SBP_CFG} | grep ^z | awk '{print $2}'`

A1_TMP=`cat ../entropy/entropy_result_final.txt | grep ^A | awk  '{print $3}'`
A1_TMP_ERR=`cat ../entropy/entropy_result_final.txt | grep ^A | awk '{print (-$4+$6)/2}'`
KMOD_B=`cat ../entropy/entropy_result_final.txt | grep ^gamma | awk  '{print $3}'`
KMOD_BERR=`cat ../entropy/entropy_result_final.txt | grep ^gamma | awk '{print (-$4+$6)/2}'`
KMOD_A=`echo ${A1_TMP} ${R500} ${KMOD_B}| awk '{print $1*(1/100)^$3}'`
KMOD_AERR=`echo ${KMOD_A} ${A1_TMP} ${A1_TMP_ERR} | awk '{print $1*$3/$2}'`
BKG=`cat ../entropy/entropy_result_final.txt | grep ^bkg | awk '{print $3}'`
KMOD_K0=`echo ${KMOD_A} ${KMOD_B} ${BKG} ${R500} | awk '{print $3}'` # '{print $3/($1*($4/100)^$2+$3)}'`
KMOD_K0ERR=`cat ../entropy/entropy_result_final.txt | grep ^bkg | awk '{print (-$4+$6)/2}'`

K200=`echo ${T200} ${Z} | awk '{print 362*$1*(sqrt(0.27*(1+$2)^3+0.73))^(-4/3)*(0.27/0.3)^(-4/3)}'`
#K200=`echo ${T200} ${Z} | awk '{print $1*(1.45*(10^(-4))*0.27/0.3*(1+$2)^3)^(-2/3)}'`
#K200=`echo ${R200} ${A1_TMP} ${GAMMA1} ${BKG} | awk '{print $2*($1/100)^$3+$4 }'`
#K200=`echo ${A1_TMP} ${GAMMA1} ${BKG} ${R200} | awk '{print $1*($4/100)^$2+$3}'`
echo "T200=${T200} K200=${K200}"
A0=`echo ${K200} ${R200} ${R500} | awk '{print 1.48*$1*(1/$2)^1.22}'`
A0_ERR=`echo ${A0} | awk '{print $1*0.05}'`
GAMMA0=1.15
GAMMA0_ERR=0.04
K0=`echo 0.1 ${K200} | awk '{print $1*$2}'`
K0_LOW=0
K0_UP=`echo 0.35 ${K200} | awk '{print $1*$2}'`
M=`echo ${M200} | awk '{print $1/1e14}'`
M_LOW=`echo ${M} | awk '{print $1*0.8}'`
M_UP=`echo ${M} | awk '{print $1*1.2}'`
#cd ${FIT_BASE_DIR}/$2
if [ -e ${CFG_FILE} ] ; then
    mv ${CFG_FILE} ${CFG_FILE}_bak
fi
echo "T0    0.00     0   0.1    " >>${CFG_FILE}
echo "N1    ${N1}   ${N1}  ${N1}    " >>${CFG_FILE}
echo "N2    ${N2}   ${N2_ERR}   " >>${CFG_FILE}
echo "rs    ${RS}   ${RS_ERR}     " >>${CFG_FILE}
echo "a0    1.0   0.5    " >>${CFG_FILE}
echo "gamma0    ${GAMMA0}   0.4    " >>${CFG_FILE}
echo "k0    ${K0}   ${K0_LOW}   ${K0_UP}    " >>${CFG_FILE}
echo "kmod_a    ${KMOD_A}   0.5    " >>${CFG_FILE}
echo "kmod_b    ${KMOD_B}   0.4    " >>${CFG_FILE}
echo "kmod_c 0.5 1" >>${CFG_FILE}
echo "kmod_k0    ${KMOD_K0}   ${KMOD_K0ERR}   " >>${CFG_FILE}
echo "a1 ${KMOD_A} ${KMOD_AERR} " >>${CFG_FILE}
echo "gamma1 ${KMOD_B} ${KMOD_BERR}" >>${CFG_FILE}
echo "k1 ${KMOD_K0} ${KMOD_K0ERR}" >>${CFG_FILE}
echo "N3    0.0    0.0      0.9  " >>${CFG_FILE}
echo "z     ${Z}    ${Z}    ${Z}    " >>${CFG_FILE}
echo "rho     ${RHO}    ${RHO_ERR}         " >>${CFG_FILE}
echo "a     0.18    0.06    T" >>${CFG_FILE}
echo "beta  0.5     0.5     0.5     " >>${CFG_FILE}
echo "nt    0.8     0.25    " >>${CFG_FILE}
echo "nm    0.2     0.0   0.5     " >>${CFG_FILE}
echo "R200  ${R200} ${R200} ${R200} " >>${CFG_FILE}
echo "R500  ${R500} ${R500} ${R500} " >>${CFG_FILE}
echo "x     1.0     1.0       1.0       " >>${CFG_FILE}

cp ${CFG_FILE} ${FIT_BASE_DIR}/$2/
