#####
#useage:
#${SCRIP_DIR}/gen_zzh_param.sh  <T_AVE> <REDSHIFT>
#generate param_zzh.txt file
########
CFG_FILE=param_zzh_for_py.txt
#SCRIP_DIR=`echo $0 | awk -F"/gen_zzh_param.sh" '{print $1}'`
#cd $1
#FIT_BASE_DIR="/home/astro/astro_data/mass_sample/entropy"
T_AVE=$1
Z=$2

H0=71
OMEGA_M=0.27
EZ=`echo ${OMEGA_M} ${Z}| awk '{print ($1*(1+$2)^3+1-$1)^0.5}'`
G=6.673e-11 
KPC=3.08562e19 #m
MSUN=2e30 #kg
PI=3.1415926
MP=1.02e-27 #kg
KEV=1.602e-16 #J
#echo  ${PI} ${G} ${MP} ${RHO} ${R500} ${MSUN} ${KPC} ${T_AVE} ${KEV} 
N1=`echo ${PI} ${G} ${MP} ${KEV} ${KPC} ${MSUN} | awk '{print 8*$1*$2*0.61*$3*$6/$4/3/$5}'` #kev*kpc/Msun

R200=`echo ${T_AVE} ${EZ} ${H0} | awk '{print 2.77*$3/70*($1/10)^0.5/$2*1000}'` #kpc
R500=`echo ${R200} | awk '{print $1*0.66}'` #kpc Pratt et al. 2010
M500=`echo ${T_AVE} ${EZ} ${H0} | awk '{print 0.236*1e14*$1^1.76/$2*70/$3}'` #Msun
M200=`echo ${M500} | awk '{print $1/0.717}'` #kpc
C200=`echo ${M200} ${H0} | awk '{print 10^(0.905-0.101*log($1/(10^12*70/$2))/log(10))}'`
RS=`echo ${R200} ${C200} | awk '{print $1/$2}'` #kpc
RS_ERR=`echo ${RS} | awk '{print $1*0.5}'` #kpc
RHO=`echo "${M200} ${RS} ${C200}" | awk '{print $1/4/3.1415926/$2^3/(log(1+$3)-$3/(1+$3))}'` #Msun/kpc^3
RHO_ERR=`echo ${RHO} | awk '{print $1*4}'`
N2=1.0
N2_ERR=0.2

T200=`echo ${R200} ${M200} ${G} ${KPC} ${MSUN} ${KEV}  ${MP}| awk '{print $3*$2*0.61*$7/2/$1*$5/$4/$6}'`


K200=`echo ${T200} ${Z} | awk '{print 362*$1*(sqrt(0.27*(1+$2)^3+0.73))^(-4/3)*(0.27/0.3)^(-4/3)}'`
#K200=`echo ${T200} ${Z} | awk '{print $1*(1.45*(10^(-4))*0.27/0.3*(1+$2)^3)^(-2/3)}'`
#K200=`echo ${R200} ${A1_TMP} ${GAMMA1} ${BKG} | awk '{print $2*($1/100)^$3+$4 }'`
#K200=`echo ${A1_TMP} ${GAMMA1} ${BKG} ${R200} | awk '{print $1*($4/100)^$2+$3}'`
echo "T200=${T200} K200=${K200}"
A0=`echo ${K200} ${R200}  | awk '{print 1.41*$1*(1/$2)^1.1}'`
echo "A0=${A0}"
A0_ERR=`echo ${A0} | awk '{print $1*0.05}'`
GAMMA0=1.1
GAMMA0_ERR=0.3
K0=`echo 0.1 ${K200} | awk '{print $1*$2}'`
K0_LOW=0
K0_UP=`echo 0.35 ${K200} | awk '{print $1*$2}'`

if [ -e ${CFG_FILE} ] ; then
    mv ${CFG_FILE} ${CFG_FILE}_bak
fi
echo "T0    0.00     0.0   0.1    " >>${CFG_FILE}
echo "N1    ${N1}   ${N1}  ${N1}    " >>${CFG_FILE}
echo "N2    ${N2}   ${N2_ERR}   " >>${CFG_FILE}
echo "rs    ${RS}   ${RS_ERR}     " >>${CFG_FILE}
echo "a0    1.0   2.0    " >>${CFG_FILE}
echo "gamma0    1.1   0.2    " >>${CFG_FILE}
echo "k0    320   0   1000    " >>${CFG_FILE}
echo "kmod_a    1.0   1.0    " >>${CFG_FILE}
echo "kmod_b    1.1   0.2    " >>${CFG_FILE}
echo "kmod_c 0.0 0.5" >>${CFG_FILE}
echo "kmod_k0    100   300   " >>${CFG_FILE}
echo "a1 1.0 1.0 " >>${CFG_FILE}
echo "gamma1 1.1 0.55" >>${CFG_FILE}
echo "k1 150 150" >>${CFG_FILE}
echo "N3    1.0    0.95      5.0  " >>${CFG_FILE}
echo "z     ${Z}    ${Z}    ${Z}    " >>${CFG_FILE}
echo "rho     ${RHO}    ${RHO_ERR}         " >>${CFG_FILE}
echo "a     0.18    0.06    T" >>${CFG_FILE}
echo "beta  0.5     0.5     0.5     " >>${CFG_FILE}
echo "nt    0.8     0.25    " >>${CFG_FILE}
echo "nm    0.2     0.0   0.5     " >>${CFG_FILE}
echo "R200  ${R200} ${R200} ${R200} " >>${CFG_FILE}
echo "R500  ${R500} ${R500} ${R500} " >>${CFG_FILE}
echo "x     1.0     1.0       1.0       " >>${CFG_FILE}
echo "T_ave ${T_AVE}" >> ${CFG_FILE}

#cp ${CFG_FILE} ${FIT_BASE_DIR}/$2/
