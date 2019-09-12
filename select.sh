SCRIP_DIR=/mnt/data/mass_sample/scrips/mass_profile_all/mass_profile
R500=`cat final_result.txt | grep ^r500 | awk '{print $2}'`
cd ../spc/profile
Z=`cat deproj*.log | grep redshift | sort -g | tail -n1 |  awk '{print $4}'`
#echo $Z
A=`cat tcl_temp_profile.txt | wc -l`
KPC_PER_PIXEL=`${SCRIP_DIR}/calc_distance $Z | grep ^kpc_per_pixel | awk '{print $2}'`
#R_MAX_PIXEL=`cat tcl_temp_profile.txt | tail -n1 | awk '{print $1}'`
cd ../../img
SBP_R_MAX_PIXEL=`cat sbprofile.txt | tail -n1 | awk '{print $1}'`
SBP_R_MAX_KPC=`echo ${SBP_R_MAX_PIXEL} ${KPC_PER_PIXEL} | awk '{print $1*$2}'`
FLAG_R=`echo "${R500}*0.93 < ${SBP_R_MAX_KPC}" | bc -l` 
cd ../mass
if [ $A -eq 6 ] ; then
   pwd #>> /mnt/data/mass_sample/entropy/sample_possible.lis
else
    if [ $A -gt 4 ] ; then 
        if [ ${FLAG_R} -eq 1 ]; then
#            TEMP_INNER=`cat tcl_temp_profile.txt | head -n1 | awk '{print $3}'`
#            TEMP_HIGH=` cat tcl_temp_profile.txt | awk '{print $3}' | sort -g | tail -n1`
#        RATE=`echo " ${TEMP_INNER} ${TEMP_HIGH}" | awk '{print $1/$2}'`
#        FLAG_T=`echo " ${RATE} < 0.8 " | bc -l`
#        if [ ${FLAG_T} -lt 10 ] ; then
            pwd >> /mnt/data/mass_sample/entropy/sample_r500.lis
#        fi
    fi
fi
fi
