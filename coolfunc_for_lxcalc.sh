#!/bin/sh
#
# unalias -a
#
###########################################################
## Task:                                                 ##
## Calc `cooling function' data according to             ##
## given `temperature profile'                           ##
##                                                       ##
## NOTE:                                                 ##
## given `tprofile': <radius> <temperature>              ##
## calc `cooling function' by invoking `XSPEC'           ##
## using model `wabs*apec'                               ##
##                                                       ##
## LIweitiaNux <liweitianux@gmail.com>                   ##
## August 17, 2012                                       ##
###########################################################

## cmdline arguments {{{
if [ $# -ne 5 ] && [ $# -ne 8 ]; then
    printf "usage:\n"
    printf "    `basename $0` <tprofile> <avg_abund> <nH> <redshift> <coolfunc_outfile> [<e_low> <e_high> <cfunc_erg_outfile>]\n"
    exit 1
fi
base_path=/home/astro/scripts/mass_profile_all/mass_profile
TPROFILE=$1
ABUND_VAL=$2
N_H=$3
REDSHIFT=$4
NORM=`$base_path/calc_distance $REDSHIFT|grep norm|awk '{print $2}'`
FLAG=1
echo $NORM


COOLFUNC_DAT=$5
COOLFUNC_DAT_RATIO=flux_cnt_ratio.txt

if [ ! -r "${TPROFILE}" ]; then
    printf "ERROR: given tprofile '${TPROFILE}' NOT accessiable\n"
    exit 2
fi
if [ ! -r "${ABUND_VAL}" ]; then
    printf "WARNING: '${ABUND_VAL}' NOT accessiable,use average abundance\n"
    FLAG=0
    elif [ ! -r "${N_H}" ]; then
        printf "WARNING: '${N_H}' NOT accessiable, something wrong\n"
fi

[ -e "${COOLFUNC_DAT}" ] && rm -f ${COOLFUNC_DAT}
[ -e "${COOLFUNC_DAT_RATIO}" ] && rm -f ${COOLFUNC_DAT_RATIO}
## arguments }}}

## specify variable name outside while loop
## otherwise the inside vars invisible
XSPEC_CF_XCM="_coolfunc_calc.xcm"
[ -e "${XSPEC_CF_XCM}" ] && rm -f ${XSPEC_CF_XCM}
if [ $FLAG -eq 1 ]; then
E_LOW=0.01
E_HIGH=100
COOLFUNC_DAT_RATIO="cnt_ratio.txt"
## generate xspec script {{{
cat >> ${XSPEC_CF_XCM} << _EOF_
## XSPEC Tcl script
## calc cooling function data
##
## generated by: `basename $0`
## date: `date`

set xs_return_results 1
set xs_echo_script 0
# set tcl_precision 12
dummyrsp .01 100 4096
## set basic data {{{
#set nh ${N_H}
set redshift 0
#set abund_val ${ABUND_VAL}
set norm ${NORM}
## basic }}}

## xspec related {{{
# debug settings {{{
chatter 0
# debug }}}
query yes
abund angr
dummyrsp 0.3 11.0 1024
# load model 'wabs*apec' to calc cooling function
model wabs*apec & 0.01 & 1.0 & 0.5  & \${redshift} & \${norm} & /*
## xspec }}}

## set input and output filename
set tpro_fn "${TPROFILE}"
set apro_fn "${ABUND_VAL}"
set nhpro_fn "${N_H}"
set cf_fn "${COOLFUNC_DAT}"
set cff_fn "${COOLFUNC_DAT_RATIO}"
if { [ file exists \${cf_fn} ] } {
    exec rm -fv \${cf_fn}
}

if { [ file exists \${cff_fn} ] } {
    exec rm -fv \${cff_fn}
}

## open files
set tpro_fd [ open \${tpro_fn} r ]
set apro_fd [ open \${apro_fn} r ]
set nhpro_fd [ open \${nhpro_fn} r ]
set cf_fd [ open \${cf_fn} w ]
set cff_fd [ open \${cff_fn} w ]

## read data from tprofile line by line
while { [ gets \${tpro_fd} tpro_line ] != -1 } {
    # gets one line
    gets \${apro_fd} apro_line
    gets \${nhpro_fd} nhpro_line
    scan \${tpro_line} "%f %f" radius temp_val
    scan \${apro_line} "%f %f" radius_a abund_val
    scan \${nhpro_line} "%f %f" radius_nh nh_val
    #puts "radius: \${radius}, temperature: \${temp_val}"
    # set temperature value
    newpar 1 0
    newpar 2 \${temp_val}
    newpar 3 \${abund_val}
    # calc flux & tclout
    flux 0.01 100
    tclout flux 1
    scan \${xspec_tclout} "%f %f %f %f" cf_data holder holder holder
    #puts "cf_data: \${cf_data}"
    puts \${cf_fd} "\${radius}    \${cf_data}"
#    flux 0.01 100
#    tclout flux 1
#    scan \${xspec_tclout} "%f %f %f %f" cff_data holder holder holder
#    puts \${cff_fd} "\${radius}   \${cff_data}"
}

## close opened files
close \${tpro_fd}
close \${cf_fd}

## exit
tclexit
_EOF_

## extract xcm }}}
elif [ ${FLAG} -eq 0 ]; then 

## generate xspec script {{{
cat >> ${XSPEC_CF_XCM} << _EOF_
## XSPEC Tcl script
## calc cooling function data
##
## generated by: `basename $0`
## date: `date`

set xs_return_results 1
set xs_echo_script 0
# set tcl_precision 12
dummyrsp .01 100 4096
## set basic data {{{
set nh 0
set redshift 0
set abund_val ${ABUND_VAL}
set norm ${NORM}
## basic }}}

## xspec related {{{
# debug settings {{{
chatter 0
# debug }}}
query yes
abund angr
dummyrsp 0.3 11.0 1024
# load model 'wabs*apec' to calc cooling function
model wabs*apec & \${nh} & 1.0 & \${abund_val} & \${redshift} & \${norm} & /*
## xspec }}}

## set input and output filename
set tpro_fn "${TPROFILE}"
set cf_fn "${COOLFUNC_DAT}"
set cff_fn "${COOLFUNC_DAT_RATIO}"
if { [ file exists \${cf_fn} ] } {
    exec rm -fv \${cf_fn}
}

if { [ file exists \${cff_fn} ] } {
    exec rm -fv \${cff_fn}
}

## open files
set tpro_fd [ open \${tpro_fn} r ]
set cf_fd [ open \${cf_fn} w ]
set cff_fd [ open \${cff_fn} w ]

## read data from tprofile line by line
while { [ gets \${tpro_fd} tpro_line ] != -1 } {
    # gets one line
    scan \${tpro_line} "%f %f" radius temp_val
    #puts "radius: \${radius}, temperature: \${temp_val}"
    # set temperature value
    newpar 2 \${temp_val}
    # calc flux & tclout
    flux 0.01 100
    tclout flux 1
    scan \${xspec_tclout} "%f %f %f %f" cf_data holder holder holder
    #puts "cf_data: \${cf_data}"
    puts \${cf_fd} "\${radius}    \${cf_data}"
#    flux 0.01 100.0
#    tclout flux 1
#    scan \${xspec_tclout} "%f %f %f %f" cff_data holder holder holder
#    puts \${cff_fd} "\${radius}   [expr \${cff_data}/\${cf_data}]"
}

## close opened files
close \${tpro_fd}
close \${cf_fd}


## exit
tclexit
_EOF_
## extract xcm }}}
fi
## invoke xspec to calc
printf "invoking XSPEC to calculate cooling function data ...\n"
# xspec - ${XSPEC_CF_XCM}
xspec - ${XSPEC_CF_XCM} > /dev/null

## clean
# if [ -e "${XSPEC_CF_XCM}" ]; then
#     rm -f ${XSPEC_CF_XCM}
# fi

exit 0

