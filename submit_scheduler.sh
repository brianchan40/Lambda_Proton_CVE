#!/bin/tcsh
#daynum = today's date

#set daynum=$1
#set cen = $1
set opt_weight = $2
set lamtype = $1
set sys_err_opt = $3
set OUTPUT="/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/output/"

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/ ) then
echo "dir: ${OUTPUT}/ already exists"
else
mkdir ${OUTPUT}/
endif

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/${lamtype}_debug ) then
echo "dir: ${OUTPUT}/${lamtype}_debug already exists"
rm -r ${OUTPUT}/${lamtype}_debug
endif
mkdir ${OUTPUT}/${lamtype}_debug

#set up data directory
if ( -e ${OUTPUT}/${lamtype} ) then
echo "dir: ${OUTPUT}/${lamtype} already exists"
rm -r ${OUTPUT}/${lamtype}
endif
mkdir ${OUTPUT}/${lamtype}

#submit the Scheduler template
star-submit-template -template ./Analysis.xml -entities output=${OUTPUT},opt_weight=${opt_weight},sys_err_opt=${sys_err_opt},lamtype=${lamtype}

echo "submitting job done"
