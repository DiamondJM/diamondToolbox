#!/bin/bash
# prepare ct_implant for registration

myCT=$1 # name of CT file (eg if filename=my_ct.nii, this should be my_ct)
myMR=$2 # name of MR file. (e.g. if filename=mr_pre.nii, this should be mr_pre)
workDir=$3 # where to do this work. ${mrName}.nii needs to be in this dir
resampOrient=${4-RAI} # orient CT like this

# navigate to temp (intermediates) directory
old_pwd=`pwd`
cd $workDir
echo "working dir: $workDir"

# reorient CT to match coords_ctpost_vox.csv
cp ${myCT}.nii noreorient_${myCT}.nii
cmd="3dresample -orient ${resampOrient} -prefix ${myCT} -inset noreorient_${myCT}.nii -overwrite"
echo '-----------------------'
echo $cmd
echo '-----------------------'
$cmd

#Shift CT to align better with anatomy of interest then resample it
cmd="@Align_Centers -base ${myMR}+orig -dset ${myCT}+orig"
echo '-----------------------'
echo $cmd
echo '-----------------------'
$cmd

cmd="3dresample -prefix ${myCT}_sh2${myMR}.rs -master ${myMR}+orig -input ${myCT}_shft+orig"
echo $cmd
$cmd

# Clean up outer perimeter of CT
# Remove "air voxels"
# Outputs a .ma file
cmd="@NoisySkullStrip -input ${myCT}_sh2${myMR}.rs+orig."
echo '-----------------------'
echo $cmd
echo '-----------------------'
$cmd

cd $old_pwd
