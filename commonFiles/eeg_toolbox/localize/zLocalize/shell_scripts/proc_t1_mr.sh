#!/bin/bash
# prepare a t1 mri for registration

# Notes:
#  - try masking out the interior of the intensity-inverted MR
#  - not sure if we need featuresize 1. Have had success without it

mr=$1                   # name of MR file. (e.g. if filename=mr_pre.nii, this should be mr_pre)
workDir=${2-`pwd`}              # where to do this work. ${mr}.nii needs to be in this dir
invert=${3-0}

# Mike commented out this stuff 5/23/17 - not sure why CT/alignment got rolled in here
#ct=$2
#workDir=$3              # where to do this work. ${mr}.nii needs to be in this dir
#resampOrient=${4-RAI}   # orient CT like this
#invert=${5-0}           # if 1, invert the MR
#cost=${6-lpc}           # cost function of alignment
#removeAir=${7-0}        # if 1, use @NoisySkullStrip to remove CT air vox


# navigate to work directory
old_dir=`pwd`
cd $workDir
echo "working dir: $workDir"

# ---------------------
# --- MR processing ---
# ---------------------
# output: _do
# get rid of obliquity, to keep life simple

    cmd="3dWarp -deoblique -prefix ${mr}_do ${mr}.nii"
    echo '-----------------------'
    echo $cmd
    echo '-----------------------'
    $cmd

# output: _do_amd
# Brain Only mask. Mike believes this is done so that 3dSkullStrip
# will be faster in the next step
if [ ! -f "${mr}_do_amd+orig.BRIK" ]; then
    cmd="3dAutomask -apply_prefix ${mr}_do_amd -dilate 3 ${mr}_do+orig"
    echo $cmd
    $cmd
fi


if [ $invert -ne 0 ]; then
    #  output: _do_amd_inv
    # Invert the (non-skull-stripped) MR
    cmd="3dcalc -a ${mr}_do_amd+orig. -expr 'step(a)*(-a)' -prefix ${mr}_do_amd_inv"
    echo '-----------------------'
    echo $cmd
    echo '-----------------------'
    3dcalc -a ${mr}_do_amd+orig. -expr 'step(a)*(-a)' -prefix ${mr}_do_amd_inv
    dset_mr=${mr}_do_amd_inv+orig
else
    # output: _do_amd.ns
    # strip skull
    cmd="3dSkullStrip -input ${mr}_do_amd+orig. -prefix ${mr}_do_amd.ns -blur_fwhm 2"
    echo '-----------------------'
    echo $cmd
    echo '-----------------------'
    $cmd
    dset_mr=${mr}_do_amd.ns+orig
fi

cd $old_dir


