#!/bin/sh
# align.sh performs an affine alignment

# Notes:
#  - try masking out the interior of the intensity-inverted MR
#    by using an inverted step function on the skull stripped brain,
#	 then aligning with lpa
#  - not sure if we need featuresize 1. Have had success without it

mr=$1                   # name of MR file. (e.g. if filename=mr_pre.nii, this should be mr_pre)
ct=$2					# name of CT file (eg if filename=ct_implant.nii, this should be ct_implant)
workDir=$3              # where to do this work. ${mr}.nii needs to be in this dir
resampOrient=${4-RAI}   # orient CT like this
cost=${5-lpc}           # cost function of alignment
invert=${6-0}           # if 1, invert the MR
removeAir=${7-0}        # if 1, use @NoisySkullStrip to remove CT air vox
run=${8-1} 				# if 0, don't run; just echo commands

# navigate to work directory
old_dir=`pwd`
cd $workDir
echo "working dir: $workDir"


# ---------------------
# --- MR processing ---
# ---------------------
# output: _do
# deoblique
cmd="3dWarp -deoblique -prefix ${mr}_do ${mr}.nii"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd

# output: _do_amd
# Brain Only mask. Mike believes this is done so that 3dSkullStrip
# will be faster in the next step
cmd="3dAutomask -apply_prefix ${mr}_do_amd -dilate 3 ${mr}_do+orig"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd

# output: _do_amd.ns
# strip skull
cmd="3dSkullStrip -input ${mr}_do_amd+orig. -prefix ${mr}_do_amd.ns -blur_fwhm 2"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd

if [ "$invert" = 1 ]; then
	# output: _do_amd_inv
	# Invert the (non-skull-stripped) MR, and mask out the (skull-stripped) brain
	# Note that 1-step(b) create a mask that is 0 on the interior of the skull, and 1 elsewhere
	cmd="3dcalc -a ${mr}_do_amd+orig -b ${mr}_do_amd.ns+orig -expr 'step(a)*(-a)*(1-step(b))' -prefix ${mr}_do_amd_inv"
	echo ""
	echo $cmd
	echo ""
	if [ "$run" = 1 ]; then
		3dcalc -a ${mr}_do_amd+orig -b ${mr}_do_amd.ns+orig -expr 'step(a)*(-a)*(1-step(b))' -prefix ${mr}_do_amd_inv
	fi

	# Note - this would be the invert command without the brain mask:
	# cmd="3dcalc -a ${mr}_do_amd+orig. -expr 'step(a)*(-a)' -prefix ${mr}_do_amd_inv"

	dset_mr=${mr}_do_amd_inv
else
	dset_mr=${mr}_do_amd.ns
fi


# ---------------------
# --- CT processing ---
# ---------------------

# reorient CT to match coords_ctpost_vox.csv
[ "$run" = 1 ] && cp ${ct}.nii noreorient_${ct}.nii
cmd="3dresample -orient ${resampOrient} -prefix ${ct} -inset noreorient_${ct}.nii -overwrite"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd

#Shift CT to align better with anatomy of interest then resample it
cmd="@Align_Centers -base ${mr}_do+orig -dset ${ct}+orig"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd
xfm_shft=${ct}_shft.1D

cmd="3dresample -prefix ${ct}_sh2${mr}_do.rs -master ${mr}_do+orig -input ${ct}_shft+orig"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd


if [ "$removeAir" = 1 ]; then
	# Clean up outer perimeter of CT
	# Remove "air voxels"
	# Outputs a .ma file
	cmd="@NoisySkullStrip -input ${ct}_sh2${mr}_do.rs+orig."
	echo ""
	echo $cmd
	echo ""
	[ "$run" = 1 ] && $cmd
	dset_ct=${ct}_sh2${mr}_do.rs.ma
else
	dset_ct=${ct}_sh2${mr}_do.rs
fi

# ------------------------
# --- Affine transform ---
# ------------------------
suffix=_XFMTO_${cost}_${mr}_do

cmd="align_epi_anat.py -dset1 ${dset_mr}+orig -dset2 ${dset_ct}+orig -dset2to1 -dset1_strip None -dset2_strip None"
cmd="$cmd -suffix $suffix -cost $cost -deoblique off -Allineate_opts \"-twopass -nomask\" -overwrite"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd

xfm_af=${dset_ct}${suffix}_mat.aff12.1D
xfm_full=full_${xfm_af}
cmd="cat_matvec -ONELINE $xfm_af $xfm_shft > $xfm_full"
echo ""
echo $cmd
echo ""
[ "$run" = 1 ] && $cmd


cd $old_dir

