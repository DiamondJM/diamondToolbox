#!/bin/tcsh -f


set cnt = 1

if ("$argv[$cnt]" == "-CT_path") then
@ cnt ++
set ct = "$argv[$cnt]"
endif

@ cnt ++
if ("$argv[$cnt]" == "-T1_path") then
@ cnt ++
set t1 = "$argv[$cnt]"
endif

# master grid added
cp $ct .
set ct = `3dinfo -prefix $ct`
echo "ct: $ct"
set root = `3dinfo -prefix_noext $ct | awk '{$1=$1};1'`
echo "ct_root = $root"


@Align_Centers -base $t1 -dset $ct -cm_no_amask
# Change note: Original didn't have the -cm_no_amask
# by hand: @Align_Centers -base temp_ANAT.nii -dset ct_implant.nii -cm_no_amask -overwrite


3dresample -input "${root}_shft.nii" -prefix ${root}_res_shft.nii  -master $t1  -dxyz 1 1 1 -rmode NN
# by hand: 3dresample -input ct_implant_shft.nii -prefix ct_implant_res_shft.nii -master temp_ANAT.nii -dxyz 1 1 1 -rmode NN -overwrite

align_epi_anat.py -dset1 $t1 -dset2  ${root}_res_shft.nii -dset1_strip None -dset2_strip None -dset2to1 -suffix _al  -feature_size 1  -overwrite -cost nmi -giant_move -rigid_body >> status.txt
# by hand: align_epi_anat.py -dset1 temp_ANAT.nii -dset2 ct_implant_res_shft.nii -dset1_strip None -dset2_strip None -dset2to1 -suffix -al -overwrite -feature_size 1 -cost nmi -giant_move -rigid_body

3dcopy  ${root}_res_shft_al+orig CT_RAI_res_al.nii

3dcopy $t1 ./temp_ANAT.nii

afni -com "SWITCH_UNDERLAY temp_ANAT.nii" -com "SWITCH_OVERLAY CT_RAI_res_al.nii"



HELP:
echo ""
echo "Usage: `basename $0` <-CT_path BASE> <-T1_path DSET> [-] "
echo "                     "
echo ""
echo "   Run inside /data/CT/coregistration"
echo "   tcsh @alignCTtoT1 -CT_path path_to_CT -T1_path path_to_anatomy"
echo "   "
echo ""





#set input = $2
#set CTdir = `dirname $input`
#set T1dir = `dirname $input`
#setenv AFNI_DECONFLICT OVERWRITE
#set name = "$input"

#read test

#afni $T1dir/ANAT*.nii  -dset $CTdir/CT_highresRAI.nii

#@Align_Centers -base $T1dir/ANAT_????.nii  -dset $CTdir/CT_highresRAI.nii

#align_epi_anat.py -dset1 $T1dir -dset2 $CTdir  -dset1_strip None -dset2_strip None -dset2to1 -suffix _al2ct_nmi_noclip  -feature_size 1  -overwrite -cost nmi -giant_move -master_dset2 CT_highresRAInoclip_shft.nii

#afni $T1dir  CT_highresRAInoclip_shft_al.nii

