# Script for preprocessing a T1-MRI
#   when a T2 is included with a base T1 to improve FreeSurfer's
#   recon-all surface construction. See the AFNI recommendations
#   for more details: 
#	https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/fatcat_prep/Prepreprocessing_I.html#fp-align-anat-pair
#
# Input:
#	path_t1		- all paths absolute or relative to working_dir
#	path_t2
#	path_t1_out - output file path
#	working_dir - script will cd here first
#
# Revision History
#   10/2017 MST - Created

# note that slice thickness should be in [0.75,3] mm

# Align T1 to T2 (via fat_proc_align_anat_pair), rigid transformation

path_t1=$1
path_t2=$2
path_t1_out=$3
working_dir=$4

cur_pwd=`pwd`
cd $working_dir

cp $path_t1 "${path_t1}.orig"

fat_proc_align_anat_pair \
	-in_t1w $path_t1 \
	-in_t2w $path_t2 \
	-prefix $path_t2_out \
	-out_t2w_grid

cd $cur_pwd
