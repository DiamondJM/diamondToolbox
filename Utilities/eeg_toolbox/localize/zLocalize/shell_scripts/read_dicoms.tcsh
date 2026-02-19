#!/bin/tcsh
# read_dicoms.tcsh
# Outputs nifti and afni versions of dicoms
#

# set directory options
set imgName=`echo $1`
set dicomDir=`echo $2`
set outDir=`echo $3`

# inVol = 1st dicom
set dicomList=`ls ${dicomDir}`
set inVol=`echo ${dicomDir}/$dicomList[1]`

# check if dicoms exist
if (${%inVol} == 0) then
	set runConvert=0
	echo 'No dicoms found'
else
	set runConvert=1
	#echo 'Dicoms found'
endif

# if dicoms exist, run the conversion
if ($runConvert == 1) then
	# out and afni ovls are the imgName
	# mkdir ${outDir}
	echo ${outDir}
	set outVol=`echo ${outDir}/${imgName}.nii`
	set afniVol=`echo ${outDir}/${imgName}`
	
	# do mri_convert to get nii
	echo '-------------------------------------------------------'
	echo '-------------------------------------------------------'
	echo 'Converting dicom to nii via mri_convert (freesurfer)...'
	echo '	Input volume: '$inVol
	echo '	Output volume: '$outVol
	echo "  mri_convert -it dicom -ot nii ${inVol} ${outVol}"
    echo '-------------------------------------------------------'
	echo '-------------------------------------------------------'
    
	mri_convert -it dicom -ot nii ${inVol} ${outVol}

	# eventually prompt to see if this itself worked
	# https://surfer.nmr.mgh.harvard.edu/fswiki/FsFastUnpackData
	# an alternative to mri_convert

	# do 3dcopy to get afni BRIK/HEAD
	echo '-------------------------------------------------------'
	echo '-------------------------------------------------------'
	echo 'Converting nii to brik via 3dcopy (afni)...'
	echo '	Input volume: '$outVol
	echo '	Output volume: '${afniVol}'+orig.BRIK/HEAD'
	echo "  3dcopy $outVol $afniVol"
    echo '-------------------------------------------------------'
	echo '-------------------------------------------------------'
    
	3dcopy $outVol $afniVol	

endif
