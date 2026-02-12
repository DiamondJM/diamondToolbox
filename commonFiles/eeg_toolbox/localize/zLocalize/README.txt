zLocalize (v1.1)

I. CONTENTS

  core: 	main functions for al parts of localization pipeline. These can be used independently
  helper: 	necessary helper functions called by core
  localizer: 	contains the driver, localizer_run.m, which runs our customized, full localization pipeline
  packages: 	external packages
  shell_scripts: non-matlab shell scripts


II. DEPENDENCIES

  The zLocalization pipeline depends upon a number of other (freely available except MATLAB) software packages:  

  MATLAB - Core functionality. zLocalize was testing on version R2017a.

  AFNI/SUMA - Imaging. Available at https://afni.nimh.nih.gov/download.
	This package is required because the MATLAB pipeline utilizes AFNI binary files to co-register the MR and CT
	scan, process DICOMS, normalize the surface mesh,and transform the CT coordinates into MR space.

  FreeSurfer - Surface-based neuroimaging. Available at https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall.
	This package is required to create a surface from the MRI.

  Please note: zLocalize has only been tested on MacOS operating systems.


III. REQUIRED RAW DATA

  NIFTI or DICOM files - read_dicoms.m requires you have either the raw DICOMS or a NIFTI representation
  			of both a high resolution pre-operative T1-weighted MRI and a post-operative CT

  element_info table - projection.m requires that you have a table (usually a CSV) which stores the
                        dimensions and “cut” properties of every piece of all subdural grid/strips

  CT voxel coordinates - xfm_leads.m requires that you have a table of RAI ijk-voxel (default) or millimeter 
			xyz-coordinates of every electrode on the CT

*See function documentation for more information
          
                             
INSTRUCTIONS

  Any of the functions in the core folder of this package may be used and incorporated into your
  own custom electrode localization pipeline. The localizer_run_template.m function is the driver that we
  have written to execute the core functions in sequence and manage inputs and outputs for
  a sample pipeline.






  
