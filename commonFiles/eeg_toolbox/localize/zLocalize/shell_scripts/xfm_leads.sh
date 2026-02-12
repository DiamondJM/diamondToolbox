#!/bin/sh
# xfm_leads.sh
# 

# inputs
workingDir=$1
ctFile=$2
xfm=$3
vox=$4
vox2mm=${5-"vox2mm.1D"}
mm=${6-"mm.txt"}
elec=${7-"elec"}

# navigate to working directory
cur_pwd=`pwd`
cd ${workingDir}

# convert from voxel (ijk) to xyz_mm (RAI) (if non-null vox file passed)
echo $vox
if [ -n "$vox" ]; then
  echo voxel
  cat_matvec ${ctFile}'::IJK_TO_DICOM_REAL' > $vox2mm
  Vecwarp -matvec $vox2mm -input ${vox} > $mm
fi

# create electrogrid from this
rm "${elec}.spec" "${elec}.gii"
@ElectroGrid -with_markers -prefix $elec -coords $mm

# Apply transform to electrodes
ConvertSurface -novolreg -sv ${ctFile} -i $elec.gii -ixmat_1D ${xfm} -o_gii $elec.gii -overwrite

cd $cur_pwd
