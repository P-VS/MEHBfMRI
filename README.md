# MEHB-fMRI
SPM tools for preprocessing fMRI series scanned using HyperBand and Multi Echo

This toolbox is based on the SPM12 slice timing correction tool. Please make sure SPM12 (freely downloadable from https://www.fil.ion.ucl.ac.uk/spm/) is properly installed and running in Matlab.

Installation
Coppy the folder MEHBfmri into the SPM12/toolbox folder.

Launch the toolbox in SPMT12 -> BATCH -> SPM -> tools -> Multi echo & HyperBand.

The toolbox contains 2 subfunctions:
1. Combine echoes for Multi-Echo fMRI
In multi echo fMRI, up to 5 echoes are measured per TR. Rather than fitting T2* per voxel (which is time consuming), the echo images can be combined per TR as S(t)=sum(wi.Si(t)) with i the echo number (1,2..), Si the ith echo image and wi the ith weighting factor. wi can be chosen as
  1. AVE: wi=1 (simple averaging)
  2. BS: wi=TEi (BOLD sensitivity)
  3. tSNR: wi=tSNRi (temporal SNR)
  4. tBS: wi=tSNRi*TEi (temporal BOLD sensitivity)

2. Slice time correction for HyperBand
Slice timing is done based on a matrix containing the slice timings for each file. These slice timings can be found in the json file. This json file is created during the conversion of your data from dicom to nifti using toolboxes such as dicm2nii (https://github.com/xiangruili/dicm2nii).
