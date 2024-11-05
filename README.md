# MEHBfMRI
SPM tools for preprocessing fMRI series scanned using Simultaneous MultiSlice (SMS) / MultiBand / HyperBand and Multi Echo

This toolbox is no longer acctively updated, but a new automatic SPM preprocessing including these scripts for single- and multi-echo fMRI with slice time correction, weighted echo combination and (ME-) ICA-based denoising can bbe found at https://github.com/P-VS/AutoSPM12_processing.

This toolbox is based on the SPM12 slice timing correction tool. Please make sure SPM12 (freely downloadable from https://www.fil.ion.ucl.ac.uk/spm/) is properly installed and running in Matlab.

Installation
Copy these files in a folder MEHBfmri subfolder into the SPM12/toolbox folder.

Launch the toolbox in SPMT12 -> BATCH -> SPM -> tools -> Multi echo & HyperBand.

The toolbox contains 2 subfunctions:
1. Combine echoes for Multi-Echo fMRI
In multi echo fMRI, up to 5 echoes are measured per TR. Rather than fitting T2* per voxel (which is time consuming), the echo images can be combined per TR as S(t)=sum(wi.Si(t))/sum(wi) with i the echo number (1,2..), Si the ith echo image and wi the ith weighting factor. wi can be chosen as
  1. AVE: wi=1 (simple averaging)
  2. BS: wi=TEi (BOLD sensitivity)
  3. tSNR: wi=tSNRi (temporal SNR)
  4. tBS: wi=tSNRi*TEi (temporal BOLD sensitivity)
  5. T2* weighted: T2* weighted based on a T2* map determined per dynamic (wi(t)=TEi * exp(-TEi/T2*(t)))
  6. T2* mapping: T2* mapping per dynamic

2. Slice time correction for HyperBand
Slice timing is done based on a matrix containing the slice timings for each file. These slice timings can be found in the json file. This json file is created during the conversion of your data from dicom to nifti using toolboxes such as dicm2nii (https://github.com/xiangruili/dicm2nii).

This toolbox is written by dr. Peter Van Schuerbeek from UZ Brussel (VUB, Belgium). 
