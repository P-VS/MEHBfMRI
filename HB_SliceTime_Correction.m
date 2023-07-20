function outfiles = HB_SliceTime_Correction(funcfile,jsonfile)

jsondat = fileread(jsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;
slicetimings = jsondat.SliceTiming;

[filepath,name,ext] = fileparts(funcfile);
outfiles = fullfile(filepath,['a' name '.nii']);

matlabbatch{1}.spm.tools.MEHBfmri.stcor.scans(1) = {funcfile};%cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
matlabbatch{1}.spm.tools.MEHBfmri.stcor.SliceT = slicetimings;
matlabbatch{1}.spm.tools.MEHBfmri.stcor.TR = tr;
matlabbatch{1}.spm.tools.MEHBfmri.stcor.refslice = 1;
matlabbatch{1}.spm.tools.MEHBfmri.stcor.prefix = 'a';

spm_jobman('run', matlabbatch);

end