function MEHBfmri = tbx_cfg_MEHBfmri

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MEHBfMRI')); end


% Echo data
%---------------------------------------
func              = cfg_files;
func.tag          = 'func';
func.name         = 'Functional data';
func.filter       = 'image';
func.ufilter      = '.*';
func.num          = [1 Inf];
func.help         = {'Select the fMRI data.'};

te             = cfg_entry;
te.tag         = 'te';
te.name        = 'TE (ms)';
te.help        = {'TE in ms'};
te.val         = {};
te.strtype     = 'r';
te.num         = [1 1];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
tedat         = cfg_branch;
tedat.tag     = 'tedat';
tedat.name    = 'Data';
tedat.val     = {func te};
tedat.help    = {'Data for each echo.'};

%--------------------------------------------------------------------------
% TE Data
%--------------------------------------------------------------------------
tefmri         = cfg_repeat;
tefmri.tag     = 'data';
tefmri.name    = 'TE data';
tefmri.help    = {'List of echoes'};
tefmri.values  = {tedat};
tefmri.num     = [1 Inf];

%--------------------------------------------------------------------------
% Method
%--------------------------------------------------------------------------
method         = cfg_menu;
method.tag     = 'method';
method.name    = 'Method';
method.help    = {
                  ['The weighting method to combine the different echo images ' ...
                  '(w1*S1 + w2*S2 + ...)/(w1 + w2 + ...)']
                  'AVE: simple averaging (wi=1)' 
                  'BS: BOLD sensitivity (wi=TEi)'
                  'tSNR: temporal SNR (wi=tSNRi)'
                  'tBS: temporal BOLD sensitivity (wi=tSNRi * TEi)'
}';
method.labels = {
                 'AVE'
                 'BS'
                 'tSNR'
                 'tBS'
}';
method.values = {0 1 2 3};
method.val    = {3};

%--------------------------------------------------------------------------
% Corregistration
%--------------------------------------------------------------------------
correg         = cfg_menu;
correg.tag     = 'correg';
correg.name    = 'Corregistration per echo?';
correg.help    = {
                  ['Should each echo be corregistered to the first echo per dynamic?']
}';
correg.labels = {
                 'Yes'
                 'No'
}';
correg.values = {1 0};
correg.val    = {0};

% Main structure ME-fMRI
%---------------------------------------
mems_fmri           = cfg_exbranch;
mems_fmri.tag       = 'mems_fmri';
mems_fmri.name      = 'Combine echoes for Muli-Echo fMRI';
mems_fmri.help      = {'This function combines the multiple echoes from a ME-fMRI experiment'};
mems_fmri.val       = {tefmri method correg};
mems_fmri.prog      = @(job)vout_memsfmri('run',job);
mems_fmri.vout      = @(job)vout_memsfmri('vout',job);


%--------------------------------------------------------------------------
% scans Session
%--------------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Session';
scans.help    = {'Select images to slice-time correct.'};
scans.filter  = 'image';
scans.ufilter = '.*'; %
scans.num     = [1 Inf];
scans.preview = @(f) spm_check_registration(char(f));


%--------------------------------------------------------------------------
% Slice time correction
%--------------------------------------------------------------------------
SliceT             = cfg_entry;
SliceT.tag         = 'SliceT';
SliceT.name        = 'Slice Timings';
SliceT.help        = {'Give the slice timings (see json file)'};
SliceT.val         = {};
SliceT.strtype     = 'r';
SliceT.num         = [1 Inf];

TR         = cfg_entry;
TR.tag     = 'TR';
TR.name    = 'TR (ms)';
TR.help    = {'Enter the TR in ms.'};
TR.strtype = 'r';
TR.num     = [1 1];

refslice         = cfg_entry;
refslice.tag     = 'refslice';
refslice.name    = 'Reference Slice';
refslice.help    = {'Enter the reference slice.'
                    ''
                    'If slice times are provided instead of slice indices in the previous item, this value should represent a reference time (in ms) instead of the slice index of the reference slice.'};
refslice.strtype = 'r';
refslice.num     = [1 1];

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the slice-time corrected image file(s).','Default prefix is ''a''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('slicetiming.prefix', val{:});


% Main structure Slice Time Correction
%---------------------------------------
stcor           = cfg_exbranch;
stcor.tag       = 'stcor';
stcor.name      = 'Slice time correction for HyperBand';
stcor.help      = {'This function is the addapted slice time correction for when using HyperBand'};
stcor.val       = {scans SliceT TR refslice prefix};
stcor.prog      = @(job)vout_stcor('run',job);
stcor.vout      = @(job)vout_stcor('vout',job);

% Main structure of the toolbox
%---------------------------------------
MEHBfmri         = cfg_choice;
MEHBfmri.tag     = 'MEHBfmri';
MEHBfmri.name    = 'Multi echo & HyperBand';
MEHBfmri.help    = {'This toolbox is meant to do the extra processing steps when using multi-echo fMRI or HypperBand'};
MEHBfmri.values  = {mems_fmri stcor};

function out = vout_memsfmri(cmd,job)

switch lower(cmd)
    case 'run'
        out.files=mems_fmri_run(job);
    case 'vout'
        out(1)           =cfg_dep;
        out(1).sname     =sprintf('MEMS fMRI files');
        out(1).src_output=substruct('.','files');
        out(1).tgt_spec  =cfg_findspec({{'filter','image','strtype','e'}});
end

function out = vout_stcor(cmd,job)

switch lower(cmd)
    case 'run'
        out.files=stcor_run(job);
    case 'vout'
        out(1)           =cfg_dep;
        out(1).sname     =sprintf('Slice time correction files');
        out(1).src_output=substruct('.','files');
        out(1).tgt_spec  =cfg_findspec({{'filter','image','strtype','e'}});
end