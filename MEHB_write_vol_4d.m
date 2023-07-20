function V = MEHB_write_vol_4d(V,Y)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_write_vol(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
% Y          - a 3D or 4D matrix containing the image voxels
% V (output) - data structure after modification for writing.
%
% Note that if there is no 'pinfo' field, then SPM will figure out the
% max and min values from the data and use these to automatically determine
% scalefactors.  If 'pinfo' exists, then the scalefactor in this is used.
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_vol.m 5731 2013-11-04 18:11:44Z guillaume $

ind  = cat(1,V.n);
if isfield(V,'private')
    N    = cat(1,V.private);
    mat = N(1).mat;
else
    mat = V(1).mat;
end

%-Input: fname
%--------------------------------------------------------------------------
fname = V(1).fname;

%-Input: RT
%--------------------------------------------------------------------------
RT = NaN;
if isnan(RT) && isfield(V,'private')
    if isfield(V(1).private,'timing') && isfield(V(1).private.timing,'tspace')
        RT = V(1).private.timing.tspace;
    end
end

dt = V(1).dt;

%-Create and write 4D volume image
%==========================================================================
spm_unlink(fname);

%-Create NifTI header
%--------------------------------------------------------------------------
ni         = nifti;
ni.dat     = file_array(fname,...
                        size(Y),...
                        dt,...
                        0,...
                        1,...
                        0);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = V(1).descrip;
if ~isnan(RT)
    ni.timing = struct('toffset',0, 'tspace',RT);
end
create(ni);

%-Write 4D data
%--------------------------------------------------------------------------
spm_progress_bar('Init',size(ni.dat,4),'Saving 4D image','Volumes Complete');
for i=1:size(Y,4)
    ni.dat(:,:,:,i) = Y(:,:,:,i);
    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

%-Remove .mat file if present and not necessary
%--------------------------------------------------------------------------
matfname = spm_file(fname,'ext','mat');
if spm_existfile(matfname)
    M = load(matfname);
    if isequal(fieldnames(M),{'mat'}) % contains only 'mat'
        if sum(sum(M.mat(:,:,1).^2))==0
            M.mat(:,:,1) = N(1).mat;
        end
        if sum(sum(diff(M.mat,1,3).^2))<1e-8
            spm_unlink(matfname);
        end
    end
end

%-Return spm_vol structure
%--------------------------------------------------------------------------
V = spm_vol(fname);