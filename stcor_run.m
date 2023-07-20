function [tout]=stcor_run(job)

tout=cell(numel(job.scans),1);

[pth,name,ext] = fileparts(job.scans{1});

Vin     = spm_vol(fullfile(pth,[name '.nii']));
nslices = Vin(1).dim(3);

refvol = spm_read_vols(Vin(1));
mask = MEHB_mask(refvol);

SliceT = job.SliceT;

TR=job.TR;
if TR>10
    TR=TR/1000;
end
    
if nslices ~= numel(SliceT)
    error('Mismatch between number of slices and length of ''Slice timings'' vector.');
end

%-Slice timing correction
%==========================================================================
if numel(job.scans)>1
    for i=1:numel(job.scans)
        Vin(i)   = spm_vol(job.scans{i});
    end
end
nimgo = numel(Vin);
nimg  = 2^(floor(log2(nimgo))+1);
if Vin(1).dim(3) ~= nslices
    error('Number of slices differ: %d vs %d.', nslices, Vin(1).dim(3));
end

% Create new header files
Vout  = Vin;
for k=1:nimgo
    Vout(k).fname  = spm_file(Vin(k).fname, 'prefix', job.prefix);
    Vout(k).descrip = 'Slice time correction';
    Vout(k).n = [k 1];
end
        
% Compute shifting amount from reference slice and slice timings
% Compute time difference between the acquisition time of the
% reference slice and the current slice by using slice times
% supplied in sliceorder vector
rtime=SliceT(job.refslice);
shiftamount = (SliceT - rtime)/TR;

% For loop to perform correction slice by slice

do_spm_hb_stc(mask,Vout,Vin,nimgo,nimg,nslices,shiftamount,job.scans) 

for p = 1:numel(job.scans)
    tout(p)=spm_file(job.scans(p),'prefix',job.prefix);
end

fprintf('%-40s: %30s\n','Completed',spm('time'))  %-#

%%-------------------------------------------------------------------------------------------

function do_spm_hb_stc(mask,Vout,Vin,nimgo,nimg,nslices,shiftamount,scans)

fprintf('Reading the data\n')

% Set up [time x voxels] matrix for holding image info
vol = zeros([Vin(1).dim(1:3) nimgo]);
nvol=vol;

if numel(scans)>1
    for m=1:nimgo
        vol(:,:,:,m) = spm_read_vols(Vin(m));
    end
else
    vol = spm_read_vols(Vin);
end

task = sprintf('Correcting acquisition delay: session %d', 1);
spm_progress_bar('Init',nslices,task,'planes complete');
fprintf('Start Slice Time correction\n')

for k=1:nslices

    slices = vol(:,:,k,:);

    slicemask = mask(:,:,k);
    tmp = find(slicemask>0);

    if numel(tmp)>0

        stack  = zeros([nimg numel(tmp)]);

        rslices = reshape(slices,[Vout(1).dim(1)*Vout(1).dim(2) nimgo]);
        mslices=rslices(tmp,:);
        
        % Set up shifting variables
        len     = size(stack,1);
        phi     = zeros(1,len);
        
        % Check if signal is odd or even -- impacts how Phi is reflected
        %  across the Nyquist frequency. Opposite to use in pvwave.
        OffSet  = 0;
        if rem(len,2) ~= 0, OffSet = 1; end
        
        % Phi represents a range of phases up to the Nyquist frequency
        % Shifted phi 1 to right.
        for f = 1:len/2
            phi(f+1) = -1*shiftamount(k)*2*pi/(len/f);
        end
        
        % Mirror phi about the center
        % 1 is added on both sides to reflect Matlab's 1 based indices
        % Offset is opposite to program in pvwave again because indices are 1 based
        phi(len/2+1+1-OffSet:len) = -fliplr(phi(1+1:len/2+OffSet));
            
        % Transform phi to the frequency domain and take the complex transpose
        shifter = [cos(phi) + sin(phi)*sqrt(-1)].';
        shifter = shifter(:,ones(size(stack,2),1)); % Tony's trick
         
        % Extract columns from slices
        stack(1:nimgo,:) = mslices';
            
        % Fill in continous function to avoid edge effects
        for g=1:size(stack,2)
            stack(nimgo+1:end,g) = linspace(stack(nimgo,g),...
            stack(1,g),nimg-nimgo)';
        end
        
        % Shift the columns
        stack = real(ifft(fft(stack,[],1).*shifter,[],1));
            
        % Re-insert shifted columns
        rslices(tmp,:) = stack(1:nimgo,:)';
        newslices = reshape(rslices,[Vin(1).dim(1:2) nimgo]);

        nvol(:,:,k,:) = newslices;
    end

    spm_progress_bar('Set',k);
    %fprintf(['Done slice timme correction for slice ' num2str(k) '\n'])
end

spm_progress_bar('Clear');
    
Vout = MEHB_write_vol_4d(Vout,nvol);