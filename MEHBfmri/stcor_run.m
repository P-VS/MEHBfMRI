function [out]=stcor_run(job)

out=[];

Vin     = spm_vol(job.scans{1});
nslices = Vin(1).dim(3);

SliceT = job.SliceT;

TR=job.TR;
if TR<10
    TR=TR*1000;
end
    
if nslices ~= numel(SliceT)
    error('Mismatch between number of slices and length of ''Slice timings'' vector.');
end

%-Slice timing correction
%==========================================================================
for i=1:numel(job.scans)
    Vin(i)   = spm_vol(job.scans{i});
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
    if isfield(Vout(k),'descrip')
        desc = [Vout(k).descrip ' '];
    else
        desc = '';
    end
    Vout(k).descrip = [desc 'acq-fix ref-slice ' num2str(job.refslice)];
end
Vout = spm_create_vol(Vout);

% Set up [time x voxels] matrix for holding image info
slices = zeros([Vout(1).dim(1:2) nimgo]);
stack  = zeros([nimg Vout(1).dim(1)]);

task = sprintf('Correcting acquisition delay: session %d', 1);
spm_progress_bar('Init',nslices,task,'planes complete');
        
% Compute shifting amount from reference slice and slice timings
% Compute time difference between the acquisition time of the
% reference slice and the current slice by using slice times
% supplied in sliceorder vector
rtime=SliceT(job.refslice);
shiftamount = (SliceT - rtime)*1000/TR;

% For loop to perform correction slice by slice
for k = 1:nslices
        
    % Read in slice data
    B  = spm_matrix([0 0 k]);
    for m=1:nimgo
        slices(:,:,m) = spm_slice_vol(Vin(m),B,Vin(1).dim(1:2),1);
    end
        
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
        
    % Loop over columns
    for i=1:Vout(1).dim(2)
            
        % Extract columns from slices
        stack(1:nimgo,:) = reshape(slices(:,i,:),[Vout(1).dim(1) nimgo])';
            
        % Fill in continous function to avoid edge effects
        for g=1:size(stack,2)
            stack(nimgo+1:end,g) = linspace(stack(nimgo,g),...
            stack(1,g),nimg-nimgo)';
        end
            
        % Shift the columns
        stack = real(ifft(fft(stack,[],1).*shifter,[],1));
            
        % Re-insert shifted columns
        slices(:,i,:) = reshape(stack(1:nimgo,:)',[Vout(1).dim(1) 1 nimgo]);
    end
        
    % Write out the slice for all volumes
    for p = 1:nimgo
        Vout(p) = spm_write_plane(Vout(p),slices(:,:,p),k);
        
        out=[out; spm_file(job.scans{p},'prefix',job.prefix)];
    end
    spm_progress_bar('Set',k);
end
spm_progress_bar('Clear');

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#