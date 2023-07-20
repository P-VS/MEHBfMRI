function [out]=mems_fmri_run(job)

warnstate = warning;
warning off;

spm_defaults;

nechoes=numel(job.tedat);
ntime=numel(job.tedat(1).func);

out=cell(ntime,1);

for ie=1:nechoes
    te(ie) = job.tedat(ie).te;

    for it=1:ntime
        Vfunc=spm_vol(job.tedat(ie).func{it});
        if ie==1 && it==1
            voldim = Vfunc.dim;
            tefuncdat = zeros(voldim(1),voldim(2),voldim(3),ntime,nechoes);
        end

        tefuncdat(:,:,:,it,ie) = spm_read_vols(Vfunc);
    end
end

funcdat = zeros(voldim(1),voldim(2),voldim(3),ntime);

switch job.method
    case 0 %Average (wi=1/nechoes)

        funcdat = sum(tefuncdat,5) ./ nechoes;

        cmethod = 'average';

    case 1 %BS (wi=TEi/sum(TEi))

        sum_weights = sum(te,'all');
        
        for ti=1:ntime
            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* te(ne);
                functidat = functidat ./ sum_weights;
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end
        end

        cmethod = 'BS';

    case 2 %tSNR (wi=tSNRi/sum(tSNRi)
        
        mask = MEHB_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);

        meantei = reshape(mean(tefuncdat,4),[voldim(1)*voldim(2)*voldim(3),nechoes]);
        stdtei = reshape(std(tefuncdat,0,4),[voldim(1)*voldim(2)*voldim(3),nechoes]);
        
        weights(mask_ind,:) = meantei(mask_ind,:) ./ stdtei(mask_ind,:);

        weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);
        
        sum_weights = sum(weights,4);
        weights_mask = find(sum_weights>0);

        for ti=1:ntime
            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* weights(:,:,:,ne);
                functidat(weights_mask) = functidat(weights_mask) ./ sum_weights(weights_mask);
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end 
        end

        cmethod = 'tSNR';

    case 3 %tBS (wi=tSNRi * TEi/sum(tSNRi * TEi))

        mask = MEHB_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);

        meantei = reshape(mean(tefuncdat,4),[voldim(1)*voldim(2)*voldim(3),nechoes]);
        stdtei = reshape(std(tefuncdat,0,4),[voldim(1)*voldim(2)*voldim(3),nechoes]);
        
        for ie=1:nechoes
            weights(mask_ind,ie) = repmat(te(ie),numel(mask_ind),1) .* meantei(mask_ind,ie) ./ stdtei(mask_ind,ie);
        end

        weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);
        
        sum_weights = sum(weights,4);
        weights_mask = find(sum_weights>0);

        for ti=1:ntime
            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* weights(:,:,:,ne);
                functidat(weights_mask) = functidat(weights_mask) ./ sum_weights(weights_mask);
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end 
        end

        cmethod = 'tBS';

    case 4 %T2* weigghted (wi(t)=TEi * exp(-TEi/T2*(t)),sum(TEi * exp(-TEi/T2*(t))))

        mask = my_spmbatch_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        %based on https://github.com/jsheunis/fMRwhy/tree/master
        for ti=1:ntime
    
            tifuncdat = reshape(tefuncdat(:,:,:,ti,:),[voldim(1),voldim(2),voldim(3),nechoes]);

            % Create "design matrix" X
            X = horzcat(ones(nechoes,1), -te(:));
        
            t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
        
            Y=[];
            for ne=1:nechoes
                temptefuncdat = reshape(tifuncdat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
                Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
            end
            Y = max(Y, 1e-11);
        
            % Estimate "beta matrix" by solving set of linear equations
            beta_hat = pinv(X) * log(Y);
             % Calculate S0 and T2star from beta estimation
            T2star_fit = beta_hat(2, :); %is R2*
        
            T2star_thresh_min = 1/1500; % arbitrarily chosen, same as tedana
            I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
            T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
        
            t2star(mask_ind) = T2star_fit;
            
            weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);
        
            for ne=1:nechoes
                weights(:,ne) = repmat(-te(ne),voldim(1)*voldim(2)*voldim(3),1) .* t2star(:,1);
                weights(:,ne) = exp(weights(:,ne));
                weights(:,ne) = repmat(te(ne),voldim(1)*voldim(2)*voldim(3),1) .* weights(:,ne);
            end
        
            weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);
            
            sum_weights = sum(weights,4);
            weights_mask = find(sum_weights>0);

            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* weights(:,:,:,ne);
                functidat(weights_mask) = functidat(weights_mask) ./ sum_weights(weights_mask);
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end 
        end

        cmethod = 'T2* weighted';

    case 5 %T2* mapping

        mask = my_spmbatch_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        %based on https://github.com/jsheunis/fMRwhy/tree/master
        for ti=1:ntime
    
            tifuncdat = reshape(tefuncdat(:,:,:,ti,:),[voldim(1),voldim(2),voldim(3),nechoes]);

            % Create "design matrix" X
            X = horzcat(ones(nechoes,1), -te(:));
        
            t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
        
            Y=[];
            for ne=1:nechoes
                temptefuncdat = reshape(tifuncdat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
                Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
            end
            Y = max(Y, 1e-11);
        
            % Estimate "beta matrix" by solving set of linear equations
            beta_hat = pinv(X) * log(Y);
             % Calculate S0 and T2star from beta estimation
            T2star_fit = beta_hat(2, :); %is R2*
        
            T2star_thresh_min = 1/1500; % arbitrarily chosen, same as tedana
            I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
            T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
        
            t2star(mask_ind) = T2star_fit;
            zeromask = (t2star>0);

            t2star(zeromask) = 1 ./ t2star(zeromask);

            funcdat(:,:,:,ti) = reshape(t2star,[voldim(1),voldim(2),voldim(3)]);
        end

        cmethod = 'T2* mapping';
end

for ti=1:ntime
    Vfunc = spm_vol(job.tedat(1).func{ti});
    
    [fpath,fname,~] = fileparts(Vfunc.fname);
    nfname = split(fname,'bold_e');

    Vout = Vfunc;

    Vout.fname = fullfile(fpath,['c' nfname{1} 'bold_' nfname{2} '.nii']);
    Vout.descrip = ['combine echoes - ' cmethod];
    Vout.pinfo = [1,0,0];
    Vout.dt = [spm_type('float32'),spm_platform('bigend')];
    Vout.n = [1 1];

    Vout = MEHB_write_vol_4d(Vout,funcdat(:,:,:,ti));

    out(ti) = {Vout.fname};
end