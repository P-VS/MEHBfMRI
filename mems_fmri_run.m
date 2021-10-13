function [out]=mems_fmri_run(job)

warnstate = warning;
warning off;

spm_defaults;

ne=numel(job.tedat);
nt=numel(job.tedat(1).func);

V=spm_vol(job.tedat(1).func{1});
im=spm_read_vols(V);

dim=size(im);

mask=zeros(dim(1),dim(2),dim(3));
tmp=find(im>0.10*max(im(:)));
mask(tmp)=1;

wf=ones(dim(1),dim(2),dim(3),ne);

out=cell(nt,1);

if job.method==2 || job.method==3
    for ei=1:ne
        tdat=zeros(dim(1),dim(2),dim(3),nt);
        for ti=1:nt
            V=spm_vol(job.tedat(ei).func{ti});
            tdat(:,:,:,ti)=spm_read_vols(V);
        end
        
        wf(:,:,:,ei)=mask.*mean(tdat,4)./std(tdat,0,4);
    end
end

if job.method==1 || job.method==3
    for ei=1:ne
        wf(:,:,:,ei)=wf(:,:,:,ei)*job.tedat(ei).te/1000;
    end    
end

spm_progress_bar('Init',nt,'Combine TE','Volumes done');

for ti=1:nt
    etdat=zeros(dim(1),dim(2),dim(3));
    endat=zeros(dim(1),dim(2),dim(3));
    
    for ei=1:ne
        if job.correg==1
            fprintf(['Time ' num2str(ti) ' echo ' num2str(ei) '\n'])
        
            if ei==1
                V=spm_vol(job.tedat(ei).func{ti});
                tdat=spm_read_vols(V);
            else
                job09.ref               = job.tedat(1).func(ti);
                job09.source            = job.tedat(ei).func(ti);
                job09.other             = job.tedat(ei).func(ti);
                job09.eoptions.cost_fun = 'nmi';
                job09.eoptions.sep      = [4 2];
                job09.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                job09.eoptions.fwhm     = [7 7];
                job09.roptions.interp   = 0;
                job09.roptions.wrap     = [0 0 0];
                job09.roptions.mask     = 0;
                job09.roptions.prefix   = 'r';

                reifmri = spm_run_coreg(job09);
            
                V=spm_vol(reifmri.rfiles{1});
                tdat=spm_read_vols(V);
            
                [path nm ext]=fileparts(reifmri.rfiles{1});
                delete(fullfile(path,[nm '.nii']));
            end
        else
            V=spm_vol(job.tedat(ei).func{ti});
            tdat=spm_read_vols(V);
        end
        
        if ~(job.method==1)
            etdat=etdat+wf(:,:,:,ei).*tdat;
            endat=endat+wf(:,:,:,ei);
        else
            etdat=etdat+wf(:,:,:,ei).*tdat.*tdat;
            endat=endat+wf(:,:,:,ei).*tdat;
        end
    end
    
    edat=etdat./endat;
    
    [path nm ext]=fileparts(job.tedat(1).func{ti});
    fname=['te' nm];
    VI=V;
    VI.fname=fullfile(path,[fname '.nii']);
    VI.descrip='TE combined fMRI data';
    VI=rmfield(VI,'pinfo');
    VI=spm_write_vol(VI,edat);
    
    out(ti)=spm_file(job.tedat(1).func(ti),'prefix','te');
    
    spm_progress_bar('Set',ti);
end

spm_progress_bar('Clear');