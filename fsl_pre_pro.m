function fsl_pre_pro(acq_par,data_par,fsl_path)
 
% Perform Motion Correction and Spatial Filtering on echo_1 and echo_2 data.
% Segmentation on T1 and registration with low res data. Transform high
% res mask to native space


%%

    disp('motion correction');
    
    eval(['!' fsl_path 'bin/mcflirt' ' -in ' data_par.processed_dir 'TE1_4D.nii.gz' ' -out ' data_par.processed_dir 'echo1_mc.nii.gz -plots -meanvol']);
    eval(['!' fsl_path 'bin/mcflirt' ' -in ' data_par.processed_dir 'TE2_4D.nii.gz' ' -out ' data_par.processed_dir 'echo2_mc.nii.gz -plots -meanvol']);
    
    %Take the motion corrected  data and calculate the mean across time 
    eval(['!' fsl_path 'bin/fslmaths ' [data_par.processed_dir, 'echo1_mc.nii.gz'] ' -Tmean ' [data_par.processed_dir, 'mean_te1.nii.gz']]);
    eval(['!' fsl_path 'bin/fslmaths ' [data_par.processed_dir, 'echo2_mc.nii.gz'] ' -Tmean ' [data_par.processed_dir, 'mean_te2.nii.gz']]);

    % Perform bet2 on the mean data, use a threshold of .5 %chops off some
    % of the inferior temporal poles but other than that OK
    disp('brain extraction');
    eval(['!' fsl_path 'bin/bet2 ' [data_par.processed_dir, 'mean_te1.nii.gz '] [data_par.processed_dir, 'mean_te1_brain'] ' -f 0.5 -m']);
    eval(['!' fsl_path 'bin/bet2 ' [data_par.processed_dir, 'mean_te2.nii.gz '] [data_par.processed_dir, 'mean_te2_brain'] ' -f 0.5 -m']);

    % Mask the motion corrected data with the mask to create the masked (bet) motion corrected data
    eval(['!' fsl_path 'bin/fslmaths ' [data_par.processed_dir, 'echo1_mc.nii.gz'] ' -mas ' [data_par.processed_dir, 'mean_te2_brain_mask.nii.gz '] [processed_dir, 'echo1_mc_bet.nii.gz']]);
    eval(['!' fsl_path 'bin/fslmaths ' [data_par.processed_dir, 'echo2_mc.nii.gz'] ' -mas ' [data_par.processed_dir, 'mean_te2_brain_mask.nii.gz '] [processed_dir, 'echo2_mc_bet.nii.gz']]);


%%
    disp('processing ASL data');
    % create mean perf data and subtraction timeseries
    echo1_data_nii=load_untouch_nii([data_par.processed_dir 'echo1_mc_bet.nii.gz']); %try not smoothed data
    echo1_data=double(echo1_data_nii.img);

    %CBF data create unscaled perfusion data
    for i=2:length(echo1_data)-1
    
        if mod(i,2) == 0 
            %even
            flow_data(:,:,:,i)=echo1_data(:,:,:,i)-((echo1_data(:,:,:,i-1)+echo1_data(:,:,:,i+1))./2); %control - average of surrounding tags
        else
            %odd
            flow_data(:,:,:,i)=-echo1_data(:,:,:,i)+((echo1_data(:,:,:,i-1)+echo1_data(:,:,:,i+1))./2); %-tag + average of surrounding control images 
        end  
    end
    
    flow_data(:,:,:,1)=zeros; %put in dummy data
    flow_data(:,:,:,end+1)=zeros;
    
    %try gauss filtering of ASL data rather than susan which seems not to have any moderate setting
    %means we can actually keep some CBF contrast!!!
    sigma_filt=(acq_par.sp_FWHM/2.355)./acq_par.in_plane_res;  
    for i=1:length(flow_data)
        flow_data(:,:,:,i)=imgaussian_asy(squeeze(flow_data(:,:,:,i)),sigma_filt,acq_par.slice_res/acq_par.in_plane_res);
    end
     

    nii_hdr=echo1_nii.hdr;
    nii_hdr.dime.datatype=64;
    nii_hdr.dime.bitpix=64;
    flow_nii.nii.hdr=nii_hdr;
    flow_nii.nii.img=flow_data;
   
    save_nii(flow_nii.nii, [data_par.processed_dir, 'ASL_diff.nii.gz']);
    
    perf_est.nii.img=squeeze(mean(flow_data,4));
    
%     nii_hdr=echo1_data_nii.hdr;
%     nii_hdr.dime.datatype=64;
%     nii_hdr.dime.bitpix=64;
%     nii_hdr.dime.dim(5)=1;
%     perf_est.nii.hdr=nii_hdr;
    
    perf_est.nii.hdr=dc_data.nii.hdr;

    save_nii(perf_est.nii, [data_par.processed_dir, 'perf_unscaled.nii.gz']);
    
    %%

    disp('SUSAN spatial filtering being performed');

    %spatial smoothing with susan

    %FEAT sets the SUSAN intensity to being 0.75 * the contrast between the median brain intensity and the background
    %calculate brightness threshold from echo2 

    mean_te2_nii = load_untouch_nii([data_par.processed_dir 'mean_te2_brain.nii.gz']);
    mean_te2=mean_te2_nii.img;
    mean_te2(le(mean_te2,0))=nan;
    bt=0.1*nanmedian(mean_te2(:));

    sigma_filt=sp_FWHM/2.355; % only smooth BOLD data
    
    eval(['!' fsl_path 'bin/susan ' [data_par.processed_dir, 'echo2_mc_bet.nii.gz '] num2str(bt) ' ' sigma_filt ' 3 1 1 ' [data_par.processed_dir, 'mean_te2.nii.gz '] num2str(bt) ' ' [data_par.processed_dir, 'te2_smoothed.nii.gz'] ]);
    eval(['!' fsl_path 'bin/susan ' [data_par.processed_dir, 'M0.nii.gz '] num2str(bt*5) ' ' sigma_filt ' 3 1 1 ' [data_par.processed_dir, 'mean_te2.nii.gz '] num2str(bt*5) ' ' [data_par.processed_dir, 'M0_smoothed.nii.gz'] ]);

    % bet and register mprage

    disp('BET mprage');
    betFrac=0.4; %be careful can get rid of too much GM (need to be more conservative)
    anat_nifti=[data_par.processed_dir, 'mprage.nii.gz'];
    eval(['!' fsl_path 'bin/bet2 ' anat_nifti [' ',data_par.processed_dir,'anat_brain -f '] num2str(betFrac)]);
     
    disp('epi to anat brain BBR registration');
    eval(['!' fsl_path 'bin/epi_reg --epi=' [data_par.processed_dir, 'M0.nii.gz'] ' --t1=' [data_par.processed_dir, 'mprage.nii.gz'] ' --t1brain=' [data_par.processed_dir, 'anat_brain.nii.gz'] ' --noclean --out=' [data_par.processed_dir, 'epi_anat_reg']]);
    eval(['!' fsl_path 'bin/convert_xfm -omat ' data_par.processed_dir 'anat2fmri.mat -inverse ' data_par.processed_dir 'epi_anat_reg.mat']);
    
    %apply transformtion matrix to anatomical data and segmentation masks    
    echo1_nifti_name=[data_par.processed_dir, 'M0.nii.gz'];
    eval(['!' fsl_path 'bin/flirt' ' -in ' data_par.processed_dir 'anat_brain.nii.gz' ' -ref ' echo1_nifti_name ' -out ' data_par.processed_dir 'anat_brain_low_res.nii.gz' ' -init ' data_par.processed_dir 'anat2fmri.mat' ' -applyxfm']);

    eval(['!' fsl_path 'bin/flirt' ' -in ' data_par.processed_dir 'epi_anat_reg_fast_pve_0.nii.gz' ' -ref ' echo1_nifti_name ' -out ' data_par.processed_dir 'anat_brain_pve_0_low_res.nii.gz' ' -init ' data_par.processed_dir 'anat2fmri.mat' ' -applyxfm']);
    eval(['!' fsl_path 'bin/flirt' ' -in ' data_par.processed_dir 'epi_anat_reg_fast_pve_1.nii.gz' ' -ref ' echo1_nifti_name ' -out ' data_par.processed_dir 'anat_brain_pve_1_low_res.nii.gz' ' -init ' data_par.processed_dir 'anat2fmri.mat' ' -applyxfm']);    
    eval(['!' fsl_path 'bin/flirt' ' -in ' data_par.processed_dir 'epi_anat_reg_fast_pve_2.nii.gz' ' -ref ' echo1_nifti_name ' -out ' data_par.processed_dir 'anat_brain_pve_2_low_res.nii.gz' ' -init ' data_par.processed_dir 'anat2fmri.mat' ' -applyxfm']);
    

end