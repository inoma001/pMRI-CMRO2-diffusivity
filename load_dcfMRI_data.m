function dc_data=load_dcfMRI_data(data_par,acq_par,dc_data)

echo1_data_nii=load_untouch_nii([data_par.processed_dir 'ASL_diff.nii.gz']);
dc_data.echo1_data=double(echo1_data_nii.img);

echo2_data_nii=load_untouch_nii([data_par.processed_dir 'te2_smoothed.nii.gz']);
dc_data.echo2_data=double(echo2_data_nii.img);

M0_nii = load_untouch_nii([data_par.processed_dir 'M0_smoothed.nii.gz']);
dc_data.M0_3D=double(M0_nii.img); 

gm_pve_nii = load_untouch_nii([data_par.processed_dir 'anat_brain_pve_1_low_res.nii.gz']);
gm_pve=double(gm_pve_nii.img);
dc_data.gm_pve_lr=gm_pve;

sigma_filt=(acq_par.sp_FWHM/2.355)./acq_par.in_plane_res; 
%aysemmtric 3D filter taking into account the slice width
dc_data.gm_pve=imgaussian_asy(gm_pve,sigma_filt,acq_par.slice_res/acq_par.in_plane_res);

D_prior=dc_data.gm_pve.*0.15; % assign D prior
D_prior(le(D_prior,0.05))=0.05; %set lower limit for D prior... needed due to problems with FAST for correctly identifiying grey matter
% more problematic for deep grey matter 

M0_mask=squeeze(dc_data.M0_3D(:,:,:,1));
M0_mask(le(M0_mask,acq_par.M0thr))=0; % needs to be set appropriately for data acquisition
M0_mask(ne(M0_mask,0))=1;
dc_data.D_prior=D_prior.*M0_mask;

dc_data.M0_mask=M0_mask;



end