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

perf_un_nii=load_untouch_nii([processed_dir, 'perf_unscaled.nii.gz']);
perf_un=double(perf_un_nii.img);

perf_vect=sort(perf_un(:),'descend');
max_perf=nanmedian(perf_vect(1:100));

perf_un=perf_un./max_perf;
D_prior=perf_un.*0.15; %assume max D is around 0.15
D_prior(le(D_prior,0))=0;
D_prior(isinf(D_prior))=0;
D_prior(isnan(D_prior))=0;


M0_mask=squeeze(dc_data.M0_3D(:,:,:,1));
M0_mask(le(M0_mask,acq_par.M0thr))=0; % needs to be set appropriately for data acquisition
M0_mask(ne(M0_mask,0))=1;
dc_data.D_prior=D_prior.*M0_mask;

dc_data.M0_mask=M0_mask;



end