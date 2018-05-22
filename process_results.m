function process_results(data_fit,data_par,acq_par,dc_data)

% save acq_par to a mat file for a record of acquisition/processing
% paramters
%save([data_par.results_dir 'par.mat'],'acq_par');

% Calculate CMRO2 map...
CaO2 = calc_CaO2(acq_par.oxic_arterial,acq_par.Hb);
CaO20 = mean(CaO2(1:20));
CMRO2(:,:,:)=data_fit.CBF0.*data_fit.OEF0.*CaO20.*39.34;

% Calculate M map
data_fit.k(le(data_fit.k,0))=0;
data_fit.k(ge(data_fit.k,0.8))=0;
beta=1.0;
M(:,:,:)=100*acq_par.TE2.*data_fit.k.*((data_fit.OEF0.*acq_par.Hb).^beta);

gm_mask=dc_data.gm_pve_lr;
gm_mask(lt(gm_mask,0.5))=0;
gm_mask(ne(gm_mask,0))=1;
gm_mask=gm_mask.*dc_data.M0_mask;

%N.B. these functions need NIFTI_TOOLS
write_out_nii_results(dc_data.nii.hdr,data_fit.OEF0,'OEF.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,data_fit.CVR,'CVR.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,data_fit.k,'k_calib.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,data_fit.CBF0,'CBF0.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,data_fit.D,'D.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,CMRO2,'CMRO2.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,M,'M.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,data_fit.resid,'resid.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,dc_data.gm_pve,'gm_pve.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,gm_mask,'gm_mask.nii.gz',data_par.results_dir);
write_out_nii_results(dc_data.nii.hdr,dc_data.D_prior,'D_prior.nii.gz',data_par.results_dir);

%summary GM results

gm_D=data_fit.D.*gm_mask;
gm_D(eq(gm_D,0))=nan;
disp(['GM oxygen diffusivity (ml/100g/mmHg/min) : ' num2str(nanmean(gm_D(:)))]);

gm_k=data_fit.k.*gm_mask;
gm_k(eq(gm_k,0))=nan;
disp(['GM k : ' num2str(nanmean(gm_k(:)))]);

gm_M=M.*gm_mask;
gm_M(eq(gm_M,0))=nan;
disp(['GM M : ' num2str(nanmean(gm_M(:)))]);

gm_resid=data_fit.resid.*gm_mask;
gm_resid(eq(gm_resid,0))=nan;
disp(['GM residuals : ' num2str(nanmean(gm_resid(:)))]);

gm_CBF=data_fit.CBF0.*gm_mask;
gm_CBF(eq(gm_CBF,0))=nan;
disp(['GM CBF (ml/100g/min) : ' num2str(nanmean(gm_CBF(:)))]);

gm_OEF=data_fit.OEF0.*gm_mask;
gm_OEF(eq(gm_OEF,0))=nan;
disp(['GM OEF : ' num2str(nanmean(gm_OEF(:)))]);

gm_CMRO2=CMRO2.*gm_mask;
gm_CMRO2(eq(gm_CMRO2,0))=nan;
disp(['GM CMRO2 (umol O2/100g/min) : ' num2str(nanmean(gm_CMRO2(:)))]);

gm_CVR=data_fit.CVR.*gm_mask;
gm_CVR(eq(gm_CVR,0))=nan;
disp(['GM CVR (%CBF/mmHg CO2) : ' num2str(nanmean(gm_CVR(:)))]);


end