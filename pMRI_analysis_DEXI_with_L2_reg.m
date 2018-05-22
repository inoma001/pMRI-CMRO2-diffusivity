function pMRI_analysis_DEXI_with_L2_reg(data_path,pre_process,do_align,PLD,Tag_Dur,TE2,Hb)

fsl_path='/cubric/software/fsl.versions/5.0.9/';

% pre_process = 0 : do none, 1: yes
% do_align = 0: no 1:yes
% PLD (seconds)
% Tag_Dur (seconds)
% TE2 (BOLD TE in ms)
% Hb (g/ml) (typical value 0.15 g/ml)

% data_path is a string specifying the folder containing the
% dcfMRI data and end-tidal traces.

% data_path must contain:
% TE1_4D.nii.gz, TE2_4D.nii.gz M0.nii.gz, mprage.nii.gz, and
% endtidal_traces_overlap.mat 

%%

data_par.processed_dir=data_path;
data_par.results_dir=([data_path 'results_adaptive_L2_O2_Diff/']);
mkdir(data_par.results_dir);

echo1_nii = load_untouch_nii([data_par.processed_dir 'TE1_4D.nii.gz']);

nii_hdr=echo1_nii.hdr;
nii_hdr.dime.datatype=64;
nii_hdr.dime.bitpix=64;
nii_hdr.dime.dim(5)=1;
dc_data.nii.hdr=nii_hdr;

acq_par.TR=echo1_nii.hdr.dime.pixdim(5);
acq_par.in_plane_res=echo1_nii.hdr.dime.pixdim(2);
acq_par.slice_res=echo1_nii.hdr.dime.pixdim(4);
acq_par.sp_FWHM=acq_par.in_plane_res*1.25; %4.3 mm spatial filter (for 3.4375mm in plane resolution)

acq_par.Hb=Hb;
acq_par.TE2=TE2; %ms
acq_par.Tag_Dur=Tag_Dur; %seconds
acq_par.PLD=PLD; %seconds
acq_par.sldelay=0.0367; %slice delay to add to PLD. (for ipat3 and TE = 10ms)
acq_par.M0thr=8000; 
acq_par.cut=420; 

%%

% if selected pre-process the dcfMRI data for analysis (needs to be run at
% least once)
if pre_process > 0
    fsl_pre_pro(acq_par,data_par,fsl_path)
end

% load processed dcfMRI data
dc_data=load_dcfMRI_data(data_par,acq_par,dc_data);

% if selected align end-tidal traces with CBF data (needs to be run at
% least once)
if do_align >0
    align_traces(dc_data,data_par,acq_par)
end

% load aligned end-tidal traces
load([data_par.processed_dir 'endtidal_traces.mat']);
acq_par.PaCO20=PaCO20;
acq_par.cap_arterial=cap_arterial;
acq_par.oxic_arterial=oxic_arterial;

%% load and interpolate O2 diffusivity lookup table
[acq_par,OEF_lookup] = load_oef_lookup(acq_par);

disp('analysing  data');
data_fit=full_model_reg_L2(dc_data,acq_par,OEF_lookup);

% calculate additional maps, save results, and report GM averages
process_results(data_fit,data_par,acq_par,dc_data)

end

