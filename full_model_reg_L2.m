function data_fit=full_model_reg_L2(dc_data,acq_par,OEF_lookup)
% Analyse data

[acq_par.x_ax,acq_par.y_ax,acq_par.z_ax,acq_par.data_points]=size(dc_data.echo1_data);

OEF_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);
CVR_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);
k_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);
flow_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);
D_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);
resid_est=zeros(acq_par.x_ax,acq_par.y_ax,acq_par.z_ax);

for selected_slice=1:acq_par.z_ax
%for selected_slice=8:8

disp(['processing slice ' int2str(selected_slice)]);
disp(' ');

% select slices
te1_slice(:,:,:)=squeeze(dc_data.echo1_data(:,:,selected_slice,:));
te2_slice(:,:,:)=squeeze(dc_data.echo2_data(:,:,selected_slice,:));
M0_slice(:,:)=squeeze(dc_data.M0_3D(:,:,selected_slice));
D_prior_slice(:,:)=squeeze(dc_data.D_prior(:,:,selected_slice));


[OEF_est(:,:,selected_slice),CVR_est(:,:,selected_slice),k_est(:,:,selected_slice),flow_est(:,:,selected_slice),D_est(:,:,selected_slice),resid_est(:,:,selected_slice)] = estimate_phys_params(acq_par,OEF_lookup,D_prior_slice,te1_slice,te2_slice,M0_slice,selected_slice);

end

data_fit.OEF0=OEF_est;
data_fit.CVR=CVR_est;
data_fit.k=k_est;
data_fit.CBF0=flow_est;
data_fit.D=D_est;
data_fit.resid=resid_est;

end

         
   

function [OEF_est,CVR_est,k_est,flow_est,D_est,resid_est] = estimate_phys_params(acq_par,OEF_lookup,D_prior_slice,te1_slice,te2_slice,M0_slice,selected_slice)


filter_length=acq_par.data_points-2;
% create HP filter
cut=acq_par.cut/acq_par.TR;
sigN2=(cut/sqrt(2))^2;
K=toeplitz(1/sqrt(2*pi*sigN2)*exp(-[0:(filter_length -1)].^2/(2*sigN2)));
K=spdiags(1./sum(K')',0,filter_length,filter_length)*K;
H=zeros(filter_length,filter_length); %smoothing matrix
X=[ones(filter_length,1) (1:filter_length)'];
for k=1:filter_length
    W=diag(K(k,:));
    Hat=X*pinv(W*X)*W;
    H(k,:)=Hat(k,:);
end

%HPfilt is the filtering matrix used to premultiply the data and the design by
HPfilt=eye(filter_length)-H;

PLD=acq_par.PLD+acq_par.sldelay*(selected_slice-1); %add slice delay (seconds)


%%

%  [ CVR     D      k    flow (ml/100g/min]

x0=[0        0       0      0  ];


%%
CVR_est=zeros(acq_par.x_ax,acq_par.y_ax);
OEF_est=zeros(acq_par.x_ax,acq_par.y_ax);
k_est=zeros(acq_par.x_ax,acq_par.y_ax);
flow_est=zeros(acq_par.x_ax,acq_par.y_ax);
D_est=zeros(acq_par.x_ax,acq_par.y_ax);
resid_est=zeros(acq_par.x_ax,acq_par.y_ax);


current_x=1;
% 
for i=1:acq_par.x_ax
    for j=1:acq_par.y_ax

% for i=16:16
%     for j=40:40

        
        if i>current_x %output progress of anlaysis to screen
            current_x=i;
            disp(['analysing row ' int2str(i)]);
        end
                
          if D_prior_slice(i,j)>0
             
                    te1_in=squeeze(te1_slice(i,j,:));
                    te2_in=squeeze(te2_slice(i,j,:));
                    M0b=M0_slice(i,j); 
                    D_prior=D_prior_slice(i,j);
                    
                    %CBF data
                    echo1_in=te1_in(2:end-1); %already processed data in main function

                    BOLD_data=((te2_in(2:acq_par.data_points-1)+((te2_in(1:acq_par.data_points-2)+te2_in(3:acq_par.data_points))./2)))./2; %average
                    echo2_in=HPfilt*BOLD_data+mean(BOLD_data(1:20));
                    
                    [estimates,fval] = nested_BOLD_ASL_solver(OEF_lookup,acq_par,D_prior,echo1_in,echo2_in,M0b,PLD,HPfilt,x0); %Uses a nested function to pass extra parameters to the minimization routine
            
                    CVR_est(i,j)=real((estimates(1).*(1.15))+4);
                    flow_est(i,j)=real((estimates(4).*(34))+60);
                    
                    D_weight=0.05;
                    D_centre=D_prior;
                    
                    D_est(i,j) = real((estimates(2)*(D_weight))+D_centre);
        
                    flow_est(ge(flow_est,170))=170;
                    flow_est(le(flow_est,1))=1;
                    D_est(ge(D_est,0.28))=0.28;
                    D_est(le(D_est,0.002))=0.002;
  
                    OEF_est(i,j) = lininterp2(OEF_lookup.CBF_HR, OEF_lookup.D_HR, OEF_lookup.OEF_2D, flow_est(i,j), D_est(i,j));

                    k_est(i,j)=real(estimates(3))+0.05;

                    %% residul calculation (excluding regularization parameters
                    
                    CaO20 = mean(calc_CaO2(acq_par.cap_arterial(1:20),acq_par.Hb));
                    SvO2=((CaO20*(1-OEF_est(i,j)))/(1.34*acq_par.Hb));
                    
                    [ echo1_ts,echo2_ts] = rebuild_DEXI_PCASL(acq_par.cap_arterial,acq_par.oxic_arterial,flow_est(i,j), M0b, SvO2, CVR_est(i,j), k_est(i,j) ,acq_par.TE2, acq_par.Tag_Dur, PLD, acq_par.Hb);
                    echo1_est=echo1_ts(2:end-1);
                    echo1_est=echo1_est';

                    %HP filter echo 2 estimate 

                    BOLD_data=((echo2_ts(2:end-1)).*mean(echo2_in(1:20)))'; %dBOLD/BOLDo * BOLDo = dBOLD
                    echo2_est=HPfilt*BOLD_data;
                    echo2_est=echo2_est+mean(echo2_in(1:20))-mean(echo2_est(1:20));

                    resid_est(i,j)=norm([0.5.*(echo1_est-echo1_in)./500; (1-0.5).*(echo2_est-echo2_in)./500]);
                    
                    
         else
                    CVR_est(i,j)=nan;
                    D_est(i,j)=nan;
                    OEF_est(i,j)=nan;        
                    k_est(i,j)=nan;
                    flow_est(i,j)=nan;
                    resid_est(i,j)=nan;
         end
                    
    end
end
        
        
              

end