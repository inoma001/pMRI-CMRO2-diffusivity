function [ echo1_ts,echo2_ts] = rebuild_DEXI_PCASL(d_PaCO2,PaO2,CBF0, M0b, SvO2, CVR, phys_K,TE2, tag_dur, PLD, Hb)


%% SIGNAL GENERATION


        
%T1 Blood vector related to PaO2
%R1a=(PaO2.*1.59E-4)+0.59;

SaO20 = mean(calc_SaO2(PaO2(1:20)));
One_Y=1-SaO20;
R1a=((PaO2.*1.527E-4)+0.1713*One_Y)+0.5848; %fitting to in-vivo data in 10.1002/mrm.23137 R2=0.992 
T1a=1./R1a;
%using 1-Ya data helps with non-linearity around lower PaO2 values.





flow_ch=(CVR.*d_PaCO2./100)+1;

%%

%calculate BOLD signal change for the second echo

% beta=0.91; %Griffeth 2011 parameters 
% alpha=0.14;

% beta=1.5; %M only makes sense with this value (if trying to convert K)
% alpha=0.18; %0.18

beta=1; %Alberto optimisation 
alpha=0.06;

fi = 1.34 ;                  % [ml]O2/[g]Hb
CaO2 = calc_CaO2(PaO2,Hb);
CaO20 = mean(CaO2(1:20));


%fractional change in deoxyhemoglobin concentration

%frac_deox=(1./flow_ch)-(1./(OEF0*Hb)).*((1./fi).*(CaO2-(1./flow_ch).*CaO20)+Hb.*((1./flow_ch)-1));
%frac_deox=(CaO20*OEF0/(fi*Hb))/(1-(1-OEF0)*CaO20/(fi*Hb)).*(1./flow_ch)+((1-CaO2./(fi*Hb))./(1-(1-OEF0).*CaO20/(fi*Hb)));


frac_deox=(1./flow_ch)-(1./((1-SvO2)*Hb)).*((1./fi).*(CaO2-(1./flow_ch).*CaO20)+Hb.*((1./flow_ch)-1));
dR2star=(phys_K*((1-SvO2)*Hb)^beta).*(flow_ch.^(alpha).*(frac_deox.^beta)-1); %alpha and beta 

%dR2star=(phys_K*(OEF0*Hb)^beta).*(flow_ch.^(alpha).*(frac_deox.^beta)-1); %alpha and beta  


flow=CBF0/6000; %convert ml/100g/min input into ml/ml/s for Woolrich equations
alpha_eff=0.85;

%'Inter-Vendor Reproducibility of Pseudo-Continuous Arterial Spin Labeling
%at 3 Tesla'
% Henri J. M. M. Mutsaerts , Rebecca M. E. Steketee, Dennis F. R. Heijtel, Joost P. A. Kuijer, Matthias J. P. van Osch, Charles B. L. M. Majoie, Marion Smits, Aart J. Nederveen
lambda=0.9;
%Sn_b=((2*alpha_eff*M0b*flow)./(lambda.*R1a)).*flow_ch.*(exp(-PLD.*R1a)-exp(-(tag_dur+PLD).*R1a))
BGS_sf=0.88; %for nominal CBF0 in Effect of background suppression on CBF quantitation in pseudo continuous arterial spin labeling 

Sn_b=((2.*alpha_eff.*BGS_sf.*flow.*flow_ch.*T1a.*M0b).*(1-exp(-tag_dur./T1a)))./(lambda.*exp(PLD./T1a));

%calculate forward model of echo 1 timeseries
echo1_ts=Sn_b;

%calculate forward model of echo 2 timeseries
echo2_ts=-TE2.*dR2star; %dBOLD/BOLDo 



end