function align_traces(dc_data,data_par,acq_par)

%% Align end-tidal traces with CBF data
load([data_par.processed_dir 'endtidal_traces_overlap.mat']);
o2_trace=o2_trace.*760./100; % convert % o2 in mmHg
o2_trace_over=o2_trace;
co2_trace_over=co2_trace;
A.data(:,1)=co2_trace;
A.data(:,2)=o2_trace;


%Use ASL data to align traces
%ASL_ts is the spatially smoothed subtraction te1 timeseries
ASL_ts=dc_data.echo1_data(:,:,:,2:end-1);
% create ASL mask using mean signal
mean_ASL=squeeze(nanmean(ASL_ts(:,:,:,1:20),4)); %take mean of first few timepoints to get baseline flow

%%

alpha=0.85;
lambda=0.9;
R1a=1/1.65;
M0b=12000;

tag_dur=acq_par.Tag_Dur;
PLD=acq_par.PLD;

disp('calculating baseline flow');
CBF0=(6000.*mean_ASL.*lambda.*R1a)./ (2.*alpha.*M0b.*(exp(-PLD.*R1a)-exp(-(tag_dur+PLD).*R1a)));

CBF_mask=CBF0;
CBF_mask(le(CBF_mask,30))=0;
CBF_mask(ge(CBF_mask,80))=0;
CBF_mask(ne(CBF_mask,0))=1;


disp('masking time series');
for i=1:length(ASL_ts)
    ASL_ts(:,:,:,i)=ASL_ts(:,:,:,i).*CBF_mask;
end
%%

%create mean timeseries
ASL_ts(le(ASL_ts,0))=nan;
te1_av1(:,:,:)=squeeze(nanmean(ASL_ts,1));
te1_av2(:,:)=squeeze(nanmean(te1_av1,1));
te1_av3=squeeze(nanmean(te1_av2,1));

%convert to perfusion ... to find resting flow
CBF=(6000.*te1_av3.*lambda.*R1a)./(2.*alpha.*M0b.*(exp(-PLD.*R1a)-exp(-(tag_dur+PLD).*R1a)));
flow0=mean(CBF(1:20));

co2_trace_long=A.data(:,1);
o2_trace_long=A.data(:,2);


%% default parameters .. change to get a better fit to the data
  %should ideally allow flow0 and CVR to be fit to the data
CVR=2.5;

% not relevant parameters but needed to use existing code
Hb=acq_par.Hb;
OEF0=0.35;
K=0.15;
TE2=acq_par.TE2;


%% Use cross-correlation with modelled time series to get lag
    
echo1_ts_long = rebuild_DEXI_PCASL(co2_trace_long,o2_trace_long,flow0, M0b, OEF0, CVR, K,TE2, tag_dur, PLD, Hb);
echo1_est_long=echo1_ts_long(2:end-1);
[acor,lag] = xcorr(echo1_est_long,te1_av3);
[~,I] = max(abs(acor));
lagDiff = lag(I)+1;


%%

% can put CVR fitting here to make user assessment easier

%%

x=lagDiff; %initial offset estimate (must be 2 or greater)
x_offset=x
p=length(te1_av3);

if x<2 %make sure data exists
    x=2;
end


%% User check of alignment


while x>0 %user enters zero to quit optimisation of trace alignment
    
    co2_trace=co2_trace_long(x-1:p+x);
    co2_trace=co2_trace-mean(co2_trace(1:20));
    o2_trace=o2_trace_long(x-1:p+x);

    x_is=x;
    
    %use forward model to predict ASL and BOLD timecourse from end-tidals
    echo1_ts = rebuild_DEXI_PCASL(co2_trace,o2_trace,flow0, M0b, OEF0, CVR, K,TE2, tag_dur, PLD, Hb);
    echo1_est=echo1_ts(2:end-1);
       
    figure;
    plot(te1_av3)
    hold on
    plot(echo1_est);
    hold off
    
    prompt = 'Please input increment in number of TRs: +ve <-  : -ve ->  (enter 0 to accept) : ';
    x_input = round(input(prompt));
    if x_input == 0
        x=0;
        close all
    else
        x_offset=x_offset+x_input
        x=x_offset;
        close all
        
  
        
        co2_trace=co2_trace_long(x-1:p+x);
        co2_trace=co2_trace-mean(co2_trace(1:20));
        o2_trace=o2_trace_long(x-1:p+x);
        %use forward model to predict ASL and BOLD timecourse from end-tidals
        echo1_ts = rebuild_DEXI_PCASL(co2_trace,o2_trace,flow0, M0b, OEF0, CVR, K,TE2, tag_dur, PLD, Hb);
        echo1_est=echo1_ts(2:end-1);
        
        
        figure;
        plot(te1_av3)
        hold on
        plot(echo1_est);
        hold off
    end
     

end

cap_arterial=co2_trace';
oxic_arterial=o2_trace';
PaCO20=mean(co2_trace_long(x_is-1:x_is+19));

save([data_par.processed_dir 'endtidal_traces.mat'],'cap_arterial','oxic_arterial','PaCO20');

end

