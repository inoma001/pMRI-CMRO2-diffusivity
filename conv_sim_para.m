function para = conv_sim_para(trace,tmax,inv_alpha,TR)
%Convolve gamma variate with respiratory trace %this speeds up fitting
%loads rather than using gampdf

alpha=1/inv_alpha;
t=(0:TR/10:TR*(5*inv_alpha+30)); %create time axis, increases in length for broader gamma functions  %higher sample rate for gamma function to 
%capture shape properly
y=(t.^alpha).*(tmax^(-alpha)).*exp(alpha).*exp((-alpha.*t)./tmax);


%convolve with trace
para=conv(y,trace)./sum(y);
para=para(1:length(trace));


end