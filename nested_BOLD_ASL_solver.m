function [x,fval] = nested_BOLD_ASL_solver(OEF_lookup,acq_par,D_prior,echo1_in,echo2_in,M0b,PLD,HPfilt,x0)
 

D_weight=0.05; %wich is std of uniform distribution from 0 to 0.18 (so of width 0.18)
D_centre=D_prior;

%TE1/TE2 objective function weighting
obj_weight=0.5; 

options=optimset('display','off'); %set bounds to stay within D lookup table 
[x,fval]=lsqnonlin(@nestedfun,x0,[-3.4,-D_centre/D_weight,-0.05/1,-60/34],[Inf (0.275-D_centre)/D_weight Inf (170-60)/34],options);

% [x,fval]=lsqnonlin(@nestedfun,x0,[-3.4,-D_centre/D_weight,-0.05/1,-60/34],[5.2 (0.275-D_centre)/D_weight Inf (170-60)/34],options);

%%

    % Nested function that computes the objective function
    function y = nestedfun(x)


        
        CVR = (x(1)*(1.15))+4;  %intermediate value good for GM/WM 1.15
        D = real((x(2)*(D_weight))+D_centre);
        K = x(3)+0.05;
        flow0=(x(4)*(34))+60; %flow estimate in ml/100g/min
        

        flow0(ge(flow0,170))=170;
        flow0(le(flow0,1))=1;
        D(ge(D,0.275))=0.275;
        D(le(D,0.0021))=0.0021;
        
        OEF0 = lininterp2(OEF_lookup.CBF_HR, OEF_lookup.D_HR, OEF_lookup.OEF_2D, flow0, D);
        
        
        CaO20 = mean(calc_CaO2(acq_par.oxic_arterial(1:20),acq_par.Hb));
        SvO2=((CaO20*(1-OEF0))/(1.34*acq_par.Hb));
        
        [ echo1_ts,echo2_ts] = rebuild_DEXI_PCASL(acq_par.cap_arterial,acq_par.oxic_arterial,flow0, M0b, SvO2, CVR, K, acq_par.TE2, acq_par.Tag_Dur, PLD, acq_par.Hb);
        echo1_est=echo1_ts(2:end-1);
        echo1_est=echo1_est';

        %HP filter echo 2 estimate 

        BOLD_data=((echo2_ts(2:end-1)).*mean(echo2_in(1:20)))'; %dBOLD/BOLDo * BOLDo = dBOLD
        echo2_est=HPfilt*BOLD_data;
        echo2_est=echo2_est+mean(echo2_in(1:20))-mean(echo2_est(1:20)); %make sure baseline is the same


%% Objective function
obj_norm=norm([obj_weight.*(echo1_est-echo1_in)./500; (1-obj_weight).*(echo2_est-echo2_in)./500]);

reg_factor=obj_norm/15; %sim optimisation 23.5 (25) from optimisation
%for minimal bias and rms error in D and OEF use both for regularisation
y=[obj_weight.*(echo1_est-echo1_in)./500; (1-obj_weight).*(echo2_est-echo2_in)./500; reg_factor*x(2); 2*reg_factor*(OEF0-0.4)];


    end

end
