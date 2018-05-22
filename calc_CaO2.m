% u CaO2 data la PaO2

% [CaO2] =  [ml]O2/[ml]bl  
% per trasformare in [CaO2] =  [umol]O2/[ml]bl devo fare CaO2*34.39[umol/ml]  

function CaO2 = calc_CaO2(PaO2,Hb)

%% fixed values
fi = 1.34 ;     % [ml]O2/[g]Hb
%Hb = 0.15 ;       % [g]Hb/[ml]blood
eps = 0.000031;   % [ml]O2/([ml]bl*[mmHg])

%%

SaO2 = calc_SaO2(PaO2);

CaO2 = fi*Hb*SaO2+PaO2*eps;