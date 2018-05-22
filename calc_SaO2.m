
% calcola SaO2 data la PaO2

function SaO2 = calc_SaO2(PaO2)

% Severinghaus'

SaO2 = 1./ ( ( 23400./(PaO2.^3+150.*PaO2) ) +1 ); 

