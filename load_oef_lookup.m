function [acq_par,OEF_lookup] = load_oef_lookup(acq_par)

pH = 6.1 + log10(24/(0.03*acq_par.PaCO20)); %Henderson:Hasselbalch equation: pH = 6.1 + log {HCO3−/(0.03 × Pco2)}. (assuming HCO3-=24 from J Lab Clin Med. 2003 Dec;142(6):414-20.)
acq_par.P50 = 221.87-26.37*pH; % J Lab Clin Med. 2003 Dec;142(6):414-20. Regulation of hemoglobin affinity for oxygen by carbonic anhydrase

load('oxygen_model_lookup_Hb_CBF_P50_D_OEF.mat');
%resample to high-res 2D lookup matrix for subject parameters for quick linear interpolation
CBF_HR=linspace(1,180,2000);
D_HR=linspace(0.002,0.28,2000);
OEF_2D=squeeze(interpn(Hb_vect,CBF_vect,P50_vect,D_vect,OEF_mat,acq_par.Hb,CBF_HR,acq_par.P50,D_HR,'spline'));

OEF_lookup.OEF_2D=OEF_2D;
OEF_lookup.CBF_HR=CBF_HR;
OEF_lookup.D_HR=D_HR;

end
