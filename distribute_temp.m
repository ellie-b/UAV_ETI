%% Calculate temperature from AWS temperature, absorbed radiation, and distributed absorbed radiation

function [dist_temp] = distribute_temp(SWaws,SWabs,T)

Tc = T-273.15;
tf = 9.7721;
Tres_aws = Tc-(SWaws*tf);

dist_temp = (SWabs*tf) + Tres_aws + 273.15;